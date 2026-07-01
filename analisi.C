//-------c++----------------//
#include <iostream>
#include <string>
#include <sstream>
#include <fstream>
#include <stdlib.h>
#include <stdio.h>
#include <vector>
#include <array>
#include <mutex>
#include <atomic>
#include <chrono>
#include <iomanip>
#include <numeric>
#include <functional>

//------ROOT----------------//
#include <TTree.h>
#include <TTreeReader.h>
#include "TTreeReaderValue.h"
#include "TTreeReaderArray.h"
#include <TBranch.h>
#include <TFile.h>
#include <TH1.h>
#include <TF1.h>
#include <TGraph.h>
#include <TROOT.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TSystem.h>
#include <TString.h>
#include "TFileMerger.h"
#include <ROOT/TTreeProcessorMT.hxx>

//------Custom headers----------------//
#include "src/Analyzer.hpp"
#include "include/ConfigFile.hpp"


void analisi(){

  const int reserve_length = 1100;

  ConfigFile cf("beta_config.ini");

  // IN-OUT FILES
  std::string Filename   = cf.Value("HEADER","original_input_filename");
  std::string outFilename = cf.Value("HEADER","output_filename");
  const char *filename        = Filename.c_str();
  const char *output_filename = outFilename.c_str();
  cout << "Analysis of file " << Filename << " started\n";
  cout << "Output file:      " << outFilename << "\n\n";

  // CHANNELS
  const int active_channels = cf.Value("CHANNELS","active_channels");
  const int ch_mcp          = cf.Value("CHANNELS","ch_mcp");

  // RANGES
  const bool  pmax_search_range = (cf.Value("RANGES","search_range") != 0);
  const float srg_min           = cf.Value("RANGES","pmax_search_range_min");
  const float srg_max           = cf.Value("RANGES","pmax_search_range_max");
  const int   small_range       = cf.Value("RANGES","small_range");
  const int   mcp_range_flag    = cf.Value("RANGES","mcp_range");
  const float search_around_pmax = cf.Value("RANGES","small_range_interval");
  const float search_around_mcp  = cf.Value("RANGES","mcp_range_interval");

  // PARAMETERS
  const int   number_points_gaus_fit = cf.Value("PARAMETERS","number_points_gaus_fit");
  const int   n_points_baseline      = cf.Value("PARAMETERS","n_points_baseline");
  const int   ADC_conversion         = cf.Value("PARAMETERS","ADC_conversion");
  const float ADC_conversion_factor  = cf.Value("PARAMETERS","ADC_conversion_factor");
  const float time_const             = cf.Value("PARAMETERS","time_scalar");
  const float voltage_const          = cf.Value("PARAMETERS","voltage_scalar");
  const unsigned int maxIndex        = cf.Value("PARAMETERS","sampling_points");
  const float mcp_delay              = cf.Value("PARAMETERS","mcp_delay");
  const int   n_threads              = cf.Value("PARAMETERS","n_threads");
  const Long64_t max_events          = cf.Value("PARAMETERS","max_events");

  // TRACKER
  const int join_txt_tracker = cf.Value("TRACKER","join_txt_tracker");
  const int tracker_offset   = cf.Value("TRACKER","tracker_offset");

  // Pre-load channel enable/invert flags once
  std::vector<int> enable_flag(active_channels + 2);
  std::vector<int> invert_flag(active_channels + 2);
  for(int ch=0; ch<active_channels; ch++){
    enable_flag[ch] = cf.Value("ACTIVE_CHANNEL", Form("ch%i", ch));
    invert_flag[ch] = cf.Value("INVERT_SIGNAL",  Form("ch%i", ch));
  }
  enable_flag[active_channels]   = cf.Value("ACTIVE_CHANNEL","trg0");
  invert_flag[active_channels]   = cf.Value("INVERT_SIGNAL", "trg0");
  enable_flag[active_channels+1] = cf.Value("ACTIVE_CHANNEL","trg1");
  invert_flag[active_channels+1] = cf.Value("INVERT_SIGNAL", "trg1");

  // Get total number of entries for progress reporting (and tracker preload sizing)
  Long64_t nEntries = 0;
  {
    TFile *tmpf = TFile::Open(filename);
    nEntries = ((TTree*)tmpf->Get("wfm"))->GetEntries();
    tmpf->Close();
  }
  std::cout << "Total events to process: " << nEntries << "\n\n";

  // Pre-load tracker CSV into a vector indexed by global TTree entry number.
  // This eliminates the serial CSV dependency and allows parallel event processing.
  struct TrackerPos { float x1=0, y1=0, x2=0, y2=0, chi2=0; };
  std::vector<TrackerPos> tracker_data;

  if(join_txt_tracker == 1){
    std::string TF = cf.Value("TRACKER","tracker_filename");
    {
      tracker_data.resize(nEntries);
    }
    std::ifstream csv(TF.c_str(), std::ios::in);
    if(!csv.is_open()){
      std::cout << "Failed to open tracker file: " << TF << std::endl;
    }else{
      std::cout << "Pre-loading tracker file: " << TF << std::endl;
      int fntrk; double x1,y1,x2,y2,chi2;
      for(Long64_t i=0; i<(Long64_t)tracker_data.size(); i++){
        if(i > tracker_offset){
          csv >> fntrk >> x1 >> y1 >> x2 >> y2 >> chi2;
          tracker_data[i] = {float(x1), float(y1), float(x2), float(y2), float(chi2)};
        }
      }
      std::cout << "Tracker data loaded (" << tracker_data.size() << " entries)\n\n";
    }
  }

  // Pre-warm TFormula/Cling JIT for "gaus" in single-threaded context.
  // Without this, the first concurrent TF1("f","gaus",...) calls from multiple
  // threads trigger simultaneous Cling JIT compilation, which is not thread-safe
  // and causes heap corruption.
  { TF1 f_pw("f_gaus_prewarm","gaus",0,1); f_pw.AddToGlobalList(false); }

  // Prevent TH1 objects from auto-registering in gROOT's global directory.
  // TH1D("pmax_histo",...) in GausFit functions would register on construction
  // (before SetDirectory(nullptr) is called), causing concurrent name collisions
  // across threads and corruption of ROOT's internal object list.
  TH1::AddDirectory(kFALSE);

  // Enable ROOT implicit multithreading (0 = use all available cores)
  if(n_threads > 0) ROOT::EnableImplicitMT(n_threads);
  else              ROOT::EnableImplicitMT();
  std::cout << "Running with " << ROOT::GetThreadPoolSize() << " threads\n\n";

  // Temp file management: TTreeProcessorMT calls the lambda once per TTree cluster,
  // potentially from different threads. Each invocation writes its own partial file.
  std::vector<std::string> temp_files;
  std::mutex               temp_files_mutex;

  // Base path for temp files (output filename without .root)
  std::string outBase = outFilename;
  if(outBase.size() > 5 && outBase.substr(outBase.size()-5) == ".root")
    outBase.resize(outBase.size()-5);

  // Progress tracking
  std::atomic<Long64_t> events_done{0};
  std::atomic<int>      last_print_pct{-1};
  auto t_start = std::chrono::steady_clock::now();

  ROOT::TTreeProcessorMT tp(filename, "wfm");

  tp.Process([&](TTreeReader &reader){

    // ----- output variables (local to this invocation) -----
    std::vector<float> Pmax1, PmaxFit, negPmax1Fit;
    std::vector<float> Tmax1, Tmax1Fit, negTmax1Fit;
    std::vector<float> Area1, Area1_new, Area_NC, Area_fixed_window;
    std::vector<float> RiseTime1Fit, rms1, chi2;
    std::vector<std::vector<double>> CFD1Fit, WIDTH1;
    float x_pos1=0, y_pos1=0, x_pos2=0, y_pos2=0, chi2_trk=0;
    Long64_t event = 0;

    Pmax1.reserve(20);            PmaxFit.reserve(20);      negPmax1Fit.reserve(20);
    Tmax1.reserve(20);            Tmax1Fit.reserve(20);     negTmax1Fit.reserve(20);
    Area1.reserve(20);            Area1_new.reserve(20);    Area_NC.reserve(20);
    Area_fixed_window.reserve(20); RiseTime1Fit.reserve(20); rms1.reserve(20);
    CFD1Fit.reserve(20);          WIDTH1.reserve(20);        chi2.reserve(20);

    // ----- lazy output file creation (only for clusters with events to write) -----
    // Clusters entirely past max_events would otherwise create thousands of empty files.
    TFile *outFile = nullptr;
    TTree *outTree = nullptr;

    auto ensure_output = [&](){
      if(outFile) return;
      std::string temp_name;
      {
        std::lock_guard<std::mutex> lk(temp_files_mutex);
        temp_name = Form("%s_part%04zu.root", outBase.c_str(), temp_files.size());
        temp_files.push_back(temp_name);
      }
      outFile = TFile::Open(temp_name.c_str(), "RECREATE");
      outTree = new TTree("Analysis","Analysis");
      outTree->Branch("event",             &event);
      outTree->Branch("pmax",              "std::vector<float>",               &Pmax1);
      outTree->Branch("pmax_fit",          "std::vector<float>",               &PmaxFit);
      outTree->Branch("negpmax_fit",       "std::vector<float>",               &negPmax1Fit);
      outTree->Branch("tmax",              "std::vector<float>",               &Tmax1);
      outTree->Branch("tmax_fit",          "std::vector<float>",               &Tmax1Fit);
      outTree->Branch("negtmax_fit",       "std::vector<float>",               &negTmax1Fit);
      outTree->Branch("area",              "std::vector<float>",               &Area1);
      outTree->Branch("area_new",          "std::vector<float>",               &Area1_new);
      outTree->Branch("area_nc",           "std::vector<float>",               &Area_NC);
      outTree->Branch("area_fixed_window", "std::vector<float>",               &Area_fixed_window);
      outTree->Branch("risetime",          "std::vector<float>",               &RiseTime1Fit);
      outTree->Branch("cfd",               "std::vector<std::vector<double>>", &CFD1Fit);
      outTree->Branch("width",             "std::vector<std::vector<double>>", &WIDTH1);
      outTree->Branch("rms",               "std::vector<float>",               &rms1);
      outTree->Branch("x_pos1",   &x_pos1);
      outTree->Branch("y_pos1",   &y_pos1);
      outTree->Branch("x_pos2",   &x_pos2);
      outTree->Branch("y_pos2",   &y_pos2);
      outTree->Branch("chi2",     &chi2);
      outTree->Branch("chi2_trk", &chi2_trk);
    };

    // ----- input readers -----
    TTreeReaderArray<float> timeReader(reader, "time");
    std::vector<TTreeReaderArray<float>> voltageReader1;
    for(int ch=0; ch<active_channels; ch++)
      voltageReader1.push_back(TTreeReaderArray<float>(reader, Form("w%i", ch)));
    voltageReader1.push_back(TTreeReaderArray<float>(reader, "trg0"));
    voltageReader1.push_back(TTreeReaderArray<float>(reader, "trg1"));

    // ----- thread-local analyzer objects -----
    Analyzer *a1      = new Analyzer();
    Analyzer *a_check = new Analyzer();

    // ----- working buffers -----
    std::vector<float> w1_check, t1_check, w1_inner, t1_inner;
    w1_check.reserve(reserve_length);  t1_check.reserve(reserve_length);
    w1_inner.reserve(reserve_length);  t1_inner.reserve(reserve_length);

    std::vector<double> cf_inner, width_inner;
    cf_inner.reserve(7);  width_inner.reserve(7);

    float search_range[2]       = {0,0};
    float search_range_global[2] = {srg_min, srg_max};
    float search_range_final[2] = {0,0};
    float max_p_check=0, max_t_check=0, baseline_correction=0;

    std::pair<float, unsigned int> tp_pair1_small{0.,0};
    std::array<float, 3>           fit_array_small = {0.,0.,0.};
    std::pair<float, float>        tp_pair1_fit_small{0.,0.};
    std::pair<float, unsigned int> tp_pair1{0.,0};
    std::array<float, 3>           fit_array = {0.,0.,0.};
    std::pair<float, float>        tp_pair1_fit{0.,0.};
    std::pair<float, unsigned int> neg_tp_pair1{0.,0};
    std::array<float, 3>           neg_fit_array = {0.,0.,0.};
    std::pair<float, float>        neg_tp_pair1_fit{0.,0.};


    Long64_t cluster_event_count = 0;

    while(reader.Next()){

      const Long64_t globalEntry = reader.GetCurrentEntry();
      if(max_events > 0 && globalEntry >= max_events) break;
      ++cluster_event_count;

      // clear output vectors
      Pmax1.clear();       PmaxFit.clear();     negPmax1Fit.clear();
      Tmax1.clear();       Tmax1Fit.clear();    negTmax1Fit.clear();
      Area1.clear();       Area1_new.clear();   Area_NC.clear();
      Area_fixed_window.clear(); RiseTime1Fit.clear(); rms1.clear();
      CFD1Fit.clear();     WIDTH1.clear();      chi2.clear();

      // Tracker position from pre-loaded vector
      if(join_txt_tracker==1 && globalEntry < (Long64_t)tracker_data.size()){
        x_pos1   = tracker_data[globalEntry].x1;
        y_pos1   = tracker_data[globalEntry].y1;
        x_pos2   = tracker_data[globalEntry].x2;
        y_pos2   = tracker_data[globalEntry].y2;
        chi2_trk = tracker_data[globalEntry].chi2;
      }else{
        x_pos1 = y_pos1 = x_pos2 = y_pos2 = chi2_trk = 0;
      }


      ///////// RANGE DETERMINATION /////////
      if(small_range==1){

        max_p_check = 0;
        max_t_check = 0;

        for(int ch=0; ch<active_channels; ch++){
          if(ch==ch_mcp || enable_flag[ch]!=1) continue;

          w1_check.clear();
          t1_check.clear();

          if(invert_flag[ch]==1){
            for(unsigned int i=0; i<voltageReader1.at(ch).GetSize(); i++){
              w1_check.push_back(ADC_conversion==1 ? float(-voltageReader1.at(ch)[i])*ADC_conversion_factor
                                                   : float(-voltageReader1.at(ch)[i]));
              t1_check.push_back(timeReader[i]);
            }
          }else{
            for(unsigned int i=0; i<voltageReader1.at(ch).GetSize(); i++){
              w1_check.push_back(ADC_conversion==1 ? float(voltageReader1.at(ch)[i])*ADC_conversion_factor
                                                   : float(voltageReader1.at(ch)[i]));
              t1_check.push_back(timeReader[i]);
            }
          }

          *a_check = Analyzer(w1_check, t1_check);
          a_check->Correct_Baseline(n_points_baseline);

          tp_pair1_small     = a_check->Find_Signal_Maximum(pmax_search_range, search_range_global);
          fit_array_small    = a_check->Pmax_with_GausFit(tp_pair1_small, maxIndex, number_points_gaus_fit);
          tp_pair1_fit_small = std::make_pair(fit_array_small[0], fit_array_small[1]);

          if(tp_pair1_fit_small.first*voltage_const > max_p_check){
            max_p_check = tp_pair1_fit_small.first*voltage_const;
            max_t_check = tp_pair1_fit_small.second*time_const;
          }
        }

        if(max_t_check>search_range_global[0] && max_t_check<search_range_global[1]){
          search_range[0] = max_t_check - search_around_pmax;
          search_range[1] = max_t_check + search_around_pmax;
        }else{
          search_range[0] = search_range_global[0];
          search_range[1] = search_range_global[1];
        }

      }else if(mcp_range_flag==1){

        if(enable_flag[ch_mcp]==1){
          w1_check.clear();
          t1_check.clear();

          if(invert_flag[ch_mcp]==1){
            for(unsigned int i=0; i<voltageReader1.at(ch_mcp).GetSize(); i++){
              w1_check.push_back(ADC_conversion==1 ? float(-voltageReader1.at(ch_mcp)[i])*ADC_conversion_factor
                                                   : float(-voltageReader1.at(ch_mcp)[i]));
              t1_check.push_back(timeReader[i]);
            }
          }else{
            for(unsigned int i=0; i<voltageReader1.at(ch_mcp).GetSize(); i++){
              w1_check.push_back(ADC_conversion==1 ? float(voltageReader1.at(ch_mcp)[i])*ADC_conversion_factor
                                                   : float(voltageReader1.at(ch_mcp)[i]));
              t1_check.push_back(timeReader[i]);
            }
          }

          *a_check = Analyzer(w1_check, t1_check);
          a_check->Correct_Baseline(n_points_baseline);
          tp_pair1_small  = a_check->Find_Signal_Maximum(pmax_search_range, search_range_global);
          fit_array_small = a_check->Pmax_with_GausFit(tp_pair1_small, maxIndex, number_points_gaus_fit);
        }

        if(fit_array_small[1]>search_range_global[0] && fit_array_small[1]<search_range_global[1]){
          search_range[0] = fit_array_small[1] - mcp_delay - search_around_mcp;
          search_range[1] = fit_array_small[1] - mcp_delay + search_around_mcp;
        }else{
          search_range[0] = search_range_global[0];
          search_range[1] = search_range_global[1];
        }

      }else{
        search_range[0] = search_range_global[0];
        search_range[1] = search_range_global[1];
      }
      ///////// END OF RANGE DETERMINATION /////////


      for(int ch_counter=0; ch_counter<(active_channels+2); ch_counter++){

        if(ch_counter==ch_mcp || ch_counter==active_channels || ch_counter==active_channels+1){
          search_range_final[0] = search_range_global[0];
          search_range_final[1] = search_range_global[1];
        }else{
          search_range_final[0] = search_range[0];
          search_range_final[1] = search_range[1];
        }

        if(enable_flag[ch_counter] != 1) continue;

        w1_inner.clear();
        t1_inner.clear();

        if(invert_flag[ch_counter]==1){
          for(unsigned int i=0; i<voltageReader1.at(ch_counter).GetSize(); i++){
            w1_inner.push_back(ADC_conversion==1 ? float(-voltageReader1.at(ch_counter)[i])*ADC_conversion_factor
                                                 : float(-voltageReader1.at(ch_counter)[i]));
            t1_inner.push_back(timeReader[i]);
          }
        }else{
          for(unsigned int i=0; i<voltageReader1.at(ch_counter).GetSize(); i++){
            w1_inner.push_back(ADC_conversion==1 ? float(voltageReader1.at(ch_counter)[i])*ADC_conversion_factor
                                                 : float(voltageReader1.at(ch_counter)[i]));
            t1_inner.push_back(timeReader[i]);
          }
        }

        if(w1_inner.size()<maxIndex || t1_inner.size()<maxIndex){
          cout << "Voltage or Time vector less than " << maxIndex
               << " entries. Skipping ch " << ch_counter << " event " << globalEntry << endl;
          continue;
        }

        if(w1_inner.size() != t1_inner.size()){
          cout << "Different number of entries in Voltage and Time vectors."
               << " Skipping ch " << ch_counter << " event " << globalEntry << endl;
          continue;
        }

        *a1 = Analyzer(w1_inner, t1_inner);
        baseline_correction = a1->Correct_Baseline(n_points_baseline);

        for(int i=0; i<int(w1_inner.size()); i++) w1_inner[i] -= baseline_correction;

        tp_pair1     = a1->Find_Signal_Maximum(pmax_search_range, search_range_final);
        fit_array    = a1->Pmax_with_GausFit(tp_pair1, maxIndex, number_points_gaus_fit);
        tp_pair1_fit = std::make_pair(fit_array[0], fit_array[1]);
        chi2.push_back(fit_array[2]);

        neg_tp_pair1     = a1->Find_Negative_Signal_Maximum(pmax_search_range, search_range_final);
        neg_fit_array    = a1->Negative_Pmax_with_GausFit(neg_tp_pair1, maxIndex, number_points_gaus_fit);
        neg_tp_pair1_fit = std::make_pair(neg_fit_array[0], neg_fit_array[1]);

        Pmax1.push_back(tp_pair1.first*voltage_const);          //mV
        Tmax1.push_back(a1->Get_Tmax(tp_pair1)*time_const);     //ns
        PmaxFit.push_back(tp_pair1_fit.first*voltage_const);    //mV
        Tmax1Fit.push_back(tp_pair1_fit.second*time_const);     //ns
        negPmax1Fit.push_back(neg_tp_pair1_fit.first*voltage_const);  //mV
        negTmax1Fit.push_back(neg_tp_pair1_fit.second*time_const);    //ns

        Area1.push_back(a1->Find_Pulse_Area(tp_pair1)*voltage_const*time_const);

        const float noise_rms = a1->Find_Noise(n_points_baseline)*voltage_const;
        Area_NC.push_back(a1->Area_NC(tp_pair1, 5, 5, noise_rms));                                                                          //mV*ns
        Area1_new.push_back(a1->New_Pulse_Area(tp_pair1_fit, tp_pair1.second, "Simpson", search_range_final)*voltage_const*time_const);     //mV*ns
        Area_fixed_window.push_back(a1->Pulse_Integration_with_Fixed_Window_Size_with_GausFit(tp_pair1_fit, tp_pair1.second, "Simpson", 1, 1)*voltage_const*time_const); //mV*ns

        RiseTime1Fit.push_back(a1->Find_Rise_Time_with_GausFit(tp_pair1_fit, tp_pair1.second, 0.1, 0.9)*time_const); //ns

        cf_inner.clear();
        width_inner.clear();
        for(int jj=0; jj<7; jj++){
          const float t_rising = a1->Rising_Edge_CFD_Time_with_GausFit(10+jj*10, tp_pair1_fit, tp_pair1.second)*time_const;
          cf_inner.push_back(double(t_rising));
          width_inner.push_back(double(a1->Falling_Edge_CFD_Time_with_GausFit(10+jj*10, tp_pair1_fit, tp_pair1.second)*time_const - t_rising));
        }
        CFD1Fit.push_back(cf_inner);
        WIDTH1.push_back(width_inner);

        rms1.push_back(noise_rms); //mV
      }

      ensure_output();
      event = globalEntry;
      outTree->Fill();

    } // end while(reader.Next())

    // Progress report at end of each cluster (thread-safe, prints every 5%)
    if(cluster_event_count > 0){
      Long64_t done = events_done.fetch_add(cluster_event_count) + cluster_event_count;
      int pct = (nEntries > 0) ? int(100.0 * done / nEntries) : 0;
      int prev = last_print_pct.exchange(pct);
      if(pct > prev && pct % 5 == 0){
        auto now  = std::chrono::steady_clock::now();
        double elapsed = std::chrono::duration<double>(now - t_start).count();
        double rate    = (elapsed > 0) ? done / elapsed : 0;
        double eta_s   = (rate > 0)    ? (nEntries - done) / rate : -1;
        std::lock_guard<std::mutex> lk(temp_files_mutex);
        std::cout << "[" << std::setw(3) << pct << "%] "
                  << done << "/" << nEntries << " ev"
                  << "  rate: " << std::fixed << std::setprecision(0) << rate << " ev/s"
                  << "  ETA: ";
        if(eta_s >= 0){
          int h = int(eta_s)/3600, m = (int(eta_s)%3600)/60, s = int(eta_s)%60;
          if(h > 0) std::cout << h << "h ";
          if(h > 0 || m > 0) std::cout << m << "m ";
          std::cout << s << "s";
        } else {
          std::cout << "?";
        }
        std::cout << "\n" << std::flush;
      }
    } // end if(cluster_event_count > 0)

    if(outFile){
      outFile->Write();
      outFile->Close();
      delete outFile;
    }

    delete a1;
    delete a_check;

  }); // end tp.Process


  // Merge all partial files into the final output
  std::cout << "\nMerging " << temp_files.size() << " partial files...\n";
  TFileMerger merger(kTRUE);
  merger.OutputFile(output_filename);
  for(auto &f : temp_files) merger.AddFile(f.c_str());
  Bool_t merge_ok = merger.Merge();

  if(merge_ok){
    for(auto &f : temp_files) gSystem->Unlink(f.c_str());
    std::cout << "Done. Output: " << outFilename << "\n";
  } else {
    std::cerr << "ERROR: merge failed! Partial files kept in:\n";
    for(auto &f : temp_files) std::cerr << "  " << f << "\n";
  }
}


int main(){
  analisi();
  return 0;
}
