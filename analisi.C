//-------c++----------------//
#include <iostream>
#include <string>
#include <sstream>
#include <fstream>
#include <stdlib.h>
#include <stdio.h>
#include <vector>
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
#include <TThread.h>
#include <TROOT.h>
#include <TStyle.h>
#include <TImage.h>
#include <TCanvas.h>
#include <TSystem.h>
#include <TString.h>


//------Custom headers----------------//
#include "src/Analyzer.hpp"
#include "include/ConfigFile.hpp"


void analisi(){

  //ROOT::EnableImplicitMT(6);
  //ROOT::EnableThreadSafety();

  const int reserve_length = 1100;

  ConfigFile cf("beta_config.ini");

  // IN-OUT FILES & TREE (HEADER)
  std::string Filename = cf.Value("HEADER","original_input_filename");
  std::cout << "Analysis of file " << Filename << " started" << endl;
  const char *filename = Filename.c_str();
  TFile *file = TFile::Open(filename);
  TTreeReader myReader("wfm", file);

  std::string outFilename = cf.Value("HEADER","output_filename");
  const char *output_filename = outFilename.c_str();
  TFile *OutputFile = new TFile(output_filename,"recreate");
  TTree *OutTree = new TTree("Analysis","Analysis");
  cout << "\n\n The output file will be: " << outFilename << " \n\n";

  // CHANNELS
  int active_channels = cf.Value("CHANNELS","active_channels");
  int ch_mcp          = cf.Value("CHANNELS","ch_mcp");

  // RANGES
  bool pmax_search_range;
  if(cf.Value("RANGES","search_range") == 0) pmax_search_range = false;
  else                                        pmax_search_range = true;
  float search_range[2]        = {0,0};
  float search_range_final[2]  = {0,0};
  float search_range_global[2] = {0,0};
  search_range_global[0]   = cf.Value("RANGES","pmax_search_range_min");
  search_range_global[1]   = cf.Value("RANGES","pmax_search_range_max");
  int   small_range        = cf.Value("RANGES","small_range");
  int   mcp_range          = cf.Value("RANGES","mcp_range");
  float search_around_pmax = cf.Value("RANGES","small_range_interval");
  float search_around_mcp  = cf.Value("RANGES","mcp_range_interval");

  // PARAMETERS
  int   number_points_gaus_fit = cf.Value("PARAMETERS","number_points_gaus_fit");
  int   n_points_baseline      = cf.Value("PARAMETERS","n_points_baseline");
  int   ADC_conversion         = cf.Value("PARAMETERS","ADC_conversion");
  float ADC_conversion_factor  = cf.Value("PARAMETERS","ADC_conversion_factor");
  const float time_const       = cf.Value("PARAMETERS","time_scalar");
  const float voltage_const    = cf.Value("PARAMETERS","voltage_scalar");
  unsigned int maxIndex        = cf.Value("PARAMETERS","sampling_points");
  float mcp_delay              = cf.Value("PARAMETERS","mcp_delay");

  // TRACKER
  int join_txt_tracker = cf.Value("TRACKER","join_txt_tracker");

  // Open tracker CSV file if needed
  std::ifstream Filein;
  if(join_txt_tracker == 1){
    std::string TF = cf.Value("TRACKER","tracker_filename");
    Filein.open(TF.c_str(), std::ios::in);
    if(!Filein.is_open()) std::cout << "Failed to open tracker file: " << TF << std::endl;
    else                  std::cout << "Opened tracker file: " << TF << std::endl;
  }

  int tracker_offset = cf.Value("TRACKER","tracker_offset");
  double fxtrk1=0, fytrk1=0, fxtrk2=0, fytrk2=0, fchi2trk1=0;
  int fntrk=0;

  // Pre-load channel enable/invert flags once (avoids per-event config lookups)
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

  // Output branches
  std::vector<float> Pmax1;
  std::vector<float> PmaxFit;
  std::vector<float> negPmax1Fit;
  std::vector<float> Tmax1;
  std::vector<float> Tmax1Fit;
  std::vector<float> negTmax1Fit;
  std::vector<float> Area1;
  std::vector<float> Area1_new;
  std::vector<float> Area_NC;
  std::vector<float> Area_fixed_window;
  std::vector<float> RiseTime1Fit;
  std::vector<std::vector<double>> CFD1Fit;
  std::vector<std::vector<double>> WIDTH1;
  std::vector<float> chi2;
  std::vector<float> rms1;
  float x_pos1, y_pos1, x_pos2, y_pos2, chi2_trk;

  Pmax1.reserve(20);
  PmaxFit.reserve(20);
  negPmax1Fit.reserve(20);
  Tmax1.reserve(20);
  Tmax1Fit.reserve(20);
  negTmax1Fit.reserve(20);
  Area1.reserve(20);
  Area1_new.reserve(20);
  Area_NC.reserve(20);
  Area_fixed_window.reserve(20);
  RiseTime1Fit.reserve(20);
  rms1.reserve(20);
  CFD1Fit.reserve(20);
  WIDTH1.reserve(20);
  chi2.reserve(20);

  Analyzer *a1      = new Analyzer();
  Analyzer *a_check = new Analyzer();

  int event;

  OutTree->Branch("event",             &event);
  OutTree->Branch("pmax",              "std::vector<float>",               &Pmax1);
  OutTree->Branch("pmax_fit",          "std::vector<float>",               &PmaxFit);
  OutTree->Branch("negpmax_fit",       "std::vector<float>",               &negPmax1Fit);
  OutTree->Branch("tmax",              "std::vector<float>",               &Tmax1);
  OutTree->Branch("tmax_fit",          "std::vector<float>",               &Tmax1Fit);
  OutTree->Branch("negtmax_fit",       "std::vector<float>",               &negTmax1Fit);
  OutTree->Branch("area",              "std::vector<float>",               &Area1);
  OutTree->Branch("area_new",          "std::vector<float>",               &Area1_new);
  OutTree->Branch("area_nc",           "std::vector<float>",               &Area_NC);
  OutTree->Branch("area_fixed_window", "std::vector<float>",               &Area_fixed_window);
  OutTree->Branch("risetime",          "std::vector<float>",               &RiseTime1Fit);
  OutTree->Branch("cfd",               "std::vector<std::vector<double>>", &CFD1Fit);
  OutTree->Branch("width",             "std::vector<std::vector<double>>", &WIDTH1);
  OutTree->Branch("rms",               "std::vector<float>",               &rms1);
  OutTree->Branch("x_pos1",   &x_pos1);
  OutTree->Branch("y_pos1",   &y_pos1);
  OutTree->Branch("x_pos2",   &x_pos2);
  OutTree->Branch("y_pos2",   &y_pos2);
  OutTree->Branch("chi2",     &chi2);
  OutTree->Branch("chi2_trk", &chi2_trk);

  int j_counter = 0;

  std::vector<TTreeReaderArray<float>> voltageReader1;
  TTreeReaderArray<float> timeReader(myReader, "time");

  for(int ch=0; ch<active_channels; ch++){
    voltageReader1.push_back(TTreeReaderArray<float>(myReader, Form("w%i", ch)));
  }
  voltageReader1.push_back(TTreeReaderArray<float>(myReader, "trg0"));
  voltageReader1.push_back(TTreeReaderArray<float>(myReader, "trg1"));

  std::vector<float> w1_check;
  std::vector<float> t1_check;
  w1_check.reserve(reserve_length);
  t1_check.reserve(reserve_length);

  int enable_channel_1 = 0;
  int invert_channel_1 = 0;

  float max_p_check_plane1 = 0;
  float max_t_check_plane1 = 0;
  float baseline_correction = 0;

  std::pair<float, unsigned int> tp_pair1_small{0.,0};
  std::array<float, 3>           fit_array_small = {0.,0.,0.};
  std::pair<float, float>        tp_pair1_fit_small{0.,0.};

  std::vector<float> w1_inner;
  std::vector<float> t1_inner;
  w1_inner.reserve(reserve_length);
  t1_inner.reserve(reserve_length);

  std::pair<float, unsigned int> tp_pair1{0.,0};
  std::array<float, 3>           fit_array = {0.,0.,0.};
  std::pair<float, float>        tp_pair1_fit{0.,0.};
  std::pair<float, unsigned int> neg_tp_pair1{0.,0};
  std::array<float, 3>           neg_fit_array = {0.,0.,0.};
  std::pair<float, float>        neg_tp_pair1_fit{0.,0.};

  std::vector<double> cf_inner;
  std::vector<double> width_inner;
  cf_inner.reserve(7);
  width_inner.reserve(7);

  int entry_count = 0;

  while(myReader.Next()){  // && j_counter < 200000

    w1_check.clear();
    t1_check.clear();
    w1_inner.clear();
    t1_inner.clear();
    Pmax1.clear();
    PmaxFit.clear();
    negPmax1Fit.clear();
    Tmax1.clear();
    Tmax1Fit.clear();
    negTmax1Fit.clear();
    Area1.clear();
    Area1_new.clear();
    Area_NC.clear();
    Area_fixed_window.clear();
    RiseTime1Fit.clear();
    rms1.clear();
    CFD1Fit.clear();
    WIDTH1.clear();
    chi2.clear();

    // Read tracker position from CSV (same offset logic as the former TrkMerge_wfm)
    if(join_txt_tracker==1 && entry_count > tracker_offset){
      Filein >> fntrk >> fxtrk1 >> fytrk1 >> fxtrk2 >> fytrk2 >> fchi2trk1;
    }
    x_pos1   = (join_txt_tracker==1) ? float(fxtrk1)    : 0;
    y_pos1   = (join_txt_tracker==1) ? float(fytrk1)    : 0;
    x_pos2   = (join_txt_tracker==1) ? float(fxtrk2)    : 0;
    y_pos2   = (join_txt_tracker==1) ? float(fytrk2)    : 0;
    chi2_trk = (join_txt_tracker==1) ? float(fchi2trk1) : 0;


    ///////// RANGE DETERMINATION /////////
    if(small_range==1){

      max_p_check_plane1 = 0;
      max_t_check_plane1 = 0;

      for(int ch_counter=0; ch_counter<active_channels; ch_counter++){

        if(ch_counter == ch_mcp) continue;

        w1_check.clear();
        t1_check.clear();

        enable_channel_1 = enable_flag[ch_counter];
        invert_channel_1 = invert_flag[ch_counter];

        if(enable_channel_1 == 1){

          if(invert_channel_1 == 1){
            for(unsigned int i=0; i<voltageReader1.at(ch_counter).GetSize(); i++){
              w1_check.push_back(ADC_conversion==1 ? float(-voltageReader1.at(ch_counter)[i])*ADC_conversion_factor
                                                   : float(-voltageReader1.at(ch_counter)[i]));
              t1_check.push_back(timeReader[i]);
            }
          }else{
            for(unsigned int i=0; i<voltageReader1.at(ch_counter).GetSize(); i++){
              w1_check.push_back(ADC_conversion==1 ? float(voltageReader1.at(ch_counter)[i])*ADC_conversion_factor
                                                   : float(voltageReader1.at(ch_counter)[i]));
              t1_check.push_back(timeReader[i]);
            }
          }

          *a_check = Analyzer(w1_check, t1_check);
          baseline_correction = a_check->Correct_Baseline(n_points_baseline);

          tp_pair1_small     = a_check->Find_Signal_Maximum(pmax_search_range, search_range_global);
          fit_array_small    = a_check->Pmax_with_GausFit(tp_pair1_small, maxIndex, number_points_gaus_fit);
          tp_pair1_fit_small = std::make_pair(fit_array_small[0], fit_array_small[1]);

          if(tp_pair1_fit_small.first*voltage_const > max_p_check_plane1){
            max_p_check_plane1 = tp_pair1_fit_small.first*voltage_const;
            max_t_check_plane1 = tp_pair1_fit_small.second*time_const;
          }
        }
      }

      if(max_t_check_plane1>search_range_global[0] && max_t_check_plane1<search_range_global[1]){
        search_range[0] = max_t_check_plane1 - search_around_pmax;
        search_range[1] = max_t_check_plane1 + search_around_pmax;
      }else{
        search_range[0] = search_range_global[0];
        search_range[1] = search_range_global[1];
      }

    }else if(mcp_range==1){

      w1_check.clear();
      t1_check.clear();

      enable_channel_1 = enable_flag[ch_mcp];
      invert_channel_1 = invert_flag[ch_mcp];

      if(enable_channel_1 == 1){

        if(invert_channel_1 == 1){
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
        baseline_correction = a_check->Correct_Baseline(n_points_baseline);

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

      // MCP and trigger channels always use the global search range
      if(ch_counter==ch_mcp || ch_counter==active_channels || ch_counter==active_channels+1){
        search_range_final[0] = search_range_global[0];
        search_range_final[1] = search_range_global[1];
      }else{
        search_range_final[0] = search_range[0];
        search_range_final[1] = search_range[1];
      }

      w1_inner.clear();
      t1_inner.clear();

      enable_channel_1 = enable_flag[ch_counter];
      invert_channel_1 = invert_flag[ch_counter];

      if(enable_channel_1 == 1){

        if(invert_channel_1 == 1){
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
          cout << "Voltage or Time vector less than " << maxIndex << " entries. Skipping whole event" << endl;
          continue;
        }

        if(w1_inner.size() != t1_inner.size()){
          cout << "Different number of entries in Voltage and Time vectors. Skipping whole event" << endl;
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
        Area_NC.push_back(a1->Area_NC(tp_pair1, 5, 5, noise_rms));                                                                         //mV*ns
        Area1_new.push_back(a1->New_Pulse_Area(tp_pair1_fit, tp_pair1.second, "Simpson", search_range_final)*voltage_const*time_const);    //mV*ns
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
    }

    event = j_counter;
    OutTree->Fill();

    if(j_counter%10000 == 0) cout << "processed events: " << j_counter << endl;
    j_counter++;
    entry_count++;
  }

  OutTree->Write();
  OutputFile->Write();
  OutputFile->Close();
}


int main(){
  analisi();
  return 0;
}
