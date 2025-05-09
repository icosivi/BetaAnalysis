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
//#include <glib-2.0/glib.h>
//#include <glib-2.0/glib/gprintf.h>
//#include <glib.h>
//#include <glib/gprintf.h>
//#include <gtk/gtk.h>
//#include <unistd.h>
//#include <dirent.h>

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


//------Custom headers----------------//
#include "src/Analyzer.hpp"
#include "include/ConfigFile.hpp"


void analisi( ){

  int n_segments = 10;
  float time_window_size = 50e-9;

  //ROOT::EnableImplicitMT(6);
  //ROOT::EnableThreadSafety();

  //TH1F *pmax_histo = new TH1F("pmax_histo","pmax_histo",100,-1,1);
  const int reserve_length=1100;

  //Config file definition
  ConfigFile cf("beta_config.ini");


  // IN-OUT FILES & TREE (HEADER)
  std::string Filename = cf.Value("HEADER","input_filename");
  std::cout << "Anaysis of file " << Filename << " started" << endl; 
  const char *filename = Filename.c_str();
  TFile *file = TFile::Open(filename);
  TTree *itree = dynamic_cast<TTree*>(file->Get("wfm"));
  TTreeReader myReader("wfm", file);
  
  std::string outFilename = cf.Value("HEADER","output_filename");
  const char *output_filename = outFilename.c_str();
  TFile *OutputFile = new TFile(output_filename,"recreate");
  TTree *OutTree = new TTree("Analysis","Analysis");
  cout<<"\n\n The output file will be: "<< outFilename <<" \n\n";


  // CHANNELS
  int active_channels = cf.Value("CHANNELS","active_channels");
  int active_channel[8] = {0,0,0,0,0,0,0,0};
  active_channel[0] = cf.Value("ACTIVE_CHANNEL", "ch1" );
  active_channel[1] = cf.Value("ACTIVE_CHANNEL", "ch2" );
  active_channel[2] = cf.Value("ACTIVE_CHANNEL", "ch3" );
  active_channel[3] = cf.Value("ACTIVE_CHANNEL", "ch4" );
  active_channel[4] = cf.Value("ACTIVE_CHANNEL", "ch5" );
  active_channel[5] = cf.Value("ACTIVE_CHANNEL", "ch6" );
  active_channel[6] = cf.Value("ACTIVE_CHANNEL", "ch7" );
  active_channel[7] = cf.Value("ACTIVE_CHANNEL", "ch8" );
  int invert_channel[8] = {0,0,0,0,0,0,0,0};
  invert_channel[0] = cf.Value("INVERT_SIGNAL", "ch1" );
  invert_channel[1] = cf.Value("INVERT_SIGNAL", "ch2" );
  invert_channel[2] = cf.Value("INVERT_SIGNAL", "ch3" );
  invert_channel[3] = cf.Value("INVERT_SIGNAL", "ch4" );
  invert_channel[4] = cf.Value("INVERT_SIGNAL", "ch5" );
  invert_channel[5] = cf.Value("INVERT_SIGNAL", "ch6" );
  invert_channel[6] = cf.Value("INVERT_SIGNAL", "ch7" );
  invert_channel[7] = cf.Value("INVERT_SIGNAL", "ch8" );

  // RANGES
  bool pmax_search_range;
  if( cf.Value("RANGES", "search_range") == 0 ) pmax_search_range = false;
  else pmax_search_range = true;
  float search_range[2] = {0,0};
  float search_range_final[2] = {0,0};
  search_range_final[0] = cf.Value("RANGES", "pmax_search_range_min" ) ;
  search_range_final[1] = cf.Value("RANGES", "pmax_search_range_max" ) ;
  

  // PARAMETERS
  int number_points_gaus_fit = cf.Value("PARAMETERS","number_points_gaus_fit");
  int n_points_baseline = cf.Value("PARAMETERS","n_points_baseline");
  const float time_const = cf.Value("PARAMETERS","time_scalar");  
  const float voltage_const = cf.Value("PARAMETERS","voltage_scalar");
  unsigned int maxIndex = cf.Value("PARAMETERS","sampling_points");
  float tot_levels[2] = { float(cf.Value("PARAMETERS","tot_rising")), float(cf.Value("PARAMETERS","tot_falling")) };

  // TRACKER
  int join_txt_tracker = cf.Value("TRACKER", "join_txt_tracker" );

  // POWER SUPPLY
  int ps_channel[4] = {0,0,0,0};
  ps_channel[0] = cf.Value("PS", "ps_channel0" ) ;
  ps_channel[1] = cf.Value("PS", "ps_channel1" ) ;
  ps_channel[2] = cf.Value("PS", "ps_channel2" ) ;
  ps_channel[3] = cf.Value("PS", "ps_channel3" ) ;

  
  std::vector<float> Pmax1;
  std::vector<float> PmaxFit;
  std::vector<float> negPmax1Fit;
  std::vector<float> Tmax1;
  std::vector<float> Tmax1Fit;
  std::vector<float> negTmax1Fit;
  std::vector<float> Area1;
  std::vector<float> UArea1;
  std::vector<float> Area1_new;
  std::vector<float> UArea1_new;
  //std::vector<float> DC_Area1;
  std::vector<float> Area_NC;
  //std::vector<float> Area_NC_pos; 
  std::vector<float> Area_fixed_window;
  std::vector<float> RiseTime1Fit;
  std::vector<float> FallTime1Fit;
  std::vector<float> dVdt1Fit;
  std::vector<float> dVdt1Fit_2080;
  std::vector<std::vector<double>> CFD1Fit;
  std::vector<std::vector<double>> WIDTH1;
  std::vector<float> chi2;
  std::vector<float> t_thr1;
  std::vector<float> tot1;
  std::vector<float> rms1;
  std::vector<std::vector<float>> w1 ; //to be commented for skipping the waveform;
  std::vector<std::vector<float>> t1 ; //to be commented for skipping the waveform;
  //float x_pos1, y_pos1,x_pos2, y_pos2, chi2_trk ;
   std::vector<float> i_current;
  std::vector<float> v_bias;
  
  Pmax1.reserve(20);
  PmaxFit.reserve(20);
  negPmax1Fit.reserve(20);
  Tmax1.reserve(20);
  Tmax1Fit.reserve(20);
  negTmax1Fit.reserve(20);
  Area1.reserve(20);
  UArea1.reserve(20);
  Area1_new.reserve(20);
  UArea1_new.reserve(20);
  //DC_Area1.reserve(20);
  Area_NC.reserve(20);
  //Area_NC_pos.reserve(20);
  Area_fixed_window.reserve(20);
  RiseTime1Fit.reserve(20);
  FallTime1Fit.reserve(20);
  dVdt1Fit.reserve(20);
  dVdt1Fit_2080.reserve(20);
  t_thr1.reserve(20);
  tot1.reserve(20);
  rms1.reserve(20);
  CFD1Fit.reserve(20);
  WIDTH1.reserve(20);
  chi2.reserve(20);
  //w1.reserve(20);//to be commented for skipping the waveform;
  //t1.reserve(20);//to be commented for skipping the waveform;
  i_current.reserve(20);
  v_bias.reserve(20);
  Analyzer *a1=new Analyzer();
  Analyzer *a_check=new Analyzer();
  
  int event;
  //int evt_delta = 0; //for tracker sync

  float timestamp;
  
  OutTree->Branch("event",&event);
  //OutTree->Branch("evt_delta",&evt_delta); //for tracker sync
  //OutTree->Branch("w", "std::vector<std::vector<float>>", &w1);
  //OutTree->Branch("t", "std::vector<std::vector<float>>" ,&t1);
  OutTree->Branch("pmax", "std::vector<float>",&Pmax1);
  OutTree->Branch("pmax_fit", "std::vector<float>",&PmaxFit);
  //OutTree->Branch("negpmax", "std::vector<float>",&negPmax1Fit);
  OutTree->Branch("tmax", "std::vector<float>",&Tmax1);
  OutTree->Branch("tmax_fit", "std::vector<float>",&Tmax1Fit);
  //OutTree->Branch("negtmax", "std::vector<float>",&negTmax1Fit);
  OutTree->Branch("area", "std::vector<float>",&Area1);
  //OutTree->Branch("uarea", "std::vector<float>",&UArea1);
  OutTree->Branch("area_new", "std::vector<float>",&Area1_new);
  //OutTree->Branch("uarea_new", "std::vector<float>",&UArea1_new);
  //OutTree->Branch("dc_area", "std::vector<float>",&DC_Area1);
  OutTree->Branch("area_nc", "std::vector<float>",&Area_NC);
  //OutTree->Branch("area_nc_pos", "std::vector<float>",&Area_NC_pos);
  OutTree->Branch("area_fixed_window", "std::vector<float>",&Area_fixed_window);
  //OutTree->Branch("risetime", "std::vector<float>",&RiseTime1Fit);
  ///OutTree->Branch("falltime", "std::vector<float>",&FallTime1Fit);
  OutTree->Branch("dvdt", "std::vector<float>",&dVdt1Fit);
  OutTree->Branch("dvdt_2080", "std::vector<float>",&dVdt1Fit_2080);
  OutTree->Branch("cfd", "std::vector<std::vector<double>>",&CFD1Fit);
  //OutTree->Branch("width", "std::vector<std::vector<double>>",&WIDTH1);
  //OutTree->Branch("t_thr", "std::vector<float>",&t_thr1);  // time at which a certain thr (in V) is passed
  //OutTree->Branch("tot", "std::vector<float>",&tot1);
  OutTree->Branch("rms", "std::vector<float>",&rms1);
  //OutTree->Branch("x_pos1", &x_pos1);
  //OutTree->Branch("y_pos1", &y_pos1);
  //OutTree->Branch("x_pos2", &x_pos2);
  //OutTree->Branch("y_pos2", &y_pos2);
  OutTree->Branch("chi2", &chi2);
  //OutTree->Branch("chi2_trk", &chi2_trk);
  OutTree->Branch("time",&timestamp);
  OutTree->Branch("I", "std::vector<float>", &i_current);
  OutTree->Branch("V", "std::vector<float>", &v_bias);
      
  int j_counter = 0;
  
  std::vector<TTreeReaderArray<Double32_t>> voltageReader1 ;
  std::vector<TTreeReaderArray<Double32_t>> timeReader1 ;
  
  //TTreeReaderValue<float> x1Reader(myReader, "xtrk1" );
  //TTreeReaderValue<float> y1Reader(myReader, "ytrk1" );
  //TTreeReaderValue<float> x2Reader(myReader, "xtrk2" );
  //TTreeReaderValue<float> y2Reader(myReader, "ytrk2" );
  //TTreeReaderValue<float> chi2Reader(myReader, "chi2trk" );
      
  for(int ch_counter=1; ch_counter<=active_channels; ch_counter++ ){

    if(active_channel[ch_counter-1]==1){

      voltageReader1.push_back(TTreeReaderArray<double>(myReader, Form("w%i",ch_counter) ));  
      timeReader1.push_back(TTreeReaderArray<double>(myReader, Form("t%i",ch_counter) ));

    }    
  }

  TTreeReaderArray<double> currentReader1(myReader,"i_current") ;
  TTreeReaderArray<double> biasReader1(myReader,"v_bias") ;
  TTreeReaderValue<double> tstampReader1(myReader,"i_timestamp") ;

  int ps_total = 0 ;
  
  for(int ps_counter=0; ps_counter<4; ps_counter++){

    if(ps_channel[ps_counter] == 1){ 

      ps_total++ ;

    }

  }

  std::vector<float> w1_check;
  std::vector<float> t1_check;
  w1_check.reserve(reserve_length);
  t1_check.reserve(reserve_length);

  int enable_channel_1 = 0;
  int invert_channel_1 = 0;

  float max_p_check_plane1 = 0;
  float max_t_check_plane1 = 0;

  float baseline_correction = 0;

  std::vector<float> w1_inner;
  std::vector<float> t1_inner;
  w1_inner.reserve(reserve_length);
  t1_inner.reserve(reserve_length);

  std::pair<float, unsigned int> tp_pair1{0.,0};
  std::array<float, 3> fit_array = {0.,0.,0.};
  std::pair<float, float> tp_pair1_fit{0.,0.};
	std::pair<float, unsigned int> neg_tp_pair1{0.,0}; 
  std::array<float, 3> neg_fit_array = {0.,0.,0.};
	std::pair<float, float> neg_tp_pair1_fit{0.,0.};

  std::vector<double> cf_inner ;
	std::vector<double> width_inner ;
  cf_inner.reserve(7);
  width_inner.reserve(7);
  
  while(myReader.Next() ){ //  && j_counter<10000

    timestamp = *tstampReader1;
    i_current.clear();
    v_bias.clear();
    
    for( int ps_counter=0; ps_counter<ps_total; ps_counter++ ){

      //i_current.push_back( *currentReader1.at(ps_counter) ) ;
      //v_bias.push_back( *biasReader1.at(ps_counter) ) ;
      i_current.push_back( currentReader1[ps_counter] ) ;
      v_bias.push_back( biasReader1[ps_counter] ) ;

    }
  
    
    /*if(join_txt_tracker==1){

      x_pos1 = *x1Reader;
      y_pos1 = *y1Reader;
      x_pos2 = *x2Reader;
      y_pos2 = *y2Reader;
      chi2_trk = *chi2Reader;
      
    }else{

      x_pos1 = 0;
      y_pos1 = 0;
      x_pos2 = 0;
      y_pos2 = 0;
      chi2_trk = 0;
    }*/

    for(int ns=0; ns<n_segments; ns++){

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
      UArea1.clear();
      Area1_new.clear();
      Area_fixed_window.clear();
      UArea1_new.clear();
      //DC_Area1.clear();
      Area_NC.clear();
      //Area_NC_pos.clear();
      RiseTime1Fit.clear();
      FallTime1Fit.clear();
      dVdt1Fit.clear();
      dVdt1Fit_2080.clear();
      t_thr1.clear();
      tot1.clear();
      rms1.clear();
      CFD1Fit.clear();
      WIDTH1.clear();
      chi2.clear();
      w1.clear();//to be commented for skipping the waveform;
      t1.clear();//to be commented for skipping the waveform;
    
     if(j_counter>500){
     int active_ch_counter = 0;
     for( int ch_counter=1; ch_counter<=active_channels; ch_counter++ ){
      //for( int ch_counter=0; ch_counter<active_channels; ch_counter++ ){
          
      w1_inner.clear();
      t1_inner.clear();

      //cout<<voltageReader1.at(active_ch_counter).GetSize()<<endl;
  
      /*if(ch_counter < active_channels ){
          
        enable_channel_1 = cf.Value("ACTIVE_CHANNEL", Form("ch%i", ch_counter) );
        invert_channel_1 = cf.Value("INVERT_SIGNAL", Form("ch%i", ch_counter) );
          
      }*/
  
 	    if( active_channel[ch_counter-1]==1 ){
       //if( enable_channel_1 == 1){
  
 	      if( invert_channel[ch_counter-1]==1 ){

          //cout<<ns*(voltageReader1.at(active_ch_counter).GetSize()/n_segments)<<"    "<<(ns+1)*(voltageReader1.at(active_ch_counter).GetSize()/n_segments)<<endl;
    
	        for(unsigned int i=ns*(voltageReader1.at(active_ch_counter).GetSize()/n_segments); i<(ns+1)*(voltageReader1.at(active_ch_counter).GetSize()/n_segments) ; i++){
  
            w1_inner.push_back(-voltageReader1.at(active_ch_counter).At(i));
	 			    t1_inner.push_back(timeReader1.at(active_ch_counter).At(i));
  
 	        }
  
 	      }else{
  
	        for(unsigned int i=0; i<voltageReader1.at(active_ch_counter).GetSize();i++){
  
	    	    w1_inner.push_back(voltageReader1.at(active_ch_counter).At(i));
	 			    t1_inner.push_back(timeReader1.at(active_ch_counter).At(i));
  
 	        }
 	      }
    
 	      if(w1_inner.size()<maxIndex || t1_inner.size()<maxIndex){
  
 	    	  cout<<"Voltage or Time vector less than 1000 entries. Skipping whole event"<<endl;
 	    	  continue;
 	    	
        }
  
 	    	if(w1_inner.size()==0 || t1_inner.size()==0){
  
 	    	  cout<<"Voltage or Time vector empty. Skipping whole event"<<endl;
 	    	  continue;
 	    	
        }
  
 	    	if(w1_inner.size()!= t1_inner.size()){
  
 	    	  cout<<"Different number of entries in Voltage and Time vectors. Skipping whole event"<<endl;
 	    		continue;

 	    	}

        search_range_final[0] = search_range_final[0]+float(ns)*time_window_size;
        search_range_final[1] = search_range_final[1]+float(ns)*time_window_size;
  
	    	*a1=Analyzer( w1_inner, t1_inner );
        baseline_correction = a1->Correct_Baseline(n_points_baseline); // we do not want signals in the first 5 ns, otherwise baseline correction is biased

        for(int i=0; i<int(w1_inner.size()); i++) w1_inner.at(i) = w1_inner.at(i) - baseline_correction ;

        w1.push_back( w1_inner );//to be commented for skipping the waveform;
        t1.push_back( t1_inner );//to be commented for skipping the waveform;
     
        tp_pair1 = a1->Find_Signal_Maximum(pmax_search_range,search_range_final); 
        fit_array = a1->Pmax_with_GausFit(tp_pair1,maxIndex,number_points_gaus_fit);
        tp_pair1_fit = std::make_pair( fit_array[0], fit_array[1] ) ;
        chi2.push_back( fit_array[2] );

	    	neg_tp_pair1 = a1->Find_Negative_Signal_Maximum(pmax_search_range,search_range_final); 
        neg_fit_array = a1->Negative_Pmax_with_GausFit(neg_tp_pair1,maxIndex,number_points_gaus_fit);
	    	neg_tp_pair1_fit = std::make_pair( neg_fit_array[0], neg_fit_array[1] ) ;

	    	
        Pmax1.push_back( tp_pair1.first*voltage_const ) ; //mV
        Tmax1.push_back(  a1->Get_Tmax(tp_pair1)*time_const ) ; //ns
        PmaxFit.push_back( tp_pair1_fit.first*voltage_const ) ; //mV
        Tmax1Fit.push_back(  tp_pair1_fit.second*time_const ) ; //ns

        negPmax1Fit.push_back(  neg_tp_pair1_fit.first*voltage_const ); //mV
	    	negTmax1Fit.push_back(  neg_tp_pair1_fit.second*time_const ) ; //ns

        Area1.push_back(  a1->Find_Pulse_Area(tp_pair1)*voltage_const*time_const ) ; 
        //DC_Area1.push_back( a1->DC_Area(baseline_correction)*voltage_const*time_const ); //mV*ns
        Area_NC.push_back( a1->Area_NC(tp_pair1, 5, 5, a1->Find_Noise(n_points_baseline)*voltage_const*time_const) ); //mV*ns
        //Area_NC_pos.push_back( a1->Area_NC_pos(tp_pair1, 5, 5, a1->Find_Noise(n_points_baseline)*voltage_const) ); //mV*ns
        Area1_new.push_back( a1->New_Pulse_Area(tp_pair1_fit,tp_pair1.second,"Simpson",search_range_final)*voltage_const*time_const ) ;//mV*ns 
        Area_fixed_window.push_back( a1->Pulse_Integration_with_Fixed_Window_Size_with_GausFit(tp_pair1_fit,tp_pair1.second,"Simpson", 1, 1)*voltage_const*time_const ); //mV*ns

        RiseTime1Fit.push_back( a1->Find_Rise_Time_with_GausFit(tp_pair1_fit, tp_pair1.second, 0.1, 0.9)*time_const ) ; //ns
        FallTime1Fit.push_back( a1->Find_Fall_Time_with_GausFit(tp_pair1_fit, tp_pair1.second, 0.1, 0.9)*time_const ) ; //ns
        dVdt1Fit.push_back( a1->Find_Dvdt_with_GausFit(20,0,tp_pair1_fit,tp_pair1.second)*(voltage_const/time_const) ) ;  //mV/ns
	    	dVdt1Fit_2080.push_back( a1->Find_Dvdt2080_with_GausFit(0,tp_pair1_fit,tp_pair1.second)*(voltage_const/time_const) );  //mV/ns

        //std::vector<float> cf_inner ;
	    	//std::vector<float> width_inner ;
        //cf_inner.reserve(7);
        //width_inner.reserve(7);
        cf_inner.clear();
        width_inner.clear();

        for(int jj=0;jj<7;jj++){
  
	    		  cf_inner.push_back( double(a1->Rising_Edge_CFD_Time_with_GausFit(10+jj*10,tp_pair1_fit,tp_pair1.second)*time_const) ) ;
	    		  width_inner.push_back( double((a1->Falling_Edge_CFD_Time_with_GausFit(10+jj*10,tp_pair1_fit,tp_pair1.second)*time_const) - 
                                 (a1->Rising_Edge_CFD_Time_with_GausFit(10+jj*10,tp_pair1_fit,tp_pair1.second)*time_const)) ) ;
  
	    	}

        CFD1Fit.push_back( cf_inner ) ;
	    	WIDTH1.push_back( width_inner ) ;

	      UArea1.push_back( a1->Find_Undershoot_Area(tp_pair1)*voltage_const*time_const ); //mV*ns
        UArea1_new.push_back( a1->New_Undershoot_Area(tp_pair1_fit,neg_tp_pair1_fit, neg_tp_pair1.second,"Simpson",search_range_final)*voltage_const*time_const ) ;//mV*ns
         
	    	tot1.push_back( a1->Find_Time_Over_Threshold(tot_levels[0],tp_pair1,tot_levels[1])*time_const ) ; //ns
	    	rms1.push_back( a1->Find_Noise(n_points_baseline)*voltage_const ) ; //mV

        active_ch_counter++;
  
	    }	
     }
    
      event=j_counter;
      OutTree->Fill();

       /*if(x_pos1>-800){
      
          OutTree->Fill();

       } */
       
      }

      if(j_counter%10000 == 0) cout<<"processed events:"<<j_counter<<endl;
      j_counter++;

    }
  
  }

  OutTree->Write();
  OutputFile->Write();
  OutputFile->Close();
  
}


int main(){

analisi();

return 0;

}
