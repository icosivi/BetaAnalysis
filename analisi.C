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

  //ROOT::EnableImplicitMT(6);
  //ROOT::EnableThreadSafety();

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
  int ch_mcp = cf.Value("CHANNELS", "ch_mcp" );

  // RANGES
  bool pmax_search_range;
  if( cf.Value("RANGES", "search_range") == 0 ) pmax_search_range = false;
  else pmax_search_range = true;
  double search_range[2] = {0,0};
  double search_range_final[2] = {0,0};
  double search_range_global[2] = {0,0};
  search_range_global[0] = cf.Value("RANGES", "pmax_search_range_min" ) ;
  search_range_global[1] = cf.Value("RANGES", "pmax_search_range_max" ) ;
  int small_range = cf.Value("RANGES", "small_range" );
  int mcp_range = cf.Value("RANGES", "mcp_range" );
  double search_around_pmax = cf.Value("RANGES", "small_range_interval" );
  double search_around_mcp = cf.Value("RANGES", "mcp_range_interval" );

  // PARAMETERS
  int number_points_gaus_fit = cf.Value("PARAMETERS","number_points_gaus_fit");
  int n_points_baseline = cf.Value("PARAMETERS","n_points_baseline");
  int ADC_conversion = cf.Value("PARAMETERS","ADC_conversion");
  double ADC_conversion_factor = cf.Value("PARAMETERS","ADC_conversion_factor");
  double temporal_bin_width = cf.Value("PARAMETERS","temporal_bin_width"); //0.2; 0.0488;
  const double time_const = cf.Value("PARAMETERS","time_scalar");  
  const double voltage_const = cf.Value("PARAMETERS","voltage_scalar");
  unsigned int maxIndex = cf.Value("PARAMETERS","sampling_points");
  double tot_levels[2] = { cf.Value("PARAMETERS","tot_rising"), cf.Value("PARAMETERS","tot_falling") };
  double mcp_delay = cf.Value("PARAMETERS","mcp_delay");

  // TRACKER
  int join_txt_tracker = cf.Value("TRACKER", "join_txt_tracker" );

  
  std::vector<double> Pmax1;
  std::vector<double> PmaxFit;
  std::vector<double> negPmax1Fit;
  std::vector<double> Tmax1;
  std::vector<double> Tmax1Fit;
  std::vector<double> negTmax1Fit;
  std::vector<double> Area1;
  std::vector<double> UArea1;
  std::vector<double> Area1_new;
  std::vector<double> UArea1_new;
  std::vector<double> DC_Area1;
  std::vector<double> Area_NC;
  std::vector<double> Area_NC_pos; 
  std::vector<double> Area_fixed_window;
  std::vector<double> RiseTime1Fit;
  std::vector<double> FallTime1Fit;
  std::vector<double> dVdt1Fit;
  std::vector<double> dVdt1Fit_2080;
  std::vector<std::vector<double>> CFD1Fit;
  std::vector<std::vector<double>> WIDTH1;
  std::vector<double> chi2;
  std::vector<double> t_thr1;
  std::vector<double> tot1;
  std::vector<double> rms1;
  std::vector<std::vector<double>> w1 ; //to be commented for skipping the waveform;
  std::vector<std::vector<double>> t1 ; //to be commented for skipping the waveform;
  double x_pos1, y_pos1,x_pos2, y_pos2, chi2_trk ;
  
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
  DC_Area1.reserve(20);
  Area_NC.reserve(20);
  Area_NC_pos.reserve(20);
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
  Analyzer *a1=new Analyzer();
  Analyzer *a_check=new Analyzer();
  
  int event;
  //int evt_delta = 0; //for tracker sync
  
  OutTree->Branch("event",&event);
  //OutTree->Branch("evt_delta",&evt_delta); //for tracker sync
  //OutTree->Branch("w", "std::vector<std::vector<double>>", &w1);
  //OutTree->Branch("t", "std::vector<std::vector<double>>" ,&t1);
  OutTree->Branch("pmax", "std::vector<double>",&Pmax1);
  OutTree->Branch("pmax_fit", "std::vector<double>",&PmaxFit);
  //OutTree->Branch("negpmax", "std::vector<double>",&negPmax1Fit);
  OutTree->Branch("tmax", "std::vector<double>",&Tmax1);
  OutTree->Branch("tmax_fit", "std::vector<double>",&Tmax1Fit);
  //OutTree->Branch("negtmax", "std::vector<double>",&negTmax1Fit);
  OutTree->Branch("area", "std::vector<double>",&Area1);
  //OutTree->Branch("uarea", "std::vector<double>",&UArea1);
  OutTree->Branch("area_new", "std::vector<double>",&Area1_new);
  //OutTree->Branch("uarea_new", "std::vector<double>",&UArea1_new);
  //OutTree->Branch("dc_area", "std::vector<double>",&DC_Area1);
  OutTree->Branch("area_nc", "std::vector<double>",&Area_NC);
  OutTree->Branch("area_nc_pos", "std::vector<double>",&Area_NC_pos);
  OutTree->Branch("area_fixed_window", "std::vector<double>",&Area_fixed_window);
  //OutTree->Branch("risetime", "std::vector<double>",&RiseTime1Fit);
  ///OutTree->Branch("falltime", "std::vector<double>",&FallTime1Fit);
  OutTree->Branch("dvdt", "std::vector<double>",&dVdt1Fit);
  OutTree->Branch("dvdt_2080", "std::vector<double>",&dVdt1Fit_2080);
  OutTree->Branch("cfd", "std::vector<std::vector<double>>",&CFD1Fit);
  //OutTree->Branch("width", "std::vector<std::vector<double>>",&WIDTH1);
  //OutTree->Branch("t_thr", "std::vector<double>",&t_thr1);  // time at which a certain thr (in V) is passed
  //OutTree->Branch("tot", "std::vector<double>",&tot1);
  OutTree->Branch("rms", "std::vector<double>",&rms1);
  OutTree->Branch("x_pos1", &x_pos1);
  OutTree->Branch("y_pos1", &y_pos1);
  OutTree->Branch("x_pos2", &x_pos2);
  OutTree->Branch("y_pos2", &y_pos2);
  OutTree->Branch("chi2", &chi2);
  OutTree->Branch("chi2_trk", &chi2_trk);
      
  int j_counter = 0;
  
  std::vector<TTreeReaderArray<Double32_t>> voltageReader1 ;
  
  TTreeReaderValue<float> x1Reader(myReader, "xtrk1" );
  TTreeReaderValue<float> y1Reader(myReader, "ytrk1" );
  TTreeReaderValue<float> x2Reader(myReader, "xtrk2" );
  TTreeReaderValue<float> y2Reader(myReader, "ytrk2" );
  TTreeReaderValue<float> chi2Reader(myReader, "chi2trk" );
      
  for(int ch_counter=0; ch_counter<active_channels; ch_counter++ ){
  
    voltageReader1.push_back(TTreeReaderArray<Double32_t>(myReader, Form("w%i",ch_counter) ));
  
  }
  
  voltageReader1.push_back(TTreeReaderArray<Double32_t>(myReader, "trg0" )); 
  voltageReader1.push_back(TTreeReaderArray<Double32_t>(myReader, "trg1" ));

  std::vector<double> w1_check;
  std::vector<double> t1_check;
  w1_check.reserve(221560);
  t1_check.reserve(221560);

  int enable_channel_1 = 0;
  int invert_channel_1 = 0;

  double max_p_check_plane1 = 0;
  double max_t_check_plane1 = 0;

  double baseline_correction = 0;

  std::pair<double, unsigned int> tp_pair1_small{0.,0}; 
  std::array<double, 3> fit_array_small = {0.,0.,0.};
  std::pair<double, double> tp_pair1_fit_small{0.,0.};

  std::vector<double> w1_inner;
  std::vector<double> t1_inner;
  w1_inner.reserve(221560);
  t1_inner.reserve(221560);

  std::pair<double, unsigned int> tp_pair1{0.,0};
  std::array<double, 3> fit_array = {0.,0.,0.};
  std::pair<double, double> tp_pair1_fit{0.,0.};
	std::pair<double, unsigned int> neg_tp_pair1{0.,0}; 
  std::array<double, 3> neg_fit_array = {0.,0.,0.};
	std::pair<double, double> neg_tp_pair1_fit{0.,0.};

  std::vector<double> cf_inner ;
	std::vector<double> width_inner ;
  cf_inner.reserve(7);
  width_inner.reserve(7);
  
  
  while(myReader.Next() && j_counter<100000 ){ //  && j_counter<10000

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
    DC_Area1.clear();
    Area_NC.clear();
    Area_NC_pos.clear();
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
  
    
    if(join_txt_tracker==1){

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
    }
    
    
    //int enable_channel_1 = 0;
    //int invert_channel_1 = 0;

    
    ///////// BEGINNING OF "SMALL-RANGE" PART /////////
    if(small_range==1){
    
      //double max_p_check_plane1 = 0;
      //double max_t_check_plane1 = 0;
      max_p_check_plane1 = 0;
      max_t_check_plane1 = 0;
  
      for( int ch_counter=0; ch_counter<active_channels; ch_counter++ ){

       if(ch_counter!=ch_mcp){
            
        //std::vector<double> w1_check;
        //std::vector<double> t1_check;
        //w1_check.reserve(221560);
        //t1_check.reserve(221560);
        w1_check.clear();
        t1_check.clear();
          
        enable_channel_1 = cf.Value("ACTIVE_CHANNEL", Form("ch%i", ch_counter) );
        invert_channel_1 = cf.Value("INVERT_SIGNAL", Form("ch%i", ch_counter) );
      
          
        if( enable_channel_1 == 1){
      
          if( invert_channel_1 == 1 ){
        
            for(unsigned int i=0; i<voltageReader1.at(ch_counter).GetSize();i++){
      
              if(ADC_conversion==1) w1_check.push_back( double(-voltageReader1.at(ch_counter).At(i))*ADC_conversion_factor );
              else w1_check.push_back( double(-voltageReader1.at(ch_counter).At(i)) );
              t1_check.push_back( double(i)*temporal_bin_width ); 
      
            }
      
          }else{
      
            for(unsigned int i=0; i<voltageReader1.at(ch_counter).GetSize();i++){
      
              if(ADC_conversion==1) w1_check.push_back( double(voltageReader1.at(ch_counter).At(i))*ADC_conversion_factor );
              else w1_check.push_back( double(voltageReader1.at(ch_counter).At(i)) );
              t1_check.push_back( double(i)*temporal_bin_width );
      
            }
          }
      
          *a_check=Analyzer( w1_check, t1_check );
          //double baseline_correction = a_check->Correct_Baseline(n_points_baseline); 
          baseline_correction = a_check->Correct_Baseline(n_points_baseline);
      
          //std::pair<double, unsigned int> tp_pair1_small = a_check->Find_Signal_Maximum(pmax_search_range,search_range_global); 
          //std::array<double, 3> fit_array_small = a_check->Pmax_with_GausFit(tp_pair1_small,maxIndex,7);  
          //std::pair<double, double> tp_pair1_fit_small{0.,0.}; 
          tp_pair1_small = a_check->Find_Signal_Maximum(pmax_search_range,search_range_global); 
          fit_array_small = a_check->Pmax_with_GausFit(tp_pair1_small,maxIndex,number_points_gaus_fit);
          tp_pair1_fit_small = std::make_pair( fit_array_small[0], fit_array_small[1] ) ;
    
          if( tp_pair1_fit_small.first*voltage_const > max_p_check_plane1){
    
            max_p_check_plane1 = tp_pair1_fit_small.first*voltage_const;
            max_t_check_plane1 = tp_pair1_fit_small.second*time_const; 
  
          }  
        } 
       }
      }
  
      if(max_t_check_plane1>search_range_global[0] && max_t_check_plane1<search_range_global[1]){
      
        search_range[0] = max_t_check_plane1 - search_around_pmax ;
        search_range[1] = max_t_check_plane1 + search_around_pmax ;
  
      }else{
  
        search_range[0] = search_range_global[0] ;
        search_range[1] = search_range_global[1] ;
  
      }

    }else if(mcp_range==1){
    
      //double max_p_check_plane1 = 0;
      //double max_t_check_plane1 = 0;
      max_p_check_plane1 = 0;
      max_t_check_plane1 = 0;
            
      //std::vector<double> w1_check;
      //std::vector<double> t1_check;
      //w1_check.reserve(221560);
      //t1_check.reserve(221560);
      w1_check.clear();
      t1_check.clear();
          
      enable_channel_1 = cf.Value("ACTIVE_CHANNEL", Form("ch%i", ch_mcp) );
      invert_channel_1 = cf.Value("INVERT_SIGNAL", Form("ch%i", ch_mcp) );
      
          
        if( enable_channel_1 == 1){
      
          if( invert_channel_1 == 1 ){
        
            for(unsigned int i=0; i<voltageReader1.at(ch_mcp).GetSize();i++){
      
              if(ADC_conversion==1) w1_check.push_back( double(-voltageReader1.at(ch_mcp).At(i))*ADC_conversion_factor );
              else w1_check.push_back( double(-voltageReader1.at(ch_mcp).At(i)) );
              t1_check.push_back( double(i)*temporal_bin_width ); 
      
            }
      
          }else{
      
            for(unsigned int i=0; i<voltageReader1.at(ch_mcp).GetSize();i++){
      
              if(ADC_conversion==1) w1_check.push_back( double(voltageReader1.at(ch_mcp).At(i))*ADC_conversion_factor );
              else w1_check.push_back( double(voltageReader1.at(ch_mcp).At(i)) );
              t1_check.push_back( double(i)*temporal_bin_width );
      
            }
          }
      
          *a_check=Analyzer( w1_check, t1_check );
          //double baseline_correction = a_check->Correct_Baseline(n_points_baseline); 
          baseline_correction = a_check->Correct_Baseline(n_points_baseline);
      
          //std::pair<double, unsigned int> tp_pair1_small = a_check->Find_Signal_Maximum(pmax_search_range,search_range_global); 
          //std::array<double, 3> fit_array_small = a_check->Pmax_with_GausFit(tp_pair1_small,maxIndex,7);  
          //std::pair<double, double> tp_pair1_fit_small{0.,0.}; 
          tp_pair1_small = a_check->Find_Signal_Maximum(pmax_search_range,search_range_global); 
          fit_array_small = a_check->Pmax_with_GausFit(tp_pair1_small,maxIndex,number_points_gaus_fit);
     
        } 
       
      
  
      if(fit_array_small[1]>search_range_global[0] && fit_array_small[1]<search_range_global[1]){
      
        search_range[0] = fit_array_small[1] - mcp_delay - search_around_mcp ;
        search_range[1] = fit_array_small[1] - mcp_delay + search_around_mcp ;
  
      }else{
  
        search_range[0] = search_range_global[0] ;
        search_range[1] = search_range_global[1] ;
  
      }

    }else{

      search_range[0] = search_range_global[0] ;
      search_range[1] = search_range_global[1] ;

    }
    ///////// END OF "MCP RANGE" PART /////////
    

 


  
  
    for( int ch_counter=0; ch_counter<(active_channels+2); ch_counter++ ){

      if(ch_counter==ch_mcp || ch_counter==16 || ch_counter==17){

        search_range_final[0] = search_range_global[0] ;
        search_range_final[1] = search_range_global[1] ;

      }else{

        search_range_final[0] = search_range[0] ;
        search_range_final[1] = search_range[1] ;

      }
          
      //std::vector<double> w1_inner;
      //std::vector<double> t1_inner;
      //w1_inner.reserve(221560);
      //t1_inner.reserve(221560);
      w1_inner.clear();
      t1_inner.clear();
  
      if(ch_counter < active_channels ){
          
        enable_channel_1 = cf.Value("ACTIVE_CHANNEL", Form("ch%i", ch_counter) );
        invert_channel_1 = cf.Value("INVERT_SIGNAL", Form("ch%i", ch_counter) );
          
      }else if(ch_counter == active_channels){
          
        enable_channel_1 = cf.Value("ACTIVE_CHANNEL", "trg0" );
        invert_channel_1 = cf.Value("INVERT_SIGNAL", "trg0" );
          
      }else if(ch_counter == active_channels+1 ){
          
        enable_channel_1 = cf.Value("ACTIVE_CHANNEL", "trg1" );
        invert_channel_1 = cf.Value("INVERT_SIGNAL", "trg1" );
          
      }

  
 	    if( enable_channel_1 == 1){
  
 	      if( invert_channel_1 == 1 ){
    
	        for(unsigned int i=0; i<voltageReader1.at(ch_counter).GetSize();i++){
  
            if(ADC_conversion==1) w1_inner.push_back( double(-voltageReader1.at(ch_counter).At(i))*ADC_conversion_factor );
            else w1_inner.push_back( double(-voltageReader1.at(ch_counter).At(i)) );
            t1_inner.push_back( double(i)*temporal_bin_width ); 
  
 	        }
  
 	      }else{
  
	        for(unsigned int i=0; i<voltageReader1.at(ch_counter).GetSize();i++){
  
	    	    if(ADC_conversion==1) w1_inner.push_back( double(voltageReader1.at(ch_counter).At(i))*ADC_conversion_factor );
            else w1_inner.push_back( double(voltageReader1.at(ch_counter).At(i)) );
            t1_inner.push_back( double(i)*temporal_bin_width );
  
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
  
	    	*a1=Analyzer( w1_inner, t1_inner );
	    	//double baseline_correction = a1->Correct_Baseline(n_points_baseline); // we do not want signals in the first 5 ns, otherwise baseline correction is biased
        baseline_correction = a1->Correct_Baseline(n_points_baseline); // we do not want signals in the first 5 ns, otherwise baseline correction is biased

        for(int i=0; i<w1_inner.size(); i++) w1_inner.at(i) = w1_inner.at(i) - baseline_correction ;

        w1.push_back( w1_inner );//to be commented for skipping the waveform;
        t1.push_back( t1_inner );//to be commented for skipping the waveform;
  
	    	//std::pair<double, unsigned int> tp_pair1 = a1->Find_Signal_Maximum(pmax_search_range,search_range_final); 
        //std::array<double, 3> fit_array = a1->Pmax_with_GausFit(tp_pair1,maxIndex,number_points_gaus_fit);
        //std::pair<double, double> tp_pair1_fit{0.,0.};
	    	//std::pair<double, unsigned int> neg_tp_pair1 = a1->Find_Negative_Signal_Maximum(pmax_search_range,search_range_final); 
	    	//std::pair<double, double> neg_tp_pair1_fit = a1->Negative_Pmax_with_GausFit(neg_tp_pair1,maxIndex);   
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
        DC_Area1.push_back( a1->DC_Area(baseline_correction)*voltage_const*time_const ); //mV*ns
        Area_NC.push_back( a1->Area_NC(tp_pair1, 5, 5, a1->Find_Noise(n_points_baseline)*voltage_const) ); //mV*ns
        Area_NC_pos.push_back( a1->Area_NC_pos(tp_pair1, 5, 5, a1->Find_Noise(n_points_baseline)*voltage_const) ); //mV*ns
        Area1_new.push_back( a1->New_Pulse_Area(tp_pair1_fit,tp_pair1.second,"Simpson",search_range_final)*voltage_const*time_const ) ;//mV*ns 
        Area_fixed_window.push_back( a1->Pulse_Integration_with_Fixed_Window_Size_with_GausFit(tp_pair1_fit,tp_pair1.second,"Simpson", 1, 1)*voltage_const*time_const ); //mV*ns

        RiseTime1Fit.push_back( a1->Find_Rise_Time_with_GausFit(tp_pair1_fit, tp_pair1.second, 0.1, 0.9)*time_const ) ; //ns
        FallTime1Fit.push_back( a1->Find_Fall_Time_with_GausFit(tp_pair1_fit, tp_pair1.second, 0.1, 0.9)*time_const ) ; //ns
        dVdt1Fit.push_back( a1->Find_Dvdt_with_GausFit(20,0,tp_pair1_fit,tp_pair1.second)*(voltage_const/time_const) ) ;  //mV/ns
	    	dVdt1Fit_2080.push_back( a1->Find_Dvdt2080_with_GausFit(0,tp_pair1_fit,tp_pair1.second)*(voltage_const/time_const) );  //mV/ns

        //std::vector<double> cf_inner ;
	    	//std::vector<double> width_inner ;
        //cf_inner.reserve(7);
        //width_inner.reserve(7);
        cf_inner.clear();
        width_inner.clear();

        for(int jj=0;jj<7;jj++){
  
	    		  cf_inner.push_back( a1->Rising_Edge_CFD_Time_with_GausFit(10+jj*10,tp_pair1_fit,tp_pair1.second)*time_const ) ;
	    		  width_inner.push_back( (a1->Falling_Edge_CFD_Time_with_GausFit(10+jj*10,tp_pair1_fit,tp_pair1.second)*time_const) - 
                                 (a1->Rising_Edge_CFD_Time_with_GausFit(10+jj*10,tp_pair1_fit,tp_pair1.second)*time_const) ) ;
  
	    	}

        CFD1Fit.push_back( cf_inner ) ;
	    	WIDTH1.push_back( width_inner ) ;

	      UArea1.push_back( a1->Find_Undershoot_Area(tp_pair1)*voltage_const*time_const ); //mV*ns
        UArea1_new.push_back( a1->New_Undershoot_Area(tp_pair1_fit,neg_tp_pair1_fit, neg_tp_pair1.second,"Simpson",search_range_final)*voltage_const*time_const ) ;//mV*ns
         
	    	tot1.push_back( a1->Find_Time_Over_Threshold(tot_levels[0],tp_pair1,tot_levels[1])*time_const ) ; //ns
	    	rms1.push_back( a1->Find_Noise(n_points_baseline)*voltage_const ) ; //mV
  
	    }	
    }
  
    event=j_counter;
    

    if(x_pos1>-800){
      
      OutTree->Fill();
      //if(j_counter%10000 == 0) cout<<x_pos1<<" "<<y_pos1<<endl;

    } 
  
    if(j_counter%10000 == 0) cout<<"processed events:"<<j_counter<<endl;
    j_counter++;
  
  }

  OutTree->Write();
  OutputFile->Write();
  OutputFile->Close();
  
}


int main(){

analisi();

return 0;

}
