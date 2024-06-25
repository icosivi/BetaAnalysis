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
#include <sys/stat.h>
#include <dirent.h>
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
#include <TROOT.h>
#include <TGraph.h>
#include <TThread.h>
#include <TROOT.h>
#include <TStyle.h>
#include <TImage.h>
#include <TCanvas.h>
#include <TSystemFile.h>
#include <TSystemDirectory.h>

//------Custom headers----------------//
#include "src/Analyzer.hpp"
#include "include/ConfigFile.hpp"


void analisi(){

  //Config file definition
  ConfigFile cf("beta_config.ini");

  //time window is the DAQ time window, that you can check on the oscilloscope. search range is the window where signals occur
  double search_range[2] = {0,0};
  double time_window[2] = {0,0}; 

  bool pmax_search_range;

  if( cf.Value("HEADER", "search_range") == 0 ) pmax_search_range = false;
  else pmax_search_range = true;

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

  int ps_channel[4] = {0,0,0,0};
  ps_channel[0] = cf.Value("HEADER", "ps_channel0" ) ;
  ps_channel[1] = cf.Value("HEADER", "ps_channel1" ) ;
  ps_channel[2] = cf.Value("HEADER", "ps_channel2" ) ;
  ps_channel[3] = cf.Value("HEADER", "ps_channel3" ) ;

  search_range[0] = cf.Value("HEADER", "pmax_search_range_min" ) ;
  search_range[1] = cf.Value("HEADER", "pmax_search_range_max" ) ;

  double tot_levels[2] = { cf.Value("HEADER","tot_rising"), cf.Value("HEADER","tot_falling") };

  int n_points_baseline = cf.Value("HEADER","n_points_baseline");

  
  std::string Filename = cf.Value("HEADER","input_filename");
  std::cout << "Anaysis of file " << Filename << " started" << endl;
  const char *filename = Filename.c_str();
  TFile *file = TFile::Open(filename);
  TTree *itree = dynamic_cast<TTree*>(file->Get("wfm"));
  TTreeReader myReader("wfm", file);


  std::string outFilename = cf.Value("HEADER","output_filename");
  cout<<" "<<endl;
  cout<<"The output file will be: "<<endl;
  cout<<outFilename<<endl;
  cout<<" "<<endl;
  const char *output_filename = outFilename.c_str();
  TFile *OutputFile = new TFile(output_filename,"recreate");
  TTree *OutTree = new TTree("Analysis","Analysis");


  /*  old way of saving files: fill config file with path to the raw and stats directories, then provide the name of the raw; the stats will be saved in the stats directory with stats+raw_filename
  std::string path = cf.Value("HEADER","filename_path");
  std::string file_in = cf.Value("HEADER","input_filename");
  std::string Filename = path+"raw/"+file_in;
  std::cout << "Anaysis of file " << Filename << " started" << endl;
  const char *filename = Filename.c_str();
  TFile *file = TFile::Open(filename);
  TTree *itree = dynamic_cast<TTree*>(file->Get("wfm"));
  TTreeReader myReader("wfm", file);

  // Output file & tree
  std::string delimiter = "Sr";
  std::string token_pre = file_in.substr(0, file_in.find(delimiter));
  std::string token_post = file_in.substr(file_in.find(delimiter));
  std::string outDir = path+"stats/"+token_pre;
  const char *outdir = outDir.c_str();
  std::string outFilename = path+"stats/"+token_pre+"stats_"+token_post;

  cout<<" "<<endl;
  cout<<"The output file will be: "<<endl;
  cout<<outFilename<<endl;
  cout<<" "<<endl;

  int check;

  struct stat st;
  if( stat( outdir, &st ) == 0){

    cout<<"output directory already exists"<<endl;

  }else{

    check = mkdir(outdir, 0777);
    if(check==0) cout<<"directory succesfully created"<<endl;
    else cout<<"something went wrong in creating the output directory..."<<endl;

  }

  const char *output_filename = outFilename.c_str();
  TFile *OutputFile = new TFile(output_filename,"recreate");
  TTree *OutTree = new TTree("Analysis","Analysis");*/


 
  // Variable declaration and Analyzer object 
  const double time_const = cf.Value("HEADER","time_scalar");  
  const double voltage_const = cf.Value("HEADER","voltage_scalar");
  unsigned int maxIndex = cf.Value("HEADER","sampling_points"); //number of sampling points, default is 1002, assuming a 50ns time window with 20 GS
  int ch_number = cf.Value("HEADER","active_channels");

  std::vector<double> Pmax1;
  std::vector<double> Pmax1Fit;
  std::vector<double> negPmax1Fit;
  std::vector<double> Tmax1;
  std::vector<double> Tmax1Fit;
  std::vector<double> negTmax1Fit;
  std::vector<double> Area1;
  std::vector<double> UArea1;
  std::vector<double> Area1_new;
  std::vector<double> UArea1_new;
  std::vector<double> DC_Area1;
  std::vector<double> RiseTime1Fit;
  std::vector<double> FallTime1Fit;
  std::vector<double> dVdt1Fit;
  std::vector<double> dVdt1Fit_2080;
  std::vector<std::vector<double>> CFD1Fit;
  std::vector<std::vector<double>> WIDTH1;
  std::vector<double> t_thr1;
  std::vector<double> tot1;
  std::vector<double> rms1;
  std::vector<std::vector<double>> w1;
  std::vector<std::vector<double>> t1;

  std::vector<double> i_current;
  std::vector<double> v_bias;
  

  Pmax1.reserve(20);
  Pmax1Fit.reserve(20);
  negPmax1Fit.reserve(20);
  Tmax1.reserve(20);
  Tmax1Fit.reserve(20);
  negTmax1Fit.reserve(20);
  Area1.reserve(20);
  UArea1.reserve(20);
  Area1_new.reserve(20);
  UArea1_new.reserve(20);
  DC_Area1.reserve(20);
  RiseTime1Fit.reserve(20);
  FallTime1Fit.reserve(20);
  dVdt1Fit.reserve(20);
  dVdt1Fit_2080.reserve(20);
  t_thr1.reserve(20);
  tot1.reserve(20);
  rms1.reserve(20);
  CFD1Fit.reserve(20);
  WIDTH1.reserve(20);
  w1.reserve(20);
  t1.reserve(20);
  Analyzer *a1=new Analyzer();

  i_current.reserve(20);
  v_bias.reserve(20);

  int event;
  double timestamp;

  OutTree->Branch("event",&event);
  OutTree->Branch("time",&timestamp);
  OutTree->Branch("w", "std::vector<std::vector<double>>", &w1);
  OutTree->Branch("t", "std::vector<std::vector<double>>" ,&t1);
  OutTree->Branch("pmax", "std::vector<double>",&Pmax1Fit);
  OutTree->Branch("negpmax", "std::vector<double>",&negPmax1Fit);
  OutTree->Branch("tmax", "std::vector<double>",&Tmax1Fit);
  OutTree->Branch("negtmax", "std::vector<double>",&negTmax1Fit);
  OutTree->Branch("area", "std::vector<double>",&Area1);
  OutTree->Branch("uarea", "std::vector<double>",&UArea1);
  OutTree->Branch("area_new", "std::vector<double>",&Area1_new);
  OutTree->Branch("uarea_new", "std::vector<double>",&UArea1_new);
  OutTree->Branch("dc_area", "std::vector<double>",&DC_Area1);
  OutTree->Branch("risetime", "std::vector<double>",&RiseTime1Fit);
  OutTree->Branch("falltime", "std::vector<double>",&FallTime1Fit);
  OutTree->Branch("dvdt", "std::vector<double>",&dVdt1Fit);
  OutTree->Branch("dvdt_2080", "std::vector<double>",&dVdt1Fit_2080);
  OutTree->Branch("cfd", "std::vector<std::vector<double>>",&CFD1Fit);
  OutTree->Branch("width", "std::vector<std::vector<double>>",&WIDTH1);
  OutTree->Branch("t_thr", "std::vector<double>",&t_thr1);  // time at which a certain thr (in V) is passed
  OutTree->Branch("tot", "std::vector<double>",&tot1);
  OutTree->Branch("rms", "std::vector<double>",&rms1);
  
  OutTree->Branch("I", "std::vector<double>", &i_current);
  OutTree->Branch("V", "std::vector<double>", &v_bias);
  
 
  int j_counter = 0;

  std::vector<TTreeReaderArray<double>> voltageReader1 ;
  std::vector<TTreeReaderArray<double>> timeReader1 ;
  TTreeReaderArray<double> currentReader1(myReader,"i_current") ;
  TTreeReaderArray<double> biasReader1(myReader,"v_bias") ;

  //int enable_channel = 0;
  //for(int ch_counter=1; ch_counter<=ch_number; ch_counter++ ){
  for(int ch_counter=1; ch_counter<=8; ch_counter++ ){

    if(active_channel[ch_counter-1]==1){

      voltageReader1.push_back(TTreeReaderArray<double>(myReader, Form("w%i",ch_counter) ));  
      timeReader1.push_back(TTreeReaderArray<double>(myReader, Form("t%i",ch_counter) ));

    }    
  }

  TTreeReaderValue<double> tstampReader1(myReader,"i_timestamp") ;

  int ps_total = 0 ;
  
  for(int ps_counter=0; ps_counter<4; ps_counter++){

    if(ps_channel[ps_counter] == 1){ 

      //currentReader1.push_back( TTreeReaderValue<double>(myReader, Form("I%i",ps_counter) ) );
      //biasReader1.push_back( TTreeReaderValue<double>(myReader, Form("V%i",ps_counter) ) );
      ps_total++ ;

    }

  }



  while(myReader.Next()){

    timestamp = *tstampReader1;
    i_current.clear();
    v_bias.clear();
    
    for( int ps_counter=0; ps_counter<ps_total; ps_counter++ ){

      //i_current.push_back( *currentReader1.at(ps_counter) ) ;
      //v_bias.push_back( *biasReader1.at(ps_counter) ) ;
      i_current.push_back( currentReader1[ps_counter] ) ;
      v_bias.push_back( biasReader1[ps_counter] ) ;
      //i_current.push_back( 0 ) ;
      //v_bias.push_back( 0 ) ;

    }


    Pmax1.clear();
    Pmax1Fit.clear();
    negPmax1Fit.clear();
    Tmax1.clear();
    Tmax1Fit.clear();
    negTmax1Fit.clear();
    Area1.clear();
    UArea1.clear();
    Area1_new.clear();
    UArea1_new.clear();
    DC_Area1.clear();
    RiseTime1Fit.clear();
    FallTime1Fit.clear();
    dVdt1Fit.clear();
    dVdt1Fit_2080.clear();
    t_thr1.clear();
    tot1.clear();
    rms1.clear();
    CFD1Fit.clear();
    WIDTH1.clear();
    w1.clear();
    t1.clear();
    
    int active_ch_counter = 0;
    for( int ch_counter=1; ch_counter<=8; ch_counter++ ){
      

      std::vector<double> w1_inner;
      std::vector<double> t1_inner;

      w1_inner.reserve(221560);
      t1_inner.reserve(221560);

      //int enable_channel_1 = cf.Value("ACTIVE_CHANNEL", Form("ch%i", ch_counter) );
      //int invert_channel_1 = cf.Value("INVERT_SIGNAL", Form("ch%i", ch_counter) );
 	    

 	    if( active_channel[ch_counter-1]==1){

        if( j_counter == 0 ){

          time_window[0] = timeReader1.at(active_ch_counter).At(0);
          time_window[1] = timeReader1.at(active_ch_counter).At(timeReader1.at(active_ch_counter).GetSize()-1);
          
          cout<<" "<<endl;
          cout<<"Time window: "<<time_window[0]<<"; "<<time_window[1]<<endl;
        
        }

 		    if( invert_channel[ch_counter-1]==1 ){

	 		    for(unsigned int i=0; i<voltageReader1.at(active_ch_counter).GetSize();i++){

	 			    w1_inner.push_back(-voltageReader1.at(active_ch_counter).At(i));
	 			    t1_inner.push_back(timeReader1.at(active_ch_counter).At(i));

 			    }

 		    }else{

	 		    for(unsigned int i=0; i<voltageReader1.at(active_ch_counter).GetSize();i++){

	 			    w1_inner.push_back(voltageReader1.at(active_ch_counter).At(i));
	 			    t1_inner.push_back(timeReader1.at(active_ch_counter).At(i));

 			    }
 		    }

 		    w1.push_back( w1_inner );
 		    t1.push_back( t1_inner );


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
		    double baseline_correction = a1->Correct_Baseline(n_points_baseline);

		    std::pair<double, unsigned int> tp_pair1 = a1->Find_Signal_Maximum(pmax_search_range,search_range); 
		    std::pair<double, double> tp_pair1_fit = a1->Pmax_with_GausFit(tp_pair1,maxIndex);
		    std::pair<double, unsigned int> neg_tp_pair1 = a1->Find_Negative_Signal_Maximum(pmax_search_range,search_range); 
		    std::pair<double, double> neg_tp_pair1_fit = a1->Negative_Pmax_with_GausFit(neg_tp_pair1,maxIndex);      
		    Pmax1.push_back( tp_pair1.first*voltage_const ) ; //mV
		    Pmax1Fit.push_back( tp_pair1_fit.first*voltage_const ) ; //mV
		    negPmax1Fit.push_back(  neg_tp_pair1_fit.first*voltage_const ); //mV
		    Tmax1.push_back(  a1->Get_Tmax(tp_pair1)*time_const ) ; //ns
		    Tmax1Fit.push_back(  tp_pair1_fit.second*time_const ) ; //ns
		    negTmax1Fit.push_back(  neg_tp_pair1_fit.second*time_const ) ; //ns
		    Area1.push_back(  a1->Find_Pulse_Area(tp_pair1)*voltage_const*time_const ) ; //mV*ns
		    UArea1.push_back( a1->Find_Undershoot_Area(tp_pair1)*voltage_const*time_const ); //mV*ns
		    Area1_new.push_back( a1->New_Pulse_Area(tp_pair1_fit,tp_pair1.second,"Simpson",search_range)*voltage_const*time_const ) ;//mV*ns  
		    UArea1_new.push_back( a1->New_Undershoot_Area(tp_pair1_fit,neg_tp_pair1_fit, neg_tp_pair1.second,"Simpson",search_range)*voltage_const*time_const ) ;//mV*ns
		    DC_Area1.push_back( a1->DC_Area(baseline_correction)*voltage_const*time_const ); //mV*ns
        RiseTime1Fit.push_back( a1->Find_Rise_Time_with_GausFit(tp_pair1_fit, tp_pair1.second, 0.1, 0.9)*time_const ) ; //ns
		    FallTime1Fit.push_back( a1->Find_Fall_Time_with_GausFit(tp_pair1_fit, tp_pair1.second, 0.1, 0.9)*time_const ) ; //ns
		    dVdt1Fit.push_back( a1->Find_Dvdt_with_GausFit(20,0,tp_pair1_fit,tp_pair1.second)*(voltage_const/time_const) ) ;  //mV/ns
		    dVdt1Fit_2080.push_back( a1->Find_Dvdt2080_with_GausFit(0,tp_pair1_fit,tp_pair1.second)*(voltage_const/time_const) );  //mV/ns
		    t_thr1.push_back( a1->Find_Time_At_Threshold_with_GausFit(tot_levels[0],tp_pair1_fit,tp_pair1.second)*time_const ); //ns  
		    tot1.push_back( a1->Find_Time_Over_Threshold(tot_levels[0],tp_pair1,tot_levels[1])*time_const ) ; //ns
		    rms1.push_back( a1->Find_Noise(n_points_baseline)*voltage_const ) ; //mV

		    std::vector<double> cf_inner ;
		    std::vector<double> width_inner ;

        cf_inner.reserve(7);
        width_inner.reserve(7);

		    for(int jj=0;jj<7;jj++){

			    cf_inner.push_back( a1->Rising_Edge_CFD_Time_with_GausFit(10+jj*10,tp_pair1_fit,tp_pair1.second)*time_const ) ;
			    width_inner.push_back( (a1->Falling_Edge_CFD_Time_with_GausFit(10+jj*10,tp_pair1_fit,tp_pair1.second)*time_const) - 
                               (a1->Rising_Edge_CFD_Time_with_GausFit(10+jj*10,tp_pair1_fit,tp_pair1.second)*time_const) ) ;

		    }

		    CFD1Fit.push_back( cf_inner ) ;
		    WIDTH1.push_back( width_inner ) ;

        active_ch_counter++;

	    }	
    }

    event=j_counter;
    
    OutTree->Fill();

    if(j_counter==0) cout<<" "<<endl;
    if(j_counter%10000 == 0) cout<<"processed events:"<<j_counter<<endl;
    j_counter++;

  }

  OutTree->Write();
  OutputFile->Write();
 
 }


int main(){

analisi();

return 0;

}
