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

  //Config file definition
  ConfigFile cf("beta_config.ini");

  const bool join_txt_tracker = true;

  //opens txt file and takes the data
  std::string line;
  std::string txtfilename=cf.Value("HEADER", "tracker_filename") ;
  std::ifstream Filein;
  std::vector<double> x_tracker_dirty, y_tracker_dirty; //have multiple tracks per event
  std::vector<int> nevent_dirty;
  std::vector<double> x_tracker, y_tracker; //all multiple tracks events are omitted
  std::vector<int> nevent;
  double x, y, z_position;
  int is_fitted, n, n_old;


  //saves tracker data on x_tracker and y_tracker
  if(join_txt_tracker){
    Filein.open(txtfilename, std::ios::in);
    if(!Filein.is_open()) std::cout << "It failed" << std::endl;
    else std::cout << "Opened file " << txtfilename << std::endl;
    while(getline(Filein, line)){
        Filein >> n >> x >> y >> z_position >> is_fitted; 
        x_tracker_dirty.push_back(x);
        y_tracker_dirty.push_back(y);
        nevent_dirty.push_back(n);
    }
    //filters out multi tracks 
    for(size_t i = 0; i < x_tracker_dirty.size(); i++){
        if(nevent_dirty[i] != nevent_dirty[i+1] || i == x_tracker_dirty.size()){
                if(nevent_dirty[i] != nevent_dirty[i-1] || i == 0){
                x_tracker.push_back(x_tracker_dirty[i]);
                y_tracker.push_back(y_tracker_dirty[i]);
                nevent.push_back(nevent_dirty[i]);
            }
        }
    }

  }

  //time window is the DAQ time window, that you can check on the oscilloscope. search range is the window where signals occur
  bool pmax_search_range;

  if( cf.Value("HEADER", "search_range") == 0 ) pmax_search_range = false;
  else pmax_search_range = true;

  double search_range[2];
  search_range[0] = cf.Value("HEADER", "pmax_search_range_min" ) ;
  search_range[1] = cf.Value("HEADER", "pmax_search_range_max" ) ;

  double tot_levels[2] = { cf.Value("HEADER","tot_rising"), cf.Value("HEADER","tot_falling") };

  int n_points_baseline = cf.Value("HEADER","n_points_baseline");

  int ADC_conversion = cf.Value("HEADER","ADC_conversion");
  double ADC_conversion_factor = cf.Value("HEADER","ADC_conversion_factor");
  double temporal_bin_width = cf.Value("HEADER","temporal_bin_width"); //0.2; 0.0488;

  //Input file
  std::string path = cf.Value("HEADER","filename_path");
  std::string file_in = cf.Value("HEADER","input_filename");
  std::string Filename = path+"raw/"+file_in;
  std::cout << "Anaysis of file " << Filename << " started" << endl; 
  const char *filename = Filename.c_str();
  TFile *file = TFile::Open(filename);
  TTree *itree = dynamic_cast<TTree*>(file->Get("Waveforms"));
  TTreeReader myReader("Waveforms", file);
  
  // Output file & tree
  std::string outFilename;
  if(join_txt_tracker) outFilename = path+"stats/stats_TRACKER_"+file_in;
  else outFilename = path+"stats/stats_"+file_in;
  const char *output_filename = outFilename.c_str();
  TFile *OutputFile = new TFile(output_filename,"recreate");
  TTree *OutTree = new TTree("Analysis","Analysis");

  cout<<" "<<endl;
  cout<<"The output file will be: "<<endl;
  cout<<outFilename<<endl;
  cout<<" "<<endl;
   
  // Variable declaration and Analyzer object 
  const double time_const = cf.Value("HEADER","time_scalar");  
  const double voltage_const = cf.Value("HEADER","voltage_scalar");
  unsigned int maxIndex = cf.Value("HEADER","sampling_points");
  int active_channels = cf.Value("HEADER","active_channels");
  
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
  double x_pos, y_pos ;
  std::string caen_label;
  
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
  
  int event;
  //int evt_delta = 0; //for tracker sync
  
  OutTree->Branch("event",&event);
  //OutTree->Branch("evt_delta",&evt_delta); //for tracker sync
  OutTree->Branch("caen_label",&caen_label);
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
  OutTree->Branch("x_pos", &x_pos);
  OutTree->Branch("y_pos", &y_pos);
      
  n = 0;
  int j_counter = 0;
  
  //std::vector<TTreeReaderArray<Double32_t>> voltageReader1 ;
  
  //TTreeReaderArray<Double32_t> posReader(myReader, "pos");
  TTreeReaderValue<string> producerReader(myReader, "producer");
  TTreeReaderValue<std::vector<std::vector<double>>> voltageReader1(myReader, "voltages" );

      
  /*for(int ch_counter=0; ch_counter<active_channels; ch_counter++ ){
  
    voltageReader1.push_back(TTreeReaderArray<Double32_t>(myReader, "voltages" ));  
  
  }*/
  
  //voltageReader1.push_back(TTreeReaderArray<Double32_t>(myReader, "trg0" )); 
  //voltageReader1.push_back(TTreeReaderArray<Double32_t>(myReader, "trg1" ));
  
  double baseline_correction = 0;
  
  while(myReader.Next() ){
    
    if(join_txt_tracker && (j_counter != nevent[n])){
    //if(join_txt_tracker && ((j_counter+1) != nevent[n]) ){
      
      j_counter++;
      continue;

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
  
    
    if(join_txt_tracker){

      x_pos = x_tracker[n];
      y_pos = y_tracker[n];

      n++;
      
    }else{

      x_pos = 0;
      y_pos = 0;

    }
    
  
    caen_label = *producerReader;

    int enable_channel_1 = 0;
    int invert_channel_1 = 0;
  
    for( int ch_counter=0; ch_counter<active_channels; ch_counter++ ){
          
      std::vector<double> w1_inner;
      std::vector<double> t1_inner;
      w1_inner.reserve(221560);
      t1_inner.reserve(221560);

      if(*producerReader=="CAEN_UZH" || *producerReader=="CAEN_IJS"){
      //if( *producerReader=="CAEN_IJS"){
  
        if(ch_counter < active_channels ){
          
          enable_channel_1 = cf.Value("ACTIVE_CHANNEL", Form("ch%i", ch_counter) );
          invert_channel_1 = cf.Value("INVERT_SIGNAL", Form("ch%i", ch_counter) );
          
        }

  
 	      if( enable_channel_1 == 1){
  
 	        if( invert_channel_1 == 1 ){
    
	          for(unsigned int i=0; i<voltageReader1->at(ch_counter).size();i++){
  
              if(ADC_conversion==1) w1_inner.push_back( double(-voltageReader1->at(ch_counter).at(i))*ADC_conversion_factor );
              else w1_inner.push_back( double(-voltageReader1->at(ch_counter).at(i)) );
              t1_inner.push_back( double(i)*temporal_bin_width ); 
  
 	          }
  
 	        }else{
  
	          for(unsigned int i=0; i<voltageReader1->at(ch_counter).size();i++){

              if(ADC_conversion==1) w1_inner.push_back( double(voltageReader1->at(ch_counter).at(i))*ADC_conversion_factor );
              else w1_inner.push_back( double(voltageReader1->at(ch_counter).at(i)) );
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
	    	  baseline_correction = a1->Correct_Baseline(n_points_baseline); // we do not want signals in the first 5 ns, otherwise baseline correction is biased

          for(int i=0; i<w1_inner.size(); i++) w1_inner.at(i) = w1_inner.at(i) - baseline_correction ;

          w1.push_back( w1_inner );
          t1.push_back( t1_inner );

      }
  
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
  
	    }	
    }
  
    event=j_counter;
        
    OutTree->Fill();
  
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



/*

/// TRACKER SYNC ///
      if( j_counter<15000 ){

        x_pos = x_tracker[n-1];
        y_pos = y_tracker[n-1];

      }else if( j_counter>=15000 && j_counter<28000 ){

        x_pos = x_tracker[n+1];
        y_pos = y_tracker[n+1];

      }else if( (j_counter>=28000 && j_counter<55000) ){

        x_pos = x_tracker[n+2];
        y_pos = y_tracker[n+2];

      }else if( j_counter>=55000 && j_counter<65000 ){

        x_pos = x_tracker[n+3];
        y_pos = y_tracker[n+3];

      }else if( j_counter>=65000 && j_counter<110000 ){

        x_pos = x_tracker[n+4];
        y_pos = y_tracker[n+4];

      }else if( j_counter>=110000 ){

        x_pos = x_tracker[n+5];
        y_pos = y_tracker[n+5];

      }


/// cross 450um ///
     x_rsd_a = 450/2*(Pmax1Fit.at(11) + Pmax1Fit.at(6) - Pmax1Fit.at(9) - Pmax1Fit.at(8)) / (Pmax1Fit.at(6) + Pmax1Fit.at(8) + Pmax1Fit.at(9) + Pmax1Fit.at(11)) ;
     y_rsd_a = 450/2*(Pmax1Fit.at(9) + Pmax1Fit.at(11) - Pmax1Fit.at(8) - Pmax1Fit.at(6)) / (Pmax1Fit.at(9) + Pmax1Fit.at(11) + Pmax1Fit.at(8) + Pmax1Fit.at(6)) ;

     x_rsd_b = 450/2*(Pmax1Fit.at(5) + Pmax1Fit.at(6) - Pmax1Fit.at(7) - Pmax1Fit.at(8)) / (Pmax1Fit.at(5) + Pmax1Fit.at(6) + Pmax1Fit.at(7) + Pmax1Fit.at(8)) ;
     y_rsd_b = 450/2*(Pmax1Fit.at(6) + Pmax1Fit.at(8) - Pmax1Fit.at(5) - Pmax1Fit.at(7)) / (Pmax1Fit.at(5) + Pmax1Fit.at(6) + Pmax1Fit.at(7) + Pmax1Fit.at(8)) ;

     x_rsd_c = 450/2*(Pmax1Fit.at(14) + Pmax1Fit.at(13) - Pmax1Fit.at(11) - Pmax1Fit.at(12)) / (Pmax1Fit.at(14) + Pmax1Fit.at(11) + Pmax1Fit.at(12) + Pmax1Fit.at(13)) ;
     y_rsd_c = 450/2*(Pmax1Fit.at(13) + Pmax1Fit.at(12) - Pmax1Fit.at(14) - Pmax1Fit.at(11)) / (Pmax1Fit.at(14) + Pmax1Fit.at(11) + Pmax1Fit.at(12) + Pmax1Fit.at(13)) ;

     x_rsd_d = 450/2*(Pmax1Fit.at(3) + Pmax1Fit.at(14) - Pmax1Fit.at(6) - Pmax1Fit.at(11)) / (Pmax1Fit.at(3) + Pmax1Fit.at(14) + Pmax1Fit.at(6) + Pmax1Fit.at(11)) ;
     y_rsd_d = 450/2*(Pmax1Fit.at(11) + Pmax1Fit.at(14) - Pmax1Fit.at(6) - Pmax1Fit.at(3)) / (Pmax1Fit.at(3) + Pmax1Fit.at(14) + Pmax1Fit.at(6) + Pmax1Fit.at(11)) ;

     x_rsd_e = 450/2*(Pmax1Fit.at(3) + Pmax1Fit.at(4) - Pmax1Fit.at(6) - Pmax1Fit.at(5)) / (Pmax1Fit.at(3) + Pmax1Fit.at(4) + Pmax1Fit.at(6) + Pmax1Fit.at(5)) ;
     y_rsd_e = 450/2*(Pmax1Fit.at(3) + Pmax1Fit.at(6) - Pmax1Fit.at(5) - Pmax1Fit.at(4)) / (Pmax1Fit.at(3) + Pmax1Fit.at(4) + Pmax1Fit.at(6) + Pmax1Fit.at(5)) ;
    
    
    double pad_x_pos[15] = {0,675,675,225,225,-225,-225,-675,-675,-675,0,-225,-225,225,225};
    double pad_y_pos[15] = {0,-225,-675,-225,-675,-675,-225,-675,-225,225,0,225,675,675,225};
    double den = 0;
    double num_x = 0;
    double num_y = 0;

    for(int i=0; i<15; i++){

      if( i!=0 || i!=1 || i!=2 || i !=10 || i!=12 || i!=13 ){

        den += Pmax1Fit.at(i);
        num_x += Pmax1Fit.at(i)*pad_x_pos[i];
        num_y += Pmax1Fit.at(i)*pad_y_pos[i];

      }
    }


    x_centroide = num_x/den - 4800. ;
    y_centroide = num_y/den - 557. ;


    double x_diff[21];
    double minv = 100000;
    int min_counter = -100;
    
    for(int i=0; i<21;i++){

      if( x_tracker[n-10+i]>-5.475 && x_tracker[n-10+i]<-4.575 && y_tracker[n-10+i]>-1.232 && y_tracker[n-10+i]<-0.332 ) x_diff[i] = fabs(x_centroide - x_tracker[n-10+i]*1000) + fabs(y_centroide - y_tracker[n-10+i]*1000);
      else x_diff[i] = 100000;

      if(x_diff[i]<minv){ 

         minv = x_diff[i];
         min_counter = i;
       }

     }

   evt_delta = min_counter-10;







    */
