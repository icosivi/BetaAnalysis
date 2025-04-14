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


void analisi_csv( ){

  int ch_mcp = 8;

  float pixel_pos[8][2] = {
    {0.150f, 0.150f},
    {0.450f, 0.150f},
    {0.150f, 0.450f},
    {0.450f, 0.450f},
    {0.150f, 0.750f},
    {0.450f, 0.750f},
    {0.150f, 1.050f},
    {0.450f, 1.050f}
  };

  int mask[8][4] = {
    {3, 5, 6, 4},
    {2, 3, 4, 1},
    {4, 6, 7, 13},
    {1, 4, 13, 0},
    {13, 7, 9, 11},
    {0, 13, 11, 15},
    {11, 9, 10, 12},
    {15, 11, 12, 14}
  };
  
   int number_samplings = 40;
   //double main_vec[4*number_samplings];

   ofstream f;
   f.open("/Users/icosivi/Desktop/DESY_TB8_DCRSD/csv/run108_40samplings_ReReco_time.csv");

   f << "x,y,t"; // Initial x and y

   for (int i = 0; i < (number_samplings*4); ++i) {
    f << ",w" << i; // Add "w" followed by the index
    }

   f << "\n";

  //ROOT::EnableImplicitMT(6);
  //ROOT::EnableThreadSafety();

  //Config file definition
  ConfigFile cf("beta_config.ini");

  bool join_txt_tracker = true;

  //time window is the DAQ time window, that you can check on the oscilloscope. search range is the window where signals occur
  bool pmax_search_range;

  if( cf.Value("HEADER", "search_range") == 0 ) pmax_search_range = false;
  else pmax_search_range = true;

  double search_range[2] = {0,0};
  double search_range_global[2] = {0,0};
  double search_range_final[2] = {0,0};
  search_range_global[0] = cf.Value("HEADER", "pmax_search_range_min" ) ;
  search_range_global[1] = cf.Value("HEADER", "pmax_search_range_max" ) ;
  
  double search_around_pmax = 1;

  double tot_levels[2] = { cf.Value("HEADER","tot_rising"), cf.Value("HEADER","tot_falling") };

  int n_points_baseline = cf.Value("HEADER","n_points_baseline");

  int ADC_conversion = cf.Value("HEADER","ADC_conversion");
  double ADC_conversion_factor = cf.Value("HEADER","ADC_conversion_factor");
  double temporal_bin_width = cf.Value("HEADER","temporal_bin_width"); //0.2; 0.0488;

  std::string Filename = cf.Value("HEADER","input_filename");
  std::cout << "Anaysis of file " << Filename << " started" << endl; 
  const char *filename = Filename.c_str();
  TFile *file = TFile::Open(filename);
  TTree *itree = dynamic_cast<TTree*>(file->Get("wfm"));
  TTreeReader myReader("wfm", file);


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
  std::vector<double> Area1_new;
  std::vector<double> rms1;
  std::vector<std::vector<double>> CFD1Fit;
  std::vector<std::vector<double>> w1 ; //to be commented for skipping the waveform;
  std::vector<std::vector<double>> t1 ; //to be commented for skipping the waveform;
  double x_pos1, y_pos1,x_pos2, y_pos2, chi2_trk ;
  
  Pmax1.reserve(20);
  Pmax1Fit.reserve(20);
  negPmax1Fit.reserve(20);
  Tmax1.reserve(20);
  Tmax1Fit.reserve(20);
  negTmax1Fit.reserve(20);
  Area1.reserve(20);
  Area1_new.reserve(20);
  rms1.reserve(20);
  CFD1Fit.reserve(20);
  
  w1.reserve(20);//to be commented for skipping the waveform;
  t1.reserve(20);//to be commented for skipping the waveform;
  Analyzer *a1=new Analyzer();
  Analyzer *a_check=new Analyzer();
  
  int event;
  //int evt_delta = 0; //for tracker sync
  
      
  //n = 0;
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
  
  
  while( myReader.Next() ){ //  && j_counter<20000
    Pmax1.clear();
    Pmax1Fit.clear();
    negPmax1Fit.clear();
    Tmax1.clear();
    Tmax1Fit.clear();
    negTmax1Fit.clear();
    Area1.clear();
    Area1_new.clear();
    rms1.clear();
    CFD1Fit.clear();
    w1.clear();//to be commented for skipping the waveform;
    t1.clear();//to be commented for skipping the waveform;
  
    
    if(join_txt_tracker){

      x_pos1 = *x1Reader;
      y_pos1 = *y1Reader;
      x_pos2 = *x2Reader;
      y_pos2 = *y2Reader;
      chi2_trk = *chi2Reader;

      //n++;
      
    }else{

      x_pos1 = 0;
      y_pos1 = 0;
      x_pos2 = 0;
      y_pos2 = 0;
      chi2_trk = 0;
    }
    
  
    int enable_channel_1 = 0;
    int invert_channel_1 = 0;


    
     ///////// BEGINNING OF "SMALL-RANGE" PART /////////
    
    
      double max_p_check_plane1 = 0;
      double max_t_check_plane1 = 0;
  
      for( int ch_counter=0; ch_counter<active_channels; ch_counter++ ){

       if(ch_counter!=ch_mcp){
            
        std::vector<double> w1_check;
        std::vector<double> t1_check;
      
        w1_check.reserve(221560);
        t1_check.reserve(221560);
          
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
          double baseline_correction = a_check->Correct_Baseline(n_points_baseline); 
      
          std::pair<double, unsigned int> tp_pair1 = a_check->Find_Signal_Maximum(pmax_search_range,search_range_global); 
          std::pair<double, double> tp_pair1_fit = a_check->Pmax_with_GausFit(tp_pair1,maxIndex);   
    
          if( tp_pair1_fit.first*voltage_const > max_p_check_plane1){
    
            max_p_check_plane1 = tp_pair1_fit.first*voltage_const;
            max_t_check_plane1 = tp_pair1_fit.second*time_const; 
  
          }  
        } 
       }
      }
  
      if(max_t_check_plane1>5 && max_t_check_plane1<195){
      
        search_range[0] = max_t_check_plane1 - search_around_pmax ;
        search_range[1] = max_t_check_plane1 + search_around_pmax ;
  
      }else{
  
        search_range[0] = search_range_global[0] ;
        search_range[1] = search_range_global[1] ;
  
      }
  
    

    ///////// END OF "SMALL-RANGE" PART /////////
    

    //search_range[0] = search_range_global[0] ;
    //search_range[1] = search_range_global[1] ;

    


  
    for( int ch_counter=0; ch_counter<(active_channels+2); ch_counter++ ){

      if(ch_counter==ch_mcp || ch_counter==16 || ch_counter==17){

        search_range_final[0] = search_range_global[0] ;
        search_range_final[1] = search_range_global[1] ;

      }else{

        search_range_final[0] = search_range[0] ;
        search_range_final[1] = search_range[1] ;

      }

      //cout<<search_range_final[0]<<"    "<<search_range_final[1]<<endl;
          
      std::vector<double> w1_inner;
      std::vector<double> t1_inner;
  
      w1_inner.reserve(221560);
      t1_inner.reserve(221560);
  
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
	    	double baseline_correction = a1->Correct_Baseline(n_points_baseline); // we do not want signals in the first 5 ns, otherwise baseline correction is biased

        for(int i=0; i<w1_inner.size(); i++) w1_inner.at(i) = w1_inner.at(i) - baseline_correction ;

        w1.push_back( w1_inner );//to be commented for skipping the waveform;
        t1.push_back( t1_inner );//to be commented for skipping the waveform;
  
	    	std::pair<double, unsigned int> tp_pair1 = a1->Find_Signal_Maximum(pmax_search_range,search_range_final); 
	    	std::pair<double, double> tp_pair1_fit = a1->Pmax_with_GausFit(tp_pair1,maxIndex);
	    	std::pair<double, unsigned int> neg_tp_pair1 = a1->Find_Negative_Signal_Maximum(pmax_search_range,search_range_final); 
	    	std::pair<double, double> neg_tp_pair1_fit = a1->Negative_Pmax_with_GausFit(neg_tp_pair1,maxIndex);      
	    	Pmax1.push_back( tp_pair1.first*voltage_const ) ; //mV
	    	Pmax1Fit.push_back( tp_pair1_fit.first*voltage_const ) ; //mV
	    	negPmax1Fit.push_back(  neg_tp_pair1_fit.first*voltage_const ); //mV
	    	Tmax1.push_back(  a1->Get_Tmax(tp_pair1)*time_const ) ; //ns
	    	Tmax1Fit.push_back(  tp_pair1_fit.second*time_const ) ; //ns
	    	negTmax1Fit.push_back(  neg_tp_pair1_fit.second*time_const ) ; //ns
	    	Area1.push_back(  a1->Find_Pulse_Area(tp_pair1)*voltage_const*time_const ) ; //mV*ns
	    	Area1_new.push_back( a1->New_Pulse_Area(tp_pair1_fit,tp_pair1.second,"Simpson",search_range_final)*voltage_const*time_const ) ;//mV*ns  
	    	rms1.push_back( a1->Find_Noise(n_points_baseline)*voltage_const ) ; //mV     

        std::vector<double> cf_inner ;
        cf_inner.reserve(7);
  
	    	for(int jj=0;jj<7;jj++) cf_inner.push_back( a1->Rising_Edge_CFD_Time_with_GausFit(10+jj*10,tp_pair1_fit,tp_pair1.second)*time_const ) ;
  
	    	CFD1Fit.push_back( cf_inner ) ;
  
	    }	
    }
  
     event=j_counter;

        double pmax_pixels[8] = {0,0,0,0,0,0,0,0};
        for( int i=0; i<8;i++) pmax_pixels[i] = Pmax1Fit.at(mask[i][0]) + Pmax1Fit.at(mask[i][1]) + Pmax1Fit.at(mask[i][2]) + Pmax1Fit.at(mask[i][3]);
  
  
        int size = sizeof(pmax_pixels) / sizeof(pmax_pixels[0]); // Calculate the size of the static array
        double max_value = pmax_pixels[0]; // Initialize with the first element
        int max_index = 0;
  
        for (int i = 1; i < size; ++i) {
          if (pmax_pixels[i] > max_value) {
              max_value = pmax_pixels[i];
              max_index = i;
        }
        }
  
  
        double p_all[15] = {Pmax1Fit.at(0),Pmax1Fit.at(1),Pmax1Fit.at(2),Pmax1Fit.at(3),Pmax1Fit.at(4),Pmax1Fit.at(5),Pmax1Fit.at(6),Pmax1Fit.at(7),Pmax1Fit.at(9),Pmax1Fit.at(10),Pmax1Fit.at(11),Pmax1Fit.at(12),Pmax1Fit.at(13),Pmax1Fit.at(14),Pmax1Fit.at(15)};
        int size_all = sizeof(p_all) / sizeof(p_all[0]); // Calculate the size of the static array
        float max_value_all = p_all[0]; // Initialize with the first element
        int max_index_all = 0;
  
        for (int i = 1; i < size_all; ++i) {
          if (p_all[i] > max_value_all) {
              max_value_all = p_all[i];
              max_index_all = i;
        }
        }

        double main_vec[4*number_samplings];
        double timestamp = 0. ;
        double timestamp2 = 0. ;
        double dt = 0. ;

        for(int j=0; j<4*number_samplings; j++)  main_vec[j]=0;


        for(int j=0; j<int(w1.at(max_index_all).size()); j++){

         if(t1.at(max_index_all).at(j)<=Tmax1Fit.at(max_index_all) && t1.at(max_index_all).at(j+1)>=Tmax1Fit.at(max_index_all)){
         //if(w1.at(max_index_all).at(j)*ADC_conversion_factor<=max_value_all && w1.at(max_index_all).at(j+1)*ADC_conversion_factor>=max_value_all){

            for(int k=0; k<number_samplings; k++){

                main_vec[k] = w1.at(mask[max_index][0]).at(j-(number_samplings/2)+1+k);
                main_vec[k+number_samplings] = w1.at(mask[max_index][1]).at(j-(number_samplings/2)+1+k);
                main_vec[k+2*number_samplings] = w1.at(mask[max_index][2]).at(j-(number_samplings/2)+1+k);
                main_vec[k+3*number_samplings] = w1.at(mask[max_index][3]).at(j-(number_samplings/2)+1+k);

            }
            timestamp = t1.at(max_index_all).at(j);
            timestamp2 = t1.at(max_index_all).at(j+1);
            break;
         }
         }

        //cout<<max_index<<"  "<<max_index_all<<endl;

       if(max_value>25 && max_index!=4 && max_index!=6 ){  
          pippo->Fill(timestamp-Tmax1Fit.at(max_index_all));
          pippo2->Fill(timestamp2-Tmax1Fit.at(max_index_all));

          if(ch_mcp<8){

            if(max_index_all<8) dt = CFD1Fit[max_index_all][4]-CFD1Fit[ch_mcp][2];
            else dt = CFD1Fit[max_index_all][4]-CFD1Fit[ch_mcp][2] + (CFD1Fit[16][4]-CFD1Fit[17][4]);

          }else{

            if(max_index_all>=8) dt = CFD1Fit[max_index_all][4]-CFD1Fit[ch_mcp][2];
            else dt = CFD1Fit[max_index_all][4]-CFD1Fit[ch_mcp][2] - (CFD1Fit[16][4]-CFD1Fit[17][4]);

          }
        
        f << x_pos1-pixel_pos[max_index][0] << "," << y_pos1-pixel_pos[max_index][1]<< "," << dt; 

         for (int i = 0; i < number_samplings*4; ++i) {
         f << "," << main_vec[i]; 
        }

        f << "\n";

       } 




  
    if(j_counter%10000 == 0) cout<<"processed events:"<<j_counter<<endl;
    j_counter++;




  
  }

   f.close();
  
}


int main(){

analisi_csv();

return 0;

}
