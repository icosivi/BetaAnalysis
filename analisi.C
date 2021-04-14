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

//------Custom headers----------------//
#include "src/Analyzer.hpp"
#include "include/ConfigFile.hpp"


void analisi(){

//CAREFUL: select search_range for signal search and oscilloscope time_window
//double search_range[2] = {-2e-9,2e-9};
//double time_window[2] = {-25e-9,25e-9}; 
//double search_range[2] = {20,30};
//double time_window[2] = {0,50};
//double search_range[2] = {-0.215e-6,-0.2e-6};
//double time_window[2] = {-0.22e-6,-0.19e-6};


//time window is the DAQ time window, that you can check on the oscilloscope. search range is the window where signals occur
double search_range[2] = {0,0};
double time_window[2] = {0,0}; 


//Config file definition
 ConfigFile cf("beta_config.ini");

//Input file
 std::string Filename = cf.Value("HEADER","input_filename");
 const char *filename = Filename.c_str();
 TFile *file = TFile::Open(filename);
 TTree *itree = dynamic_cast<TTree*>(file->Get("wfm"));
 TTreeReader myReader("wfm", file);
 
//Definition of TTreeReaders
 int ch_number = cf.Value("HEADER","active_channels");
 char *branch[ch_number];
 char *t_branch[ch_number];
 int ch_cont[4];
 int cont;

 for(int j=0;j<4;j++){

	if(itree->GetListOfBranches()->FindObject(Form("w%i",j+1)) != 0){

		ch_cont[j]=1;
		cont=j+1;
	 	
	}else{

		ch_cont[j]=0;
	}

 }


 

 for(int j=0;j<4;j++){

	if(ch_cont[j]==1){

	 	branch[j] = Form("w%i",j+1);
	 	t_branch[j] = Form("t%i",j+1);
	 		
	}else{

	 	branch[j] = Form("w%i",cont);
	 	t_branch[j] = Form("t%i",cont);
	  		
	}

 }

 TTreeReaderArray<double> voltageReader1(myReader, branch[0]);  
 TTreeReaderArray<double> timeReader1(myReader, t_branch[0]);
 TTreeReaderArray<double> voltageReader2(myReader, branch[1]);  
 TTreeReaderArray<double> timeReader2(myReader, t_branch[1]);
 TTreeReaderArray<double> voltageReader3(myReader, branch[2]);  
 TTreeReaderArray<double> timeReader3(myReader, t_branch[2]);
 TTreeReaderArray<double> voltageReader4(myReader, branch[3]);  
 TTreeReaderArray<double> timeReader4(myReader, t_branch[3]);

 //TTreeReaderValue<double> xReader(myReader,"x_pos");
 //TTreeReaderValue<double> yReader(myReader,"y_pos");
// End of TTreeReaders definition

// Output file & tree
 std::string outFilename = cf.Value("HEADER","output_filename");
 const char *output_filename = outFilename.c_str();
 TFile *OutputFile = new TFile(output_filename,"recreate");
 TTree *OutTree = new TTree("Analysis","Analysis");
 
// Variable declaration and Analyzer objects 
// Raw data have Voltage in V and Time in ns
 const double time_const = cf.Value("HEADER","time_scalar");  //multiply by this const to pass from s to ns
 const double voltage_const = cf.Value("HEADER","voltage_scalar"); // pass from V to mV
 const double conversion_factor = cf.Value("HEADER","conversion_factor_area_charge"); //pass from Area to Charge
 unsigned int maxIndex = cf.Value("HEADER","sampling_points"); //number of sampling points, default is 1002, assuming a 50ns time window with 20 GS

 double Pmax1=-1000;
 double Pmax1Fit=-1000;
 double negPmax1Fit=-1000;
 double Tmax1=-1000;
 double Tmax1Fit=-1000;
 double negTmax1Fit=-1000;
 double Area1=-1000;
 double UArea1=-1000;
 double Area1_new=-1000;
 double UArea1_new=-1000;
 double RiseTime1Fit=-1000;
 double FallTime1Fit=-1000;
 double dVdt1Fit=-1000;
 double dVdt1Fit_2080=-1000;
 double CFD1Fit[7]={-1000,-1000,-1000,-1000,-1000,-1000,-1000};
 double WIDTH1[7]={-1000,-1000,-1000,-1000,-1000,-1000,-1000};
 double t_thr1=-1000;
 double rms1=-1000;
 std::vector<double> w1;
 std::vector<double> t1;
 w1.reserve(221560);
 t1.reserve(221560);
 Analyzer *a1=new Analyzer();

 double Pmax2=-1000;
 double Pmax2Fit=-1000;
 double negPmax2Fit=-1000;
 double Tmax2=-1000;
 double Tmax2Fit=-1000;
 double negTmax2Fit=-1000;
 double Area2=-1000;
 double UArea2=-1000;
 double Area2_new=-1000;
 double UArea2_new=-1000;
 double RiseTime2Fit=-1000;
 double FallTime2Fit=-1000;
 double dVdt2Fit=-1000;
 double dVdt2Fit_2080=-1000;
 double CFD2Fit[7]={-1000,-1000,-1000,-1000,-1000,-1000,-1000};
 double WIDTH2[7]={-1000,-1000,-1000,-1000,-1000,-1000,-1000};
 double t_thr2=-1000;
 double rms2=-1000;
 std::vector<double> w2;
 std::vector<double> t2;
 w2.reserve(221560);
 t2.reserve(221560);
 Analyzer *a2=new Analyzer();

 double Pmax3=-1000;
 double Pmax3Fit=-1000;
 double negPmax3Fit=-1000;
 double Tmax3=-1000;
 double Tmax3Fit=-1000;
 double negTmax3Fit=-1000;
 double Area3=-1000;
 double UArea3=-1000;
 double Area3_new=-1000;
 double UArea3_new=-1000;
 double RiseTime3Fit=-1000;
 double FallTime3Fit=-1000;
 double dVdt3Fit=-1000;
 double dVdt3Fit_2080=-1000;
 double CFD3Fit[7]={-1000,-1000,-1000,-1000,-1000,-1000,-1000};
 double WIDTH3[7]={-1000,-1000,-1000,-1000,-1000,-1000,-1000};
 double t_thr3=-1000;
 double rms3=-1000;
 std::vector<double> w3;
 std::vector<double> t3;
 w3.reserve(221560);
 t3.reserve(221560);
 Analyzer *a3=new Analyzer();

 double Pmax4=-1000;
 double Pmax4Fit=-1000;
 double negPmax4Fit=-1000;
 double Tmax4=-1000;
 double Tmax4Fit=-1000;
 double negTmax4Fit=-1000;
 double Area4=-1000;
 double UArea4=-1000;
 double Area4_new=-1000;
 double UArea4_new=-1000;
 double RiseTime4Fit=-1000;
 double FallTime4Fit=-1000;
 double dVdt4Fit=-1000;
 double dVdt4Fit_2080=-1000;
 double CFD4Fit[7]={-1000,-1000,-1000,-1000,-1000,-1000,-1000};
 double WIDTH4[7]={-1000,-1000,-1000,-1000,-1000,-1000,-1000};
 double t_thr4=-1000;
 double rms4=-1000;
 std::vector<double> w4;
 std::vector<double> t4;
 w4.reserve(221560);
 t4.reserve(221560);
 Analyzer *a4=new Analyzer();

 int ntrig,event;
 //double x_pos,y_pos;

 OutTree->Branch("event",&event);
 OutTree->Branch("ntrig",&ntrig);

 //OutTree->Branch("x_pos",&x_pos);
 //OutTree->Branch("y_pos",&y_pos);

// Definition of Active channels
 int enable_channel_1 = cf.Value("ACTIVE_CHANNEL","ch1");
 int enable_channel_2 = cf.Value("ACTIVE_CHANNEL","ch2");
 int enable_channel_3 = cf.Value("ACTIVE_CHANNEL","ch3");
 int enable_channel_4 = cf.Value("ACTIVE_CHANNEL","ch4");

 int invert_channel_1 = cf.Value("INVERT_SIGNAL","ch1");
 int invert_channel_2 = cf.Value("INVERT_SIGNAL","ch2");
 int invert_channel_3 = cf.Value("INVERT_SIGNAL","ch3");
 int invert_channel_4 = cf.Value("INVERT_SIGNAL","ch4");

 int counter = 0;

//Output branches declaration 	 
 if( enable_channel_1 == 1){

	 OutTree->Branch("w1",&w1);
	 OutTree->Branch("t1",&t1);
	 //OutTree->Branch("Pmax1",&Pmax1,"Pmax1/D");
	 OutTree->Branch("pmax1",&Pmax1Fit,"Pmax1Fit/D");
	 OutTree->Branch("negpmax1",&negPmax1Fit,"negPmax1Fit/D");
	 //OutTree->Branch("Tmax1",&Tmax1,"Tmax1/D");
	 OutTree->Branch("tmax1",&Tmax1Fit,"Tmax1Fit/D");
	 OutTree->Branch("negtmax1",&negTmax1Fit,"negTmax1Fit/D");
	 OutTree->Branch("area1",&Area1,"Area1/D");
	 OutTree->Branch("Uarea1",&UArea1,"UArea1/D");
	 OutTree->Branch("area1_new",&Area1_new,"Area1_new/D");
	 OutTree->Branch("uarea1_new",&UArea1_new,"UArea1_new/D");
	 OutTree->Branch("risetime1",&RiseTime1Fit,"RiseTime1Fit/D");
	 OutTree->Branch("falltime1",&FallTime1Fit,"FallTime1Fit/D");
	 OutTree->Branch("dvdt1",&dVdt1Fit,"dVdt1Fit/D");
	 OutTree->Branch("dvdt1_2080",&dVdt1Fit_2080,"dVdt1Fit_2080/D");
	 OutTree->Branch("cfd1",&CFD1Fit,"CFD1Fit[7]/D");
	 OutTree->Branch("width1",&WIDTH1,"WIDTH1[7]/D");
	 OutTree->Branch("t_thr1",&t_thr1,"t_thr1/D");  // time at which a certain thr (in V) is passed
	 OutTree->Branch("rms1",&rms1,"rms1/D");

	 counter++;
 
 	}


 if( enable_channel_2 == 1){

	 OutTree->Branch("w2",&w2);
	 OutTree->Branch("t2",&t2);
	 //OutTree->Branch("Pmax2",&Pmax2,"Pmax2/D");
	 OutTree->Branch("pmax2",&Pmax2Fit,"Pmax2Fit/D");
	 OutTree->Branch("negpmax2",&negPmax2Fit,"negPmax2Fit/D");
	 //OutTree->Branch("Tmax2",&Tmax2,"Tmax2/D");
	 OutTree->Branch("tmax2",&Tmax2Fit,"Tmax2Fit/D");
	 OutTree->Branch("negtmax2",&negTmax2Fit,"negTmax2Fit/D");
	 OutTree->Branch("area2",&Area2,"Area2/D");
	 OutTree->Branch("uarea2",&UArea2,"UArea2/D");
	 OutTree->Branch("area2_new",&Area2_new,"Area2_new/D");
	 OutTree->Branch("uarea2_new",&UArea2_new,"UArea2_new/D");
	 OutTree->Branch("risetime2",&RiseTime2Fit,"RiseTime2Fit/D");
	 OutTree->Branch("falltime2",&FallTime2Fit,"FallTime2Fit/D");
	 OutTree->Branch("dvdt2",&dVdt2Fit,"dVdt2Fit/D");
	 OutTree->Branch("dvdt2_2080",&dVdt2Fit_2080,"dVdt2Fit_2080/D");
	 OutTree->Branch("cfd2",&CFD2Fit,"CFD2Fit[7]/D");
	 OutTree->Branch("width2",&WIDTH2,"WIDTH2[7]/D");
	 OutTree->Branch("t_thr2",&t_thr2,"t_thr2/D");  // time at which a certain thr (in V) is passed
	 OutTree->Branch("rms2",&rms2,"rms2/D");

	 counter++;

 	}


 if( enable_channel_3 == 1){

	 OutTree->Branch("w3",&w3);
	 OutTree->Branch("t3",&t3);
	 //OutTree->Branch("Pmax3",&Pmax3,"Pmax3/D");
	 OutTree->Branch("pmax3",&Pmax3Fit,"Pmax3Fit/D");
	 OutTree->Branch("negpmax3",&negPmax3Fit,"negPmax3Fit/D");
	 //OutTree->Branch("Tmax3",&Tmax3,"Tmax3/D");
	 OutTree->Branch("tmax3",&Tmax3Fit,"Tmax3Fit/D");
	 OutTree->Branch("negtmax3",&negTmax3Fit,"negTmax3Fit/D");
	 OutTree->Branch("area3",&Area3,"Area3/D");
	 OutTree->Branch("uarea3",&UArea3,"UArea3/D");
	 OutTree->Branch("area3_new",&Area3_new,"Area3_new/D");
	 OutTree->Branch("uarea3_new",&UArea3_new,"UArea3_new/D");
	 OutTree->Branch("risetime3",&RiseTime3Fit,"RiseTime3Fit/D");
	 OutTree->Branch("falltime3",&FallTime3Fit,"FallTime3Fit/D");
	 OutTree->Branch("dvdt3",&dVdt3Fit,"dVdt3Fit/D");
	 OutTree->Branch("dvdt3_2080",&dVdt3Fit_2080,"dVdt3Fit_2080/D");
	 OutTree->Branch("cfd3",&CFD3Fit,"CFD3Fit[7]/D");
	 OutTree->Branch("width3",&WIDTH3,"WIDTH3[7]/D");
	 OutTree->Branch("t_thr3",&t_thr3,"t_thr3/D");  // time at which a certain thr (in V) is passed
	 OutTree->Branch("rms3",&rms3,"rms3/D");

	 counter++;

 	}


 if( enable_channel_4 == 1){

	 OutTree->Branch("w4",&w4);
	 OutTree->Branch("t4",&t4);
	 //OutTree->Branch("Pmax4",&Pmax4,"Pmax4/D");
	 OutTree->Branch("pmax4",&Pmax4Fit,"Pmax4Fit/D");
	 OutTree->Branch("negpmax4",&negPmax4Fit,"negPmax4Fit/D");
	 //OutTree->Branch("Tmax4",&Tmax4,"Tmax4/D");
	 OutTree->Branch("tmax4",&Tmax4Fit,"Tmax4Fit/D");
	 OutTree->Branch("negtmax4",&negTmax4Fit,"negTmax4Fit/D");
	 OutTree->Branch("area4",&Area4,"Area4/D");
	 OutTree->Branch("uarea4",&UArea4,"UArea4/D");
	 OutTree->Branch("area4_new",&Area4_new,"Area4_new/D");
	 OutTree->Branch("uarea4_new",&UArea4_new,"UArea4_new/D");
	 OutTree->Branch("risetime4",&RiseTime4Fit,"RiseTime4Fit/D");
	 OutTree->Branch("falltime4",&FallTime4Fit,"FallTime4Fit/D");
	 OutTree->Branch("dvdt4",&dVdt4Fit,"dVdt4Fit/D");
	 OutTree->Branch("dvdt4_2080",&dVdt4Fit_2080,"dVdt4Fit_2080/D");
	 OutTree->Branch("cfd4",&CFD4Fit,"CFD4Fit[7]/D");
	 OutTree->Branch("width4",&WIDTH4,"WIDTH4[7]/D");
	 OutTree->Branch("t_thr4",&t_thr4,"t_thr4/D");  // time at which a certain thr (in V) is passed
	 OutTree->Branch("rms4",&rms4,"rms4/D");

	 counter++;

 	}
	


 int j = 0;

//Loop over raw data and fill output branches. Analyzer class manages output variables calculation 
 while(myReader.Next()){

 	if(j==0){

 		bool act_ch = false;

 		if(enable_channel_1 == 1){

 			time_window[0] = timeReader1.At(0);
 			time_window[1] = timeReader1.At(timeReader1.GetSize()-1);
 			search_range[0] = (time_window[0]+time_window[1])/2. - 100*(timeReader1.At(1)-timeReader1.At(0)) ;
 			search_range[1] = (time_window[0]+time_window[1])/2. + 100*(timeReader1.At(1)-timeReader1.At(0)) ;
 			act_ch = true;

 			cout<<"Time window: "<<time_window[0]<<"; "<<time_window[1]<<endl;
 			cout<<"Search window: "<<search_range[0]<<"; "<<search_range[1]<<endl;

 		}

 		if(enable_channel_2 == 1 && !act_ch){

 			time_window[0] = timeReader2.At(0);
 			time_window[1] = timeReader2.At(timeReader2.GetSize()-1);
 			search_range[0] = (time_window[0]+time_window[1])/2. - 100*(timeReader2.At(1)-timeReader2.At(0)) ;
 			search_range[1] = (time_window[0]+time_window[1])/2. + 100*(timeReader2.At(1)-timeReader2.At(0)) ;
 			act_ch = true;

 			cout<<"Time window: "<<time_window[0]<<"; "<<time_window[1]<<endl;
 			cout<<"Search window: "<<search_range[0]<<"; "<<search_range[1]<<endl;

 		}

 		if(enable_channel_3 == 1 && !act_ch){

 			time_window[0] = timeReader3.At(0);
 			time_window[1] = timeReader3.At(timeReader3.GetSize()-1);
 			search_range[0] = (time_window[0]+time_window[1])/2. - 100*(timeReader3.At(1)-timeReader3.At(0)) ;
 			search_range[1] = (time_window[0]+time_window[1])/2. + 100*(timeReader3.At(1)-timeReader3.At(0)) ;
 			act_ch = true;

 			cout<<"Time window: "<<time_window[0]<<"; "<<time_window[1]<<endl;
 			cout<<"Search window: "<<search_range[0]<<"; "<<search_range[1]<<endl;

 		}

 		if(enable_channel_4 == 1 && !act_ch){

 			time_window[0] = timeReader4.At(0);
 			time_window[1] = timeReader4.At(timeReader4.GetSize()-1);
 			search_range[0] = (time_window[0]+time_window[1])/2. - 100*(timeReader4.At(1)-timeReader4.At(0)) ;
 			search_range[1] = (time_window[0]+time_window[1])/2. + 100*(timeReader4.At(1)-timeReader4.At(0)) ;
 			act_ch = true;

 			cout<<"Time window: "<<time_window[0]<<"; "<<time_window[1]<<endl;
 			cout<<"Search window: "<<search_range[0]<<"; "<<search_range[1]<<endl;

 		}


 	}

 	int counter2 = 0;

 	if( enable_channel_1 == 1){

		w1.clear();
	 	t1.clear();

 		if( invert_channel_1 == 1 ){

	 		for(unsigned int i=0; i<voltageReader1.GetSize();i++){

	 			w1.push_back(-voltageReader1.At(i));
	 			t1.push_back(timeReader1.At(i));

 			}

 		}else{

	 		for(unsigned int i=0; i<voltageReader1.GetSize();i++){

	 			w1.push_back(voltageReader1.At(i));
	 			t1.push_back(timeReader1.At(i));

 			}

 			}

 		if(w1.size()<maxIndex || t1.size()<maxIndex){

 			cout<<"Voltage or Time vector less than 1000 entries. Skipping whole event"<<endl;
 			continue;
 		}

 		if(w1.size()==0 || t1.size()==0){

 			cout<<"Voltage or Time vector empty. Skipping whole event"<<endl;
 			continue;
 		}

 		if(w1.size()!= t1.size()){

 			cout<<"Different number of entries in Voltage and Time vectors. Skipping whole event"<<endl;
 			continue;
 		}


		*a1=Analyzer(w1,t1);  //to be checked: is this a correct usage of memory?
		a1->Correct_Baseline(100);

		std::pair<double, unsigned int> tp_pair1 = a1->Find_Signal_Maximum(false,search_range); 
		std::pair<double, double> tp_pair1_fit = a1->Pmax_with_GausFit(tp_pair1,maxIndex);
		std::pair<double, unsigned int> neg_tp_pair1 = a1->Find_Negative_Signal_Maximum(false,search_range); 
		std::pair<double, double> neg_tp_pair1_fit = a1->Negative_Pmax_with_GausFit(neg_tp_pair1,maxIndex);      
		Pmax1 = tp_pair1.first*voltage_const; //mV
		Pmax1Fit = tp_pair1_fit.first*voltage_const; //mV
		negPmax1Fit = neg_tp_pair1_fit.first*voltage_const; //mV
		Tmax1 = a1->Get_Tmax(tp_pair1)*time_const; //ns
		Tmax1Fit = tp_pair1_fit.second*time_const; //ns
		negTmax1Fit = neg_tp_pair1_fit.second*time_const; //ns
		Area1 = a1->Find_Pulse_Area(tp_pair1)*voltage_const*time_const; //mV*ns
		UArea1 = a1->Find_Undershoot_Area(tp_pair1)*voltage_const*time_const; //mV*ns
		Area1_new = a1->New_Pulse_Area(tp_pair1_fit,tp_pair1.second,"Simpson",search_range, time_window[0], time_window[1])*voltage_const*time_const;//mV*ns  Similar to "Area1" method, but start/end times of the pulse are defined analitically. Should be the most accurate definition
		UArea1_new = a1->New_Undershoot_Area(tp_pair1_fit,neg_tp_pair1_fit, neg_tp_pair1.second,"Simpson",search_range)*voltage_const*time_const;//mV*ns
		RiseTime1Fit = a1->Find_Rise_Time_with_GausFit(tp_pair1_fit, tp_pair1.second, 0.1, 0.9)*time_const; //ns
		FallTime1Fit = a1->Find_Fall_Time_with_GausFit(tp_pair1_fit, tp_pair1.second, 0.1, 0.9)*time_const; //ns
		dVdt1Fit = a1->Find_Dvdt_with_GausFit(20,0,tp_pair1_fit,tp_pair1.second)*(voltage_const/time_const);  //mV/ns
		dVdt1Fit_2080 = a1->Find_Dvdt2080_with_GausFit(0,tp_pair1_fit,tp_pair1.second)*(voltage_const/time_const);  //mV/ns
		t_thr1 = a1->Find_Time_At_Threshold(20,tp_pair1)*time_const; //ns    [time at which a thr (in mV !!!) is passed]
		rms1 = a1->Find_Noise(100)*voltage_const; //mV

		for(int j=0;j<7;j++){

			CFD1Fit[j]=a1->Rising_Edge_CFD_Time_with_GausFit(10+j*10,tp_pair1_fit,tp_pair1.second)*time_const;
			WIDTH1[j]=(a1->Falling_Edge_CFD_Time_with_GausFit(10+j*10,tp_pair1_fit,tp_pair1.second)*time_const) - CFD1Fit[j];

		}

		counter2++;

	}


	if( enable_channel_2 == 1){

		w2.clear();
	 	t2.clear();

 		if( invert_channel_2 == 1){

	 		for(unsigned int i=0; i<voltageReader2.GetSize();i++){

	 			w2.push_back(-voltageReader2.At(i));
	 			t2.push_back(timeReader2.At(i));

 			}

 		}else{

	 		for(unsigned int i=0; i<voltageReader2.GetSize();i++){

	 			w2.push_back(voltageReader2.At(i));
	 			t2.push_back(timeReader2.At(i));

 			}

 			}

 		if(w2.size()<maxIndex || t2.size()<maxIndex){

 			cout<<"Voltage or Time vector less than 1000 entries. Skipping whole event"<<endl;
 			continue;
 		}

 		if(w2.size()==0 || t2.size()==0){

 			cout<<"Voltage or Time vector empty. Skipping whole event"<<endl;
 			continue;
 		}

 		if(w2.size()!= t2.size()){

 			cout<<"Different number of entries in Voltage and Time vectors. Skipping whole event"<<endl;
 			continue;
 		}

		*a2=Analyzer(w2,t2);  
		a2->Correct_Baseline(100);

		std::pair<double, unsigned int> tp_pair2 = a2->Find_Signal_Maximum(false,search_range); 
		std::pair<double, double> tp_pair2_fit = a2->Pmax_with_GausFit(tp_pair2,maxIndex);
		std::pair<double, unsigned int> neg_tp_pair2 = a2->Find_Negative_Signal_Maximum(false,search_range); 
		std::pair<double, double> neg_tp_pair2_fit = a2->Negative_Pmax_with_GausFit(neg_tp_pair2,maxIndex);      
		Pmax2 = tp_pair2.first*voltage_const; //mV
		Pmax2Fit = tp_pair2_fit.first*voltage_const; //mV
		negPmax2Fit = neg_tp_pair2_fit.first*voltage_const; //mV
		Tmax2 = a2->Get_Tmax(tp_pair2)*time_const; //ns
		Tmax2Fit = tp_pair2_fit.second*time_const; //ns
		negTmax2Fit = neg_tp_pair2_fit.second*time_const; //ns
		Area2 = a2->Find_Pulse_Area(tp_pair2)*voltage_const*time_const;
		UArea2 = a2->Find_Undershoot_Area(tp_pair2)*voltage_const*time_const;
		Area2_new = a2->New_Pulse_Area(tp_pair2_fit,tp_pair2.second,"Simpson",search_range, time_window[0], time_window[1])*voltage_const*time_const;
		UArea2_new = a2->New_Undershoot_Area(tp_pair2_fit,neg_tp_pair2_fit, neg_tp_pair2.second,"Simpson",search_range)*voltage_const*time_const;//mV*ns
		RiseTime2Fit = a2->Find_Rise_Time_with_GausFit(tp_pair2_fit, tp_pair2.second, 0.1, 0.9)*time_const; //ns
		FallTime2Fit = a2->Find_Fall_Time_with_GausFit(tp_pair2_fit, tp_pair2.second, 0.1, 0.9)*time_const; //ns
		dVdt2Fit = a2->Find_Dvdt_with_GausFit(20,0,tp_pair2_fit,tp_pair2.second)*(voltage_const/time_const);
		dVdt2Fit_2080 = a2->Find_Dvdt2080_with_GausFit(0,tp_pair2_fit,tp_pair2.second)*(voltage_const/time_const);  //mV/ns
		t_thr2 = a2->Find_Time_At_Threshold(0.02,tp_pair2)*time_const;
		rms2 = a2->Find_Noise(100)*voltage_const;

		for(int j=0;j<7;j++){

			CFD2Fit[j]=a2->Rising_Edge_CFD_Time(10+j*10,tp_pair2)*time_const;
			WIDTH2[j]=a2->Falling_Edge_CFD_Time_with_GausFit(10+j*10,tp_pair2_fit,tp_pair2.second)*time_const - CFD2Fit[j];


		}

		

		counter2++;

	}


	if( enable_channel_3 == 1){

		w3.clear();
	 	t3.clear();

 		if( invert_channel_3 == 1){

	 		for(unsigned int i=0; i<voltageReader3.GetSize();i++){

	 			w3.push_back(-voltageReader3.At(i));
	 			t3.push_back(timeReader3.At(i));

 			}

 		}else{

	 		for(unsigned int i=0; i<voltageReader3.GetSize();i++){

	 			w3.push_back(voltageReader3.At(i));
	 			t3.push_back(timeReader3.At(i));

 			}

 			}

 		if(w3.size()<maxIndex || t3.size()<maxIndex){

 			cout<<"Voltage or Time vector less than 1000 entries. Skipping whole event"<<endl;
 			continue;
 		}

 		if(w3.size()==0 || t3.size()==0){

 			cout<<"Voltage or Time vector empty. Skipping whole event"<<endl;
 			continue;
 		}

 		if(w3.size()!= t3.size()){

 			cout<<"Different number of entries in Voltage and Time vectors. Skipping whole event"<<endl;
 			continue;
 		}

		*a3 = Analyzer(w3,t3);  
		a3->Correct_Baseline(100);

		std::pair<double, unsigned int> tp_pair3 = a3->Find_Signal_Maximum(false,search_range); 
		std::pair<double, double> tp_pair3_fit = a3->Pmax_with_GausFit(tp_pair3,maxIndex);
		std::pair<double, unsigned int> neg_tp_pair3 = a3->Find_Negative_Signal_Maximum(false,search_range); 
		std::pair<double, double> neg_tp_pair3_fit = a3->Negative_Pmax_with_GausFit(neg_tp_pair3,maxIndex);      
		Pmax3 = tp_pair3.first*voltage_const; //mV
		Pmax3Fit = tp_pair3_fit.first*voltage_const; //mV
		negPmax3Fit = neg_tp_pair3_fit.first*voltage_const; //mV
		Tmax3 = a3->Get_Tmax(tp_pair3)*time_const; //ns
		Tmax3Fit = tp_pair3_fit.second*time_const; //ns
		negTmax3Fit = neg_tp_pair3_fit.second*time_const; //ns
		Area3 = a3->Find_Pulse_Area(tp_pair3)*voltage_const*time_const;
		UArea3 = a3->Find_Undershoot_Area(tp_pair3)*voltage_const*time_const;
		Area3_new = a3->New_Pulse_Area(tp_pair3_fit,tp_pair3.second,"Simpson",search_range, time_window[0], time_window[1])*voltage_const*time_const;
		UArea3_new = a3->New_Undershoot_Area(tp_pair3_fit,neg_tp_pair3_fit, neg_tp_pair3.second,"Simpson",search_range)*voltage_const*time_const;//mV*ns
		RiseTime3Fit = a3->Find_Rise_Time_with_GausFit(tp_pair3_fit, tp_pair3.second, 0.1, 0.9)*time_const; //ns
		FallTime3Fit = a3->Find_Fall_Time_with_GausFit(tp_pair3_fit, tp_pair3.second, 0.1, 0.9)*time_const; //ns
		dVdt3Fit = a3->Find_Dvdt_with_GausFit(20,0,tp_pair3_fit,tp_pair3.second)*(voltage_const/time_const);
		dVdt3Fit_2080 = a3->Find_Dvdt2080_with_GausFit(0,tp_pair3_fit,tp_pair3.second)*(voltage_const/time_const);  //mV/ns
		t_thr3 = a3->Find_Time_At_Threshold(0.02,tp_pair3)*time_const;
		rms3 = a3->Find_Noise(100)*voltage_const;

		for(int j=0;j<7;j++){

			CFD3Fit[j]=a3->Rising_Edge_CFD_Time(10+j*10,tp_pair3)*time_const;
			WIDTH3[j]=a3->Falling_Edge_CFD_Time_with_GausFit(10+j*10,tp_pair3_fit,tp_pair3.second)*time_const - CFD3Fit[j];

		}

		counter2++;

	}


	if( enable_channel_4 == 1){

		w4.clear();
	 	t4.clear();

 		if( invert_channel_4 == 1){

	 		for(unsigned int i=0; i<voltageReader4.GetSize();i++){

	 			w4.push_back(-voltageReader4.At(i));
	 			t4.push_back(timeReader4.At(i));

 			}

 		}else{

	 		for(unsigned int i=0; i<voltageReader4.GetSize();i++){

	 			w4.push_back(voltageReader4.At(i));
	 			t4.push_back(timeReader4.At(i));

 			}

 			}

 		if(w4.size()<maxIndex || t4.size()<maxIndex){

 			cout<<"Voltage or Time vector less than 1000 entries. Skipping whole event"<<endl;
 			continue;
 		}

 		if(w4.size()==0 || t4.size()==0){

 			cout<<"Voltage or Time vector empty. Skipping whole event"<<endl;
 			continue;
 		}

 		if(w4.size()!= t4.size()){

 			cout<<"Different number of entries in Voltage and Time vectors. Skipping whole event"<<endl;
 			continue;
 		}

		*a4 = Analyzer(w4,t4);  
		a4->Correct_Baseline(100);

		std::pair<double, unsigned int> tp_pair4 = a4->Find_Signal_Maximum(false,search_range); 
		std::pair<double, double> tp_pair4_fit = a4->Pmax_with_GausFit(tp_pair4,maxIndex);
		std::pair<double, unsigned int> neg_tp_pair4 = a4->Find_Negative_Signal_Maximum(false,search_range); 
		std::pair<double, double> neg_tp_pair4_fit = a4->Negative_Pmax_with_GausFit(neg_tp_pair4,maxIndex);      
		Pmax4 = tp_pair4.first*voltage_const; //mV
		Pmax4Fit = tp_pair4_fit.first*voltage_const; //mV
		negPmax4Fit = neg_tp_pair4_fit.first*voltage_const; //mV
		Tmax4 = a4->Get_Tmax(tp_pair4)*time_const; //ns
		Tmax4Fit = tp_pair4_fit.second*time_const; //ns
		negTmax4Fit = neg_tp_pair4_fit.second*time_const; //ns
		Area4 = a4->Find_Pulse_Area(tp_pair4)*voltage_const*time_const;
		UArea4 = a4->Find_Undershoot_Area(tp_pair4)*voltage_const*time_const;
		Area4_new = a4->New_Pulse_Area(tp_pair4_fit,tp_pair4.second,"Simpson",search_range, time_window[0], time_window[1])*voltage_const*time_const;
		UArea4_new = a4->New_Undershoot_Area(tp_pair4_fit,neg_tp_pair4_fit, neg_tp_pair4.second,"Simpson",search_range)*voltage_const*time_const;//mV*ns
		RiseTime4Fit = a4->Find_Rise_Time_with_GausFit(tp_pair4_fit, tp_pair4.second, 0.1, 0.9)*time_const; //ns
		FallTime4Fit = a4->Find_Fall_Time_with_GausFit(tp_pair4_fit, tp_pair4.second, 0.1, 0.9)*time_const; //ns
		dVdt4Fit = a4->Find_Dvdt_with_GausFit(20,0,tp_pair4_fit,tp_pair4.second)*(voltage_const/time_const); 
		dVdt4Fit_2080 = a4->Find_Dvdt2080_with_GausFit(0,tp_pair4_fit,tp_pair4.second)*(voltage_const/time_const);  //mV/ns
		t_thr4 = a4->Find_Time_At_Threshold(0.02,tp_pair4)*time_const;
		rms4 = a4->Find_Noise(100)*voltage_const;

		for(int j=0;j<7;j++){

			CFD4Fit[j]=a4->Rising_Edge_CFD_Time(10+j*10,tp_pair4)*time_const;
			WIDTH1[j]=a4->Falling_Edge_CFD_Time_with_GausFit(10+j*10,tp_pair4_fit,tp_pair4.second)*time_const - CFD4Fit[j];


		}

	}

	ntrig=j;
	event=j;

	//x_pos = ((*xReader*1000-17438.5)*0.999951+(*yReader*1000-23197.11)*(9/195));
 	//y_pos = (0.998934*(*yReader*1000-23197.11)+(*xReader*1000-17438.5)*(2/201));

 	//x_pos = *xReader;
 	//y_pos = *yReader;
		
	OutTree->Fill();

	if(j%10000 == 0) cout<<"processed events:"<<j<<endl;
	//cout<<"processed events:"<<j<<endl;

	j++;
		

	}

	
 OutTree->Write();
 OutputFile->Write();

}


int main(){

analisi();

return 0;

}