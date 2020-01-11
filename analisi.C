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
#include "Analyzer.hpp"
#include "include/ConfigFile.hpp"


void analisi(char *filename, char *output_filename){

 //gStyle->SetOptFit(0000);


//Config file definition
 ConfigFile cf("beta_config.ini");

//Input file
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
// End of TTreeReaders definition

// Output file & tree
 TFile *OutputFile = new TFile(output_filename,"recreate");
 TTree *OutTree = new TTree("Analysis","Analysis");
 
// Variable declaration and Analyzer objects 
// Raw data have Voltage in V and Time in ns
 const double time_const = cf.Value("HEADER","time_scalar");  //multiply by this const to pass from s to ns
 const double voltage_const = cf.Value("HEADER","voltage_scalar"); // pass from V to mV

 double Pmax1=-1000;
 double Pmax1Fit=-1000;
 double negPmax1=-1000;
 double negPmax1Fit=-1000;
 double Tmax1=-1000;
 double negTmax1=-1000;
 double Tmax1Fit=-1000;
 double negTmax1Fit=-1000;
 double Area1=-1000;
 double UArea1=-1000;
 double AreaFixed1=-1000;
 double AreaFixed1Fit=-1000;
 double AreaInterpolated1=-1000;
 double AreaInterpolated1Fit=-1000;
 double Area1_new=-1000;
 double RiseTime1=-1000;
 double RiseTime1Fit=-1000;
 double dVdt1=-1000;
 double dVdt1Fit=-1000;
 double cfd1=-1000;
 double cfd1Fit=-1000;
 double CFD1Fit[7]={-1000,-1000,-1000,-1000,-1000,-1000,-1000};
 double t_thr1=-1000;
 double rms1=-1000;
 std::vector<double> w1o;
 std::vector<double> t1o;
 w1o.reserve(221560);
 t1o.reserve(221560);
 Analyzer *a1=new Analyzer();

 double Pmax2=-1000;
 double Pmax2Fit=-1000;
 double negPmax2=-1000;
 double negPmax2Fit=-1000;
 double Tmax2=-1000;
 double negTmax2=-1000;
 double Tmax2Fit=-1000;
 double negTmax2Fit=-1000;
 double Area2=-1000;
 double UArea2=-1000;
 double AreaFixed2=-1000;
 double AreaFixed2Fit=-1000;
 double AreaInterpolated2=-1000;
 double AreaInterpolated2Fit=-1000;
 double Area2_new=-1000;
 double RiseTime2=-1000;
 double RiseTime2Fit=-1000;
 double dVdt2=-1000;
 double dVdt2Fit=-1000;
 double cfd2=-1000;
 double cfd2Fit=-1000;
 double CFD2Fit[7]={-1000,-1000,-1000,-1000,-1000,-1000,-1000};
 double t_thr2=-1000;
 double rms2=-1000;
 std::vector<double> w2o;
 std::vector<double> t2o;
 w2o.reserve(221560);
 t2o.reserve(221560);
 Analyzer *a2=new Analyzer();

 

 

 double Pmax3=-1000;
 double Pmax3Fit=-1000;
 double negPmax3=-1000;
 double negPmax3Fit=-1000;
 double Tmax3=-1000;
 double negTmax3=-1000;
 double Tmax3Fit=-1000;
 double negTmax3Fit=-1000;
 double Area3=-1000;
 double UArea3=-1000;
 double AreaFixed3=-1000;
 double AreaFixed3Fit=-1000;
 double AreaInterpolated3=-1000;
 double AreaInterpolated3Fit=-1000;
 double Area3_new=-1000;
 double RiseTime3=-1000;
 double RiseTime3Fit=-1000;
 double dVdt3=-1000;
 double dVdt3Fit=-1000;
 double cfd3=-1000;
 double cfd3Fit=-1000;
 double CFD3Fit[7]={-1000,-1000,-1000,-1000,-1000,-1000,-1000};
 double t_thr3=-1000;
 double rms3=-1000;
 std::vector<double> w3o;
 std::vector<double> t3o;
 w3o.reserve(221560);
 t3o.reserve(221560);
 Analyzer *a3=new Analyzer();

 double Pmax4=-1000;
 double Pmax4Fit=-1000;
 double negPmax4=-1000;
 double negPmax4Fit=-1000;
 double Tmax4=-1000;
 double negTmax4=-1000;
 double Tmax4Fit=-1000;
 double negTmax4Fit=-1000;
 double Area4=-1000;
 double UArea4=-1000;
 double AreaFixed4=-1000;
 double AreaFixed4Fit=-1000;
 double AreaInterpolated4=-1000;
 double AreaInterpolated4Fit=-1000;
 double Area4_new=-1000;
 double RiseTime4=-1000;
 double RiseTime4Fit=-1000;
 double dVdt4=-1000;
 double dVdt4Fit=-1000;
 double cfd4=-1000;
 double cfd4Fit=-1000;
 double CFD4Fit[7]={-1000,-1000,-1000,-1000,-1000,-1000,-1000};
 double t_thr4=-1000;
 double rms4=-1000;
 std::vector<double> w4o;
 std::vector<double> t4o;
 w4o.reserve(221560);
 t4o.reserve(221560);
 Analyzer *a4=new Analyzer();

 int ntrig,event;
 OutTree->Branch("event",&event);
 OutTree->Branch("ntrig",&ntrig);

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

	 OutTree->Branch("w1o",&w1o);
	 OutTree->Branch("t1o",&t1o);
	 OutTree->Branch("Pmax1",&Pmax1,"Pmax1/D");
	 OutTree->Branch("Pmax1Fit",&Pmax1Fit,"Pmax1Fit/D");
	 OutTree->Branch("negPmax1",&negPmax1,"negPmax1/D");
	 OutTree->Branch("negPmax1Fit",&negPmax1Fit,"negPmax1Fit/D");
	 OutTree->Branch("Tmax1",&Tmax1,"Tmax1/D");
	 OutTree->Branch("negTmax1",&negTmax1,"negTmax1/D");
	 OutTree->Branch("Tmax1Fit",&Tmax1Fit,"Tmax1Fit/D");
	 OutTree->Branch("negTmax1Fit",&negTmax1Fit,"negTmax1Fit/D");
	 OutTree->Branch("Area1",&Area1,"Area1/D");
	 OutTree->Branch("UArea1",&UArea1,"UArea1/D");
	 OutTree->Branch("AreaFixed1",&AreaFixed1,"AreaFixed1/D");
	 OutTree->Branch("AreaFixed1Fit",&AreaFixed1Fit,"AreaFixed1Fit/D");
	 OutTree->Branch("AreaInterpolated1",&AreaInterpolated1,"AreaInterpolated1/D");
	 OutTree->Branch("AreaInterpolated1Fit",&AreaInterpolated1Fit,"AreaInterpolated1Fit/D");
	 OutTree->Branch("Area1_new",&Area1_new,"Area1_new/D");
	 OutTree->Branch("RiseTime1",&RiseTime1,"RiseTime1/D");
	 OutTree->Branch("RiseTime1Fit",&RiseTime1Fit,"RiseTime1Fit/D");
	 OutTree->Branch("dVdt1",&dVdt1,"dVdt1/D");
	 OutTree->Branch("dVdt1Fit",&dVdt1Fit,"dVdt1Fit/D");
	 OutTree->Branch("cfd1",&cfd1,"cfd1/D");
	 OutTree->Branch("cfd1Fit",&cfd1Fit,"cfd1Fit/D");
	 OutTree->Branch("CFD1Fit",&CFD1Fit,"CFD1Fit[7]/D");
	 OutTree->Branch("t_thr1",&t_thr1,"t_thr1/D");  // time at which a certain thr (in V) is passed
	 OutTree->Branch("rms1",&rms1,"rms1/D");

	 counter++;
 
 	}


 if( enable_channel_2 == 1){

	 OutTree->Branch("w2o",&w2o);
	 OutTree->Branch("t2o",&t2o);
	 OutTree->Branch("Pmax2",&Pmax2,"Pmax2/D");
	 OutTree->Branch("Pmax2Fit",&Pmax2Fit,"Pmax2Fit/D");
	 OutTree->Branch("Pmax2",&negPmax2,"negPmax2/D");
	 OutTree->Branch("negPmax2Fit",&negPmax2Fit,"negPmax2Fit/D");
	 OutTree->Branch("Tmax2",&Tmax2,"Tmax2/D");
	 OutTree->Branch("negTmax2",&negTmax2,"negTmax2/D");
	 OutTree->Branch("Tmax2Fit",&Tmax2Fit,"Tmax2Fit/D");
	 OutTree->Branch("negTmax2Fit",&negTmax2Fit,"negTmax2Fit/D");
	 OutTree->Branch("Area2",&Area2,"Area2/D");
	 OutTree->Branch("UArea2",&UArea2,"UArea2/D");
	 OutTree->Branch("AreaFixed2",&AreaFixed2,"AreaFixed2/D");
	 OutTree->Branch("AreaFixed2Fit",&AreaFixed2Fit,"AreaFixed2Fit/D");
	 OutTree->Branch("AreaInterpolated2",&AreaInterpolated2,"AreaInterpolated2/D");
	 OutTree->Branch("AreaInterpolated2Fit",&AreaInterpolated2Fit,"AreaInterpolated2Fit/D");
	 OutTree->Branch("Area2_new",&Area2_new,"Area2_new/D");
	 OutTree->Branch("RiseTime2",&RiseTime2,"RiseTime2/D");
	 OutTree->Branch("RiseTime2Fit",&RiseTime2Fit,"RiseTime2Fit/D");
	 OutTree->Branch("dVdt2",&dVdt2,"dVdt2/D");
	 OutTree->Branch("dVdt2Fit",&dVdt2Fit,"dVdt2Fit/D");
	 OutTree->Branch("cfd2",&cfd2,"cfd2/D");
	 OutTree->Branch("cfd2Fit",&cfd2Fit,"cfd2Fit/D");
	 OutTree->Branch("CFD2Fit",&CFD2Fit,"CFD2Fit[7]/D");
	 OutTree->Branch("t_thr2",&t_thr2,"t_thr2/D");  // time at which a certain thr (in V) is passed
	 OutTree->Branch("rms2",&rms2,"rms2/D");

	 counter++;

 	}


 if( enable_channel_3 == 1){

	 OutTree->Branch("w3o",&w3o);
	 OutTree->Branch("t3o",&t3o);
	 OutTree->Branch("Pmax3",&Pmax3,"Pmax3/D");
	 OutTree->Branch("Pmax3Fit",&Pmax3Fit,"Pmax3Fit/D");
	 OutTree->Branch("negPmax3",&negPmax3,"negPmax3/D");
	 OutTree->Branch("negPmax3Fit",&negPmax3Fit,"negPmax3Fit/D");
	 OutTree->Branch("Tmax3",&Tmax3,"Tmax3/D");
	 OutTree->Branch("negTmax3",&negTmax3,"negTmax3/D");
	 OutTree->Branch("Tmax3Fit",&Tmax3Fit,"Tmax3Fit/D");
	 OutTree->Branch("negTmax3Fit",&negTmax3Fit,"negTmax3Fit/D");
	 OutTree->Branch("Area3",&Area3,"Area3/D");
	 OutTree->Branch("UArea3",&UArea3,"UArea3/D");
	 OutTree->Branch("AreaFixed3",&AreaFixed3,"AreaFixed3/D");
	 OutTree->Branch("AreaFixed3Fit",&AreaFixed3Fit,"AreaFixed3Fit/D");
	 OutTree->Branch("AreaInterpolated3",&AreaInterpolated3,"AreaInterpolated3/D");
	 OutTree->Branch("AreaInterpolated3Fit",&AreaInterpolated3Fit,"AreaInterpolated3Fit/D");
	 OutTree->Branch("Area3_new",&Area3_new,"Area3_new/D");
	 OutTree->Branch("RiseTime3",&RiseTime3,"RiseTime3/D");
	 OutTree->Branch("RiseTime3Fit",&RiseTime3Fit,"RiseTime3Fit/D");
	 OutTree->Branch("dVdt3",&dVdt3,"dVdt3/D");
	 OutTree->Branch("dVdt3Fit",&dVdt3Fit,"dVdt3Fit/D");
	 OutTree->Branch("cfd3",&cfd3,"cfd3/D");
	 OutTree->Branch("cfd3Fit",&cfd3Fit,"cfd3Fit/D");
	 OutTree->Branch("CFD3Fit",&CFD3Fit,"CFD3Fit[7]/D");
	 OutTree->Branch("t_thr3",&t_thr3,"t_thr3/D");  // time at which a certain thr (in V) is passed
	 OutTree->Branch("rms3",&rms3,"rms3/D");

	 counter++;

 	}


 if( enable_channel_4 == 1){

	 OutTree->Branch("w4o",&w4o);
	 OutTree->Branch("t4o",&t4o);
	 OutTree->Branch("Pmax4",&Pmax4,"Pmax4/D");
	 OutTree->Branch("Pmax4Fit",&Pmax4Fit,"Pmax4Fit/D");
	 OutTree->Branch("negPmax4",&negPmax4,"negPmax4/D");
	 OutTree->Branch("negPmax4Fit",&negPmax4Fit,"negPmax4Fit/D");
	 OutTree->Branch("Tmax4",&Tmax4,"Tmax4/D");
	 OutTree->Branch("negTmax4",&negTmax4,"negTmax4/D");
	 OutTree->Branch("Tmax4Fit",&Tmax4Fit,"Tmax4Fit/D");
	 OutTree->Branch("negTmax4Fit",&negTmax4Fit,"negTmax4Fit/D");
	 OutTree->Branch("Area4",&Area4,"Area4/D");
	 OutTree->Branch("UArea4",&UArea4,"UArea4/D");
	 OutTree->Branch("AreaFixed4",&AreaFixed4,"AreaFixed4/D");
	 OutTree->Branch("AreaFixed4Fit",&AreaFixed4Fit,"AreaFixed4Fit/D");
	 OutTree->Branch("AreaInterpolated4",&AreaInterpolated4,"AreaInterpolated4/D");
	 OutTree->Branch("AreaInterpolated4Fit",&AreaInterpolated4Fit,"AreaInterpolated4Fit/D");
	 OutTree->Branch("Area4_new",&Area4_new,"Area4_new/D");
	 OutTree->Branch("RiseTime4",&RiseTime4,"RiseTime4/D");
	 OutTree->Branch("RiseTime4Fit",&RiseTime4Fit,"RiseTime4Fit/D");
	 OutTree->Branch("dVdt4",&dVdt4,"dVdt4/D");
	 OutTree->Branch("dVdt4Fit",&dVdt4Fit,"dVdt4Fit/D");
	 OutTree->Branch("cfd4",&cfd4,"cfd4/D");
	 OutTree->Branch("cfd4Fit",&cfd4Fit,"cfd4Fit/D");
	 OutTree->Branch("CFD4Fit",&CFD4Fit,"CFD4Fit[7]/D");
	 OutTree->Branch("t_thr4",&t_thr4,"t_thr4/D");  // time at which a certain thr (in V) is passed
	 OutTree->Branch("rms4",&rms4,"rms4/D");

	 counter++;

 	}
	

	

 double search_range[2] = {-2e-9,2e-9};
 double time_window[2] = {-25e-9,25e-9}; // Can be included in beta_config
 int j = 0;

//Loop over raw data and fill output branches. Analyzer class manages output variables calculation 
 while(myReader.Next()){

 	int counter2 = 0;

 	if( enable_channel_1 == 1){

		w1o.clear();
	 	t1o.clear();

 		if( invert_channel_1 == 1 ){

	 		for(unsigned int i=0; i<voltageReader1.GetSize();i++){

	 			w1o.push_back(-voltageReader1.At(i));
	 			t1o.push_back(timeReader1.At(i));

 			}

 		}else{

	 		for(unsigned int i=0; i<voltageReader1.GetSize();i++){

	 			w1o.push_back(voltageReader1.At(i));
	 			t1o.push_back(timeReader1.At(i));

 			}

 			}

 		if(w1o.size()<1000 || t1o.size()<1000){

 			cout<<"Voltage or Time vector less than 1000 entries. Skipping whole event"<<endl;
 			continue;
 		}

 		if(w1o.size()==0 || t1o.size()==0){

 			cout<<"Voltage or Time vector empty. Skipping whole event"<<endl;
 			continue;
 		}

 		if(w1o.size()!= t1o.size()){

 			cout<<"Different number of entries in Voltage and Time vectors. Skipping whole event"<<endl;
 			continue;
 		}


		*a1=Analyzer(w1o,t1o);  //to be checked: is this a correct usage of memory?
		a1->Correct_Baseline(100);

		std::pair<double, unsigned int> tp_pair1 = a1->Find_Signal_Maximum(false,search_range); 
		std::pair<double, double> tp_pair1_fit = a1->Pmax_with_GausFit(tp_pair1);
		std::pair<double, unsigned int> neg_tp_pair1 = a1->Find_Negative_Signal_Maximum(false,search_range); 
		std::pair<double, double> neg_tp_pair1_fit = a1->Negative_Pmax_with_GausFit(neg_tp_pair1);      
		Pmax1 = tp_pair1.first*voltage_const; //mV
		negPmax1 = neg_tp_pair1.first*voltage_const; //mV
		Pmax1Fit = tp_pair1_fit.first*voltage_const; //mV
		negPmax1Fit = neg_tp_pair1_fit.first*voltage_const; //mV
		Tmax1 = a1->Get_Tmax(tp_pair1)*time_const; //ns
		negTmax1 = a1->Get_Negative_Tmax(neg_tp_pair1)*time_const; //ns
		Tmax1Fit = tp_pair1_fit.second*time_const; //ns
		negTmax1Fit = neg_tp_pair1_fit.second*time_const; //ns
		Area1 = a1->Find_Pulse_Area(tp_pair1)*voltage_const*time_const; //mV*ns
		UArea1 = a1->Find_Undershoot_Area(tp_pair1)*voltage_const*time_const; //mV*ns
		AreaFixed1 = a1->Pulse_Integration_with_Fixed_Window_Size(tp_pair1,"Simpson",1,1)*voltage_const*time_const; //mV*ns  ["Simpson" or "Rectangluar" are the 2 integration methods] 3rd and 4th arguments define the time window around Tmax (in ns!!!)
		AreaFixed1Fit = a1->Pulse_Integration_with_Fixed_Window_Size_with_GausFit(tp_pair1_fit,tp_pair1.second,"Simpson",1,1)*voltage_const*time_const; //mV*ns
		AreaInterpolated1 = a1->Pulse_Area_With_Linear_Interpolate_Edge(tp_pair1,"Simpson",true,1)*voltage_const*time_const; //mV*ns  4th argument is integration stop time. If 3rd arg is TRUE, stop time is calculated relative to the maximum, while it is absolute time in FALSE case.  
		AreaInterpolated1Fit = a1->Pulse_Area_With_Linear_Interpolate_Edge_with_GausFit(tp_pair1_fit,tp_pair1.second,"Simpson",true,1)*voltage_const*time_const;
		Area1_new = a1->New_Pulse_Area(tp_pair1_fit,tp_pair1.second,"Simpson",search_range, time_window[0], time_window[1])*voltage_const*time_const;//mV*ns  Similar to "Area1" method, but start/end times of the pulse are defined analitically. Should be the most accurate definition
		RiseTime1 = a1->Find_Rise_Time(tp_pair1,0.1,0.9)*time_const; //ns
		RiseTime1Fit = a1->Find_Rise_Time_with_GausFit(tp_pair1_fit, tp_pair1.second, 0.1, 0.9)*time_const; //ns
		dVdt1 = a1->Find_Dvdt(20,0,tp_pair1)*(voltage_const/time_const);  //mV/ns     [20 is the % at which the signal derivative is computed]
		dVdt1Fit = a1->Find_Dvdt_with_GausFit(20,0,tp_pair1_fit,tp_pair1.second)*(voltage_const/time_const);  //mV/ns
		cfd1 = a1->Rising_Edge_CFD_Time(20,tp_pair1)*time_const;  //ns   [standard CFD, 1st argument is the CFD %]
		cfd1Fit = a1->Rising_Edge_CFD_Time_with_GausFit(20,tp_pair1_fit,tp_pair1.second)*time_const;  //ns 
		t_thr1 = a1->Find_Time_At_Threshold(20,tp_pair1)*time_const; //ns    [time at which a thr (in mV !!!) is passed]
		rms1 = a1->Find_Noise(100)*voltage_const; //mV

		for(int j=0;j<7;j++){

			CFD1Fit[j]=a1->Rising_Edge_CFD_Time(10+j*10,tp_pair1)*time_const;

		}

		counter2++;

	}


	if( enable_channel_2 == 1){

		w2o.clear();
	 	t2o.clear();

 		if( invert_channel_2 == 1){

	 		for(unsigned int i=0; i<voltageReader2.GetSize();i++){

	 			w2o.push_back(-voltageReader2.At(i));
	 			t2o.push_back(timeReader2.At(i));

 			}

 		}else{

	 		for(unsigned int i=0; i<voltageReader2.GetSize();i++){

	 			w2o.push_back(voltageReader2.At(i));
	 			t2o.push_back(timeReader2.At(i));

 			}

 			}

 		if(w2o.size()<1000 || t2o.size()<1000){

 			cout<<"Voltage or Time vector less than 1000 entries. Skipping whole event"<<endl;
 			continue;
 		}

 		if(w2o.size()==0 || t2o.size()==0){

 			cout<<"Voltage or Time vector empty. Skipping whole event"<<endl;
 			continue;
 		}

 		if(w2o.size()!= t2o.size()){

 			cout<<"Different number of entries in Voltage and Time vectors. Skipping whole event"<<endl;
 			continue;
 		}

		*a2=Analyzer(w2o,t2o);  
		a2->Correct_Baseline(100);

		std::pair<double, unsigned int> tp_pair2 = a2->Find_Signal_Maximum(false,search_range); 
		std::pair<double, double> tp_pair2_fit = a2->Pmax_with_GausFit(tp_pair2);
		std::pair<double, unsigned int> neg_tp_pair2 = a2->Find_Negative_Signal_Maximum(false,search_range); 
		std::pair<double, double> neg_tp_pair2_fit = a2->Negative_Pmax_with_GausFit(neg_tp_pair2);      
		Pmax2 = tp_pair2.first*voltage_const; //mV
		negPmax2 = neg_tp_pair2.first*voltage_const; //mV
		Pmax2Fit = tp_pair2_fit.first*voltage_const; //mV
		negPmax2Fit = neg_tp_pair2_fit.first*voltage_const; //mV
		Tmax2 = a2->Get_Tmax(tp_pair2)*time_const; //ns
		negTmax2 = a2->Get_Negative_Tmax(neg_tp_pair2)*time_const; //ns
		Tmax2Fit = tp_pair2_fit.second*time_const; //ns
		negTmax2Fit = neg_tp_pair2_fit.second*time_const; //ns
		Area2 = a2->Find_Pulse_Area(tp_pair2)*voltage_const*time_const;
		UArea2 = a2->Find_Undershoot_Area(tp_pair2)*voltage_const*time_const;
		AreaFixed2 = a2->Pulse_Integration_with_Fixed_Window_Size(tp_pair2,"Simpson",1,1)*voltage_const*time_const;
		AreaFixed2Fit = a2->Pulse_Integration_with_Fixed_Window_Size_with_GausFit(tp_pair2_fit,tp_pair2.second,"Simpson",1,1)*voltage_const*time_const;
		AreaInterpolated2 = a2->Pulse_Area_With_Linear_Interpolate_Edge(tp_pair2,"Simpson",true,1)*voltage_const*time_const; 
		AreaInterpolated2Fit = a2->Pulse_Area_With_Linear_Interpolate_Edge_with_GausFit(tp_pair2_fit,tp_pair2.second,"Simpson",true,1)*voltage_const*time_const;
		Area2_new = a2->New_Pulse_Area(tp_pair2_fit,tp_pair2.second,"Simpson",search_range, time_window[0], time_window[1])*voltage_const*time_const;
		RiseTime2 = a2->Find_Rise_Time(tp_pair2,0.1,0.9)*time_const;
		RiseTime2Fit = a2->Find_Rise_Time_with_GausFit(tp_pair2_fit, tp_pair2.second, 0.1, 0.9)*time_const; //ns
		dVdt2 = a2->Find_Dvdt(20,0,tp_pair2)*(voltage_const/time_const);  
		dVdt2Fit = a2->Find_Dvdt_with_GausFit(20,0,tp_pair2_fit,tp_pair2.second)*(voltage_const/time_const);
		cfd2 = a2->Rising_Edge_CFD_Time(20,tp_pair2)*time_const;
		cfd2Fit = a2->Rising_Edge_CFD_Time_with_GausFit(20,tp_pair2_fit,tp_pair2.second)*time_const;  
		t_thr2 = a2->Find_Time_At_Threshold(0.02,tp_pair2)*time_const;
		rms2 = a2->Find_Noise(100)*voltage_const;

		for(int j=0;j<7;j++){

			CFD2Fit[j]=a2->Rising_Edge_CFD_Time(10+j*10,tp_pair2)*time_const;

		}

		

		counter2++;

	}


	if( enable_channel_3 == 1){

		w3o.clear();
	 	t3o.clear();

 		if( invert_channel_3 == 1){

	 		for(unsigned int i=0; i<voltageReader3.GetSize();i++){

	 			w3o.push_back(-voltageReader3.At(i));
	 			t3o.push_back(timeReader3.At(i));

 			}

 		}else{

	 		for(unsigned int i=0; i<voltageReader3.GetSize();i++){

	 			w3o.push_back(voltageReader3.At(i));
	 			t3o.push_back(timeReader3.At(i));

 			}

 			}

 		if(w3o.size()<1000 || t3o.size()<1000){

 			cout<<"Voltage or Time vector less than 1000 entries. Skipping whole event"<<endl;
 			continue;
 		}

 		if(w3o.size()==0 || t3o.size()==0){

 			cout<<"Voltage or Time vector empty. Skipping whole event"<<endl;
 			continue;
 		}

 		if(w3o.size()!= t3o.size()){

 			cout<<"Different number of entries in Voltage and Time vectors. Skipping whole event"<<endl;
 			continue;
 		}

		*a3 = Analyzer(w3o,t3o);  
		a3->Correct_Baseline(100);

		std::pair<double, unsigned int> tp_pair3 = a3->Find_Signal_Maximum(false,search_range); 
		std::pair<double, double> tp_pair3_fit = a3->Pmax_with_GausFit(tp_pair3);
		std::pair<double, unsigned int> neg_tp_pair3 = a3->Find_Negative_Signal_Maximum(false,search_range); 
		std::pair<double, double> neg_tp_pair3_fit = a3->Negative_Pmax_with_GausFit(neg_tp_pair3);      
		Pmax3 = tp_pair3.first*voltage_const; //mV
		negPmax3 = neg_tp_pair3.first*voltage_const; //mV
		Pmax3Fit = tp_pair3_fit.first*voltage_const; //mV
		negPmax3Fit = neg_tp_pair3_fit.first*voltage_const; //mV
		Tmax3 = a3->Get_Tmax(tp_pair3)*time_const; //ns
		negTmax3 = a3->Get_Negative_Tmax(neg_tp_pair3)*time_const; //ns
		Tmax3Fit = tp_pair3_fit.second*time_const; //ns
		negTmax3Fit = neg_tp_pair3_fit.second*time_const; //ns
		Area3 = a3->Find_Pulse_Area(tp_pair3)*voltage_const*time_const;
		UArea3 = a3->Find_Undershoot_Area(tp_pair3)*voltage_const*time_const;
		AreaFixed3 = a3->Pulse_Integration_with_Fixed_Window_Size(tp_pair3,"Simpson",1,1)*voltage_const*time_const;
		AreaFixed3Fit = a3->Pulse_Integration_with_Fixed_Window_Size_with_GausFit(tp_pair3_fit,tp_pair3.second,"Simpson",1,1)*voltage_const*time_const;
		AreaInterpolated3 = a3->Pulse_Area_With_Linear_Interpolate_Edge(tp_pair3,"Simpson",true,1)*voltage_const*time_const;   
		AreaInterpolated3Fit = a3->Pulse_Area_With_Linear_Interpolate_Edge_with_GausFit(tp_pair3_fit,tp_pair3.second,"Simpson",true,1)*voltage_const*time_const;
		Area3_new = a3->New_Pulse_Area(tp_pair3_fit,tp_pair3.second,"Simpson",search_range, time_window[0], time_window[1])*voltage_const*time_const;
		RiseTime3 = a3->Find_Rise_Time(tp_pair3,0.1,0.9)*time_const;
		RiseTime3Fit = a3->Find_Rise_Time_with_GausFit(tp_pair3_fit, tp_pair3.second, 0.1, 0.9)*time_const; //ns
		dVdt3 = a3->Find_Dvdt(20,0,tp_pair3)*(voltage_const/time_const);  
		dVdt3Fit = a3->Find_Dvdt_with_GausFit(20,0,tp_pair3_fit,tp_pair3.second)*(voltage_const/time_const);
		cfd3 = a3->Rising_Edge_CFD_Time(20,tp_pair3)*time_const; 
		cfd3Fit = a3->Rising_Edge_CFD_Time_with_GausFit(20,tp_pair3_fit,tp_pair3.second)*time_const; 
		t_thr3 = a3->Find_Time_At_Threshold(0.02,tp_pair3)*time_const;
		rms3 = a3->Find_Noise(100)*voltage_const;

		for(int j=0;j<7;j++){

			CFD3Fit[j]=a3->Rising_Edge_CFD_Time(10+j*10,tp_pair3)*time_const;

		}

		counter2++;

	}


	if( enable_channel_4 == 1){

		w4o.clear();
	 	t4o.clear();

 		if( invert_channel_4 == 1){

	 		for(unsigned int i=0; i<voltageReader4.GetSize();i++){

	 			w4o.push_back(-voltageReader4.At(i));
	 			t4o.push_back(timeReader4.At(i));

 			}

 		}else{

	 		for(unsigned int i=0; i<voltageReader4.GetSize();i++){

	 			w4o.push_back(voltageReader4.At(i));
	 			t4o.push_back(timeReader4.At(i));

 			}

 			}

 		if(w4o.size()<1000 || t4o.size()<1000){

 			cout<<"Voltage or Time vector less than 1000 entries. Skipping whole event"<<endl;
 			continue;
 		}

 		if(w4o.size()==0 || t4o.size()==0){

 			cout<<"Voltage or Time vector empty. Skipping whole event"<<endl;
 			continue;
 		}

 		if(w4o.size()!= t4o.size()){

 			cout<<"Different number of entries in Voltage and Time vectors. Skipping whole event"<<endl;
 			continue;
 		}

		*a4 = Analyzer(w4o,t4o);  
		a4->Correct_Baseline(100);

		std::pair<double, unsigned int> tp_pair4 = a4->Find_Signal_Maximum(false,search_range); 
		std::pair<double, double> tp_pair4_fit = a4->Pmax_with_GausFit(tp_pair4);
		std::pair<double, unsigned int> neg_tp_pair4 = a4->Find_Negative_Signal_Maximum(false,search_range); 
		std::pair<double, double> neg_tp_pair4_fit = a4->Negative_Pmax_with_GausFit(neg_tp_pair4);      
		Pmax4 = tp_pair4.first*voltage_const; //mV
		negPmax4 = neg_tp_pair4.first*voltage_const; //mV
		Pmax4Fit = tp_pair4_fit.first*voltage_const; //mV
		negPmax4Fit = neg_tp_pair4_fit.first*voltage_const; //mV
		Tmax4 = a4->Get_Tmax(tp_pair4)*time_const; //ns
		negTmax4 = a4->Get_Negative_Tmax(neg_tp_pair4)*time_const; //ns
		Tmax4Fit = tp_pair4_fit.second*time_const; //ns
		negTmax4Fit = neg_tp_pair4_fit.second*time_const; //ns
		Area4 = a4->Find_Pulse_Area(tp_pair4)*voltage_const*time_const;
		UArea4 = a4->Find_Undershoot_Area(tp_pair4)*voltage_const*time_const;
		AreaFixed4 = a4->Pulse_Integration_with_Fixed_Window_Size(tp_pair4,"Simpson",1,1)*voltage_const*time_const;
		AreaFixed4Fit = a4->Pulse_Integration_with_Fixed_Window_Size_with_GausFit(tp_pair4_fit,tp_pair4.second,"Simpson",1,1)*voltage_const*time_const;
		AreaInterpolated4 = a4->Pulse_Area_With_Linear_Interpolate_Edge(tp_pair4,"Simpson",true,1)*voltage_const*time_const; 
		AreaInterpolated4Fit = a4->Pulse_Area_With_Linear_Interpolate_Edge_with_GausFit(tp_pair4_fit,tp_pair4.second,"Simpson",true,1)*voltage_const*time_const;
		Area4_new = a4->New_Pulse_Area(tp_pair4_fit,tp_pair4.second,"Simpson",search_range, time_window[0], time_window[1])*voltage_const*time_const;
		RiseTime4 = a4->Find_Rise_Time(tp_pair4,0.1,0.9)*time_const;
		RiseTime4Fit = a4->Find_Rise_Time_with_GausFit(tp_pair4_fit, tp_pair4.second, 0.1, 0.9)*time_const; //ns
		dVdt4 = a4->Find_Dvdt(20,0,tp_pair4)*(voltage_const/time_const); 
		dVdt4Fit = a4->Find_Dvdt_with_GausFit(20,0,tp_pair4_fit,tp_pair4.second)*(voltage_const/time_const); 
		cfd4 = a4->Rising_Edge_CFD_Time(20,tp_pair4)*time_const; 
		cfd4Fit = a4->Rising_Edge_CFD_Time_with_GausFit(20,tp_pair4_fit,tp_pair4.second)*time_const; 
		t_thr4 = a4->Find_Time_At_Threshold(0.02,tp_pair4)*time_const;
		rms4 = a4->Find_Noise(100)*voltage_const;

		for(int j=0;j<7;j++){

			CFD4Fit[j]=a4->Rising_Edge_CFD_Time(10+j*10,tp_pair4)*time_const;

		}

	}

	ntrig=j;
	event=j;
		
	OutTree->Fill();

	if(j%10000 == 0) cout<<"processed events:"<<j<<endl;
	//cout<<"processed events:"<<j<<endl;

	j++;
		

	}

	
 OutTree->Write();
 OutputFile->Write();

}