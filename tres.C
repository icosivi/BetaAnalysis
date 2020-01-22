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


void tres(char *filename, int trigger_ch, int dut_ch, int cfd){


if(cfd!=10 && cfd!=20 && cfd!=30 && cfd!=40 && cfd!=50 && cfd!=60 && cfd!=70 ){

	cout<<"Available CFDs:10,20,30,40,50,60,70. You inserted a different value, please try again "<<endl;

}else{

	int index=cfd/10 -1;

	TFile *file = TFile::Open(filename);
	TTree *itree = dynamic_cast<TTree*>(file->Get("Analysis"));
	TTreeReader myReader("Analysis", file);

	char *trigger_branch = Form("CFD%iFit",trigger_ch);
	TTreeReaderArray<double> trigger(myReader, trigger_branch); 

	char *dut_branch = Form("CFD%iFit",dut_ch);
	TTreeReaderArray<double> dut(myReader, dut_branch);

	TH1D *hist=new TH1D("hist","hist",2000,-10,10);


	while(myReader.Next()){

		hist->Fill(dut.At(index)-trigger.At(index));

	}

	hist->Draw();

 }

}

// absolute positions and propagations peed need to be added
void tres_rsd(char *filename, int trigger_ch, int cfd){

 if(cfd!=10 && cfd!=20 && cfd!=30 && cfd!=40 && cfd!=50 && cfd!=60 && cfd!=70 ){

	cout<<"Available CFDs:10,20,30,40,50,60,70. You inserted a different value, please try again "<<endl;

}else{

	 int index=cfd/10 -1;

	 TFile *file = TFile::Open(filename);
	 TTree *itree = dynamic_cast<TTree*>(file->Get("Analysis"));
	 TTreeReader myReader("Analysis", file);

	 char *trigger_branch = Form("CFD%iFit",trigger_ch);
	 TTreeReaderArray<double> trigger(myReader, trigger_branch); 

	 TTreeReaderArray<double> dut1(myReader,"CFD1Fit");
	 TTreeReaderValue<double> pmax1(myReader,"Pmax1");

	 TTreeReaderArray<double> dut2(myReader,"CFD2Fit");
	 TTreeReaderValue<double> pmax2(myReader,"Pmax2");

	 TTreeReaderArray<double> dut3(myReader,"CFD3Fit");
	 TTreeReaderValue<double> pmax3(myReader,"Pmax3");

	 TTreeReaderArray<double> dut4(myReader,"CFD4Fit");
	 TTreeReaderValue<double> pmax4(myReader,"Pmax4");

	 TH1D *hist=new TH1D("hist","hist",2000,-10,10);

	 // Absolute position of the pads. Used to calculate time resolution. In um !!!!!!
	 double x_pos1;
	 double x_pos2;
	 double x_pos3;
	 double x_pos4;

	 // propagation velocity.
	 const double speed;  //um/ns
	 const double slope = 0.00074;  // ns/um   DA ESTRARRE PARTENDO DA t vs distance dove t è il tempo non corretto

	 //hit position reconstructed with centroid method
	 double x_centroid;

	 while(myReader.Next()){

	 	if(trigger_ch!=1){

	 		x_centroid += *pmax1*x_pos1;
	 		amp_sum += *pmax1;

	 	}

	 	if(trigger_ch!=2){

	 		x_centroid += *pmax2*x_pos2;
	 		amp_sum += *pmax2;

	 	}

	 	if(trigger_ch!=3){

	 		x_centroid += *pmax3*x_pos3;
	 		amp_sum += *pmax3;

	 	}

	 	if(trigger_ch!=4){

	 		x_centroid += *pmax4*x_pos4;
	 		amp_sum += *pmax4;

	 	}

	 }

	 
	 while(myReader.Next()){

	 	double rsd_time=0;
	 	double amp_sum=0;

	 	if(trigger_ch!=1){

	 		double x_i1 = x_centroid-x_pos1;
	 		double t_offset1 = slope*x_i1;
	 		rsd_time += *pmax1*(dut1.At(index)-x_i1/speed-t_offset1);
	 		amp_sum += *pmax1;

	 	}

	 	if(trigger_ch!=2){

	 		double x_i2 = x_centroid-x_pos2;
	 		double t_offset2 = slope*x_i2;
	 		rsd_time += *pmax2*(dut2.At(index)-x_i2/speed-t_offset2);
	 		amp_sum += *pmax2;

	 	}

	 	if(trigger_ch!=3){

	 		double x_i3 = x_centroid-x_pos3;
	 		double t_offset3 = slope*x_i3;
	 		rsd_time += *pmax3*(dut3.At(index)-x_i3/speed-t_offset3);
	 		amp_sum += *pmax3;

	 	}

	 	if(trigger_ch!=4){

	 		double x_i4 = x_centroid-x_pos4;
	 		double t_offset4 = slope*x_i4;
	 		rsd_time += *pmax4*(dut4.At(index)-x_i4/speed-t_offset4);
	 		amp_sum += *pmax4;

	 	}

	 	hist->Fill(rsd_time/amp_sum-trigger.At(index));

	 }

	 hist->Draw();

 }

}