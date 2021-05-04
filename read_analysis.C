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
#include <TROOT.h>
#include <TMath.h>
#include <TTree.h>
#include <TTreeReader.h>
#include "TTreeReaderValue.h"
#include "TTreeReaderArray.h"
#include <TBranch.h>
#include <TFile.h>
#include <TH1.h>
#include <TF1.h>
#include <TGraph.h>
#include <TString.h>
#include <TThread.h>
#include <TROOT.h>
#include <TStyle.h>
#include <TImage.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TF1Convolution.h>
#include <TProfile.h>

//------Custom headers----------------//
#include "include/general.hpp"
#include "include/Chameleon.h"
#include "include/ConfigFile.hpp"


void readAnalysis( int channel ){

  gStyle->SetTitleOffset( 1.3, "x" );
  gStyle->SetTitleOffset( 1.3, "y" );
  gStyle->SetHistLineWidth(2);
  gStyle->SetOptStat(1111);
  //gStyle->SetPadTickX(1);  // To get tick marks on the opposite side of the frame
  //gStyle->SetPadTickY(1);
 
  ConfigFile cf("beta_config.ini");

  std::string Filename = cf.Value("HEADER","output_filename");
  const char *filename = Filename.c_str();
  TFile *file = TFile::Open(filename);
  TTree *itree = dynamic_cast<TTree*>(file->Get("Analysis"));
  TTreeReader myReader("Analysis", file);

  int ch = channel - 1 ; //channel starts from 1 but array indeces from 0 

  bool save_plots = false ;
  int best_cfd_index = 1 ;
  bool FAST = false ;
  int ntrig,event;
  int j = 0;

  TH1F *pmax_hist = new TH1F("pmax_hist",Form("Pmax%i", channel ),350,0,350);
  pmax_hist->GetXaxis()->SetTitle("amplitude  [mV]");
 
  TH1F *pmax_hist_wCuts = new TH1F("pmax_hist_wCuts",Form("Pmax%i", channel ),350,0,350);
  pmax_hist_wCuts->GetXaxis()->SetTitle("amplitude  [mV]");


  TH1F *area_hist = new TH1F("area_hist",Form("Area%i", channel ),200,0,200);
  area_hist->GetXaxis()->SetTitle("area  [pWb]");
 
  TH1F *area_hist_wCuts = new TH1F("area_hist_wCuts",Form("Area%i", channel ),200,0,200);
  area_hist_wCuts->GetXaxis()->SetTitle("area  [pWb]");


  TH1F *dvdt_hist = new TH1F("dvdt_hist",Form("dVdt%i", channel ),100,0,400);
  dvdt_hist->GetXaxis()->SetTitle("dVdt  [mV/ns]");


  TH1F *risetime_hist = new TH1F("risetime_hist",Form("Risetime%i", channel ),200,0,2);
  risetime_hist->GetXaxis()->SetTitle("Risetime  [ns]");


  TH1F *rms_hist = new TH1F("rms_hist",Form("RMS %i", channel ),50,0,5);
  rms_hist->GetXaxis()->SetTitle("RMS  [mV]");


  TH1F *tot_hist = new TH1F("tot_hist",Form("ToT%i", channel ),300,0,3);
  tot_hist->GetXaxis()->SetTitle("ToT  [ns]");


  TH1F *timeres = new TH1F("timeres",Form("ToA_{Ch%i} - ToA_{Trigger}", channel ),100, -1, 1);
  timeres->GetXaxis()->SetTitle("[ns]");


  
  TH1F *timeres_corrected = new TH1F("timeres_corrected",Form("ToA_{Ch%i}^{'} - ToA_{Trigger}", channel),100, -1, 1);
  timeres_corrected->GetXaxis()->SetTitle("[ns]");
  



  double trigger_pmax_low = 50. ;
  double trigger_pmax_up = 250. ;
  double trigger_area_low = 55.;
  double trigger_area_up = 250. ;

  double negpmax = -50;
  double pmax_low ;
  double pmax_up = 300 ; 
  double tmax_range[2] = {-2,2};
  double area_low ;
  double area_up = 300 ;


  itree->Project("pmax_hist", Form("pmax[%i]",ch), Form("negpmax[%i]>%f && tmax[%i]>%f && tmax[%i]<%f", ch, negpmax, ch, tmax_range[0], ch, tmax_range[1]) );
  TF1 *pmax_landau_noCuts=new TF1( "pmax_landau_noCuts", "landau", 15, 300 );
  pmax_hist->Fit("pmax_landau_noCuts","R");
  pmax_low = pmax_landau_noCuts->GetParameter(1)-pmax_landau_noCuts->GetParameter(2);
  if(pmax_low<10) pmax_low = 10;
  cout<<"Pmax low cut: "<<pmax_low<<endl;

  
  itree->Project("area_hist", Form("area[%i]",ch), Form("negpmax[%i]>%f && tmax[%i]>%f && tmax[%i]<%f", ch, negpmax, ch, tmax_range[0], ch, tmax_range[1]) );
  TF1 *area_landau_noCuts=new TF1( "area_landau_noCuts", "landau", 5, 200 );
  area_hist->Fit("area_landau_noCuts","R");
  area_low = area_landau_noCuts->GetParameter(1) - area_landau_noCuts->GetParameter(2);
  if(area_low<10) area_low = 10;
  cout<<"Area low cut: "<<area_low<<endl;

  
  itree->Project("pmax_hist_wCuts", Form("pmax[%i]",ch), Form("negpmax[%i]>%f && pmax[%i]>%f && pmax[%i]<300", ch, negpmax, ch, pmax_low, ch) );
  TF1 *pmax_landau=new TF1( "pmax_landau", "landau", 15, 300 );
  pmax_hist_wCuts->Fit("pmax_landau","R");
  double PMAX = pmax_landau->GetParameter(1);


  itree->Project("area_hist_wCuts", Form("area[%i]",ch), Form("negpmax[%i]>%f && area[%i]>%f", ch, negpmax, ch, area_low) );
  TF1 *area_landau=new TF1( "area_landau", "landau", 5, 300 );
  area_hist_wCuts->Fit("area_landau","R");
  double AREA = area_landau->GetParameter(1);


  itree->Project("rms_hist", Form("rms[%i]",ch), Form("negpmax[%i]>%f && pmax[%i]>%f && pmax[%i]<300 && area[%i]>%f", 
                  ch, negpmax, ch, pmax_low, ch, ch, area_low) );
  TF1 *rms_gaus=new TF1( "rms_gaus", "gaus", 0, 10 );
  rms_hist->Fit("rms_gaus","RQ");
  double RMS = rms_gaus->GetParameter(1);


  itree->Project("dvdt_hist", Form("dvdt[%i]",ch), Form("negpmax[%i]>%f && pmax[%i]>%f && pmax[%i]<300 && area[%i]>%f", 
                  ch, negpmax, ch, pmax_low, ch, ch, area_low) );
  TF1 *dvdt_landau=new TF1( "dvdt_landau", "landau", 0, 400 );
  dvdt_hist->Fit("dvdt_landau","RQ");
  double dVdt = dvdt_landau->GetParameter(1);


  itree->Project("risetime_hist", Form("risetime[%i]",ch), Form("negpmax[%i]>%f && pmax[%i]>%f && pmax[%i]<300 && area[%i]>%f", 
                  ch, negpmax, ch, pmax_low, ch, ch, area_low) );
  TF1 *risetime_landau=new TF1( "risetime_landau", "landau", 0, 2);
  risetime_hist->Fit("risetime_landau","RQ");
  double RISETIME = risetime_landau->GetParameter(1);



  float TimeRes[7] ;
  float x_cfd[7] = {10, 20, 30, 40, 50, 60, 70};

  TF1 *timeres_gaus=new TF1( "timeres_gaus", "gaus", -1, 1 );

  for(int cfd_index=0; cfd_index<7; cfd_index++){
  
    itree->Project("timeres", Form("cfd[%i][%i]-cfd[1][1]", ch, cfd_index), 
              Form("negpmax[%i]>%f && pmax[%i]>%f && pmax[%i]<300 && area[%i]>%f && negpmax[1]>-40 && pmax[1]>%f && pmax[1]<%f && area[1]>%f && area[1]<%f", 
                          ch, negpmax, ch, pmax_low, ch, ch, area_low, trigger_pmax_low, trigger_pmax_up, trigger_area_low, trigger_area_up) );
  
    timeres->Fit("timeres_gaus","RQ");
    TimeRes[cfd_index] = timeres_gaus->GetParameter(2)*1000;

  }

  double RESOLUTION = TimeRes[best_cfd_index] ;


  TGraph *gr_time_res = new TGraph(7,x_cfd,TimeRes);
  gr_time_res->SetTitle( Form("Time resolution vs CFD - Channel%i", ch) );
  gr_time_res->GetXaxis()->SetTitle("CFD  [%]");
  gr_time_res->GetYaxis()->SetTitle("Time resolution  [ps]");

  
  if(FAST){

    TProfile *tot_prof1 = new TProfile("tot_prof1", Form("ToA vs ToT profile %i",ch), 600, 2, 8);
    itree->Project("tot_prof1", Form("t_thr[%i]:tot[%i]", ch, ch), 
                    Form("negpmax[%i]>%f && pmax[%i]>%f && area[%i]>%f && negpmax[1]>-40 && pmax[1]>%f && pmax[1]<%f && area[1]>%f && area[1]<%f", 
                          ch, negpmax, ch, pmax_low, ch, area_low, trigger_pmax_low, trigger_pmax_up, trigger_area_low, trigger_area_up) );

    TF1 *tot_prof_fit1=new TF1( "tot_prof_fit1", "pol2", 2, 8 );
    tot_prof1->Fit("tot_prof_fit1","RQ");
    //tot_prof1->Draw();

    TProfile *tot_prof2 = new TProfile("tot_prof2","ToA vs ToT profile 1", 600, 2, 8);
    itree->Project("tot_prof2", Form("(t_thr[%i]-%f-%f*tot[%i]-%f*tot[%i]*tot[%i]):tot[%i]", 
                    ch, tot_prof_fit1->GetParameter(0), tot_prof_fit1->GetParameter(1), ch, tot_prof_fit1->GetParameter(2), ch, ch, ch), 
                    Form("negpmax[%i]>%f && pmax[%i]>%f && area[%i]>%f && negpmax[1]>-40 && pmax[1]>%f && pmax[1]<%f && area[1]>%f && area[1]<%f", 
                         ch, negpmax, ch, pmax_low, ch, area_low, trigger_pmax_low, trigger_pmax_up, trigger_area_low, trigger_area_up) );

    TF1 *tot_prof_fit1_lin=new TF1( "tot_prof_fit1_lin", "pol1", 2, 8 );
    TF1 *tot_prof_fit1_const=new TF1( "tot_prof_fit1_const", "pol0", 2, 8 );
    tot_prof2->Fit("tot_prof_fit1_lin","RQ");
    tot_prof2->Fit("tot_prof_fit1_const","RQ");
    //tot_prof2->Draw();

    TF1 *timeres_corrected_gaus=new TF1( "timeres_corrected_gaus", "gaus", -1, 1 );

    if( tot_prof_fit1_const->GetChisquare() >= tot_prof_fit1_lin->GetChisquare()){

      itree->Project("timeres_corrected", Form("(t_thr[%i]-%f-%f*tot[%i]-%f*tot[%i]*tot[%i]-%f-%f*tot[%i]) - cfd[1][1]", 
                     ch, tot_prof_fit1->GetParameter(0), tot_prof_fit1->GetParameter(1), ch, tot_prof_fit1->GetParameter(2), 
                     ch, ch, tot_prof_fit1_lin->GetParameter(0), tot_prof_fit1_lin->GetParameter(1), ch ), 
                     Form("negpmax[%i]>%f && pmax[%i]>%f && area[%i]>%f && negpmax[1]>-40 && pmax[1]>%f && pmax[1]<%f && area[1]>%f && area[1]<%f", 
                          ch, negpmax, ch, pmax_low, ch, area_low, trigger_pmax_low, trigger_pmax_up, trigger_area_low, trigger_area_up) );

      timeres_corrected->Fit("timeres_corrected_gaus","RQ");

    }else{

      itree->Project("timeres_corrected", Form("(t_thr[%i]-%f-%f*tot[%i]-%f*tot[%i]*tot[%i]) - cfd[1][1]", 
                      ch, tot_prof_fit1->GetParameter(0), tot_prof_fit1->GetParameter(1), ch, tot_prof_fit1->GetParameter(2), ch, ch), 
                      Form("negpmax[%i]>%f && pmax[%i]>%f && area[%i]>%f && negpmax[1]>-40 && pmax[1]>%f && pmax[1]<%f && area[1]>%f && area[1]<%f", 
                            ch, negpmax, ch, pmax_low, ch, area_low, trigger_pmax_low, trigger_pmax_up, trigger_area_low, trigger_area_up) );

      timeres_corrected->Fit("timeres_corrected_gaus","RQ");
    
    }


    double TimeRes_corr = timeres_corrected_gaus->GetParameter(2)*1000;
    //timeres_corrected->Draw();
    cout<<"Time res corrected: "<<TimeRes_corr<<endl;

  }


  

  if(save_plots){

    std::string s = filename;
    std::string delimiter = ".root";
    std::string token = s.substr(0, s.find(delimiter)); 

    TString string_out = token;
    string_out += Form("_PLOTS_Channel%i.root", ch) ;
    TFile f(string_out,"recreate");

    f.cd();
    pmax_hist_wCuts->Write();
    area_hist_wCuts->Write();
    dvdt_hist->Write();
    rms_hist->Write();
    risetime_hist->Write();
    timeres->Write();
    gr_time_res->Write();
    if(FAST) timeres_corrected->Write();

    f.Write();
    f.Close();

  }


}


/*int main(){

readAnalysis();

return 0;

}*/
