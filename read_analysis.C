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

  int ch = channel - 1 ; //channel starts from 1 but arrays index from 0 
 
  ConfigFile cf("beta_config.ini");

  std::string Filename = cf.Value("HEADER","output_filename");
  const char *filename = Filename.c_str();
  TFile *file = TFile::Open(filename);
  TTree *itree = dynamic_cast<TTree*>(file->Get("Analysis"));
  TTreeReader myReader("Analysis", file);

  
  int trigger_channel = cf.Value("HEADER","trigger_ch");
  int trigger_ch = trigger_channel - 1 ;

  const double conversion_factor = cf.Value("HEADER","conversion_factor_area_charge") ;
  const double time_const = cf.Value("HEADER","time_scalar");  
  const double voltage_const = cf.Value("HEADER","voltage_scalar");

  int save_plots = cf.Value("ANALYSIS","save_plots") ;
  int best_cfd_index = cf.Value("ANALYSIS","cfd_index") ; 

  double negpmax = cf.Value("ANALYSIS","negpmax") ; 
  double pmax_low = 10;
  double pmax_up = cf.Value("ANALYSIS","pmax_up") ; 
  double area_low = 10;
  double area_up = cf.Value("ANALYSIS","area_up") ;



  TH1F *pmax_hist = new TH1F("pmax_hist",Form("Pmax%i", channel ), pmax_up, 0, pmax_up);
  pmax_hist->GetXaxis()->SetTitle("amplitude  [mV]");
 
  TH1F *pmax_hist_wCuts = new TH1F("pmax_hist_wCuts",Form("Pmax%i", channel ), pmax_up, 0, pmax_up);
  pmax_hist_wCuts->GetXaxis()->SetTitle("amplitude  [mV]");


  TH1F *area_hist = new TH1F("area_hist",Form("Area%i", channel ), area_up, 0, area_up);
  area_hist->GetXaxis()->SetTitle("area  [pWb]");
 
  TH1F *area_hist_wCuts = new TH1F("area_hist_wCuts",Form("Area%i", channel ), area_up, 0, area_up);
  area_hist_wCuts->GetXaxis()->SetTitle("area  [pWb]");


  TH1F *dvdt_hist = new TH1F("dvdt_hist",Form("dVdt%i", channel ),100,0,400);
  dvdt_hist->GetXaxis()->SetTitle("dVdt  [mV/ns]");


  TH1F *risetime_hist = new TH1F("risetime_hist",Form("Risetime%i", channel ),200,0,2);
  risetime_hist->GetXaxis()->SetTitle("Risetime  [ns]");


  TH1F *rms_hist = new TH1F("rms_hist",Form("RMS %i", channel ),50,0,5);
  rms_hist->GetXaxis()->SetTitle("RMS  [mV]");


  TH1F *tot_hist = new TH1F("tot_hist",Form("ToT%i", channel ),300,0,3);
  tot_hist->GetXaxis()->SetTitle("ToT  [ns]");


  TH1F *tmax_hist = new TH1F("tmax_hist",Form("Tmax%i", channel ),600,-30,30);
  tmax_hist->GetXaxis()->SetTitle("Tmax  [ns]");
  itree->Project("tmax_hist", Form("tmax[%i]",ch), Form("pmax[%i]>20", ch) );
  double tmax_search_range_mean = tmax_hist->GetBinCenter( tmax_hist->GetMaximumBin() );
  double tmax_range[2] = { tmax_search_range_mean-0.5, tmax_search_range_mean+0.5 };
  //double tmax_range[2] = { -10, -8 };
 

  double timeres_search_range[2] = { cf.Value("ANALYSIS", "timeres_search_min" )*time_const, 
                                     cf.Value("ANALYSIS", "timeres_search_max")*time_const };
  
  
  int timeres_bin = ((timeres_search_range[1] - timeres_search_range[0])/2.)*100 ;
  TH1F *timeres = new TH1F("timeres",Form("ToA_{Ch%i} - ToA_{Trigger}", channel ), timeres_bin, timeres_search_range[0], timeres_search_range[1]);
  timeres->GetXaxis()->SetTitle("[ns]");
  


  double trigger_pmax_low = cf.Value("ANALYSIS","trigger_pmax_low") ;
  double trigger_pmax_up = cf.Value("ANALYSIS","trigger_pmax_up") ;
  double trigger_area_low = cf.Value("ANALYSIS","trigger_area_low") ;
  double trigger_area_up = cf.Value("ANALYSIS","trigger_area_up") ;


  
  itree->Project("pmax_hist", Form("pmax[%i]",ch), Form("negpmax[%i]>%f && tmax[%i]>%f && tmax[%i]<%f", ch, negpmax, ch, tmax_range[0], ch, tmax_range[1]) );
  TF1 *pmax_landau_noCuts=new TF1( "pmax_landau_noCuts", "landau", pmax_low, pmax_up );
  pmax_hist->Fit("pmax_landau_noCuts","RQ");
  pmax_low = pmax_landau_noCuts->GetParameter(1)-pmax_landau_noCuts->GetParameter(2);
  if(pmax_low<10) pmax_low = 10;
  cout<<" "<<endl;
  cout<<" "<<endl;
  cout<< Form( "////////// Beginning Channel %i ///////////", channel) <<endl;
  cout<<"Pmax low cut: "<<pmax_low<<endl;

  
  itree->Project("area_hist", Form("area[%i]",ch), Form("negpmax[%i]>%f && tmax[%i]>%f && tmax[%i]<%f", ch, negpmax, ch, tmax_range[0], ch, tmax_range[1]) );
  TF1 *area_landau_noCuts=new TF1( "area_landau_noCuts", "landau", area_low, area_up );
  area_hist->Fit("area_landau_noCuts","RQ");
  area_low = area_landau_noCuts->GetParameter(1) - area_landau_noCuts->GetParameter(2);
  if(area_low<10) area_low = 10;
  cout<<"Area low cut: "<<area_low<<endl;
  cout<<" "<<endl;
  cout<<"//////////////////////////////////////////////// "<<endl;
  cout<<" "<<endl;

  
  itree->Project("pmax_hist_wCuts", Form("pmax[%i]",ch), Form("negpmax[%i]>%f && pmax[%i]>%f && pmax[%i]<%f", ch, negpmax, ch, pmax_low, ch, pmax_up) );
  TF1 *pmax_landau=new TF1( "pmax_landau", "landau", 0, pmax_up );
  pmax_hist_wCuts->Fit("pmax_landau","RQ");
  double PMAX = pmax_landau->GetParameter(1);
  cout<<Form("Pmax MPV ch %i: ", channel)<<PMAX<<endl;


  itree->Project("area_hist_wCuts", Form("area[%i]",ch), Form("negpmax[%i]>%f && area[%i]>%f && area[%i]<%f", ch, negpmax, ch, area_low, ch, area_up) );
  TF1 *area_landau=new TF1( "area_landau", "landau", 0, area_up );
  area_hist_wCuts->Fit("area_landau","RQ");
  double AREA = area_landau->GetParameter(1);
  double CHARGE = AREA/conversion_factor; 
  cout<<Form("Charge MPV ch %i: ", channel)<<CHARGE<<endl;


  itree->Project("rms_hist", Form("rms[%i]",ch), Form("negpmax[%i]>%f && pmax[%i]>%f && pmax[%i]<%f && area[%i]>%f", 
                  ch, negpmax, ch, pmax_low, ch, pmax_up, ch, area_low) );
  TF1 *rms_gaus=new TF1( "rms_gaus", "gaus", 0, 10 );
  rms_hist->Fit("rms_gaus","RQ");
  double RMS = rms_gaus->GetParameter(1);
  cout<<Form("RMS noise ch %i: ", channel)<<RMS<<endl;


  itree->Project("dvdt_hist", Form("dvdt[%i]",ch), Form("negpmax[%i]>%f && pmax[%i]>%f && pmax[%i]<%f && area[%i]>%f && area[%i]<%f", 
                  ch, negpmax, ch, pmax_low, ch, pmax_up, ch, area_low, ch, area_up) );
  TF1 *dvdt_landau=new TF1( "dvdt_landau", "landau", 0, 400 );
  dvdt_hist->Fit("dvdt_landau","RQ");
  double dVdt = dvdt_landau->GetParameter(1);
  cout<<Form("dVdt ch %i: ", channel)<<dVdt<<endl;


  itree->Project("risetime_hist", Form("risetime[%i]",ch), Form("negpmax[%i]>%f && pmax[%i]>%f && pmax[%i]<%f && area[%i]>%f && area[%i]<%f", 
                  ch, negpmax, ch, pmax_low, ch, pmax_up, ch, area_low, ch, area_up) );
  TF1 *risetime_landau=new TF1( "risetime_landau", "landau", 0, 2);
  risetime_hist->Fit("risetime_landau","RQ");
  double RISETIME = risetime_landau->GetParameter(1);
  cout<<Form("Risetime ch %i: ", channel)<<RISETIME<<endl;



  float TimeRes[7] ;
  float x_cfd[7] = {10, 20, 30, 40, 50, 60, 70};
  TF1 *timeres_gaus=new TF1( "timeres_gaus", "gaus", timeres_search_range[0], timeres_search_range[1] );

  for(int cfd_index=0; cfd_index<7; cfd_index++){
  
    itree->Project("timeres", Form("cfd[%i][%i]-cfd[%i][1]", ch, cfd_index, trigger_ch ), 
          Form("negpmax[%i]>%f && pmax[%i]>%f && pmax[%i]<%f && area[%i]>%f && area[%i]<%f && negpmax[%i]>-40 && pmax[%i]>%f && pmax[%i]<%f && area[%i]>%f && area[%i]<%f", 
          ch, negpmax, ch, pmax_low, ch, pmax_up, ch, area_low, ch, area_up,
          trigger_ch, trigger_ch, trigger_pmax_low, trigger_ch, trigger_pmax_up, trigger_ch, trigger_area_low, trigger_ch, trigger_area_up) );
  
    timeres->Fit("timeres_gaus","RQ");
    TimeRes[cfd_index] = timeres_gaus->GetParameter(2);

  }

  double RESOLUTION = TimeRes[best_cfd_index] ;
  cout<<Form("Time resolution ch %i: ", channel)<<RESOLUTION<<endl;


  TGraph *gr_time_res = new TGraph(7,x_cfd,TimeRes);
  gr_time_res->SetTitle( Form("Time resolution vs CFD - Channel%i", ch) );
  gr_time_res->GetXaxis()->SetTitle("CFD  [%]");
  gr_time_res->GetYaxis()->SetTitle("Time resolution  [ps]");

  
  if( save_plots==1 ){

    std::string s = filename;
    std::string delimiter = ".root";
    std::string token = s.substr(0, s.find(delimiter)); 

    TString string_out = token;
    string_out += Form("_PLOTS_Channel%i.root", channel ) ;
    TFile f(string_out,"recreate");

    f.cd();
    pmax_hist_wCuts->Write();
    area_hist_wCuts->Write();
    dvdt_hist->Write();
    rms_hist->Write();
    risetime_hist->Write();
    timeres->Write();
    gr_time_res->Write();

    f.Write();
    f.Close();

  }


}


/*int main(){

readAnalysis();

return 0;

}*/
