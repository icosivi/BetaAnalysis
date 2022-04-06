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
#include <TROOT.h>
#include <TStyle.h>
#include <TImage.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TProfile.h>

//------Custom headers----------------//
//#include "include/general.hpp"
//#include "include/Chameleon.h"
#include "include/ConfigFile.hpp"


void readAnalysis( std::string Filename, int channel ){

  gStyle->SetTitleOffset( 1.3, "x" );
  gStyle->SetTitleOffset( 1.3, "y" );
  gStyle->SetHistLineWidth(2);
  gStyle->SetOptStat(1111);
  //gStyle->SetPadTickX(1);  // To get tick marks on the opposite side of the frame
  //gStyle->SetPadTickY(1);
 
  ConfigFile cf("beta_config.ini");

  //std::string Filename = cf.Value("HEADER","output_filename");
  const char *filename = Filename.c_str();
  TFile *file = TFile::Open(filename);
  TTree *itree = dynamic_cast<TTree*>(file->Get("Analysis"));
  TTreeReader myReader("Analysis", file);
 
  int save_plots = cf.Value("ANALYSIS","save_plots") ;

  double pmax_low = cf.Value("ANALYSIS","pmax_low") ;
  double pmax_up = cf.Value("ANALYSIS","pmax_up") ; 
  double area_low = cf.Value("ANALYSIS","area_low") ;
  double area_up = cf.Value("ANALYSIS","area_up") ;
  double negpmax_low = cf.Value("ANALYSIS","negpmax_low") ;

  double trigger_pmax_low = cf.Value("ANALYSIS","trigger_pmax_low") ;
  double trigger_pmax_up = cf.Value("ANALYSIS","trigger_pmax_up") ;
  double trigger_area_low = cf.Value("ANALYSIS","trigger_area_low") ;
  double trigger_area_up = cf.Value("ANALYSIS","trigger_area_up") ;
  double trigger_negpmax_low = cf.Value("ANALYSIS","trigger_negpmax_low") ;


  TH1F *pmax_hist = new TH1F("pmax_hist","Signal amplitude", pmax_up, 0, pmax_up);
  pmax_hist->GetXaxis()->SetTitle("amplitude  [mV]");

  TH1F *area_hist = new TH1F("area_hist","Signal area", area_up, 0, area_up);
  area_hist->GetXaxis()->SetTitle("area  [pWb]");
  
  TH1F *dvdt_hist = new TH1F("dvdt_hist","Slew rate", 200, 0, 1000);
  dvdt_hist->GetXaxis()->SetTitle("dV/dt  [mV/ns]");

  TH1F *risetime_hist = new TH1F("risetime_hist","Risetime", 200, 0, 2);
  risetime_hist->GetXaxis()->SetTitle("Risetime  [ns]");

  TH1F *rms_hist = new TH1F("rms_hist","RMS noise",100,0,10);
  rms_hist->GetXaxis()->SetTitle("RMS  [mV]");

  double timeres_range_low = cf.Value("ANALYSIS","timeres_range_low") ;
  double timeres_range_up = cf.Value("ANALYSIS","timeres_range_up") ;
  TH1F *timeres = new TH1F("timeres", "ToA_{DUT} - ToA_{Trigger}", 200, timeres_range_low, timeres_range_up);
  timeres->GetXaxis()->SetTitle("[ns]");

  

  
  itree->Project("pmax_hist", "pmax[0]", Form("negpmax[0]>%f && negpmax[1]>%f && pmax[1]>%f && pmax[1]<%f && area_new[1]>%f && area_new[1]<%f && pmax[1]>%f && pmax[1]<%f && area_new[1]>%f && area_new[1]<%f", 
                  negpmax_low, trigger_negpmax_low, trigger_pmax_low, trigger_pmax_up, trigger_area_low, trigger_area_up, pmax_low, pmax_up, area_low, area_up ) );
  TF1 *pmax_landau = new TF1( "pmax_landau", "landau", pmax_low, pmax_up );
  pmax_hist->Fit("pmax_landau","RQ");
  double PMAX = pmax_landau->GetParameter(1);
  cout<<"Pmax MPV: "<<PMAX<<endl;

  
  itree->Project("area_hist", "area[0]", Form("negpmax[0]>%f && negpmax[1]>%f && pmax[1]>%f && pmax[1]<%f && area_new[1]>%f && area_new[1]<%f && pmax[1]>%f && pmax[1]<%f && area_new[1]>%f && area_new[1]<%f", 
                  negpmax_low, trigger_negpmax_low, trigger_pmax_low, trigger_pmax_up, trigger_area_low, trigger_area_up, pmax_low, pmax_up, area_low, area_up ) );
  TF1 *area_landau = new TF1( "area_landau", "landau", area_low, area_up );
  area_hist->Fit("area_landau","RQ");
  double AREA = area_landau->GetParameter(1);
  cout<<"Area MPV: "<<AREA<<endl;


  itree->Project("dvdt_hist", "dvdt[0]", Form("area[%i]",ch), Form("negpmax[0]>%f && negpmax[1]>%f && pmax[1]>%f && pmax[1]<%f && area_new[1]>%f && area_new[1]<%f && pmax[1]>%f && pmax[1]<%f && area_new[1]>%f && area_new[1]<%f", 
                  negpmax_low, trigger_negpmax_low, trigger_pmax_low, trigger_pmax_up, trigger_area_low, trigger_area_up, pmax_low, pmax_up, area_low, area_up ) );
  TF1 *dvdt_landau=new TF1( "dvdt_landau", "landau", 0, 1000 );
  dvdt_hist->Fit("dvdt_landau","RQ");
  double dVdt = dvdt_landau->GetParameter(1);
  cout<<"dVdt: "<<dVdt<<endl;


  itree->Project("risetime_hist", "area[0]", Form("negpmax[0]>%f && negpmax[1]>%f && pmax[1]>%f && pmax[1]<%f && area_new[1]>%f && area_new[1]<%f && pmax[1]>%f && pmax[1]<%f && area_new[1]>%f && area_new[1]<%f", 
                  negpmax_low, trigger_negpmax_low, trigger_pmax_low, trigger_pmax_up, trigger_area_low, trigger_area_up, pmax_low, pmax_up, area_low, area_up ) );
  TF1 *risetime_landau=new TF1( "risetime_landau", "gaus", 0, 2 );
  risetime_hist->Fit("risetime_landau","RQ");
  double RISETIME = risetime_landau->GetParameter(1);
  cout<<"Risetime: "<<RISETIME<<endl;


  itree->Project("rms_hist", "rms[0]", Form("negpmax[0]>%f && negpmax[1]>%f && pmax[1]>%f && pmax[1]<%f && area_new[1]>%f && area_new[1]<%f && pmax[1]>%f && pmax[1]<%f && area_new[1]>%f && area_new[1]<%f", 
                  negpmax_low, trigger_negpmax_low, trigger_pmax_low, trigger_pmax_up, trigger_area_low, trigger_area_up, pmax_low, pmax_up, area_low, area_up ) );
  TF1 *rms_gaus=new TF1( "rms_gaus", "gaus", 0, 10 );
  rms_hist->Fit("rms_gaus","RQ");
  double RMS = rms_gaus->GetParameter(1);
  cout<<"RMS noise: "<<RMS<<endl;
  


  TF1 *timeres_gaus=new TF1( "timeres_gaus", "gaus", timeres_range_low, timeres_range_up );
  itree->Project("timeres", Form("cfd[0][%i]-cfd[1][1]", cfd_index), 
          Form("negpmax[0]>%f && negpmax[1]>%f && pmax[1]>%f && pmax[1]<%f && area_new[1]>%f && area_new[1]<%f && pmax[1]>%f && pmax[1]<%f && area_new[1]>%f && area_new[1]<%f", 
                  negpmax_low, trigger_negpmax_low, trigger_pmax_low, trigger_pmax_up, trigger_area_low, trigger_area_up, pmax_low, pmax_up, area_low, area_up ) );
  
  timeres->Fit("timeres_gaus","RQ");
  double RESOLUTION = timeres_gaus->GetParameter(2);
  cout<<Form("Time resolution ch %i: ", channel)<<RESOLUTION<<endl;




  
  if( save_plots==1 ){

    std::string s = filename;
    std::string delimiter = ".root";
    std::string token = s.substr(0, s.find(delimiter)); 

    TString string_out = token;
    string_out += "_PLOTS.root" ;
    TFile f(string_out,"recreate");

    f.cd();
    pmax_hist->Write();
    area_hist->Write();
    dvdt_hist->Write();
    rms_hist->Write();
    risetime_hist->Write();
    timeres->Write();

    f.Write();
    f.Close();

  }


}