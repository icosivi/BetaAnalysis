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


void FAST( int channel ){

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


  TH1F *tot_hist = new TH1F("tot_hist",Form("ToT%i", channel ),1000,0,20);
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


  itree->Project("tot_hist", Form("tot[%i]",ch), Form("negpmax[%i]>%f && pmax[%i]>%f && area[%i]>%f", 
                  ch, negpmax, ch, pmax_low, ch, area_low) );
  TF1 *tot_landau=new TF1( "tot_landau", "landau", 0, 20 );
  tot_hist->Fit("tot_landau","RQ");
  double tot_range[2] = { tot_landau->GetParameter(1)-3*tot_landau->GetParameter(2), tot_landau->GetParameter(1)+3*tot_landau->GetParameter(2) };
  //double tot_range[2] = { 3, 10 };


  TF1 *timeres_gaus=new TF1( "timeres_gaus", "gaus", timeres_search_range[0], timeres_search_range[1] );

  itree->Project("timeres", Form("t_thr[%i]-cfd[%i][1]", ch, trigger_ch ), 
          Form("negpmax[%i]>%f && pmax[%i]>%f && area[%i]>%f && negpmax[%i]>-40 && pmax[%i]>%f && pmax[%i]<%f && area[%i]>%f && area[%i]<%f", 
          ch, negpmax, ch, pmax_low, ch, area_low,
          trigger_ch, trigger_ch, trigger_pmax_low, trigger_ch, trigger_pmax_up, trigger_ch, trigger_area_low, trigger_ch, trigger_area_up) );
  
  timeres->Fit("timeres_gaus","RQ");
  double RESOLUTION = timeres_gaus->GetParameter(2);
  cout<<Form("Time resolution ch %i: ", channel)<<RESOLUTION<<endl;


  int tot_bin = ((tot_range[1]-tot_range[0])/2.)*100 ;
  TProfile *tot_prof = new TProfile("tot_prof", Form("ToA vs ToT profile %i",channel), tot_bin, tot_range[0], tot_range[1] );

  itree->Project("tot_prof", Form("t_thr[%i]:tot[%i]", ch, ch), 
                  Form("negpmax[%i]>%f && pmax[%i]>%f && area[%i]>%f && negpmax[%i]>-40 && pmax[%i]>%f && pmax[%i]<%f && area[%i]>%f && area[%i]<%f", 
                        ch, negpmax, ch, pmax_low, ch, area_low, trigger_ch, trigger_ch, trigger_pmax_low, trigger_ch, trigger_pmax_up,
                        trigger_ch, trigger_area_low, trigger_ch, trigger_area_up) );

  //itree->Project("tot_prof", Form("t_thr[%i]:tot[%i]", ch, ch) );

  double chi2[3]={0,0,0};

  TF1 *tot_prof_fit0=new TF1( "tot_prof_fit0", "[0]/x", tot_range[0], tot_range[1] );
  tot_prof_fit0->SetParameter(0, tmax_hist->GetMean()*tot_hist->GetMean() );
  tot_prof->Fit("tot_prof_fit0","RQ");
  chi2[0] = tot_prof_fit0->GetChisquare() ;


  TF1 *tot_prof_fit1=new TF1( "tot_prof_fit1", "pol1", tot_range[0], tot_range[1] );
  tot_prof->Fit("tot_prof_fit1","RQ");
  chi2[1] = tot_prof_fit1->GetChisquare() ;


  TF1 *tot_prof_fit2=new TF1( "tot_prof_fit2", "pol2", tot_range[0], tot_range[1] );
  tot_prof->Fit("tot_prof_fit2","RQ");
  chi2[2] = tot_prof_fit2->GetChisquare() ;


  int lowest_chi2 = 0;
  double l_chi2 = chi2[0] ;

  for( int k=0; k<3; k++){

    if( chi2[k] < l_chi2 ){

      l_chi2 = chi2[k] ;
      lowest_chi2 = k ;

    }
  }

  //TProfile *tot_prof_corrected = new TProfile("tot_prof_corrected","ToA vs ToT profile 1", tot_bin, tot_range[0], tot_range[1]);
  //TF1 *tot_prof_fit_lin=new TF1( "tot_prof_fit_lin", "pol1", tot_range[0], tot_range[1] );
  //TF1 *tot_prof_fit_const=new TF1( "tot_prof_fit_const", "pol0", tot_range[0], tot_range[1] );

  TF1 *timeres_corrected_gaus=new TF1( "timeres_corrected_gaus", "gaus", -2, 2 );
  TH1F *timeres_corrected = new TH1F("timeres_corrected",Form("ToA_{Ch%i} - ToA_{Trigger}", channel ),
                                      400, -2, 2);
  timeres_corrected->GetXaxis()->SetTitle("[ns]");


  TProfile *tot_prof_corrected = new TProfile("tot_prof_corrected", Form("ToA vs ToT profile %i",channel), tot_bin, tot_range[0], tot_range[1] );

  if( lowest_chi2 == 0 ){

    itree->Project("timeres_corrected", Form("(t_thr[%i]-%f/tot[%i]) - cfd[%i][1]", 
                  ch, tot_prof_fit0->GetParameter(0), ch, trigger_ch), 
                  Form("negpmax[%i]>%f && pmax[%i]>%f && area[%i]>%f && negpmax[%i]>-40 && pmax[%i]>%f && pmax[%i]<%f && area[%i]>%f && area[%i]<%f", 
                        ch, negpmax, ch, pmax_low, ch, area_low, trigger_ch, trigger_ch, trigger_pmax_low, trigger_ch, trigger_pmax_up,
                        trigger_ch, trigger_area_low, trigger_ch, trigger_area_up) );

    timeres_corrected->Fit("timeres_corrected_gaus","RQ");

    //tot_prof_corrected->Fit("tot_prof_fit_lin","RQ");
    //tot_prof_corrected->Fit("tot_prof_fit_const","RQ");
  

  }else if( lowest_chi2 == 1){

    itree->Project("timeres_corrected", Form("(t_thr[%i]-%f-%f*tot[%i]) - cfd[%i][1]", 
                  ch, tot_prof_fit1->GetParameter(0), tot_prof_fit1->GetParameter(1), ch, trigger_ch), 
                  Form("negpmax[%i]>%f && pmax[%i]>%f && area[%i]>%f && negpmax[%i]>-40 && pmax[%i]>%f && pmax[%i]<%f && area[%i]>%f && area[%i]<%f", 
                        ch, negpmax, ch, pmax_low, ch, area_low, trigger_ch, trigger_ch, trigger_pmax_low, trigger_ch, trigger_pmax_up,
                        trigger_ch, trigger_area_low, trigger_ch, trigger_area_up) );

    timeres_corrected->Fit("timeres_corrected_gaus","RQ");

    //tot_prof_corrected->Fit("tot_prof_fit_lin","RQ");
    //tot_prof_corrected->Fit("tot_prof_fit_const","RQ");

     itree->Project("tot_prof_corrected", Form("(t_thr[%i]-%f*tot[%i]):tot[0]", 
                  ch, tot_prof_fit1->GetParameter(1), ch), 
                  Form("negpmax[%i]>%f && pmax[%i]>%f && area[%i]>%f && negpmax[%i]>-40 && pmax[%i]>%f && pmax[%i]<%f && area[%i]>%f && area[%i]<%f", 
                        ch, negpmax, ch, pmax_low, ch, area_low, trigger_ch, trigger_ch, trigger_pmax_low, trigger_ch, trigger_pmax_up,
                        trigger_ch, trigger_area_low, trigger_ch, trigger_area_up) );
  
  }else if( lowest_chi2 == 2 ){

    itree->Project("timeres_corrected", Form("(t_thr[%i]-%f-%f*tot[%i]-%f*tot[%i]*tot[%i]) - cfd[%i][1]", 
                  ch, tot_prof_fit2->GetParameter(0), tot_prof_fit2->GetParameter(1), ch, tot_prof_fit2->GetParameter(2), ch, ch, trigger_ch), 
                  Form("negpmax[%i]>%f && pmax[%i]>%f && area[%i]>%f && negpmax[%i]>-40 && pmax[%i]>%f && pmax[%i]<%f && area[%i]>%f && area[%i]<%f", 
                        ch, negpmax, ch, pmax_low, ch, area_low, trigger_ch, trigger_ch, trigger_pmax_low, trigger_ch, trigger_pmax_up,
                        trigger_ch, trigger_area_low, trigger_ch, trigger_area_up) );

    timeres_corrected->Fit("timeres_corrected_gaus","RQ");

    //tot_prof_corrected->Fit("tot_prof_fit_lin","RQ");
    //tot_prof_corrected->Fit("tot_prof_fit_const","RQ");

    itree->Project("tot_prof_corrected", Form("(t_thr[%i]-%f*tot[%i]-%f*tot[%i]*tot[%i]):tot[0]", 
                  ch, tot_prof_fit2->GetParameter(1), ch, tot_prof_fit2->GetParameter(2), ch, ch), 
                  Form("negpmax[%i]>%f && pmax[%i]>%f && area[%i]>%f && negpmax[%i]>-40 && pmax[%i]>%f && pmax[%i]<%f && area[%i]>%f && area[%i]<%f", 
                        ch, negpmax, ch, pmax_low, ch, area_low, trigger_ch, trigger_ch, trigger_pmax_low, trigger_ch, trigger_pmax_up,
                        trigger_ch, trigger_area_low, trigger_ch, trigger_area_up) );

 }


  double TimeRes_corr = timeres_corrected_gaus->GetParameter(2);
  //timeres->Draw();
  timeres_corrected->Draw("sames");
  //tot_prof->Draw();
  //tot_prof_corrected->Draw("sames");
  cout<<"Time res corrected: "<<TimeRes_corr<<endl;
  cout<<" "<<endl;
  cout<< Form( "////////// Stopping Channel %i ///////////", channel) <<endl;
  cout<<" "<<endl;
  cout<<" "<<endl;
  


  

  if( save_plots==1 ){

    std::string s = filename;
    std::string delimiter = ".root";
    std::string token = s.substr(0, s.find(delimiter)); 

    TString string_out = token;
    string_out += Form("_FAST_Channel%i.root", channel ) ;
    TFile f(string_out,"recreate");

    f.cd();
    tot_hist->Write();
    tot_prof->Write();
    timeres->Write();
    timeres_corrected->Write();

    f.Write();
    f.Close();

  }


}


/*int main(){

readAnalysis();

return 0;

}*/
