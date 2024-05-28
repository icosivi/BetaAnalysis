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
#include <TProfile2D.h>
#include <TH1.h>
#include <TH3.h>
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
#include <TChain.h>

//------Custom headers----------------//
#include "include/general.hpp"
#include "include/Chameleon.h"
#include "include/ConfigFile.hpp"


void readAnalysis( ){

  gStyle->SetTitleOffset( 1.3, "x" );
  gStyle->SetTitleOffset( 1.3, "y" );
  gStyle->SetHistLineWidth(2);
  gStyle->SetOptStat(0);
  //gStyle->SetPadTickX(1);  // To get tick marks on the opposite side of the frame
  //gStyle->SetPadTickY(1);
 
  ConfigFile cf("beta_config.ini");

  TChain *itree = new TChain("Analysis","");
  //itree->Add("/media/daq/UFSD-Disk2/RSD2/stats/stats_croci15boxes-RSD2_Run5_x_0_800_y_0_800_100waveforms_200V_250V_300V_RenameThis.root");
  itree->Add("/media/daq/MartaBackup/RSD2/stats/stats_fiocchi14-RSD2_Run7_x_0_800_y_0_800_100waveforms_330V_train.root");
  //itree->Add("/media/daq/UFSD-Disk2/RSD2/stats/stats_croci15boxes-RSD2_Run5_x_0_800_y_0_800_100waveforms_200V_250V_300V_RenameThis_1.root");
  //itree->Add("/media/daq/UFSD-Disk1/RSD2/stats_croci450RSD2-x_230_1700_10-y_230_1700_10.root");
  //itree->Add("/media/daq/UFSD-Disk1/RSD2/stats_croci450RSD2-x_230_1700_10-y_230_1700_10_1.root");
  //itree->Add("/media/daq/UFSD-Disk1/RSD2/stats_croci450RSD2-x_230_1700_10-y_230_1700_10_2.root");

  //itree->Add("/media/daq/UFSD-Disk2/RSD2/stats/stats_croci14fiocchi-RSD2_Run3_x_0_900_y_0_900.root");
  
  int trigger_channel = cf.Value("HEADER","trigger_ch");
  int trigger_ch = trigger_channel - 1 ;

  const double conversion_factor = cf.Value("HEADER","conversion_factor_area_charge") ;
  const double time_const = cf.Value("HEADER","time_scalar");  
  const double voltage_const = cf.Value("HEADER","voltage_scalar");

  int save_plots = cf.Value("ANALYSIS","save_plots") ;
  int best_cfd_index = cf.Value("ANALYSIS","cfd_index") ; 

  double negpmax = cf.Value("ANALYSIS","negpmax") ; 
  double pmax_low = 0;
  double pmax_up = cf.Value("ANALYSIS","pmax_up") ; 
  double area_low = 0;
  double area_up = cf.Value("ANALYSIS","area_up") ;

  TH3F *pmax_3d[16];
  TProfile2D *pmax_2d[16];
  TCanvas *cp[16];

  TH3F *negpmax_3d[16];
  TProfile2D *negpmax_2d[16];
  TCanvas *cnp[16];

  TH3F *area_3d[16];
  TProfile2D *area_2d[16];
  TCanvas *ca[16];

  TH3F *tmax_3d[16];
  TProfile2D *tmax_2d[16];
  TCanvas *ct[16];

  TH3F *negtmax_3d[16];
  TProfile2D *negtmax_2d[16];
  TCanvas *cnt[16];

  TH3F *cfd_3d[16];
  TProfile2D *cfd_2d[16];
  TCanvas *ccfd[16];

  TH3F *rise_3d[16];
  TProfile2D *rise_2d[16];
  TCanvas *crise[16];

  for(int k=0; k<15; k++){  
  
    pmax_3d[k] = new TH3F(Form("pmax_3d_map_ch%i", k), Form("Pmax ch%i", k),90,0,900,90,0,900,300,0,300);
  
    itree->Project(Form("pmax_3d_map_ch%i", k), Form("pmax[%i]:x_pos:y_pos", k), 
                   Form("pmax[%i]>%f && pmax[%i]<%f && area_new[%i]>%f && cfd[%i][1]>0 && risetime[%i]>0 && risetime[%i]<2 && negpmax[%i]>-350 && negpmax[%i]<0 ", 
                   k, pmax_low, k, pmax_up, k, area_low, k, k, k, k, k) );

    pmax_2d[k] = pmax_3d[k]->Project3DProfile("xy");
    pmax_2d[k]->GetXaxis()->SetTitle("x  [um]");
    pmax_2d[k]->GetYaxis()->SetTitle("y  [um]");
    pmax_2d[k]->GetYaxis()->SetTitleOffset(1.5);
    pmax_2d[k]->GetZaxis()->SetTitle("Pmax  [mV]");
    pmax_2d[k]->GetZaxis()->SetRangeUser(0,50);


    /*negpmax_3d[k] = new TH3F(Form("negpmax_3d_map_ch%i", k), Form("negPmax ch%i", k),90,0,900,90,0,900,350,-350,0);
  
    itree->Project(Form("negpmax_3d_map_ch%i", k), Form("negpmax[%i]:x_pos:y_pos", k), 
                   Form("pmax[%i]>%f && pmax[%i]<%f && area_new[%i]>%f && cfd[%i][1]>0 && risetime[%i]>0 && risetime[%i]<2 && negpmax[%i]>-350 && negpmax[%i]<0",
                   k, pmax_low, k, pmax_up, k, area_low, k, k, k, k, k) );

    negpmax_2d[k] = negpmax_3d[k]->Project3DProfile("xy");
    negpmax_2d[k]->GetXaxis()->SetTitle("x  [um]");
    negpmax_2d[k]->GetYaxis()->SetTitle("y  [um]");
    negpmax_2d[k]->GetYaxis()->SetTitleOffset(1.5);
    negpmax_2d[k]->GetZaxis()->SetTitle("negPmax  [mV]");
    negpmax_2d[k]->GetZaxis()->SetRangeUser(-350,0);


    area_3d[k] = new TH3F(Form("area_3d_map_ch%i", k), Form("Signal area ch%i", k),90,0,900,90,0,900,300,0,300);
  
    itree->Project(Form("area_3d_map_ch%i", k), Form("area_new[%i]:x_pos:y_pos", k), 
                   Form("pmax[%i]>%f && pmax[%i]<%f && area_new[%i]>%f && cfd[%i][1]>0 && risetime[%i]>0 && risetime[%i]<2 && negpmax[%i]>-350 && negpmax[%i]<0",
                   k, pmax_low, k, pmax_up, k, area_low, k, k, k, k, k) );

    area_2d[k] = area_3d[k]->Project3DProfile("xy");
    area_2d[k]->GetXaxis()->SetTitle("x  [um]");
    area_2d[k]->GetYaxis()->SetTitle("y  [um]");
    area_2d[k]->GetYaxis()->SetTitleOffset(1.5);
    area_2d[k]->GetZaxis()->SetTitle("Area  [pWb]");
    area_2d[k]->GetZaxis()->SetRangeUser(0,40);


    tmax_3d[k] = new TH3F(Form("tmax_3d_map_ch%i", k), Form("Tmax ch%i", k),90,0,900,90,0,900,5000,0,50);
  
    itree->Project(Form("tmax_3d_map_ch%i", k), Form("tmax[%i]:x_pos:y_pos", k), 
                   Form("pmax[%i]>%f && pmax[%i]<%f && area_new[%i]>%f && cfd[%i][1]>0 && risetime[%i]>0 && risetime[%i]<2 && negpmax[%i]>-350 && negpmax[%i]<0", 
                   k, pmax_low, k, pmax_up, k, area_low, k, k, k, k, k) );

    tmax_2d[k] = tmax_3d[k]->Project3DProfile("xy");
    tmax_2d[k]->GetXaxis()->SetTitle("x  [um]");
    tmax_2d[k]->GetYaxis()->SetTitle("y  [um]");
    tmax_2d[k]->GetYaxis()->SetTitleOffset(1.5);
    tmax_2d[k]->GetZaxis()->SetTitle("Tmax  [ns]");
    tmax_2d[k]->GetZaxis()->SetRangeUser(0,50);


    negtmax_3d[k] = new TH3F(Form("negtmax_3d_map_ch%i", k), Form("negTmax ch%i", k),90,0,900,90,0,900,5000,0,50);
  
    itree->Project(Form("negtmax_3d_map_ch%i", k), Form("negtmax[%i]:x_pos:y_pos", k), 
                   Form("pmax[%i]>%f && pmax[%i]<%f && area_new[%i]>%f && cfd[%i][1]>0 && risetime[%i]>0 && risetime[%i]<2 && negpmax[%i]>-350 && negpmax[%i]<0", 
                   k, pmax_low, k, pmax_up, k, area_low, k, k, k, k, k) );

    negtmax_2d[k] = negtmax_3d[k]->Project3DProfile("xy");
    negtmax_2d[k]->GetXaxis()->SetTitle("x  [um]");
    negtmax_2d[k]->GetYaxis()->SetTitle("y  [um]");
    negtmax_2d[k]->GetYaxis()->SetTitleOffset(1.5);
    negtmax_2d[k]->GetZaxis()->SetTitle("negTmax  [ns]");
    negtmax_2d[k]->GetZaxis()->SetRangeUser(0,50);


    cfd_3d[k] = new TH3F(Form("cfd_3d_map_ch%i", k), Form("CFD20 ch%i", k),90,0,900,90,0,900,5500,0,55);
  
    itree->Project(Form("cfd_3d_map_ch%i", k), Form("cfd[%i][1]:x_pos:y_pos", k), 
                   Form("pmax[%i]>%f && pmax[%i]<%f && area_new[%i]>%f && cfd[%i][1]>0 && risetime[%i]>0 && risetime[%i]<2 && negpmax[%i]>-350 && negpmax[%i]<0", 
                   k, pmax_low, k, pmax_up, k, area_low, k, k, k, k, k) );

    cfd_2d[k] = cfd_3d[k]->Project3DProfile("xy");
    cfd_2d[k]->GetXaxis()->SetTitle("x  [um]");
    cfd_2d[k]->GetYaxis()->SetTitle("y  [um]");
    cfd_2d[k]->GetYaxis()->SetTitleOffset(1.5);
    cfd_2d[k]->GetZaxis()->SetTitle("CFD20  [ns]");
    cfd_2d[k]->GetZaxis()->SetRangeUser(0,55);


    rise_3d[k] = new TH3F(Form("rise_3d_map_ch%i", k), Form("Risetime ch%i", k),90,0,900,90,0,900,700,0,0.7);
  
    itree->Project(Form("rise_3d_map_ch%i", k), Form("risetime[%i]:x_pos:y_pos", k), 
                   Form("pmax[%i]>%f && pmax[%i]<%f && area_new[%i]>%f && cfd[%i][1]>0 && risetime[%i]>0 && risetime[%i]<2 && negpmax[%i]>-350 && negpmax[%i]<0", 
                   k, pmax_low, k, pmax_up, k, area_low, k, k, k, k, k) );

    rise_2d[k] = rise_3d[k]->Project3DProfile("xy");
    rise_2d[k]->GetXaxis()->SetTitle("x  [um]");
    rise_2d[k]->GetYaxis()->SetTitle("y  [um]");
    rise_2d[k]->GetYaxis()->SetTitleOffset(1.5);
    rise_2d[k]->GetZaxis()->SetTitle("Risetime  [ns]");
    rise_2d[k]->GetZaxis()->SetRangeUser(0,0.7);*/

  }

  for(int j=0; j<15; j++){

     cp[j] = new TCanvas(Form("cp%i",j),Form("cp%i",j),1200,1200);
     cp[j]->cd();
     pmax_2d[j]->Draw("colz");
     pmax_2d[j]->SaveAs( Form("/media/daq/MartaBackup/RSD2/documents_pics/root_files/fiocchi14_Run7_330V/pmax_2d-map_ch%i.root", j) );

     /*cnp[j] = new TCanvas(Form("cnp%i",j),Form("cnp%i",j),1200,1200);
     cnp[j]->cd();
     negpmax_2d[j]->Draw("colz");
     negpmax_2d[j]->SaveAs( Form("/media/daq/UFSD-Disk2/RSD2/documents_pics/root_files/negpmax_2d-map_ch%i.root", j) );

     ca[j] = new TCanvas(Form("ca%i",j),Form("ca%i",j),1200,1200);
     ca[j]->cd();
     area_2d[j]->Draw("colz");
     area_2d[j]->SaveAs( Form("/media/daq/UFSD-Disk2/RSD2/documents_pics/root_files/area_2d-map_ch%i.root", j) );

     ct[j] = new TCanvas(Form("ct%i",j),Form("ct%i",j),1200,1200);
     ct[j]->cd();
     tmax_2d[j]->Draw("colz");
     tmax_2d[j]->SaveAs( Form("/media/daq/UFSD-Disk2/RSD2/documents_pics/root_files/tmax_2d-map_ch%i.root", j) );

     cnt[j] = new TCanvas(Form("cnt%i",j),Form("cnt%i",j),1200,1200);
     cnt[j]->cd();
     negtmax_2d[j]->Draw("colz");
     negtmax_2d[j]->SaveAs( Form("/media/daq/UFSD-Disk2/RSD2/documents_pics/root_files/negtmax_2d-map_ch%i.root", j) );

     ccfd[j] = new TCanvas(Form("ccfd%i",j),Form("ccfd%i",j),1200,1200);
     ccfd[j]->cd();
     cfd_2d[j]->Draw("colz");
     cfd_2d[j]->SaveAs( Form("/media/daq/UFSD-Disk2/RSD2/documents_pics/root_files/cfd20_2d-map_ch%i.root", j) );

     crise[j] = new TCanvas(Form("crise%i",j),Form("crise%i",j),1200,1200);
     crise[j]->cd();
     rise_2d[j]->Draw("colz");
     rise_2d[j]->SaveAs( Form("/media/daq/UFSD-Disk2/RSD2/documents_pics/root_files/rise_2d-map_ch%i.root", j) );*/

  }


  /*TH1F *pmax_distrib[7];
  TH1F *area_distrib[7];
  TH1F *negpmax_distrib[7];
  TH1F *tmax_distrib[7];
  TH1F *rise_distrib[7];
  //float xx_pos[7]={350,500,700,350,450,600,700};
  //float yy_pos[7]={925,925,925,1100,1050,800,750};
  float xx_pos[8]={100,200,300,350,450,550,650,750};
  float yy_pos[8]={300,300,300,300,300,300,300,300};

  TCanvas *c_pmax_distrib[7];
  TCanvas *c_area_distrib[7];
  TCanvas *c_negpmax_distrib[7];
  TCanvas *c_tmax_distrib[7];
  TCanvas *c_rise_distrib[7];

  int pads[12]={1,2,3,4,5,6,8,9,10,11,13,14};

  for(int l=0; l<8; l++){
 
    pmax_distrib[l] = new TH1F( Form("pmax_pos_%i",l), Form("Pmax distribution - x=%f y=%f", xx_pos[l], yy_pos[l]), 16, 0, 16 );
    area_distrib[l] = new TH1F( Form("area_pos_%i",l), Form("Area distribution - x=%f y=%f", xx_pos[l], yy_pos[l]), 16, 0, 16 );
    negpmax_distrib[l] = new TH1F( Form("negpmax_pos_%i",l), Form("negPmax distribution - x=%f y=%f", xx_pos[l], yy_pos[l]), 16, 0, 16 );
    tmax_distrib[l] = new TH1F( Form("tmax_pos_%i",l), Form("Tmax distribution - x=%f y=%f", xx_pos[l], yy_pos[l]), 16, 0, 16 );
    rise_distrib[l] = new TH1F( Form("rise_pos_%i",l), Form("Risetime distribution - x=%f y=%f", xx_pos[l], yy_pos[l]), 16, 0, 16 );

    for(int k=0; k<12; k++){

      TH1F *pmax_distribution = new TH1F("pmax_distribution","pmax_distribution",100,0,100);
      itree->Project("pmax_distribution", Form("pmax[%i]", pads[k]), 
                     Form("pmax[%i]>%f && pmax[%i]<%f && area_new[%i]>%f && cfd[%i][1]>0 && risetime[%i]>0 && risetime[%i]<2 && negpmax[%i]>-350 && negpmax[%i]<0 && x_pos>%f && x_pos<%f && y_pos>%f && y_pos<%f", 
                           pads[k], pmax_low, pads[k], pmax_up, pads[k], area_low, pads[k], pads[k], pads[k], pads[k], pads[k], (xx_pos[l]-10), (xx_pos[l]+10), (yy_pos[l]-10), (yy_pos[l]+10) ) );

      pmax_distrib[l]->Fill( pads[k], pmax_distribution->GetMean() );
      delete pmax_distribution;


      TH1F *area_distribution = new TH1F("area_distribution","area_distribution",50,0,50);
      itree->Project("area_distribution", Form("area[%i]", pads[k]), 
                     Form("pmax[%i]>%f && pmax[%i]<%f && area_new[%i]>%f && cfd[%i][1]>0 && risetime[%i]>0 && risetime[%i]<2 && negpmax[%i]>-350 && negpmax[%i]<0 && x_pos>%f && x_pos<%f && y_pos>%f && y_pos<%f", 
                           pads[k], pmax_low, pads[k], pmax_up, pads[k], area_low, pads[k], pads[k], pads[k], pads[k], pads[k], (xx_pos[l]-10), (xx_pos[l]+10), (yy_pos[l]-10), (yy_pos[l]+10) ) );

      area_distrib[l]->Fill( pads[k], area_distribution->GetMean() );
      delete area_distribution;


      TH1F *negpmax_distribution = new TH1F("negpmax_distribution","negpmax_distribution",300,-300,0);
      itree->Project("negpmax_distribution", Form("negpmax[%i]", pads[k]), 
                     Form("pmax[%i]>%f && pmax[%i]<%f && area_new[%i]>%f && cfd[%i][1]>0 && risetime[%i]>0 && risetime[%i]<2 && negpmax[%i]>-350 && negpmax[%i]<0 && x_pos>%f && x_pos<%f && y_pos>%f && y_pos<%f", 
                           pads[k], pmax_low, pads[k], pmax_up, pads[k], area_low, pads[k], pads[k], pads[k], pads[k], pads[k], (xx_pos[l]-10), (xx_pos[l]+10), (yy_pos[l]-10), (yy_pos[l]+10) ) );

      negpmax_distrib[l]->Fill( pads[k], negpmax_distribution->GetMean() );
      delete negpmax_distribution;


      TH1F *tmax_distribution = new TH1F("tmax_distribution","Tmax_distribution",500,0,50);
      itree->Project("tmax_distribution", Form("tmax[%i]", pads[k]), 
                     Form("pmax[%i]>%f && pmax[%i]<%f && area_new[%i]>%f && cfd[%i][1]>0 && risetime[%i]>0 && risetime[%i]<2 && negpmax[%i]>-350 && negpmax[%i]<0 && x_pos>%f && x_pos<%f && y_pos>%f && y_pos<%f", 
                           pads[k], pmax_low, pads[k], pmax_up, pads[k], area_low, pads[k], pads[k], pads[k], pads[k], pads[k], (xx_pos[l]-10), (xx_pos[l]+10), (yy_pos[l]-10), (yy_pos[l]+10) ) );

      tmax_distrib[l]->Fill( pads[k], tmax_distribution->GetMean() );
      delete tmax_distribution;


      TH1F *rise_distribution = new TH1F("rise_distribution","Risetime_distribution",500,0,2);
      itree->Project("rise_distribution", Form("risetime[%i]", pads[k]), 
                     Form("pmax[%i]>%f && pmax[%i]<%f && area_new[%i]>%f && cfd[%i][1]>0 && risetime[%i]>0 && risetime[%i]<2 && negpmax[%i]>-350 && negpmax[%i]<0 && x_pos>%f && x_pos<%f && y_pos>%f && y_pos<%f", 
                           pads[k], pmax_low, pads[k], pmax_up, pads[k], area_low, pads[k], pads[k], pads[k], pads[k], pads[k], (xx_pos[l]-10), (xx_pos[l]+10), (yy_pos[l]-10), (yy_pos[l]+10) ) );

      rise_distrib[l]->Fill( pads[k], rise_distribution->GetMean() );
      delete rise_distribution;

    }

    }


  for(int l=0; l<8; l++){

    c_pmax_distrib[l] = new TCanvas(Form("c_pmax_distrib_pos%i",l),Form("c_pmax_distrib_pos%i",l),1200,1200);
    c_pmax_distrib[l]->cd();
    pmax_distrib[l]->Draw("hist");
    pmax_distrib[l]->SaveAs( Form("/media/daq/UFSD-Disk2/RSD2/documents_pics/root_files/pmax_distribution_pos%i.root", l) );

    c_area_distrib[l] = new TCanvas(Form("c_area_distrib_pos%i",l),Form("c_area_distrib_pos%i",l),1200,1200);
    c_area_distrib[l]->cd();
    area_distrib[l]->Draw("hist");
    area_distrib[l]->SaveAs( Form("/media/daq/UFSD-Disk2/RSD2/documents_pics/root_files/area_distribution_pos%i.root", l) );

    c_negpmax_distrib[l] = new TCanvas(Form("c_negpmax_distrib_pos%i",l),Form("c_negpmax_distrib_pos%i",l),1200,1200);
    c_negpmax_distrib[l]->cd();
    negpmax_distrib[l]->Draw("hist");
    negpmax_distrib[l]->SaveAs( Form("/media/daq/UFSD-Disk2/RSD2/documents_pics/root_files/negpmax_distribution_pos%i.root", l) );

    c_tmax_distrib[l] = new TCanvas(Form("c_tmax_distrib_pos%i",l),Form("c_tmax_distrib_pos%i",l),1200,1200);
    c_tmax_distrib[l]->cd();
    tmax_distrib[l]->Draw("hist");
    tmax_distrib[l]->SaveAs( Form("/media/daq/UFSD-Disk2/RSD2/documents_pics/root_files/tmax_distribution_pos%i.root", l) );

    c_rise_distrib[l] = new TCanvas(Form("c_rise_distrib_pos%i",l),Form("c_rise_distrib_pos%i",l),1200,1200);
    c_rise_distrib[l]->cd();
    rise_distrib[l]->Draw("hist");
    rise_distrib[l]->SaveAs( Form("/media/daq/UFSD-Disk2/RSD2/documents_pics/root_files/rise_distribution_pos%i.root", l) );
    
  }*/



}
