//{
#include <Riostream.h>
#include <TROOT.h>
#include <TObject.h>
#include <TSystem.h>
#include "TStopwatch.h"
#include "TMath.h"
#include "TAxis.h"
#include "TFile.h"
#include "TH1.h"
#include "TH1D.h"
#include "TCanvas.h"
#include <TTree.h>
#include <TBranch.h>
#include <TSystemDirectory.h>
#include <TSystemFile.h>
#include <TList.h>


#include "include/general.hpp"
#include "include/Chameleon.h"
#include "include/ConfigFile.hpp"
#include "src/Analyzer.hpp"

void Compilatore(){


  gSystem->CompileMacro("src/general.cpp","kg");
  gSystem->CompileMacro("src/Chameleon.cpp","kg");
  gSystem->CompileMacro("src/ConfigFile.cpp","kg");
  gSystem->CompileMacro("src/Analyzer.cpp","kg");
  gSystem->CompileMacro("analisi_TB.C","kg");
  gSystem->CompileMacro("read_analysis_digi.C","kg");
  //gSystem->CompileMacro("FAST.C","kg");


  /*const char *ext=".root"
  const char *dirname="/home/federico/Scrivania/RSD_Digitizer/wafer2-100-200-root/wafer2-100-200-shifted" ;

  TSystemDirectory dir(dirname, dirname); 
  TList *files = dir.GetListOfFiles(); 
  if (files) { 

    TSystemFile *file; 
    TString fname; 
    TIter next(files);

    while ((file=(TSystemFile*)next())) { 
      
      fname = file->GetName(); 
      if (!file->IsDirectory() && fname.EndsWith(ext)) {  

      ConfigFile cf("beta_config.ini");
      int ch_number = cf.Value("HEADER","active_channels");
      int trigger_channel = cf.Value("HEADER","trigger_ch");

      gROOT->ProcessLine(Form( "analisi( %s )", fname.Data() ) );

      } 
    } 
  }*/ 
  
  gROOT->ProcessLine( "analisi_TB( )" );
  //gROOT->ProcessLine("readAnalysis()");
  /*for(int j=1; j<=ch_number; j++){

    if( j != trigger_channel ){

     gROOT->ProcessLine( Form("readAnalysis(%i)", j) );
     //gROOT->ProcessLine( Form("FAST(%i)", j) );

    }

  }*/
  //gROOT->ProcessLine( "readAnalysis(1)" );
  //gROOT->ProcessLine( "FAST(1)" );

}

