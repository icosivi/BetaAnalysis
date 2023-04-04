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
#include <TSystemFile.h>
#include <TSystemDirectory.h>

#include "include/general.hpp"
#include "include/Chameleon.h"
#include "include/ConfigFile.hpp"
#include "src/Analyzer.hpp"

void read_analysis_compiler(){


  gSystem->CompileMacro("src/general.cpp","kg");
  gSystem->CompileMacro("src/Chameleon.cpp","kg");
  gSystem->CompileMacro("src/ConfigFile.cpp","kg");
  gSystem->CompileMacro("src/Analyzer.cpp","kg");
  gSystem->CompileMacro("analisi_new.C","kg");
  gSystem->CompileMacro("read_analysis_new.C","kg");

  gROOT->ProcessLine("read_analysis_new()");
  
}

