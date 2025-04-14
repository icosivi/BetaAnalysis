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

void Compilatore_csv(){

  gSystem->CompileMacro("src/general.cpp","kg");
  gSystem->CompileMacro("src/Chameleon.cpp","kg");
  gSystem->CompileMacro("src/ConfigFile.cpp","kg");
  gSystem->CompileMacro("src/Analyzer.cpp","kg");
  gSystem->CompileMacro("analisi_csv_300um.C","kg");
  //gSystem->CompileMacro("analisi_csv_500um.C","kg");

  gROOT->ProcessLine("analisi_csv()");

}

