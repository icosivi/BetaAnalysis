//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Fri Feb 18 15:10:04 2022 by ROOT version 6.22/09
// from TTree Analysis/Analysis
// found on file: stats_croci14fiocchi-RSD2_Run3_x_0_900_y_0_900.root
//////////////////////////////////////////////////////////

#ifndef n_pads_h
#define n_pads_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.
#include "vector"
#include "vector"

class n_pads {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   Int_t           event;
   Int_t           ntrig;
   vector<vector<double> > *w;
   vector<vector<double> > *t;
   vector<double>  *pmax;
   vector<double>  *negpmax;
   vector<double>  *tmax;
   vector<double>  *negtmax;
   vector<double>  *area;
   vector<double>  *uarea;
   vector<double>  *area_new;
   vector<double>  *uarea_new;
   vector<double>  *risetime;
   vector<double>  *falltime;
   vector<double>  *dvdt;
   vector<double>  *dvdt_2080;
   vector<vector<double> > *cfd;
   vector<vector<double> > *width;
   vector<double>  *t_thr;
   vector<double>  *tot;
   vector<double>  *rms;
   Double_t        x_pos;
   Double_t        y_pos;

   // List of branches
   TBranch        *b_event;   //!
   TBranch        *b_ntrig;   //!
   TBranch        *b_w;   //!
   TBranch        *b_t;   //!
   TBranch        *b_pmax;   //!
   TBranch        *b_negpmax;   //!
   TBranch        *b_tmax;   //!
   TBranch        *b_negtmax;   //!
   TBranch        *b_area;   //!
   TBranch        *b_uarea;   //!
   TBranch        *b_area_new;   //!
   TBranch        *b_uarea_new;   //!
   TBranch        *b_risetime;   //!
   TBranch        *b_falltime;   //!
   TBranch        *b_dvdt;   //!
   TBranch        *b_dvdt_2080;   //!
   TBranch        *b_cfd;   //!
   TBranch        *b_width;   //!
   TBranch        *b_t_thr;   //!
   TBranch        *b_tot;   //!
   TBranch        *b_rms;   //!
   TBranch        *b_x_pos;   //!
   TBranch        *b_y_pos;   //!

   n_pads(TTree *tree=0);
   virtual ~n_pads();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef n_pads_cxx
n_pads::n_pads(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("/media/daq/MartaBackup/RSD2/stats/stats_fiocchi14-RSD2_Run7_x_0_800_y_0_800_100waveforms_330V_train.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("/media/daq/MartaBackup/RSD2/stats/stats_fiocchi14-RSD2_Run7_x_0_800_y_0_800_100waveforms_330V_train.root");
      }
      f->GetObject("Analysis",tree);

   }
   Init(tree);
}

n_pads::~n_pads()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t n_pads::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t n_pads::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void n_pads::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set object pointer
   w = 0;
   t = 0;
   pmax = 0;
   negpmax = 0;
   tmax = 0;
   negtmax = 0;
   area = 0;
   uarea = 0;
   area_new = 0;
   uarea_new = 0;
   risetime = 0;
   falltime = 0;
   dvdt = 0;
   dvdt_2080 = 0;
   cfd = 0;
   width = 0;
   t_thr = 0;
   tot = 0;
   rms = 0;
   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("event", &event, &b_event);
   fChain->SetBranchAddress("ntrig", &ntrig, &b_ntrig);
   fChain->SetBranchAddress("w", &w, &b_w);
   fChain->SetBranchAddress("t", &t, &b_t);
   fChain->SetBranchAddress("pmax", &pmax, &b_pmax);
   fChain->SetBranchAddress("negpmax", &negpmax, &b_negpmax);
   fChain->SetBranchAddress("tmax", &tmax, &b_tmax);
   fChain->SetBranchAddress("negtmax", &negtmax, &b_negtmax);
   fChain->SetBranchAddress("area", &area, &b_area);
   fChain->SetBranchAddress("uarea", &uarea, &b_uarea);
   fChain->SetBranchAddress("area_new", &area_new, &b_area_new);
   fChain->SetBranchAddress("uarea_new", &uarea_new, &b_uarea_new);
   fChain->SetBranchAddress("risetime", &risetime, &b_risetime);
   fChain->SetBranchAddress("falltime", &falltime, &b_falltime);
   fChain->SetBranchAddress("dvdt", &dvdt, &b_dvdt);
   fChain->SetBranchAddress("dvdt_2080", &dvdt_2080, &b_dvdt_2080);
   fChain->SetBranchAddress("cfd", &cfd, &b_cfd);
   fChain->SetBranchAddress("width", &width, &b_width);
   fChain->SetBranchAddress("t_thr", &t_thr, &b_t_thr);
   fChain->SetBranchAddress("tot", &tot, &b_tot);
   fChain->SetBranchAddress("rms", &rms, &b_rms);
   fChain->SetBranchAddress("x_pos", &x_pos, &b_x_pos);
   fChain->SetBranchAddress("y_pos", &y_pos, &b_y_pos);
   Notify();
}

Bool_t n_pads::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void n_pads::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t n_pads::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef n_pads_cxx
