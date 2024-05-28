#define n_pads_cxx
#include "n_pads.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

void n_pads::Loop()
{

   TH2F *pad = new TH2F("pad","# active pads with thr 20mV",80,0,800,80,0,800);
   TH2F *pad_bis = new TH2F("pad2","pad2",80,0,800,80,0,800);
   //TH1F *pad = new TH1F("pad","pad",5,0,5);

   float thr = 20. ;

   if (fChain == 0) return;

   Long64_t nentries = fChain->GetEntriesFast();

   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   

      float counter=0;

      if(pmax->at(1)>thr) counter++ ;
      if(pmax->at(2)>thr) counter++ ;
      if(pmax->at(3)>thr) counter++ ;
      if(pmax->at(4)>thr) counter++ ;
      if(pmax->at(5)>thr) counter++ ;
      if(pmax->at(6)>thr) counter++ ;
      if(pmax->at(8)>thr) counter++ ;
      if(pmax->at(9)>thr) counter++ ;
      if(pmax->at(10)>thr) counter++ ;
      if(pmax->at(11)>thr) counter++ ;
      if(pmax->at(13)>thr) counter++ ;
      if(pmax->at(14)>thr) counter++ ;

      pad->Fill(x_pos, y_pos, counter);
      pad_bis->Fill(x_pos,y_pos);
      //pad->Fill(counter);

      nbytes += nb;
   }

   pad->Divide(pad_bis);
   pad->Draw("colz");
   //pad->Draw();

}
