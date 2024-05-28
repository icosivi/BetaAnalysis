#define csv_maker_cxx
#include "csv_maker.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <Riostream.h>

void csv_maker::Loop()
{

   TH1F *sample = new TH1F("sample","sample",2048,-1024,1024);

   if (fChain == 0) return;

   Long64_t nentries = fChain->GetEntriesFast();

   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry); 

      if(pmax->at(10)>20){

      int counter=0 ;
      float max=-100000. ;

      for( int i=0; i<1000; i++ ){

        if( w->at(10).at(i)>max ){

          max = w->at(10).at(i) ;
          counter = i ;

        }

      }

      int counter_laser=0 ;
      float max_laser=-100000. ;

      for( int i=0; i<1000; i++ ){

        if( w->at(16).at(i)>max_laser ){

          max_laser = w->at(16).at(i) ;
          counter_laser = i ;

        }

      }

      sample->Fill(counter-counter_laser);
      

      }


      nbytes += nb;
   }

   sample->Draw();

   //f.close();
}
