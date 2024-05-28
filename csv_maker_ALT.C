#define csv_maker_cxx
#include "csv_maker.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <Riostream.h>

void csv_maker::Loop()
{

   if (fChain == 0) return;

   ofstream f;
   f.open("/media/daq/Elements/RSD2/csv/stats_fiocchi14-RSD2_Run7_x_5_805_y_5_805_10waveforms_330V_test_Samplings.csv");

   f << "x" << "," << "y" << ","
     << "pmax0_0" << "," << "pmax0_1" << "," << "pmax0_2" << "," << "pmax0_3" << "," << "pmax0_4" << "," 
     << "pmax1_0" << "," << "pmax1_1" << "," << "pmax1_2" << "," << "pmax1_3" << "," << "pmax1_4" << ","  
     << "pmax2_0" << "," << "pmax2_1" << "," << "pmax2_2" << "," << "pmax2_3" << "," << "pmax2_4" << "," 
     << "pmax3_0" << "," << "pmax3_1" << "," << "pmax3_2" << "," << "pmax3_3" << "," << "pmax3_4" << "," 
     << "pmax4_0" << "," << "pmax4_1" << "," << "pmax4_2" << "," << "pmax4_3" << "," << "pmax4_4" << "," 
     << "pmax5_0" << "," << "pmax5_1" << "," << "pmax5_2" << "," << "pmax5_3" << "," << "pmax5_4" << "," 
     << "pmax6_0" << "," << "pmax6_1" << "," << "pmax6_2" << "," << "pmax6_3" << "," << "pmax6_4" << "," 
     << "pmax7_0" << "," << "pmax7_1" << "," << "pmax7_2" << "," << "pmax7_3" << "," << "pmax7_4" << ","
     << "pmax8_0" << "," << "pmax8_1" << "," << "pmax8_2" << "," << "pmax8_3" << "," << "pmax8_4" << "," 
     << "pmax9_0" << "," << "pmax9_1" << "," << "pmax9_2" << "," << "pmax9_3" << "," << "pmax9_4" << "," 
     << "pmax10_0" << "," << "pmax10_1" << "," << "pmax10_2" << "," << "pmax10_3" << "," << "pmax10_4" << "," 
     << "pmax11_0" << "," << "pmax11_1" << "," << "pmax11_2" << "," << "pmax11_3" << "," << "pmax11_4" << "," 
     << "pmax12_0" << "," << "pmax12_1" << "," << "pmax12_2" << "," << "pmax12_3" << "," << "pmax12_4" << "," 
     << "pmax13_0" << "," << "pmax13_1" << "," << "pmax13_2" << "," << "pmax13_3" << "," << "pmax13_4" << "," 
     << "pmax14_0" << "," << "pmax14_1" << "," << "pmax14_2" << "," << "pmax14_3" << "," << "pmax14_4" << "," 
     << "pmax15_0" << "," << "pmax15_1" << "," << "pmax15_2" << "," << "pmax15_3" << "," << "pmax15_4" << "\n";

   Long64_t nentries = fChain->GetEntriesFast();

   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry); 


      int counter_laser=0 ;
      float max_laser=-100000. ;

      for( int i=0; i<1000; i++ ){

        if( w->at(16).at(i)>max_laser ){

          max_laser = w->at(16).at(i) ;
          counter_laser = i ;

        }

      }
      
    
      if(     (pmax->at(1)>=0 && pmax->at(1)<300) && (pmax->at(2)>=0 && pmax->at(2)<300) 
           && (pmax->at(3)>=0 && pmax->at(3)<300) && (pmax->at(4)>=0 && pmax->at(4)<300) && (pmax->at(5)>=0 && pmax->at(5)<300) 
           && (pmax->at(6)>=0 && pmax->at(6)<300) && (pmax->at(8)>=0 && pmax->at(8)<300) 
           && (pmax->at(9)>=0 && pmax->at(9)<300) && (pmax->at(10)>=0 && pmax->at(10)<300) && (pmax->at(11)>=0 && pmax->at(11)<300)
           && (pmax->at(13)>=0 && pmax->at(13)<300) && (pmax->at(14)>=0 && pmax->at(14)<300) 
           && (negpmax->at(1)<=0 && negpmax->at(1)>-300) && (negpmax->at(2)<=0 && negpmax->at(2)>-300) 
           && (negpmax->at(3)<=0 && negpmax->at(3)>-300) && (negpmax->at(4)<=0 && negpmax->at(4)>-300) && (negpmax->at(5)<=0 && negpmax->at(5)>-300) 
           && (negpmax->at(6)<=0 && negpmax->at(6)>-300) && (negpmax->at(8)<=0 && negpmax->at(8)>-300) 
           && (negpmax->at(9)<=0 && negpmax->at(9)>-300) && (negpmax->at(10)<=0 && negpmax->at(10)>-300) && (negpmax->at(11)<=0 && negpmax->at(11)>-300)
           && (negpmax->at(13)<=0 && negpmax->at(13)>-300) && (negpmax->at(14)<=0 && negpmax->at(14)>-300) ){

      

      f<< x_pos << "," << y_pos << ","  
      << 0 << "," << 0 << "," << 0 << "," << 0 << "," << 0 << ","
      << w->at(1).at(counter_laser-2) << "," << w->at(1).at(counter_laser-1) << "," << w->at(1).at(counter_laser) << "," << w->at(1).at(counter_laser+1) << "," << w->at(1).at(counter_laser+2) << ","
      << w->at(2).at(counter_laser-2) << "," << w->at(2).at(counter_laser-1) << "," << w->at(2).at(counter_laser) << "," << w->at(2).at(counter_laser+1) << "," << w->at(2).at(counter_laser+2) << ","
      << w->at(3).at(counter_laser-2) << "," << w->at(3).at(counter_laser-1) << "," << w->at(3).at(counter_laser) << "," << w->at(3).at(counter_laser+1) << "," << w->at(3).at(counter_laser+2) << ","
      << w->at(4).at(counter_laser-2) << "," << w->at(4).at(counter_laser-1) << "," << w->at(4).at(counter_laser) << "," << w->at(4).at(counter_laser+1) << "," << w->at(4).at(counter_laser+2) << ","
      << w->at(5).at(counter_laser-2) << "," << w->at(5).at(counter_laser-1) << "," << w->at(5).at(counter_laser) << "," << w->at(5).at(counter_laser+1) << "," << w->at(5).at(counter_laser+2) << ","
      << w->at(6).at(counter_laser-2) << "," << w->at(6).at(counter_laser-1) << "," << w->at(6).at(counter_laser) << "," << w->at(6).at(counter_laser+1) << "," << w->at(6).at(counter_laser+2) << ","
      << 0 << "," << 0 << "," << 0 << "," << 0 << "," << 0 << ","
      << w->at(8).at(counter_laser-2) << "," << w->at(8).at(counter_laser-1) << "," << w->at(8).at(counter_laser) << "," << w->at(8).at(counter_laser+1) << "," << w->at(8).at(counter_laser+2) << ","
      << w->at(9).at(counter_laser-2) << "," << w->at(9).at(counter_laser-1) << "," << w->at(9).at(counter_laser) << "," << w->at(9).at(counter_laser+1) << "," << w->at(9).at(counter_laser+2) << ","
      << w->at(10).at(counter_laser-2) << "," << w->at(10).at(counter_laser-1) << "," << w->at(10).at(counter_laser) << "," << w->at(10).at(counter_laser+1) << "," << w->at(10).at(counter_laser+2) << ","
      << w->at(11).at(counter_laser-2) << "," << w->at(11).at(counter_laser-1) << "," << w->at(11).at(counter_laser) << "," << w->at(11).at(counter_laser+1) << "," << w->at(11).at(counter_laser+2) << ","
      << 0 << "," << 0 << "," << 0 << "," << 0 << "," << 0 << ","
      << w->at(13).at(counter_laser-2) << "," << w->at(13).at(counter_laser-1) << "," << w->at(13).at(counter_laser) << "," << w->at(13).at(counter_laser+1) << "," << w->at(13).at(counter_laser+2) << ","
      << w->at(14).at(counter_laser-2) << "," << w->at(14).at(counter_laser-1) << "," << w->at(14).at(counter_laser) << "," << w->at(14).at(counter_laser+1) << "," << w->at(14).at(counter_laser+2) << "," 
      << 0 << "," << 0 << "," << 0 << "," << 0 << "," << 0 << "\n";


     }


      nbytes += nb;
   }

   f.close();
}
