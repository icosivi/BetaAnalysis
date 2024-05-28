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
   f.open("/home/daq/hdd8TB/RSD2/csv/stats_boxes15-RSD2_Run26_x_0_500_y_0_500_100waveforms_240V_10umstep_49.6pc_W3.csv");

   /*f << "x" << "," << "y" << ","
     << "pmax0" << "," << "pmax1" << "," << "pmax2" << "," << "pmax3" << "," << "pmax4" << ","
     << "negpmax0" << "," << "negpmax1" << "," << "negpmax2" << "," << "negpmax3" << "," << "negpmax4" << "\n";*/

   f << "x" << "," << "y" << ","
     << "pmax0" << "," << "pmax1" << "," << "pmax2" << "," << "pmax3" << "," << "pmax4" << "," << "pmax5" << "," << "pmax6" << "," << "pmax7" << ","
     << "pmax8" << "," << "pmax9" << "," << "pmax10" << "," << "pmax11" << "," << "pmax12" << ","
     << "negpmax0" << "," << "negpmax1" << "," << "negpmax2" << "," << "negpmax3" << "," << "negpmax4" << "," << "negpmax5" << "," << "negpmax6" << "," << "negpmax7" << ","
     << "negpmax8" << "," << "negpmax9" << "," << "negpmax10" << "," << "negpmax11" << "negpmax12" << "\n";

   /*f << "x" << "," << "y" << ","
     << "pmax0" << "," << "pmax1" << "," << "pmax2" << "," << "pmax3" << "," << "pmax4" << "," << "pmax5" << "," << "pmax6" << "," << "pmax7" << ","
     << "pmax8" << "," << "pmax9" << "," << "pmax10" << "," << "pmax11" << "," << "pmax12" << "," << "pmax13" << "," << "pmax14" << "," << "pmax15" << ","
     << "negpmax0" << "," << "negpmax1" << "," << "negpmax2" << "," << "negpmax3" << "," << "negpmax4" << "," << "negpmax5" << "," << "negpmax6" << "," << "negpmax7" << ","
     << "negpmax8" << "," << "negpmax9" << "," << "negpmax10" << "," << "negpmax11" << "," << "negpmax12" << "," << "negpmax13" << "," << "negpmax14" << "," << "negpmax15" << "\n";*/

   Long64_t nentries = fChain->GetEntriesFast();

   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry); 


      /*if(     (pmax->at(0)>=0 && pmax->at(0)<300) && (pmax->at(1)>=0 && pmax->at(1)<300) && (pmax->at(2)>=0 && pmax->at(2)<300) 
           && (pmax->at(3)>=0 && pmax->at(3)<300) && (pmax->at(4)>=0 && pmax->at(4)<300) && (pmax->at(5)>=0 && pmax->at(5)<300) 
           && (pmax->at(6)>=0 && pmax->at(6)<300) && (pmax->at(7)>=0 && pmax->at(7)<300) && (pmax->at(8)>=0 && pmax->at(8)<300) 
           && (pmax->at(9)>=0 && pmax->at(9)<300) && (pmax->at(10)>=0 && pmax->at(10)<300) && (pmax->at(11)>=0 && pmax->at(11)<300) 
           && (pmax->at(12)>=0 && pmax->at(12)<300) && (pmax->at(13)>=0 && pmax->at(13)<300) && (pmax->at(14)>=0 && pmax->at(14)<300)
           && (pmax->at(15)>=0 && pmax->at(15)<300) 
           && (negpmax->at(0)<=0 && negpmax->at(0)>-300) && (negpmax->at(1)<=0 && negpmax->at(1)>-300) && (negpmax->at(2)<=0 && negpmax->at(2)>-300) 
           && (negpmax->at(3)<=0 && negpmax->at(3)>-300) && (negpmax->at(4)<=0 && negpmax->at(4)>-300) && (negpmax->at(5)<=0 && negpmax->at(5)>-300) 
           && (negpmax->at(6)<=0 && negpmax->at(6)>-300) && (negpmax->at(7)<=0 && negpmax->at(7)>-300) && (negpmax->at(8)<=0 && negpmax->at(8)>-300) 
           && (negpmax->at(9)<=0 && negpmax->at(9)>-300) && (negpmax->at(10)<=0 && negpmax->at(10)>-300) && (negpmax->at(11)<=0 && negpmax->at(11)>-300)
           && (negpmax->at(12)<=0 && negpmax->at(12)>-300) && (negpmax->at(13)<=0 && negpmax->at(13)>-300) && (negpmax->at(14)<=0 && negpmax->at(14)>-300)
           && (negpmax->at(15)<=0 && negpmax->at(15)>-300) ){ */ 

      if(     (pmax->at(0)>=0 && pmax->at(0)<300) && (pmax->at(1)>=0 && pmax->at(1)<300) && (pmax->at(2)>=0 && pmax->at(2)<300) 
           && (pmax->at(3)>=0 && pmax->at(3)<300) && (pmax->at(4)>=0 && pmax->at(4)<300) && (pmax->at(5)>=0 && pmax->at(5)<300) 
           && (pmax->at(6)>=0 && pmax->at(6)<300) && (pmax->at(7)>=0 && pmax->at(7)<300) && (pmax->at(8)>=0 && pmax->at(8)<300) 
           && (pmax->at(9)>=0 && pmax->at(9)<300) && (pmax->at(10)>=0 && pmax->at(10)<300) && (pmax->at(11)>=0 && pmax->at(11)<300) && (pmax->at(12)>=0 && pmax->at(12)<300)
           && (negpmax->at(0)<=0 && negpmax->at(0)>-300) && (negpmax->at(1)<=0 && negpmax->at(1)>-300) && (negpmax->at(2)<=0 && negpmax->at(2)>-300) 
           && (negpmax->at(3)<=0 && negpmax->at(3)>-300) && (negpmax->at(4)<=0 && negpmax->at(4)>-300) && (negpmax->at(5)<=0 && negpmax->at(5)>-300) 
           && (negpmax->at(6)<=0 && negpmax->at(6)>-300) && (negpmax->at(7)<=0 && negpmax->at(7)>-300) && (negpmax->at(8)<=0 && negpmax->at(8)>-300) 
           && (negpmax->at(9)<=0 && negpmax->at(9)>-300) && (negpmax->at(10)<=0 && negpmax->at(10)>-300) && (negpmax->at(11)<=0 && negpmax->at(11)>-300)
           && (negpmax->at(12)<=0 && negpmax->at(12)>-300) ){ 

    /*if(       (pmax->at(0)>=0 && pmax->at(0)<300) && (pmax->at(1)>=0 && pmax->at(1)<300) && (pmax->at(2)>=0 && pmax->at(2)<300) 
           && (pmax->at(3)>=0 && pmax->at(3)<300) && (pmax->at(4)>=0 && pmax->at(4)<300) 
           && (negpmax->at(0)<=0 && negpmax->at(0)>-300) && (negpmax->at(1)<=0 && negpmax->at(1)>-300) && (negpmax->at(2)<=0 && negpmax->at(2)>-300) 
           && (negpmax->at(3)<=0 && negpmax->at(3)>-300) && (negpmax->at(4)<=0 && negpmax->at(4)>-300) ){*/

      /*f<< x_pos << "," << y_pos << ","  
      << pmax->at(0) << "," << pmax->at(1) << "," << pmax->at(2) << "," << pmax->at(3) << "," << pmax->at(4) << "," 
      << negpmax->at(0) << "," << negpmax->at(1) << "," << negpmax->at(2) << "," << negpmax->at(3) << ","  << negpmax->at(4) << "\n";*/

      f<< x_pos << "," << y_pos << ","  
      << pmax->at(0) << "," << pmax->at(1) << "," << pmax->at(2) << "," << pmax->at(3) << ","
      << pmax->at(4) << "," << pmax->at(5) << "," << pmax->at(6) << "," << pmax->at(7) << ","
      << pmax->at(8) << "," << pmax->at(9) << "," << pmax->at(10) << "," << pmax->at(11) << ","  << pmax->at(12) << ","
      << negpmax->at(0) << "," << negpmax->at(1) << "," << negpmax->at(2) << "," << negpmax->at(3) << ","
      << negpmax->at(4) << "," << negpmax->at(5) << "," << negpmax->at(6) << "," << negpmax->at(7) << ","
      << negpmax->at(8) << "," << negpmax->at(9) << "," << negpmax->at(10) << "," << negpmax->at(11) << "," << negpmax->at(12)
      << "\n";

      /*f<< x_pos << "," << y_pos << ","  
      << pmax->at(0) << "," << pmax->at(1) << "," << pmax->at(2) << "," << pmax->at(3) << "," 
      << pmax->at(4) << "," << pmax->at(5) << "," << pmax->at(6) << "," << pmax->at(7) << ","  
      << pmax->at(8) << "," << pmax->at(9) << "," << pmax->at(10) << "," << pmax->at(11) << "," 
      << pmax->at(12) << "," << pmax->at(13) << "," << pmax->at(14) << "," << pmax->at(15) << ","
      << negpmax->at(0) << "," << negpmax->at(1) << "," << negpmax->at(2) << "," << negpmax->at(3) << "," 
      << negpmax->at(4) << "," << negpmax->at(5) << "," << negpmax->at(6) << "," << negpmax->at(7) << ","  
      << negpmax->at(8) << "," << negpmax->at(9) << "," << negpmax->at(10) << "," << negpmax->at(11) << "," 
      << negpmax->at(12) << "," << negpmax->at(13) << "," << negpmax->at(14) << "," << negpmax->at(15) << "\n";*/

     }


      nbytes += nb;
   }

   f.close();
}


/*

<< area->at(0) << "," << area->at(1) << "," << area->at(2) << "," << area->at(3) << "," 
      << area->at(4) << "," << area->at(5) << "," << area->at(6) << "," << area->at(7) << ","
      << dvdt->at(0) << "," << dvdt->at(1) << "," << dvdt->at(2) << "," << dvdt->at(3) << "," 
      << dvdt->at(4) << "," << dvdt->at(5) << "," << dvdt->at(6) << "," << dvdt->at(7) << ","
      << width->at(0).at(1) << "," << width->at(1).at(1) << "," << width->at(2).at(1) << "," << width->at(3).at(1) << "," 
      << width->at(4).at(1) << "," << width->at(5).at(1) << "," << width->at(6).at(1) << "," << width->at(7).at(1) << "," 


*/

/*<< "area0" << "," << "area1" << "," << "area2" << "," << "area3" << "," << "area4" << "," << "area5" << "," << "area6" << "," << "area7" << ","
     << "area8" << "," << "area9" << "," << "area10" << "," << "area11" << "," << "area12" << "," << "area13" << "," << "area14" << "," << "area15" << ","
     << "cfd0" << "," << "cfd1" << "," << "cfd2" << "," << "cfd3" << "," << "cfd4" << "," << "cfd5" << "," << "cfd6" << "," << "cfd7" << ","
     << "cfd8" << "," << "cfd9" << "," << "cfd10" << "," << "cfd11" << "," << "cfd12" << "," << "cfd13" << "," << "cfd14" << "," << "cfd15" << ","
     << "tmax0" << "," << "tmax1" << "," << "tmax2" << "," << "tmax3" << "," << "tmax4" << "," << "tmax5" << "," << "tmax6" << "," << "tmax7" << ","
     << "tmax8" << "," << "tmax9" << "," << "tmax10" << "," << "tmax11" << "," << "tmax12" << "," << "tmax13" << "," << "tmax14" << "," << "tmax15" << */



/*<< area_new->at(0) << "," << area_new->at(1) << "," << area_new->at(2) << "," << area_new->at(3) << "," 
      << area_new->at(4) << "," << area_new->at(5) << "," << area_new->at(6) << "," << area_new->at(7) << ","  
      << area_new->at(8) << "," << area_new->at(9) << "," << area_new->at(10) << "," << area_new->at(11) << "," 
      << area_new->at(12) << "," << area_new->at(13) << "," << area_new->at(14) << "," << area_new->at(15) << ","
      << cfd->at(0).at(1) << "," << cfd->at(1).at(1) << "," << cfd->at(2).at(1) << "," << cfd->at(3).at(1) << "," 
      << cfd->at(4).at(1) << "," << cfd->at(5).at(1) << "," << cfd->at(6).at(1) << "," << cfd->at(7).at(1) << ","
      << cfd->at(8).at(1) << "," << cfd->at(9).at(1) << "," << cfd->at(10).at(1) << "," << cfd->at(11).at(1) << "," 
      << cfd->at(12).at(1) << "," << cfd->at(13).at(1) << "," << cfd->at(14).at(1) << "," << cfd->at(15).at(1) << ","
      << tmax->at(0) << "," << tmax->at(1) << "," << tmax->at(2) << "," << tmax->at(3) << "," 
      << tmax->at(4) << "," << tmax->at(5) << "," << tmax->at(6) << "," << tmax->at(7) << ","
      << tmax->at(8) << "," << tmax->at(9) << "," << tmax->at(10) << "," << tmax->at(11) << "," 
      << tmax->at(12) << "," << tmax->at(13) << "," << tmax->at(14) << "," << tmax->at(15) << */
