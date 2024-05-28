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
   f.open("/home/daq/hdd8TB/RSD2/csv/time_croci1300-RSD2_Run23_x_0_1500_y_0_1500_100waveforms_370V_10umstep_49.3pc_W15_2.csv");

   f << "x" << "," << "y" << ","<< "t" << ","
     << "pmax0" << "," << "pmax1" << "," << "pmax2" << "," << "pmax3" << "," << "pmax4" << ","
     << "negpmax0" << "," << "negpmax1" << "," << "negpmax2" << "," << "negpmax3" << "," << "negpmax4" << "," 
     << "tmax0" << "," << "tmax1" << "," << "tmax2" << "," << "tmax3" << "," << "tmax4" << "\n";

   /*f << "x" << "," << "y" << ","
     << "pmax0" << "," << "pmax1" << "," << "pmax2" << "," << "pmax3" << "," << "pmax4" << "," << "pmax5" << "," << "pmax6" << "," << "pmax7" << ","
     << "pmax8" << "," << "pmax9" << "," << "pmax10" << "," << "pmax11" << "," << "pmax12" << ","
     << "negpmax0" << "," << "negpmax1" << "," << "negpmax2" << "," << "negpmax3" << "," << "negpmax4" << "," << "negpmax5" << "," << "negpmax6" << "," << "negpmax7" << ","
:q     << "negpmax8" << "," << "negpmax9" << "," << "negpmax10" << "," << "negpmax11" << "\n";*/

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
		cout<<x_pos<<endl;

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
           && (negpmax->at(15)<=0 && negpmax->at(15)>-300) ){  */
      /*if(     (pmax->at(0)>=0 && pmax->at(0)<300) && (pmax->at(1)>=0 && pmax->at(1)<300) && (pmax->at(2)>=0 && pmax->at(2)<300) 
           && (pmax->at(3)>=0 && pmax->at(3)<300) && (pmax->at(4)>=0 && pmax->at(4)<300) 
           && (negpmax->at(0)<=0 && negpmax->at(0)>-300) && (negpmax->at(1)<=0 && negpmax->at(1)>-300) && (negpmax->at(2)<=0 && negpmax->at(2)>-300) 
           && (negpmax->at(3)<=0 && negpmax->at(3)>-300) && (negpmax->at(4)<=0 && negpmax->at(4)>-300) 
            ){*/

        // && (pmax->at(1)+pmax->at(2)+pmax->at(3)+pmax->at(4))>100

          if(x_pos>750 && y_pos>750) cout<<"pippo"<<endl;
          //float pmax_v[5] = {0,0,0,0,0};

          //for(int i=0; i<5;i++) if( pmax->at(i)>10 ) pmax_v[i] =  pmax->at(i);

          float t_reco = ( (pmax->at(1)*pmax->at(1)*(tmax->at(1)-0.14))+(pmax->at(2)*pmax->at(2)*tmax->at(2))
                  +(pmax->at(3)*pmax->at(3)*(tmax->at(3)-0.47))+(pmax->at(4)*pmax->at(4)*tmax->at(4)) )/
                  ( pmax->at(1)*pmax->at(1)+pmax->at(2)*pmax->at(2)+pmax->at(3)*pmax->at(3)+pmax->at(4)*pmax->at(4) ) ;
                  
          f<< x_pos << "," << y_pos << "," << t_reco-tmax->at(5) << ","  
          << pmax->at(0) << "," << pmax->at(1) << "," << pmax->at(2) << "," << pmax->at(3) << "," << pmax->at(4) << "," 
          << negpmax->at(0) << "," << negpmax->at(1) << "," << negpmax->at(2) << "," << negpmax->at(3) << ","  << negpmax->at(4) << "," 
          << tmax->at(0) << "," << tmax->at(1)-0.14 << "," << tmax->at(2) << "," << tmax->at(3)-0.47 << "," << tmax->at(4) << "\n";


          /*f<< t_reco-tmax->at(5) << ","  
          << pmax->at(0) << "," << pmax->at(1) << "," << pmax->at(2) << "," << pmax->at(3) << "," << pmax->at(4) << "," 
          << negpmax->at(0) << "," << negpmax->at(1) << "," << negpmax->at(2) << "," << negpmax->at(3) << ","  << negpmax->at(4) << "," 
          << tmax->at(0) << "," << tmax->at(1)-0.14 << "," << tmax->at(2) << "," << tmax->at(3)-0.47 << "," << tmax->at(4) << "\n";*/

       /* f<< x_pos << "," << y_pos << ","  
      << pmax->at(0) << "," << pmax->at(1) << "," << pmax->at(2) << "," << pmax->at(3) << "," 
      << pmax->at(4) << "," << pmax->at(5) << "," << pmax->at(6) << "," << pmax->at(7) << ","  
      << pmax->at(8) << "," << pmax->at(9) << "," << pmax->at(10) << "," << pmax->at(11) << "," 
      << pmax->at(12) << ","
      << negpmax->at(0) << "," << negpmax->at(1) << "," << negpmax->at(2) << "," << negpmax->at(3) << "," 
      << negpmax->at(4) << "," << negpmax->at(5) << "," << negpmax->at(6) << "," << negpmax->at(7) << ","  
      << negpmax->at(8) << "," << negpmax->at(9) << "," << negpmax->at(10) << "," << negpmax->at(11) << "," 
      << negpmax->at(12) << "\n";*/

      /*f<< x_pos << "," << y_pos << ","  
      << pmax->at(0) << "," << pmax->at(1) << "," << pmax->at(2) << "," << pmax->at(3) << "," 
      << pmax->at(4) << "," << pmax->at(5) << "," << pmax->at(6) << "," << pmax->at(7) << ","  
      << pmax->at(8) << "," << pmax->at(9) << "," << pmax->at(10) << "," << pmax->at(11) << "," 
      << pmax->at(12) << "," << pmax->at(13) << "," << pmax->at(14) << "," << pmax->at(15) << ","
      << negpmax->at(0) << "," << negpmax->at(1) << "," << negpmax->at(2) << "," << negpmax->at(3) << "," 
      << negpmax->at(4) << "," << negpmax->at(5) << "," << negpmax->at(6) << "," << negpmax->at(7) << ","  
      << negpmax->at(8) << "," << negpmax->at(9) << "," << negpmax->at(10) << "," << negpmax->at(11) << "," 
      << negpmax->at(12) << "," << negpmax->at(13) << "," << negpmax->at(14) << "," << negpmax->at(15) << "\n";*/

     //}


      nbytes += nb;
   }

   f.close();
}
