//////////////////////////////////////
////                                //
////  Dvdt implementation           //
////                                //
//////////////////////////////////////


//==============================================================================
// Headers

#include "../include/Waveform_Analysis.hpp"
#include "../include/general.hpp"

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
//------ROOT----------------//
#include <TH1.h>
#include <TF1.h>
#include <TGraph.h>
#include <TThread.h>
#include <TROOT.h>
#include <TStyle.h>
#include <TImage.h>
#include <TCanvas.h>


//===========================================================================================================================================

double Find_Dvdt(
  const int           fraction,
  const int           ndif,
  std::vector<double> voltageVec,
  std::vector<double> timeVec,
  const std::pair<double,unsigned int> Pmax
)
{
    // function to dv/dt at a given constant fraction value.

    double time_difference = 0.0;
    double dvdt            = 0.0;

    int ifraction       = 0;

    time_difference = timeVec.at(1) - timeVec.at(0);

    double pmax = Pmax.first;
    unsigned int imax = Pmax.second;

    for( int j = imax; j>-1; j--)
    {
      if( voltageVec.at(j) <= pmax*double(fraction)/100)
      {
        ifraction = j;

        break;
      }
    }//find index of first point before constant fraction of pulse

    if(ifraction == voltageVec.size()-1)ifraction--;

    if(ndif == 0)
    {
      dvdt = (voltageVec.at(ifraction+1) - voltageVec.at(ifraction))/time_difference;
    }

    else
    {
      dvdt = (voltageVec.at(ifraction+ndif) - voltageVec.at(ifraction-ndif))/(time_difference*(ndif*2));
    }

    return dvdt;
}
