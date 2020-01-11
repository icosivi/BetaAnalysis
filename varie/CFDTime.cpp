//////////////////////////////////////
////                                //
////  CFD Time implementation       //
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


//==============================================================================


/*==============================================================================
Find the time at the given CFD percentage.
  param npoints     := number of sampled points in a waveform.
  param fraction    := CFD percentage
  param w           := waveform
  param t           := time trace of the waveform
  param pmax_holder := the Pmax tuple that contains the index

  return the time at the given CFD
==============================================================================*/
double Rising_Edge_CFD_Time(
  const double fraction,
  std::vector<double> voltageVec,
  std::vector<double> timeVec,
  const std::pair<double,unsigned int> Pmax
)
{
    // function to calculate index of constant fraction - not truly a constant fraction discriminator

    double pmax = Pmax.first;
    unsigned int imax = Pmax.second;

    double time_fraction = 0.0;
    unsigned int ifraction = 0;

    bool failure = true;

    for( int j = imax; j>-1; j-- )
    {
      if( voltageVec.at(j) <= pmax*fraction/100.0)
      {
        ifraction     = j;              //find index of first point before constant fraction of pulse
        time_fraction = timeVec.at(j);
        failure = false;
        break;
      }
    }
    if(ifraction == voltageVec.size()-1)ifraction--;

    if( failure ) time_fraction = timeVec.at(0);
    else time_fraction = time_fraction + (timeVec.at(ifraction+1) - timeVec.at(ifraction))* (pmax*fraction/100.0 - voltageVec.at(ifraction)) /(voltageVec.at(ifraction+1) - voltageVec.at(ifraction));

    return time_fraction;
}


//===========================================================================================================================================

double Falling_Edge_CFD_Time(
  const int           fraction,
  std::vector<double> voltageVec,
  std::vector<double> timeVec,
  const std::pair<double,unsigned int> Pmax
)
{
    // function to calculate index of constant fraction - not truly a constant fraction discriminator

    double pmax = 0.0, time_fraction = 0.0;

    pmax = Pmax.first;

    const unsigned int imax = Pmax.second;

    int ifraction = 0;
    //std::cout<<"Start loop\n";
    for( std::size_t j = imax, npoints = voltageVec.size(); j<npoints; j++)
    {
        if( voltageVec.at(j) <= pmax*double(fraction)/100.0 )
        {
            ifraction = j;              //find index of first point before constant fraction of pulse
            time_fraction = timeVec.at(j);
            break;
        }
    }
    if( ifraction == 0 ){ifraction = voltageVec.size()-1; }
    //std::cout<<"finish loop\n";
    time_fraction = time_fraction + (timeVec.at(ifraction-1) - timeVec.at(ifraction))* (pmax*double(fraction)/100.0 - voltageVec.at(ifraction)) /(voltageVec.at(ifraction-1) - voltageVec.at(ifraction));
    //time_fraction = time_fraction + (timeVec[ifraction-1] - timeVec[ifraction])* (pmax*double(fraction)/100.0 - voltageVec[ifraction]) /(voltageVec[ifraction-1] - voltageVec[ifraction]);
    //std::cout<<"before return\n";
    return time_fraction;
}
