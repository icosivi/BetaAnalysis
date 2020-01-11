///////////////////////////////////////////
////                                     //
////  Time at threshoold implementation  //
////                                    //
/////////////////////////////////////////


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
#include <math.h>



/*==============================================================================
Find the time at where the voltage crosses the threshold level.

  param thresholdLevel := the threshold level
  param voltage        := the vector of voltage value of the waveform
  param time           := the vector of time value of the waveform
  param pmax_holder    := a std pair of the Pmax value and its index in the wave

  return : the time at the threshold.

/=============================================================================*/
double Find_Time_At_Threshold(
  const double                thresholdLevel,
  std::vector<double>         voltageVec,
  std::vector<double>         timeVec,
  const std::pair<double,unsigned int> Pmax
)
{
  double timeAtThreshold = 0.0, timeBelowThreshold = 0.0;

  unsigned int timeBelowThreshold_index = 0;

  unsigned int pmax_index = Pmax.second;
  double pmax = Pmax.first;
  std::size_t npoints = voltageVec.size();

  if( pmax_index == npoints-1 ) pmax_index = pmax_index - 1;//preventing out of range

    if( pmax < thresholdLevel ){ return 9999.0;}
    else
    {
      for( int i = pmax_index; i > -1; i-- )
      {
        if( voltageVec.at(i) <= thresholdLevel )
        {
          timeBelowThreshold_index = i;

          timeBelowThreshold       = timeVec.at(i);

          break;
        }
      }

      timeAtThreshold = xlinearInter( timeBelowThreshold, voltageVec.at(timeBelowThreshold_index), timeVec.at(timeBelowThreshold_index+1), voltageVec.at(timeBelowThreshold_index+1), thresholdLevel );

      return timeAtThreshold;
    }
}


/*==============================================================================
Find the time at where the voltage crosses the threshold level. (2)
Alternative implementation. A simpler version, it dose not require the Pmax.
Best for simple step function trigger.

  param npoints        := number of points in a waveform.
  param thresholdLevel := the threshold level
  param voltage        := the vector of voltage value of the waveform
  param time           := the vector of time value of the waveform

  return : the time at the threshold.

/=============================================================================*/
double Find_Time_At_Threshold(
  const double                thresholdLevel,
  std::vector<double>         voltage_v,
  std::vector<double>         time_v
)
{
  for( std::size_t i = 1, max = voltage_v.size(); i < max; i++ ) //starting at 1 to avoid boundary underflow problem.
  {
    if( voltage_v.at(i) >= thresholdLevel )
    {
      return xlinearInter( time_v.at(i-1), voltage_v.at(i-1), time_v.at(i), voltage_v.at(i), thresholdLevel );
    }
  }
  return time_v.at(0); //return the first point of the time vector if not finnding anything.
}



//==============================================================================
//==============================================================================
//find the time that across the threshold.

void Get_TimeAcrossThreshold(
  const double        thresholdLevel,
  std::vector<double> voltageVec,
  std::vector<double> timeVec,
  std::vector<double> &time_at_threshold_v,
  const unsigned int  expect_count
)
{
  bool rise_edge = true;
  bool fall_edge = false;

  for( std::size_t i = 1, npoints = voltageVec.size(); i < npoints; i++ )
  {
    if( time_at_threshold_v.size() == expect_count ) break;

    if( rise_edge )
    {
      if( voltageVec.at(i) >= thresholdLevel )
      {
        time_at_threshold_v.push_back( xlinearInter( timeVec.at(i-1), voltageVec.at(i-1), timeVec.at(i), voltageVec.at(i), thresholdLevel ) );
        rise_edge = false;
        fall_edge = true;
      }
    }
    else if( fall_edge )
    {
      if( voltageVec.at(i) <= thresholdLevel )
      {
        time_at_threshold_v.push_back( xlinearInter( timeVec.at(i-1), voltageVec.at(i-1), timeVec.at(i), voltageVec.at(i), thresholdLevel ) );
        rise_edge = true;
        fall_edge = false;
      }
    }
    else{}
  }

  if( time_at_threshold_v.size() < expect_count )
  {
    time_at_threshold_v.push_back( 1.0e6 );
  }
}
