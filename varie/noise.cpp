//////////////////////////////////////
////                                //
////  Noise implementation          //
////                                //
//////////////////////////////////////


//==============================================================================
// Headers

#include "Waveform_Analysis.hpp"
#include "general.hpp"

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
#include <TTreeReader.h>
#include "TTreeReaderValue.h"
#include "TTreeReaderArray.h"

/*==============================================================================
Find the noise of a signal baseline.
  param npoints := number of sampled points in a waveform.
  param w       := waveform
  param inoise  := number of points for the noise calculation.

  return noise
==============================================================================*/
double Find_Noise(
  std::vector<double> voltageVec,
  const unsigned int inoise
)
{
  double rms = 0.0, mean = 0.0, var = 0.0;

  for( unsigned int j = 0; j < inoise; j++)
  {
    rms  += voltageVec.at(j)*voltageVec.at(j);
    mean += voltageVec.at(j);
  }

  mean = mean / inoise;
  rms  = rms/inoise;
  var  = rms - mean * mean;
  rms  = pow(var, 0.5);

  return rms;
}


/*==============================================================================
Find the noise of a signal baseline. (Alternate implementation of calculating noise)
  param w              := waveform
  param fractional_pts := fraction of points of the waveform for the noise calculation.

  return noise
==============================================================================*/
double Find_Noise2(
  std::vector<double> voltageVec,
  const double fractional_pts
)
{
  double rms = 0.0, mean = 0.0, var = 0.0;

  int inoise = fractional_pts * voltageVec.size();

  for( std::size_t j = 0, max = inoise; j < max; j++ )
  {
    rms  += voltageVec.at(j)*voltageVec.at(j);
    mean += voltageVec.at(j);
  }

  mean = mean/inoise;
  rms  = rms/inoise;
  var  = rms - mean * mean;
  rms  = pow(var, 0.5);

  return rms;
}


/*==============================================================================
Find the noise of a signal on the back baseline.
  param npoints    := number of sampled points in a waveform.
  param w          := waveform
  param t          := time trace of the waveform
  param start_time := staring time on the time trace for the noise calculation
  param end+time   := ending time on the time trace for the noise calculation

  return noise on the back baseline
==============================================================================*/
double Find_Noise_On_Back_Baseline(
  std::vector<double> voltageVec,
  std::vector<double> timeVec,
  double start_time,
  double end_time
)
{
  double rms = 0.0, mean = 0.0, var = 0.0;
  int counter = 0;
  std::size_t npoints = voltageVec.size();

  for( unsigned int j = 0; j < npoints; j++)
  {
    if( start_time <= timeVec.at(j) && timeVec.at(j) <= end_time )
    {
      rms  += voltageVec.at(j)*voltageVec.at(j);
      mean += voltageVec.at(j);
      counter++;
    }
  }

  mean = mean / counter;
  rms  = rms/counter;
  var  = rms - mean * mean;
  rms  = pow(var, 0.5);
  return rms;
}
