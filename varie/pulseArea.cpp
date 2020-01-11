//////////////////////////////////////
////                                //
////  Pulse Area implementation     //
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
//------ROOT----------------//
#include <TH1.h>
#include <TF1.h>
#include <TGraph.h>
#include <TThread.h>
#include <TROOT.h>
#include <TStyle.h>
#include <TImage.h>
#include <TCanvas.h>
#include <TTreeReader.h>
#include "TTreeReaderValue.h"
#include "TTreeReaderArray.h"

//==============================================================================


/*==============================================================================
Finding the pulse area of a signal. Stop when it crosses zero (legacy)
  std::vector<double> voltageVec := waveform
  std::vector<double> timeVec := time trace of the waveform
  const std::pair<double,unsigned int> Pmax := tuple that contains pmax and its index

  return: (double) value of the pulse area
==============================================================================*/
double Find_Pulse_Area(
  std::vector<double> voltageVec,
  std::vector<double> timeVec,
  const std::pair<double,unsigned int> Pmax
)
{
    double pulse_area = 0.0;
    const double time_difference = timeVec.at(1) - timeVec.at(0);

    const unsigned int imax = Pmax.second;
    unsigned int istart = 0;
    unsigned int iend;
    std::size_t npoints = voltageVec.size();

    for( int j = imax; j>-1; j-- ) // find index of start of pulse
    {
      if(voltageVec.at(j) <= 0) //stop after crossing zero
      {
        istart = j;
        break;
      }
    }
    for( unsigned int j = imax; j< npoints; j++ ) // find index of end of pulse
    {
      if(voltageVec.at(j) <= 0)
      {
        iend = j;
        break;
      }
      if( j == npoints-1 )
      {
        iend = j;
      }
    }

    for( unsigned int j = istart; j < iend; j++ )
    {
      pulse_area = pulse_area + voltageVec.at(j)/1000.0;
    }

    pulse_area = pulse_area * time_difference /1.0E12;

    return pulse_area; // collected pulse area, assuming voltage is in volts, time is in seconds
}


/*==============================================================================
Finding the undershoot area of a signal. Stop when it crosses zero
  std::vector<double> voltageVec := waveform
  std::vector<double> timeVec := time trace of the waveform
  const std::pair<double,unsigned int> Pmax := tuple that contains pmax and its index

  return: (double) value of the pulse area
==============================================================================*/
double Find_Udershoot_Area(
  std::vector<double> voltageVec,
  std::vector<double> timeVec,
  const std::pair<double, unsigned int> Pmax
)
{
    double undershoot_area = 0.0;
    const double time_difference = timeVec.at(1) - timeVec.at(0);

    const unsigned int imax = Pmax.second;
    unsigned int istart = 0;
    unsigned int iend;
    std::size_t npoints = voltageVec.size();

    for( unsigned int j = imax; j < npoints; j++ ) // find index of start of pulse
    {
      if( voltageVec.at(j) <= 0) //stop after crossing zero
      {
        istart = j;
        break;
      }
    }
    for( unsigned int j = istart; j< npoints; j++ ) // find index of end of pulse
    {
      if( voltageVec.at(j) >= 0)
      {
        iend = j;
        break;
      }
      if( j == npoints-1 )
      {
        iend = j;
      }
    }
    for( unsigned int j = istart; j < iend; j++ )
    {
      undershoot_area = undershoot_area + voltageVec.at(j)/1000.0;
    }

    undershoot_area = undershoot_area * time_difference /1.0E12;

    return undershoot_area; // collected undershoot area, assuming voltage is in volts, time is in seconds
}


/*==============================================================================
Finding the undershoot area of a signal with extended time window.
  std::vector<double> voltageVec := waveform
  std::vector<double> timeVec := time trace of the waveform
  const std::pair<double,unsigned int> Pmax := the Pmax tuple that contains the index
  double extendedIntRange := exteneded time for pulse integration, assuming in unit of ps.

  return: (double) undershoot area.
==============================================================================*/
double Find_Udershoot_Area(
  std::vector<double> voltageVec,
  std::vector<double> timeVec,
  const std::pair<double,unsigned int> Pmax,
  double extendedIntRange
)
{
    double undershoot_area = 0.0;
    const double time_difference = timeVec.at(1) - timeVec.at(0);

    const unsigned int imax = Pmax.second;
    unsigned int istart = 0;
    unsigned int iend;
    std::size_t npoints = voltageVec.size();

    double timeOfMaximum = timeVec.at(imax);

    for( unsigned int j = imax; j < npoints; j++ ) // find index of start of pulse
    {
      if( voltageVec.at(j) <= 0 ) //stop after crossing zero
      {
        istart = j;
        break;
      }
    }
    for( unsigned int j = istart; j< npoints; j++ )
    {
      if( j == npoints - 1 ) iend = j;
      else
      {
        if( (timeVec.at(j)-timeOfMaximum) >= extendedIntRange )
        {
          iend = j;
          break;
        }
      }
    }

    for( unsigned int j = istart; j < iend; j++ )
    {
      undershoot_area = undershoot_area + voltageVec.at(j)/1000.0;
    }

    undershoot_area = undershoot_area * time_difference /1.0E12;

    return undershoot_area; // collected undershoot area, assuming voltage is in volts, time is in seconds
}


/*==============================================================================
find the pulse area with interpolated edge.

==============================================================================*/
double Pulse_Area_With_Linear_Interpolate_Edge(
  std::vector<double> voltageVec,
  std::vector<double> timeVec,
  const std::pair<double,unsigned int> Pmax,
  std::string integration_option
)
{
  double pulse_area = 0.0;
  const double time_difference = timeVec.at(1) - timeVec.at(0);

  unsigned int imax = Pmax.second;
  std::size_t npoints = voltageVec.size();

  if( imax == npoints-1 ) imax = imax - 1;//preventing out of range.

  const double _20pmax = Pmax.first * 0.20;
  const double _10pmax = Pmax.first * 0.10;
  double _20pmax_time = 0.0;
  double _10pmax_time = 0.0;
  unsigned int istart = 0;
  unsigned int iend;
  double start_time = 0.0;
  bool found_20pmax = false;
  bool found_10pmax = false;

  for( int j = imax; j>-1; j-- ) // find index of start of pulse
  {
    if( !found_20pmax )
    {
      if( voltageVec.at(j) <= _20pmax ) //stop after crossing zero
      {
        _20pmax_time = xlinearInter( timeVec.at(j), voltageVec.at(j), timeVec.at(j+1), voltageVec.at(j+1), _20pmax );
        found_20pmax = true;
      }
    }

    if( !found_10pmax )
    {
      if( voltageVec.at(j) <= _10pmax ) //stop after crossing zero
      {
        _10pmax_time = xlinearInter( timeVec.at(j), voltageVec.at(j), timeVec.at(j+1), voltageVec.at(j+1), _10pmax );
        found_10pmax = true;
      }
    }

    if( found_10pmax && found_20pmax ) break;
  }

  start_time = xlinearInter( _10pmax_time, _10pmax, _20pmax_time, _20pmax, 0.0 );

  for ( unsigned int j = imax; j< npoints; j++ ) // find index of end of pulse
  {
    if( voltageVec.at(j) <= 0)
    {
      iend = j;
      break;
    }
    if( j == npoints-1 )
    {
      iend = j;
    }
  }

  for ( unsigned int j = 0; j < iend; j++ )
  {
    if( timeVec.at(j) >= start_time )
    {
      istart = j;
      pulse_area = pulse_area + ((timeVec.at(j)-start_time)/1.0E12) * voltageVec.at(j)/1000.0;
      break;
    }
  }

  //===========Simpson's rule=======
  if( integration_option.compare("Simpson") == 0 )
  {
    std::vector<double> integration_y;
    for ( unsigned int j = istart; j < iend; j++ )
    {
      integration_y.push_back( voltageVec.at(j) );
    }
    for( std::size_t i = 0, max = integration_y.size(); i < max; i++ )
    {
      if( i == 0 ) pulse_area = pulse_area + (time_difference/3.0/1.0E12) * (integration_y.at(i)/1000.0);
      else if( i == integration_y.size()-1 ) pulse_area = pulse_area + (time_difference/3.0/1.0E12) * (integration_y.at(i)/1000.0);
      else if( i % 2 == 0 ) pulse_area = pulse_area + 2 * (time_difference/3.0/1.0E12) * (integration_y.at(i)/1000.0);
      else pulse_area = pulse_area + 4 * (time_difference/3.0/1.0E12) * (integration_y.at(i)/1000.0);
    }
  }

  //===========Rectrangluar=========
  else
  {
    for( unsigned int j = istart; j < iend; j++ )
    {
      pulse_area = pulse_area + (time_difference/1.0E12) * voltageVec.at(j)/1000.0;
    }
  }

  return pulse_area; // collected pulse area, assuming voltage is in volts, time is in seconds
}


/*==============================================================================
find the pulse area with interpolated edge, extend to fixed time window.

==============================================================================*/
double Pulse_Area_With_Linear_Interpolate_Edge(
  std::vector<double> voltageVec, //voltage
  std::vector<double> timeVec, //time
  const std::pair<double,unsigned int> Pmax,
  std::string integration_option,
  bool relativeTimeWindow,
  double stopTime
)
{
  double pulse_area = 0.0;
  const double time_difference = timeVec.at(1) - timeVec.at(0);

  unsigned int imax = Pmax.second;

  double timeOfMaximum = timeVec.at(imax);
  std::size_t npoints = voltageVec.size();

  if( imax == npoints-1 ) imax = imax - 1;//preventing out of range.

  const double _20pmax = Pmax.first * 0.20;
  const double _10pmax = Pmax.first * 0.10;
  double _20pmax_time = 0.0;
  double _10pmax_time = 0.0;
  unsigned int istart = 0;
  unsigned int iend;
  double start_time = 0.0;
  bool found_20pmax = false;
  bool found_10pmax = false;

  for( int j = imax; j>-1; j-- ) // find index of start of pulse
  {
    if( !found_20pmax )
    {
      if( voltageVec.at(j) <= _20pmax ) //stop after crossing zero
      {
        _20pmax_time = xlinearInter( timeVec.at(j), voltageVec.at(j), timeVec.at(j+1), voltageVec.at(j+1), _20pmax );
        found_20pmax = true;
      }
    }

    if( !found_10pmax )
    {
      if( voltageVec.at(j) <= _10pmax ) //stop after crossing zero
      {
        _10pmax_time = xlinearInter( timeVec.at(j), voltageVec.at(j), timeVec.at(j+1), voltageVec.at(j+1), _10pmax );
        found_10pmax = true;
      }
    }

    if( found_10pmax && found_20pmax ) break;
  }

  start_time = xlinearInter( _10pmax_time, _10pmax, _20pmax_time, _20pmax, 0.0 );


  for( unsigned int j = imax; j< npoints; j++ ) // find index of end of pulse
  {
    if( j == npoints - 1 ) iend = j;
    else
    {
      if( relativeTimeWindow )
      {
        if( (timeVec.at(j)-timeOfMaximum) >= stopTime )
        {
          iend = j;
          break;
        }
      }
      else
      {
        if( timeVec.at(j) >= stopTime )
        {
          iend = j;
          break;
        }
      }
    }
  }

  for( unsigned int j = 0; j < iend; j++ )
  {
    if( timeVec.at(j) >= start_time )
    {
      istart = j;
      pulse_area = pulse_area + ((timeVec.at(j)-start_time)/1.0E12) * voltageVec.at(j)/1000.0;
      break;
    }
  }

  //===========Simpson's rule=======
  if( integration_option.compare("Simpson") == 0 )
  {
    std::vector<double> integration_y;
    for ( unsigned int j = istart; j < iend; j++ )
    {
      integration_y.push_back( voltageVec.at(j) );
    }
    for( std::size_t i = 0, max = integration_y.size(); i < max; i++ )
    {
      if( i == 0 ) pulse_area = pulse_area + (time_difference/3.0/1.0E12) * (integration_y.at(i)/1000.0);
      else if( i == integration_y.size()-1 ) pulse_area = pulse_area + (time_difference/3.0/1.0E12) * (integration_y.at(i)/1000.0);
      else if( i % 2 == 0 ) pulse_area = pulse_area + 2 * (time_difference/3.0/1.0E12) * (integration_y.at(i)/1000.0);
      else pulse_area = pulse_area + 4 * (time_difference/3.0/1.0E12) * (integration_y.at(i)/1000.0);
    }
  }

  //===========Rectrangluar=========
  else
  {
    for( unsigned int j = istart; j < iend; j++ )
    {
      pulse_area = pulse_area + (time_difference/1.0E12) * voltageVec.at(j)/1000.0;
    }
  }

  return pulse_area; // collected pulse area, assuming voltage is in volts, time is in seconds
}


/*==============================================================================
find the pulse area with fixed time window size.
1ns before signal and 3ns after signal

==============================================================================*/
double Pulse_Integration_with_Fixed_Window_Size(
  std::vector<double> voltageVec,
  std::vector<double> timeVec,
  const std::pair<double,unsigned int> Pmax,
  std::string integration_option,
  const double t_beforeSignal,
  const double t_afterSignal
)
{
  double pulse_area = 0.0;
  const double time_difference = timeVec.at(1) - timeVec.at(0);
  double tRange[2] = {t_beforeSignal, t_afterSignal};

  unsigned int imax = Pmax.second;

  double timeOfMaximum = timeVec.at(imax);
  std::size_t npoints = voltageVec.size();

  if( imax == npoints-1 ) imax = imax - 1;//preventing out of range.

  const double _20pmax = Pmax.first * 0.20;
  const double _10pmax = Pmax.first * 0.10;
  double _20pmax_time = timeVec.at(0);
  double _10pmax_time = -timeVec.at(0);
  double start_time = -timeVec.at(0);
  bool found_20pmax = false;
  bool found_10pmax = false;

  for( int j = imax; j>-1; j-- ) // find index of start of pulse
  {
    if( !found_20pmax )
    {
      if( voltageVec.at(j) <= _20pmax ) //stop after crossing zero
      {
        _20pmax_time = xlinearInter( timeVec.at(j), voltageVec.at(j), timeVec.at(j+1), voltageVec.at(j+1), _20pmax );
        found_20pmax = true;
      }
    }

    if( !found_10pmax )
    {
      if( voltageVec.at(j) <= _10pmax ) //stop after crossing zero
      {
        _10pmax_time = xlinearInter( timeVec.at(j), voltageVec.at(j), timeVec.at(j+1), voltageVec.at(j+1), _10pmax );
        found_10pmax = true;
      }
    }

    if( found_10pmax && found_20pmax ) break;
  }

  start_time = xlinearInter( _10pmax_time, _10pmax, _20pmax_time, _20pmax, 0.0 );

  /*double timeOfFirstZeroCrossing;
  for ( unsigned int j = imax; j< npoints; j++ ) // find index of end of pulse
  {
    if( j == npoints - 1 ){ timeOfFirstZeroCrossing = timeVec.at(j); }
    else
    {
      if( voltageVec.at(j) <= 0.0 )
      {
        timeOfFirstZeroCrossing = timeVec.at(j);
        break;
      }
    }
  }*/

  std::vector<double> integration_voltage_vector;
  for( unsigned int i = 0; i < npoints; i++)
  {
    if( timeVec.at(i) >= (start_time-tRange[0]) && timeVec.at(i) <= timeOfMaximum + tRange[1])
    //if( t.at(i) >= (start_time) && t.at(i) <= 3000.0)
    {
      integration_voltage_vector.push_back( voltageVec.at(i) );
    }
  }

  //===========Simpson's rule=======
  if( integration_option.compare("Simpson") == 0 )
  {
    for( std::size_t i = 0, max = integration_voltage_vector.size(); i < max; i++ )
    {
      if( i == 0 ) pulse_area = pulse_area + (time_difference/3.0/1.0E12) * (integration_voltage_vector.at(i)/1000.0);
      else if( i == integration_voltage_vector.size()-1 ) pulse_area = pulse_area + (time_difference/3.0/1.0E12) * (integration_voltage_vector.at(i)/1000.0);
      else if( i % 2 == 0 ) pulse_area = pulse_area + 2 * (time_difference/3.0/1.0E12) * (integration_voltage_vector.at(i)/1000.0);
      else pulse_area = pulse_area + 4 * (time_difference/3.0/1.0E12) * (integration_voltage_vector.at(i)/1000.0);
    }
  }
  //================================
  //===========Rectrangluar=========
  else
  {
    for ( std::size_t j = 0, max = integration_voltage_vector.size(); j < max; j++ )
    {
      pulse_area = pulse_area + (time_difference/1.0E12) * integration_voltage_vector.at(j)/1000.0;
    }
  }
  //================================

  return pulse_area; // collected pulse area, assuming voltage is in volts, time is in seconds
}
