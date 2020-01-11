////////////////////////////////////////////
////                                      //
////  Multiple Signal Peak implementation //
////                                      //
///////////////////////////////////////////


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



/*==============================================================================
Finding Multiple Signal Peaks
  param w          := waveform
  param t          := time trace of the waveform
  param StartIndex := index of the waveform vector to start the next search

  return : Pmax value and its index in the waveform vector.
==============================================================================*/
std::pair <double, unsigned int> Find_Identical_Peak(
  std::vector<double> voltageVec,
  std::vector<double> timeVec,
  unsigned int StartIndex,
  bool limitSearchRegion,
  double min_search_range,
  double max_search_range
)
{
    double          pmax        = 0.0;
    unsigned int    pmax_index  = StartIndex;
    bool   StrangePeak = true;
    bool   FirstPoint  = true;
    std::size_t npoints = voltageVec.size();

    //std::cout << "StartIndex: " << StartIndex << std::endl;
    if(StartIndex)
    {
      if( limitSearchRegion )
      {
        for( unsigned int j = StartIndex; j < npoints; j++)
        {
          if( j == StartIndex )
          {
              pmax = voltageVec.at(j);
              pmax_index = j;
              FirstPoint = false;
          }
          if( min_search_range <= timeVec.at(j) && timeVec.at(j) <= max_search_range ) //zoom in to find the Pmax
          {
            if( j != StartIndex &&  voltageVec.at(j) == pmax )
            {
              pmax = voltageVec.at(j);
              pmax_index = j;
              StrangePeak = false;
            }
          }
        }
        if( StrangePeak && FirstPoint) //previous sick point will be sent to here.
        {
          pmax = 10000.0;
          pmax_index = npoints;
        }
      }
    else
    {
      for( unsigned int j = StartIndex; j < npoints; j++)
      {
        if(j == StartIndex)
        {
          pmax = voltageVec.at(j);
          pmax_index = j;
        }
        if(j != StartIndex && voltageVec.at(j) == pmax)
        {
            pmax = voltageVec.at(j);
            pmax_index = j;
        }
      }
    }
  }
  else
  {
    if( limitSearchRegion )
    {
      for( unsigned int j = 0; j < npoints; j++)
      {
        if( min_search_range <= timeVec.at(j) && timeVec.at(j) <= max_search_range ) //zoom in to find the Pmax
        {
          if( FirstPoint && voltageVec.at(j) > 0 )
          {
            pmax = voltageVec.at(j);
            pmax_index = j;
            FirstPoint = false;
            StrangePeak = false;
          }
          if( voltageVec.at(j) > pmax )
          {
            pmax = voltageVec.at(j);
            pmax_index = j;
          }
        }
      }
      if( StrangePeak ) //if it cannot find any positive point within the range, return this number
      {
        pmax = 10000.0;
        pmax_index = npoints;
      }
    }
    else
    {
      for( unsigned int j = 0; j < npoints; j++)
      {
        //std::cout << "here increment " << j << std::endl;
        if (j == 0)
        {
          pmax = 0.0;
          pmax_index = 0;
        }
        if (j != 0 && voltageVec.at(j) > pmax)
        {
          pmax = voltageVec.at(j);
          pmax_index = j;
        }
      }
    }
  }
    //std::cout << "pmax_index: " << pmax_index << std::endl;

    return std::make_pair(pmax, pmax_index);
}



//==============================================================================
//==============================================================================
//find all the pmax for multiple signals

void Get_PmaxTmax_Of_Multiple_Singal(
  const double        assist_threshold,
  std::vector<double> voltageVec,
  std::vector<double> timeVec,
  std::vector<double> &multiple_singal_pmax_v,
  std::vector<double> &multiple_singal_tmax_v,
  std::vector<int>    &indexing_v
)
{
  double pmax = 0.0;
  int pmax_index = 0;
  bool candidate_signal = false;
  bool noisy_event = true;
  std::size_t npoints = voltageVec.size();

  for( unsigned int i = 0; i < npoints; i++ )
  {
    if( !candidate_signal )
    {
      if( voltageVec.at(i) >= assist_threshold )
      {
        if( voltageVec.at(i) > pmax )
        {
          pmax = voltageVec.at(i);
          pmax_index = i;
          candidate_signal = true;
          if( noisy_event ) noisy_event = false;
        }
      }
    }
    else
    {
      if( voltageVec.at(i) > pmax )
      {
        pmax = voltageVec.at(i);
        pmax_index = i;
      }
      else if( (voltageVec.at(i) < assist_threshold) && (assist_threshold - voltageVec.at(i)) <= (assist_threshold/2.0) )
      {
        multiple_singal_pmax_v.push_back( pmax );
        multiple_singal_tmax_v.push_back( timeVec.at(pmax_index) );
        indexing_v.push_back( pmax_index );
        pmax = voltageVec.at(i);
        pmax_index = i;
        candidate_signal = false;
      }
      else{}
    }
  }

  if( noisy_event )
  {
    multiple_singal_pmax_v.push_back( 1.23456e12 );
    multiple_singal_tmax_v.push_back( 1.23456e12 );
    indexing_v.push_back( 0 );
    std::cout<<"Noisy" <<std::endl;
  }
}



/*==============================================================================
Multiple signal peaks counter.
  param

  return number of signal peaks.
===============================================================================*/
int Signal_Peak_Counter(
  std::vector<double> &voltageVec,
  std::vector<double> timeVec,
  double assisting_threshold
)
{
  std::vector<double> multiple_singal_pmax_v = {};
  multiple_singal_pmax_v.reserve(timeVec.size());

  double pmax = 0.0;
  //unsigned int pmax_index = 0;

  bool candidate_signal = false;
  bool noisy_event = true;

  std::size_t npoints = voltageVec.size();

  for( unsigned int i = 0; i < npoints; i++ )
  {
    if( !candidate_signal )
    {
      if( voltageVec.at(i) >= assisting_threshold )
      {
        if( voltageVec.at(i) > pmax )
        {
          pmax = voltageVec.at(i);
          //pmax_index = i;
          candidate_signal = true;
        }
      }
    }
    else
    {
      if( voltageVec.at(i) > pmax )
      {
        pmax = voltageVec.at(i);
        //pmax_index = i;
      }
      else if( (voltageVec.at(i) < assisting_threshold) && (assisting_threshold - voltageVec.at(i)) <= (assisting_threshold) )
      {
        multiple_singal_pmax_v.push_back( pmax );
        pmax = voltageVec.at(i);
        //pmax_index = i;
        candidate_signal = false;
        if( noisy_event ) noisy_event = false;
      }
      else{}
    }
  }

  if( !noisy_event )
  {
    return multiple_singal_pmax_v.size();
  }
  else return 0;
}
