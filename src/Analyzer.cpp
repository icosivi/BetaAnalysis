#include "../include/general.hpp"
#include "Analyzer.hpp"

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
#include <Riostream.h>
#include <TTree.h>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>
#include <TTreeReaderArray.h>
#include <TBranch.h>
#include <TFile.h>
#include <TH1.h>
#include <TF1.h>
#include <TGraph.h>
#include <TThread.h>
#include <TROOT.h>
#include <TStyle.h>
#include <TImage.h>
#include <TCanvas.h>






ClassImp(Analyzer)

Analyzer::Analyzer(std::vector<double> voltage, std::vector<double> time):TObject(),
pvoltage(voltage),
ptime(time){

}


Analyzer::Analyzer():TObject(),
pvoltage(0),
ptime(0){

}


Analyzer::Analyzer(const Analyzer &a):TObject(),
pvoltage(a.pvoltage),
ptime(a.ptime){

}


Analyzer::~Analyzer(){

//delete pvoltage;
//delete ptime;

}


double Analyzer::Correct_Baseline( int ptN ){

  double mean =0;

  for(std::size_t j = 0, max = ptN; j < max; j++){mean += this->pvoltage.at(j);}

  mean = mean/ptN;

  for(std::size_t j = 0, max = this->pvoltage.size(); j < max; j++){this->pvoltage.at(j) = this->pvoltage.at(j)- mean;}

  return mean;

}


double Analyzer::Find_Noise( const unsigned int inoise){

  double rms = 0.0, mean = 0.0, var = 0.0;

  for( unsigned int j = 0; j < inoise; j++)
  {
    rms  += this->pvoltage.at(j)*this->pvoltage.at(j);
    mean += this->pvoltage.at(j);
  }

  mean = mean / inoise;
  rms  = rms/inoise;
  var  = rms - mean * mean;
  rms  = pow(var, 0.5);

  return rms;
}


std::pair<double, unsigned int> Analyzer::Find_Signal_Maximum(bool confineSearchRegion, double searchRange[2]){

	  double          pmax       = 0.0;
    unsigned int    pmaxIndex  = 0;
    bool   firstPoint  = true;
    std::size_t npoints = this->pvoltage.size();

    if( confineSearchRegion )
    {
      for( std::size_t j = 0; j < npoints; j++)
      {
        if( this->ptime.at(j) >= searchRange[0] && this->ptime.at(j) <= searchRange[1] ) //zoom in to find the Pmax
        {
            if( firstPoint ){ pmaxIndex = j; firstPoint = false; }
            if( this->pvoltage.at(j) > pmax )
            {
              pmax      = this->pvoltage.at(j);
              pmaxIndex = j;
            }
        }
      }
    }
    else
    {
      for( std::size_t j = 0; j < npoints; j++ )
      {
          if (j == 0)
          {
            pmax = this->pvoltage.at(j);
            pmaxIndex = j;
          }
          if (j != 0 && this->pvoltage.at(j) > pmax)
          {
            pmax = this->pvoltage.at(j);
            pmaxIndex = j;
        }
      }
    }

    return std::make_pair( pmax, pmaxIndex);



}


std::pair<double, double> Analyzer::Pmax_with_GausFit(const std::pair<double, unsigned int> Pmax, unsigned int maxIndex){

  double pmax, tmax;
  unsigned int pmaxIndex = Pmax.second;
  double time_bin = this->ptime.at(1)-this->ptime.at(0);

  if( pmaxIndex > 5 && pmaxIndex < maxIndex-5 ){

    double time_min = this->ptime.at(pmaxIndex-3);
    double time_max = this->ptime.at(pmaxIndex+3);
    TH1D pmax_histo("pmax_histo","pmax_histo",7,time_min,time_max);

    bool good_fit = true;

    for(int i=0; i<7; i++){

      if(Pmax.first*this->pvoltage.at(pmaxIndex-3+i)>0) pmax_histo.Fill( this->ptime.at(pmaxIndex-3+i) , this->pvoltage.at(pmaxIndex-3+i) );
      else{

        good_fit = false;
        break;

      }

     }

  if(good_fit){

    TF1 f("f","gaus",time_min,time_max);
    f.SetParameter(0,Pmax.first);
    f.SetParameter(1,this->ptime.at(Pmax.second));  //pmaxIndex*time_bin
    f.SetParameter(2,7*time_bin);
    pmax_histo.Fit("f","RN0Q");
    pmax = f.GetParameter(0);
    tmax = f.GetParameter(1);

  } else {

   pmax = Pmax.first;
   tmax = this->ptime.at(Pmax.second);

  }

  } else {

   pmax = Pmax.first;
   tmax = this->ptime.at(Pmax.second);

  }


  if( fabs(pmax-Pmax.first)<0.2*fabs(Pmax.first) ) return std::make_pair( pmax, tmax);
  else return std::make_pair( Pmax.first, this->ptime.at(Pmax.second) );

}


std::pair<double, unsigned int> Analyzer::Find_Negative_Signal_Maximum( bool confineSearchRegion, double searchRange[2]){

    double pmax = 0.0;
    unsigned int pmaxIndex  = 0;
    bool firstPoint  = true;
    std::size_t npoints = this->pvoltage.size();

    if( confineSearchRegion )
    {
      for( std::size_t j = 0; j < npoints; j++)
      {
        if( this->ptime.at(j) >= searchRange[0] && this->ptime.at(j) <= searchRange[1] ) //zoom in to find the Pmax
        {
            if( firstPoint ){ pmaxIndex = j; firstPoint = false; }
            if( this->pvoltage.at(j) < pmax )
            {
              pmax = this->pvoltage.at(j);
              pmaxIndex = j;
            }
        }
      }
    }
    else
    {
      for( std::size_t j = 0, max = npoints; j < max; j++ )
      {
          if (j == 0)
          {
            pmax = this->pvoltage.at(j);
            pmaxIndex = j;
          }
          if (j != 0 && this->pvoltage.at(j) < pmax)
          {
            pmax = this->pvoltage.at(j);
            pmaxIndex = j;
        }
      }
    }

    return std::make_pair( pmax, pmaxIndex);
}


std::pair<double, double> Analyzer::Negative_Pmax_with_GausFit(const std::pair<double, unsigned int> NegPmax, unsigned int maxIndex){

  double pmax, tmax;
  unsigned int pmaxIndex = NegPmax.second;
  double time_bin = this->ptime.at(1)-this->ptime.at(0);

  if( pmaxIndex > 5 && pmaxIndex < maxIndex-5 ){

    double time_min = this->ptime.at(pmaxIndex-3);
    double time_max = this->ptime.at(pmaxIndex+3);
    TH1D pmax_histo("pmax_histo","pmax_histo",7,time_min,time_max);

    bool good_fit = true;

    for(int i=0; i<7; i++){

      if(NegPmax.first*this->pvoltage.at(pmaxIndex-3+i)>0) pmax_histo.Fill( this->ptime.at(pmaxIndex-3+i) , -this->pvoltage.at(pmaxIndex-3+i) );
      else{

        good_fit = false;
        break;

      }

     }

  if(good_fit){

    TF1 f("f","gaus",time_min,time_max);
    f.SetParameter(0,-NegPmax.first);
    f.SetParameter(1,this->ptime.at(NegPmax.second));  //pmaxIndex*time_bin
    f.SetParameter(2,7*time_bin);
    pmax_histo.Fit("f","RN0Q");
    pmax = -f.GetParameter(0);
    tmax = f.GetParameter(1);

  } else {

   pmax = NegPmax.first;
   tmax = this->ptime.at(NegPmax.second);

  }

  } else {

   pmax = NegPmax.first;
   tmax = this->ptime.at(NegPmax.second);

  }

  return std::make_pair( pmax, tmax);

}


double Analyzer::Get_Tmax(const std::pair<double, unsigned int> Pmax){

  double tmax = this->ptime.at(Pmax.second);
  return tmax;


}


double Analyzer::Get_Negative_Tmax(const std::pair<double, unsigned int> NegPmax){

  double tmax = this->ptime.at(NegPmax.second);
  return tmax;


}


double Analyzer::DC_Area(double baseline_correction){ //


  double dc_area = 0;
  double time_difference = this->ptime.at(1) - this->ptime.at(0);

  for(int j=5; j<(this->pvoltage.size()-10); j++ ){ 

    dc_area += (this->pvoltage.at(j))*time_difference ;  //-baseline_correction
    //std::cout<<(this->pvoltage.at(j))*time_difference<<std::endl;

  }

  return dc_area ;

}


double Analyzer::Find_Pulse_Area(const std::pair<double, unsigned int> Pmax){

    double pulse_area = 0.0;
    const double time_difference = this->ptime.at(1) - this->ptime.at(0);

    const unsigned int imax = Pmax.second;
    unsigned int istart = 0;
    unsigned int iend = 0;
    std::size_t npoints = this->pvoltage.size();

    for( int j = imax; j>-1; j-- ) // find index of start of pulse
    {
      if(this->pvoltage.at(j) <= 0) //stop after crossing zero
      {
        istart = j;
        break;
      }
    }
    for( unsigned int j = imax; j< npoints; j++ ) // find index of end of pulse
    {
      if(this->pvoltage.at(j) <= 0)
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
      pulse_area = pulse_area + this->pvoltage.at(j);    ///1000.0;
    }

    pulse_area = pulse_area * time_difference;    ///1.0E12;

    return pulse_area; // collected pulse area, assuming voltage is in V, time is in s

}


double Analyzer::Find_Undershoot_Area(const std::pair<double, unsigned int> Pmax){

    double undershoot_area = 0.0;
    const double time_difference = this->ptime.at(1) - this->ptime.at(0);

    const unsigned int imax = Pmax.second;
    unsigned int istart = 0;
    unsigned int iend = 0;
    std::size_t npoints = this->pvoltage.size();

    for( unsigned int j = imax; j < npoints; j++ ) // find index of start of pulse
    {
      if( this->pvoltage.at(j) <= 0) //stop after crossing zero
      {
        istart = j;
        break;
      }
    }
    for( unsigned int j = istart; j< npoints; j++ ) // find index of end of pulse
    {
      if( this->pvoltage.at(j) >= 0)
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
      undershoot_area = undershoot_area + this->pvoltage.at(j);   ///1000.0;
    }

    undershoot_area = undershoot_area * time_difference;    ///1.0E12;

    return undershoot_area; // collected undershoot area, assuming voltage is in V, time is in s


}


double Analyzer::Pulse_Integration_with_Fixed_Window_Size(const std::pair<double,unsigned int> Pmax, std::string integration_option, 
                                                          double t_beforeSignal, double t_afterSignal){
  
  double pulse_area = 0.0;
  const double time_difference = this->ptime.at(1) - this->ptime.at(0);
  double tRange[2] = {t_beforeSignal*10e-9, t_afterSignal*10e-9};

  unsigned int imax = Pmax.second;

  double timeOfMaximum = this->ptime.at(imax);
  std::size_t npoints = this->pvoltage.size();

  if( imax == npoints-1 ) imax = imax - 1;//preventing out of range.

  const double _20pmax = Pmax.first * 0.20;
  const double _10pmax = Pmax.first * 0.10;
  double _20pmax_time = this->ptime.at(0);
  double _10pmax_time = -this->ptime.at(0);
  double start_time = -this->ptime.at(0);
  bool found_20pmax = false;
  bool found_10pmax = false;

  for( int j = imax; j>-1; j-- ) // find index of start of pulse
  {
    if( !found_20pmax )
    {
      if( this->pvoltage.at(j) <= _20pmax ) //stop after crossing zero
      {
        _20pmax_time = xlinearInter( this->ptime.at(j), this->pvoltage.at(j), this->ptime.at(j+1), this->pvoltage.at(j+1), _20pmax );
        found_20pmax = true;
      }
    }

    if( !found_10pmax )
    {
      if( this->pvoltage.at(j) <= _10pmax ) //stop after crossing zero
      {
        _10pmax_time = xlinearInter( this->ptime.at(j), this->pvoltage.at(j), this->ptime.at(j+1), this->pvoltage.at(j+1), _10pmax );
        found_10pmax = true;
      }
    }

    if( found_10pmax && found_20pmax ) break;
  }

  start_time = xlinearInter( _10pmax_time, _10pmax, _20pmax_time, _20pmax, 0.0 );

  std::vector<double> integration_voltage_vector;
  for( unsigned int i = 0; i < npoints; i++)
  {
    if( this->ptime.at(i) >= (start_time-tRange[0]) && this->ptime.at(i) <= timeOfMaximum + tRange[1])
    //if( t.at(i) >= (start_time) && t.at(i) <= 3000.0)
    {
      integration_voltage_vector.push_back( this->pvoltage.at(i) );
    }
  }

  //===========Simpson's rule=======
  if( integration_option.compare("Simpson") == 0 )
  {
    for( std::size_t i = 0, max = integration_voltage_vector.size(); i < max; i++ )
    {
      if( i == 0 ) pulse_area = pulse_area + (time_difference/3.0) * (integration_voltage_vector.at(i));
      else if( i == integration_voltage_vector.size()-1 ) pulse_area = pulse_area + (time_difference/3.0) * (integration_voltage_vector.at(i));
      else if( i % 2 == 0 ) pulse_area = pulse_area + 2 * (time_difference/3.0) * (integration_voltage_vector.at(i));
      else pulse_area = pulse_area + 4 * (time_difference/3.0) * (integration_voltage_vector.at(i));
    }
  }
  //================================
  //===========Rectangluar=========
  else
  {
    for ( std::size_t j = 0, max = integration_voltage_vector.size(); j < max; j++ )
    {
      pulse_area = pulse_area + (time_difference) * integration_voltage_vector.at(j);
    }
  }
  //================================

  return pulse_area; // collected pulse area, assuming voltage is in volts, time is in seconds
}


double Analyzer::Pulse_Integration_with_Fixed_Window_Size_with_GausFit(const std::pair<double,double> Pmax, unsigned int imax, 
                                                                       std::string integration_option, double t_beforeSignal, double t_afterSignal){
  
  double pulse_area = 0.0;
  const double time_difference = this->ptime.at(1) - this->ptime.at(0);
  double tRange[2] = {t_beforeSignal*10e-9, t_afterSignal*10e-9};

  //unsigned int imax = Pmax.second;

  double timeOfMaximum = Pmax.second;
  std::size_t npoints = this->pvoltage.size();

  if( imax == npoints-1 ) imax = imax - 1;//preventing out of range.

  const double _20pmax = Pmax.first * 0.20;
  const double _10pmax = Pmax.first * 0.10;
  double _20pmax_time = this->ptime.at(0);
  double _10pmax_time = -this->ptime.at(0);
  double start_time = -this->ptime.at(0);
  bool found_20pmax = false;
  bool found_10pmax = false;

  for( int j = imax; j>-1; j-- ) // find index of start of pulse
  {
    if( !found_20pmax )
    {
      if( this->pvoltage.at(j) <= _20pmax ) //stop after crossing zero
      {
        _20pmax_time = xlinearInter( this->ptime.at(j), this->pvoltage.at(j), this->ptime.at(j+1), this->pvoltage.at(j+1), _20pmax );
        found_20pmax = true;
      }
    }

    if( !found_10pmax )
    {
      if( this->pvoltage.at(j) <= _10pmax ) //stop after crossing zero
      {
        _10pmax_time = xlinearInter( this->ptime.at(j), this->pvoltage.at(j), this->ptime.at(j+1), this->pvoltage.at(j+1), _10pmax );
        found_10pmax = true;
      }
    }

    if( found_10pmax && found_20pmax ) break;
  }

  start_time = xlinearInter( _10pmax_time, _10pmax, _20pmax_time, _20pmax, 0.0 );

  std::vector<double> integration_voltage_vector;
  for( unsigned int i = 0; i < npoints; i++)
  {
    if( this->ptime.at(i) >= (start_time-tRange[0]) && this->ptime.at(i) <= timeOfMaximum + tRange[1])
    //if( t.at(i) >= (start_time) && t.at(i) <= 3000.0)
    {
      integration_voltage_vector.push_back( this->pvoltage.at(i) );
    }
  }

  //===========Simpson's rule=======
  if( integration_option.compare("Simpson") == 0 )
  {
    for( std::size_t i = 0, max = integration_voltage_vector.size(); i < max; i++ )
    {
      if( i == 0 ) pulse_area = pulse_area + (time_difference/3.0) * (integration_voltage_vector.at(i));
      else if( i == integration_voltage_vector.size()-1 ) pulse_area = pulse_area + (time_difference/3.0) * (integration_voltage_vector.at(i));
      else if( i % 2 == 0 ) pulse_area = pulse_area + 2 * (time_difference/3.0) * (integration_voltage_vector.at(i));
      else pulse_area = pulse_area + 4 * (time_difference/3.0) * (integration_voltage_vector.at(i));
    }
  }
  //================================
  //===========Rectangluar=========
  else
  {
    for ( std::size_t j = 0, max = integration_voltage_vector.size(); j < max; j++ )
    {
      pulse_area = pulse_area + (time_difference) * integration_voltage_vector.at(j);
    }
  }
  //================================

  return pulse_area; // collected pulse area, assuming voltage is in volts, time is in seconds
}


double Analyzer::Pulse_Area_With_Linear_Interpolate_Edge( const std::pair<double,unsigned int> Pmax, std::string integration_option, 
                                                          bool relativeTimeWindow, double StopTime ){

  double pulse_area = 0.0;
  const double time_difference = this->ptime.at(1) - this->ptime.at(0);

  unsigned int imax = Pmax.second;

  double timeOfMaximum = this->ptime.at(imax);
  std::size_t npoints = this->pvoltage.size();

  double stopTime = StopTime*10e-9;

  if( imax == npoints-1 ) imax = imax - 1;//preventing out of range.

  const double _20pmax = Pmax.first * 0.20;
  const double _10pmax = Pmax.first * 0.10;
  double _20pmax_time = 0.0;
  double _10pmax_time = 0.0;
  unsigned int istart = 0;
  unsigned int iend = 0;
  double start_time = 0.0;
  bool found_20pmax = false;
  bool found_10pmax = false;

  for( int j = imax; j>-1; j-- ) // find index of start of pulse
  {
    if( !found_20pmax )
    {
      if( this->pvoltage.at(j) <= _20pmax ) //stop after crossing zero
      {
        _20pmax_time = xlinearInter( this->ptime.at(j), this->pvoltage.at(j), this->ptime.at(j+1), this->pvoltage.at(j+1), _20pmax );
        found_20pmax = true;
      }
    }

    if( !found_10pmax )
    {
      if( this->pvoltage.at(j) <= _10pmax ) //stop after crossing zero
      {
        _10pmax_time = xlinearInter( this->ptime.at(j), this->pvoltage.at(j), this->ptime.at(j+1), this->pvoltage.at(j+1), _10pmax );
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
        if( (this->ptime.at(j)-timeOfMaximum) >= stopTime )
        {
          iend = j;
          break;
        }
      }
      else
      {
        if( this->ptime.at(j) >= stopTime )
        {
          iend = j;
          break;
        }
      }
    }
  }

  for( unsigned int j = 0; j < iend; j++ )
  {
    if( this->ptime.at(j) >= start_time )
    {
      istart = j;
      pulse_area = pulse_area + ((this->ptime.at(j)-start_time)) * this->pvoltage.at(j);
      break;
    }
  }

  //===========Simpson's rule=======
  if( integration_option.compare("Simpson") == 0 )
  {
    std::vector<double> integration_y;
    for ( unsigned int j = istart; j < iend; j++ )
    {
      integration_y.push_back( this->pvoltage.at(j) );
    }
    for( std::size_t i = 0, max = integration_y.size(); i < max; i++ )
    {
      if( i == 0 ) pulse_area = pulse_area + (time_difference/3.0) * (integration_y.at(i));
      else if( i == integration_y.size()-1 ) pulse_area = pulse_area + (time_difference/3.0) * (integration_y.at(i));
      else if( i % 2 == 0 ) pulse_area = pulse_area + 2 * (time_difference/3.0) * (integration_y.at(i));
      else pulse_area = pulse_area + 4 * (time_difference/3.0) * (integration_y.at(i));
    }
  }

  //===========Rectangular=========
  else
  {
    for( unsigned int j = istart; j < iend; j++ )
    {
      pulse_area = pulse_area + (time_difference) * this->pvoltage.at(j);
    }
  }

  return pulse_area; // collected pulse area, assuming voltage is in volts, time is in seconds
}


double Analyzer::Pulse_Area_With_Linear_Interpolate_Edge_with_GausFit( const std::pair<double,double> Pmax, unsigned int imax, 
                                                                       std::string integration_option, bool relativeTimeWindow, double StopTime ){

  double pulse_area = 0.0;
  const double time_difference = this->ptime.at(1) - this->ptime.at(0);

  //unsigned int imax = Pmax.second;

  double timeOfMaximum = Pmax.second;
  std::size_t npoints = this->pvoltage.size();

  double stopTime = StopTime*10e-9;

  if( imax == npoints-1 ) imax = imax - 1;//preventing out of range.

  const double _20pmax = Pmax.first * 0.20;
  const double _10pmax = Pmax.first * 0.10;
  double _20pmax_time = 0.0;
  double _10pmax_time = 0.0;
  unsigned int istart = 0;
  unsigned int iend = 0;
  double start_time = 0.0;
  bool found_20pmax = false;
  bool found_10pmax = false;

  for( int j = imax; j>-1; j-- ) // find index of start of pulse
  {
    if( !found_20pmax )
    {
      if( this->pvoltage.at(j) <= _20pmax ) //stop after crossing zero
      {
        _20pmax_time = xlinearInter( this->ptime.at(j), this->pvoltage.at(j), this->ptime.at(j+1), this->pvoltage.at(j+1), _20pmax );
        found_20pmax = true;
      }
    }

    if( !found_10pmax )
    {
      if( this->pvoltage.at(j) <= _10pmax ) //stop after crossing zero
      {
        _10pmax_time = xlinearInter( this->ptime.at(j), this->pvoltage.at(j), this->ptime.at(j+1), this->pvoltage.at(j+1), _10pmax );
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
        if( (this->ptime.at(j)-timeOfMaximum) >= stopTime )
        {
          iend = j;
          break;
        }
      }
      else
      {
        if( this->ptime.at(j) >= stopTime )
        {
          iend = j;
          break;
        }
      }
    }
  }

  for( unsigned int j = 0; j < iend; j++ )
  {
    if( this->ptime.at(j) >= start_time )
    {
      istart = j;
      pulse_area = pulse_area + ((this->ptime.at(j)-start_time)) * this->pvoltage.at(j);
      break;
    }
  }

  //===========Simpson's rule=======
  if( integration_option.compare("Simpson") == 0 )
  {
    std::vector<double> integration_y;
    for ( unsigned int j = istart; j < iend; j++ )
    {
      integration_y.push_back( this->pvoltage.at(j) );
    }
    for( std::size_t i = 0, max = integration_y.size(); i < max; i++ )
    {
      if( i == 0 ) pulse_area = pulse_area + (time_difference/3.0) * (integration_y.at(i));
      else if( i == integration_y.size()-1 ) pulse_area = pulse_area + (time_difference/3.0) * (integration_y.at(i));
      else if( i % 2 == 0 ) pulse_area = pulse_area + 2 * (time_difference/3.0) * (integration_y.at(i));
      else pulse_area = pulse_area + 4 * (time_difference/3.0) * (integration_y.at(i));
    }
  }

  //===========Rectangular=========
  else
  {
    for( unsigned int j = istart; j < iend; j++ )
    {
      pulse_area = pulse_area + (time_difference) * this->pvoltage.at(j);
    }
  }

  return pulse_area; // collected pulse area, assuming voltage is in volts, time is in seconds
}


double Analyzer::Find_Rise_Time(const std::pair<double, unsigned int> Pmax, double bottom , double top){


double rise = 0.0;

  unsigned int itop = this->pvoltage.size()-2, ibottom = 0;

  bool ten = true, ninety = true;

  unsigned int imax = Pmax.second;
  double pmax = Pmax.first;

  double lowerval = pmax * bottom;
  double upperval = pmax * top;

  for( int j = imax; j > -1; j--)
  {
    if( ninety && this->pvoltage.at(j) < upperval)
    {
      itop    = j;     //find the index right below 90%
      ninety  = false;
    }
    if( ten && this->pvoltage.at(j) < lowerval)
    {
      ibottom = j;      //find the index right below 10%
      ten     = false;
    }
    if( !ten && !ninety ){ break; }
  }
  if(ibottom == this->pvoltage.size()-1){ibottom--;}
  if(itop == this->pvoltage.size()-1){itop--;}
  //std::cout<<itop<<std::endl;
  //std::cout<<ibottom<<std::endl;
  double tb = this->ptime.at(ibottom);
  double pb = this->pvoltage.at(ibottom);
  double tb_1 =  this->ptime.at(ibottom + 1);
  double pb_1 = this->pvoltage.at(ibottom + 1);

  double tt = this->ptime.at(itop);
  double pt = this->pvoltage.at(itop);
  double tt_1 =  this->ptime.at(itop + 1);
  double pt_1 = this->pvoltage.at(itop + 1);

  double tbottom = xlinearInter( tb, pb, tb_1, pb_1, lowerval);
  double ttop    = xlinearInter( tt, pt, tt_1, pt_1, upperval);

  rise  = ttop - tbottom; // rise
  return rise;


}


double Analyzer::Find_Rise_Time_with_GausFit(const std::pair<double, double> Pmax, unsigned int imax, double bottom , double top){


double rise = 0.0;

  //unsigned int itop = this->pvoltage.size()-2, ibottom = 0;
  unsigned int itop = 500, ibottom = 500;

  bool ten = true, ninety = true;

  //unsigned int imax = Pmax.second;
  double pmax = Pmax.first;

  double lowerval = pmax * bottom;
  double upperval = pmax * top;

  for( int j = imax; j > 0; j--)
  {
    if( ninety && this->pvoltage.at(j) < upperval)
    {
      itop    = j;     //find the index right below 90%
      ninety  = false;
    }
    if( ten && this->pvoltage.at(j) < lowerval)
    {
      ibottom = j;      //find the index right below 10%
      ten     = false;
    }
    if( !ten && !ninety ){ break; }
  }
  if(ibottom == this->pvoltage.size()-1){ibottom--;}
  if(itop == this->pvoltage.size()-1){itop--;}
  //std::cout<<itop<<std::endl;
  //std::cout<<ibottom<<std::endl;
  double tb = this->ptime.at(ibottom);
  double pb = this->pvoltage.at(ibottom);
  double tb_1 =  this->ptime.at(ibottom + 1);
  double pb_1 = this->pvoltage.at(ibottom + 1);

  double tt = this->ptime.at(itop);
  double pt = this->pvoltage.at(itop);
  double tt_1 =  this->ptime.at(itop + 1);
  double pt_1 = this->pvoltage.at(itop + 1);

  double tbottom = xlinearInter( tb, pb, tb_1, pb_1, lowerval);
  double ttop    = xlinearInter( tt, pt, tt_1, pt_1, upperval);

  rise  = ttop - tbottom; // rise
  return rise;


}


double Analyzer::Find_Fall_Time_with_GausFit(const std::pair<double, double> Pmax, unsigned int imax, double bottom , double top){


double rise = 0.0;

  //unsigned int itop = this->pvoltage.size()-2, ibottom = 0;
  unsigned int itop = 500, ibottom = 500;

  bool ten = true, ninety = true;
  std::size_t npoints = this->pvoltage.size()-1;

  //unsigned int imax = Pmax.second;
  double pmax = Pmax.first;

  double lowerval = pmax * bottom;
  double upperval = pmax * top;

  for( unsigned int j = imax; j < npoints; j++)
  {
    if( ninety && this->pvoltage.at(j) < upperval)
    {
      itop    = j;     //find the index right below 90%
      ninety  = false;
      //std::cout<<"pippo"<<std::endl;
    }
    if( ten && this->pvoltage.at(j) < lowerval)
    {
      ibottom = j;      //find the index right below 10%
      ten     = false;
      //std::cout<<"pippo2"<<std::endl;
    }
    if( !ten && !ninety ){ break; }
  }
  if(ibottom == this->pvoltage.size()-1){ibottom--;}
  if(itop == this->pvoltage.size()-1){itop--;}
  //std::cout<<itop<<std::endl;
  //std::cout<<ibottom<<std::endl;
  //std::cout<<" "<<std::endl;
  double tb = this->ptime.at(ibottom);
  double pb = this->pvoltage.at(ibottom);
  double tb_1 =  this->ptime.at(ibottom - 1);
  double pb_1 = this->pvoltage.at(ibottom - 1);

  double tt = this->ptime.at(itop);
  double pt = this->pvoltage.at(itop);
  double tt_1 =  this->ptime.at(itop - 1);
  double pt_1 = this->pvoltage.at(itop - 1);

  double tbottom = xlinearInter( tb, pb, tb_1, pb_1, lowerval);
  double ttop    = xlinearInter( tt, pt, tt_1, pt_1, upperval);

  rise  = ttop - tbottom; // rise
  return rise;


}


double Analyzer::Find_Dvdt(const int fraction, const int ndif, const std::pair<double,unsigned int> Pmax){

    double time_difference = 0.0;
    double dvdt = 0.0;
    unsigned int ifraction = 0;

    time_difference = this->ptime.at(1) - this->ptime.at(0);

    double pmax = Pmax.first;
    unsigned int imax = Pmax.second;

    for( int j = imax; j>-1; j--)
    {
      if( this->pvoltage.at(j) <= pmax*double(fraction)/100)
      {
        ifraction = j;

        break;
      }
    }//find index of first point before constant fraction of pulse

    if(ifraction == this->pvoltage.size()-1) ifraction--;

    if(ndif == 0)
    {
      dvdt = (this->pvoltage.at(ifraction+1) - this->pvoltage.at(ifraction))/time_difference;
    }

    else
    {
      dvdt = (this->pvoltage.at(ifraction+ndif) - this->pvoltage.at(ifraction-ndif))/(time_difference*(ndif*2));
    }

    return dvdt;

}


double Analyzer::Find_Dvdt_with_GausFit(const int fraction, const int ndif, const std::pair<double,double> Pmax, unsigned int imax){

    double time_difference = 0.0;
    double dvdt = 0.0;
    //unsigned int ifraction = 0;
    int ifraction = 0;

    time_difference = this->ptime.at(1) - this->ptime.at(0);

    double pmax = Pmax.first;
    //unsigned int imax = Pmax.second;

    for( int j = imax; j>-1; j--)
    {
      if( this->pvoltage.at(j) <= pmax*double(fraction)/100)
      {
        ifraction = j;

        break;
      }
    }//find index of first point before constant fraction of pulse

    if(ifraction == int(this->pvoltage.size())-1) ifraction--;

    if(ndif == 0 || imax<100 || imax>(this->pvoltage.size()-100))
    {
      dvdt = (this->pvoltage.at(ifraction+1) - this->pvoltage.at(ifraction))/time_difference;
    }

    else
    {
      dvdt = (this->pvoltage.at(ifraction+ndif) - this->pvoltage.at(ifraction-ndif))/(time_difference*(ndif*2));
    }

    return dvdt;

}


double Analyzer::Find_Dvdt2080_with_GausFit(const int ndif, const std::pair<double,double> Pmax, unsigned int imax){

    double time_difference = 0.0;
    double dvdt = 0.0;
    //unsigned int ifraction = 0;
    int ifraction = 0;
    int ifraction2 = 0;

    double _20pmax_time = 0;
    double _80pmax_time = 0;

    time_difference = this->ptime.at(1) - this->ptime.at(0);

    double pmax = Pmax.first;
    //unsigned int imax = Pmax.second;

    for( int j = imax; j>-1; j--)
    {
      if( this->pvoltage.at(j) <= pmax*double(20)/100)
      {
        ifraction = j;

        break;
      }
    }//find index of first point before constant fraction of pulse

    for( int j = imax; j>-1; j--)
    {
      if( this->pvoltage.at(j) <= pmax*double(80)/100)
      {
        ifraction2 = j;

        break;
      }
    }//find index of first point before constant fraction of pulse

    if(ifraction == int(this->pvoltage.size())-1) ifraction--;
    if(ifraction2 == int(this->pvoltage.size())-1) ifraction2--;

    _20pmax_time = xlinearInter( this->ptime.at(ifraction), this->pvoltage.at(ifraction), this->ptime.at(ifraction+1), 
                                 this->pvoltage.at(ifraction+1), pmax*0.2 );

    _80pmax_time = xlinearInter( this->ptime.at(ifraction2), this->pvoltage.at(ifraction2), this->ptime.at(ifraction2+1), 
                                 this->pvoltage.at(ifraction2+1), pmax*0.8 );

    dvdt = 0.6*pmax/(_80pmax_time - _20pmax_time);

    return dvdt;

}


double Analyzer::Rising_Edge_CFD_Time(const double fraction, const std::pair<double,unsigned int> Pmax){

    double pmax = Pmax.first;
    unsigned int imax = Pmax.second;

    double time_fraction = 0.0;
    unsigned int ifraction = 0;

    bool failure = true;

    for( int j = imax; j>-1; j-- )
    {
      if( this->pvoltage.at(j) <= pmax*fraction/100.0)
      {
        ifraction     = j;              //find index of first point before constant fraction of pulse
        time_fraction = this->ptime.at(j);
        failure = false;
        break;
      }
    }
    if(ifraction == this->pvoltage.size()-1) ifraction--;

    if( failure ) time_fraction = this->ptime.at(0);
    else time_fraction = time_fraction + (this->ptime.at(ifraction+1) - this->ptime.at(ifraction))* 
                        (pmax*fraction/100.0 - this->pvoltage.at(ifraction)) /(this->pvoltage.at(ifraction+1) - this->pvoltage.at(ifraction));

    return time_fraction;

}


double Analyzer::Rising_Edge_CFD_Time_with_GausFit(const double fraction, const std::pair<double, double> Pmax, unsigned int imax){

    double pmax = Pmax.first;
    //unsigned int imax = Pmax.second;

    double time_fraction = 0.0;
    unsigned int ifraction = 0;

    bool failure = true;

    for( int j = imax; j>0; j-- )
    {
      if( this->pvoltage.at(j) <= pmax*fraction/100.0)
      {
        ifraction     = j;              //find index of first point before constant fraction of pulse
        time_fraction = this->ptime.at(j);
        failure = false;
        break;
      }
    }
    if(ifraction == this->pvoltage.size()-1) ifraction--;

    if( failure ) time_fraction = this->ptime.at(0);
    else time_fraction = time_fraction + (this->ptime.at(ifraction+1) - this->ptime.at(ifraction))* 
                         (pmax*fraction/100.0 - this->pvoltage.at(ifraction)) /(this->pvoltage.at(ifraction+1) - this->pvoltage.at(ifraction));

    //cout<< this->ptime.at(ifraction) << " " << this->pvoltage.at(ifraction) << endl;
    //cout<< this->ptime.at(ifraction+1) << " " << this->pvoltage.at(ifraction+1) << endl;
    //cout<< "time at 20: " << time_fraction << " voltage at 20: " << pmax*fraction/100.0 << endl;
    //cout<<" "<<endl;
    return time_fraction;

}


double Analyzer::Falling_Edge_CFD_Time_with_GausFit(const double fraction, const std::pair<double, double> Pmax, unsigned int imax){

    double pmax = Pmax.first;
    //unsigned int imax = Pmax.second;

    double time_fraction = 0.0;
    unsigned int ifraction = 1;
    std::size_t npoints = this->pvoltage.size()-1;

    bool failure = true;

    for( unsigned int j = imax; j<npoints; j++ )
    {
      if( this->pvoltage.at(j) <= pmax*fraction/100.0)
      {
        ifraction     = j;              //find index of first point before constant fraction of pulse
        time_fraction = this->ptime.at(j);
        failure = false;
        break;
      }
    }
    if(ifraction == 0) ifraction = 10;

    if( failure ) time_fraction = this->ptime.at(npoints-1);
    else time_fraction = time_fraction + (this->ptime.at(ifraction-1) - this->ptime.at(ifraction))* 
                         (pmax*fraction/100.0 - this->pvoltage.at(ifraction)) /(this->pvoltage.at(ifraction-1) - this->pvoltage.at(ifraction));

    return time_fraction;

}



double Analyzer::Find_Time_At_Threshold_with_GausFit(const double thresholdLevel, const std::pair<double,double> Pmax, unsigned int imax){

  double thr = thresholdLevel;

  double timeAtThreshold = 0.0, timeBelowThreshold = 0.0;

  unsigned int timeBelowThreshold_index = 0;

  unsigned int pmax_index = imax;
  double pmax = Pmax.first;
  std::size_t npoints = this->pvoltage.size();

  if( pmax_index == npoints-1 ) pmax_index = pmax_index - 1;//preventing out of range

    if( pmax < thr ){ return -1000.0;}
    else
    {
      for( int i = pmax_index; i > -1; i-- )
      {
        if( this->pvoltage.at(i) <= thr )
        {
          timeBelowThreshold_index = i;
          timeBelowThreshold = this->ptime.at(i);

          break;
        }
      }

      timeAtThreshold = xlinearInter( timeBelowThreshold, this->pvoltage.at(timeBelowThreshold_index), this->ptime.at(timeBelowThreshold_index+1), 
                                      this->pvoltage.at(timeBelowThreshold_index+1), thr );

      return timeAtThreshold;
    }

}


double Analyzer::Find_Time_At_Threshold_Falling_Edge_with_GausFit(const double thresholdLevel, const std::pair<double,double> Pmax, unsigned int imax){

  double thr = thresholdLevel;

  double timeAtThreshold = 0.0, timeBelowThreshold = 0.0;

  unsigned int timeBelowThreshold_index = 0;

  unsigned int pmax_index = imax;
  double pmax = Pmax.first;
  std::size_t npoints = this->pvoltage.size();

  bool failure = true;

  if( pmax_index == npoints-1 ) pmax_index = pmax_index - 1;//preventing out of range

    if( pmax < thr ){ return -1000.0;}
    else
    {
      for( int i = pmax_index; i < npoints; i++ )
      {
        if( this->pvoltage.at(i) <= thr )
        {
          timeBelowThreshold_index = i;
          timeBelowThreshold = this->ptime.at(i);
          failure = false;

          break;
        }
      }

      if( failure ){

        timeBelowThreshold_index = npoints - 1;
        timeBelowThreshold = this->ptime.at( npoints-1 );

      }

      timeAtThreshold = xlinearInter( timeBelowThreshold, this->pvoltage.at(timeBelowThreshold_index), this->ptime.at(timeBelowThreshold_index-1), 
                                      this->pvoltage.at(timeBelowThreshold_index-1), thr );

      return timeAtThreshold;
    }

}


//Self explainatory
double Analyzer::Find_Time_Over_Threshold(const double first_thresholdLevel, const std::pair<double,unsigned int> Pmax, const double second_thresholdLevel){

  double thr1 = first_thresholdLevel;
  double thr2 = second_thresholdLevel;

  double timeAtThreshold1 = 0.0, timeBelowThreshold1 = 0.0;
  double timeAtThreshold2 = 0.0, timeBelowThreshold2 = 0.0;

  unsigned int timeBelowThreshold1_index = 0;
  unsigned int timeBelowThreshold2_index = 0;

  unsigned int pmax_index = Pmax.second;
  double pmax = Pmax.first;
  std::size_t npoints = this->pvoltage.size();

  if( pmax_index == npoints-1 ) pmax_index = pmax_index - 1;//preventing out of range

//finds time at first threshold
    if( pmax < thr1 || pmax < thr2 ){ return -1000.0;}
    else
    {
      for( int i = pmax_index; i > -1; i-- )
      {
        if( this->pvoltage.at(i) <= thr1 )
        {
          timeBelowThreshold1_index = i;

          timeBelowThreshold1 = this->ptime.at(i);

          break;
        }
      }

      timeAtThreshold1 = xlinearInter( timeBelowThreshold1, this->pvoltage.at(timeBelowThreshold1_index), this->ptime.at(timeBelowThreshold1_index+1), 
                                       this->pvoltage.at(timeBelowThreshold1_index+1), thr1 );

      //finds time at first threshold
      for( int i = pmax_index; i < (int) npoints-1 ; i++)
      {
        if( this->pvoltage.at(i) <= thr2 )
        {
          timeBelowThreshold2_index = i-1;

          timeBelowThreshold2       = this -> ptime.at(i-1);

          break;
        }
      }

      timeAtThreshold2 = xlinearInter( timeBelowThreshold2, this->pvoltage.at(timeBelowThreshold2_index), this->ptime.at(timeBelowThreshold2_index+1), 
                                       this->pvoltage.at(timeBelowThreshold2_index+1), thr2 );

    }
  return TMath::Abs(timeAtThreshold1 - timeAtThreshold2);
}


// Similar to Find_Pulse_Area but start/end times of the pulse are defined analitically. Further checks are useful, there might be bugs. There inputs args not needed!
double Analyzer::New_Pulse_Area( const std::pair<double,double> Pmax, unsigned int imax, std::string integration_option, double range[2] ){

  if(Pmax.second > range[0] && Pmax.second < range[1]){

    double pulse_area = 0.0;
    const double time_difference = this->ptime.at(1) - this->ptime.at(0);

    //unsigned int imax = Pmax.second;

    double timeOfMaximum = Pmax.second;
    std::size_t npoints = this->pvoltage.size();

    //if( imax == npoints-1 ) imax = imax - 1;//preventing out of range.

    const double _20pmax = Pmax.first * 0.20;
    const double _10pmax = Pmax.first * 0.10;
    double _20pmax_time = 0.0;
    double _10pmax_time = 0.0;
    double _20pmax_time_2 = 0.0;
    double _10pmax_time_2 = 0.0;
    unsigned int istart = 0;
    unsigned int iend = 0;
    double start_time = 0.0;
    double end_time = 0.0;
    bool found_20pmax = false;
    bool found_10pmax = false;
    bool found_20pmax_2 = false;
    bool found_10pmax_2 = false;

    for( int j = imax; j>0; j-- ) // find index of start of pulse
    {
      if( !found_20pmax )
      {
        if( this->pvoltage.at(j) <= _20pmax ) //stop after crossing zero
        {
          _20pmax_time = xlinearInter( this->ptime.at(j), this->pvoltage.at(j), this->ptime.at(j+1), this->pvoltage.at(j+1), _20pmax );
          found_20pmax = true;
        }
      }

      if( !found_10pmax )
      {
        if( this->pvoltage.at(j) <= _10pmax ) //stop after crossing zero
        {
          _10pmax_time = xlinearInter( this->ptime.at(j), this->pvoltage.at(j), this->ptime.at(j+1), this->pvoltage.at(j+1), _10pmax );
          found_10pmax = true;
        }
      }

      if( found_10pmax && found_20pmax ) break;
    }

    start_time = xlinearInter( _10pmax_time, _10pmax, _20pmax_time, _20pmax, 0.0 );


    for( unsigned int j = imax; j< npoints-1; j++ ) // find index of end of pulse
    {
      if( !found_20pmax_2 )
      {
        if( this->pvoltage.at(j) <= _20pmax ) //stop after crossing zero
        {
          _20pmax_time_2 = xlinearInter( this->ptime.at(j), this->pvoltage.at(j), this->ptime.at(j-1), this->pvoltage.at(j-1), _20pmax );
          found_20pmax_2 = true;
        }
      }

      if( !found_10pmax_2 )
      {
        if( this->pvoltage.at(j) <= _10pmax ) //stop after crossing zero
        {
          _10pmax_time_2 = xlinearInter( this->ptime.at(j), this->pvoltage.at(j), this->ptime.at(j-1), this->pvoltage.at(j-1), _10pmax );
          found_10pmax_2 = true;
        }
      }

      if( found_10pmax_2 && found_20pmax_2 ) break;
    }

    end_time = xlinearInter( _10pmax_time_2, _10pmax, _20pmax_time_2, _20pmax, 0.0 );
    //cout<<end_time<<endl;
    //end_time = Pmax.second+3e-9;
    
    //if(end_time > end_window || isinf(end_time)) end_time = Pmax.first+5e-9;
    //if(start_time > start_window || isinf(start_time)) start_time = 0;
    //if(start_time > start_window && end_time < end_window && !isinf(start_time) && !isinf(end_time)){

      for( unsigned int j = 0; j < npoints-1; j++ )
      {

        if( this->ptime.at(j) >= start_time )
        {
          istart = j;
          pulse_area = pulse_area + ((this->ptime.at(j)-start_time)) * this->pvoltage.at(j);
          break;
        }

      }

      for(unsigned int j = imax; j < npoints-1; j++)
      {

        if( this->ptime.at(j) >= end_time )
        {
          iend=j-1;
          break;
        }

      }

      //===========Simpson's rule=======
      if( integration_option.compare("Simpson") == 0 )
      {
        std::vector<double> integration_y;
        for ( unsigned int j = istart; j < iend; j++ )
        {
          integration_y.push_back( this->pvoltage.at(j) );
        }
        for( std::size_t i = 0, max = integration_y.size(); i < max; i++ )
        {
          if( i == 0 ) pulse_area = pulse_area + (time_difference/3.0) * (integration_y.at(i));
          else if( i == integration_y.size()-1 ) pulse_area = pulse_area + (time_difference/3.0) * (integration_y.at(i));
          else if( i % 2 == 0 ) pulse_area = pulse_area + 2 * (time_difference/3.0) * (integration_y.at(i));
          else pulse_area = pulse_area + 4 * (time_difference/3.0) * (integration_y.at(i));
        }
      }

      //===========Rectangular=========
      else
      {
        for( unsigned int j = istart; j < iend; j++ )
        {
          pulse_area = pulse_area + (time_difference) * this->pvoltage.at(j);
        }
      }

      return pulse_area; // collected pulse area, assuming voltage is in volts, time is in seconds

    //}else return -1000;

  }else return -1;

}


// Similar to Find_Pulse_Area but start/end times of the pulse are defined analitically. Further checks are useful, there might be bugs. There are inputs args not needed!
double Analyzer::New_Undershoot_Area( const std::pair<double,double> Pmax, const std::pair<double,double> Pmin, unsigned int imin, 
                                      std::string integration_option, double range[2]){

  if(Pmax.second > range[0] && Pmax.second < range[1]){

    double pulse_area = 0.0;
    const double time_difference = this->ptime.at(1) - this->ptime.at(0);

    //unsigned int imin = Pmin.second;
    std::size_t npoints = this->pvoltage.size();

    //if( imin == npoints-1 ) imin = imin - 1;//preventing out of range.

    const double _20pmin = Pmin.first * 0.20;
    const double _10pmin = Pmin.first * 0.10;
    double _20pmin_time2 = 0.0;
    double _10pmin_time2 = 0.0;
    double _20pmin_time = 0.0;
    double _10pmin_time = 0.0;
    unsigned int istart = 0;
    unsigned int iend = 0;
    double start_time = 0.0;
    double end_time = 0.0;
    bool found_20pmin2 = false;
    bool found_10pmin2 = false;
    bool found_20pmin = false;
    bool found_10pmin = false;

    for( int j = imin; j>0; j-- ) // find index of start of pulse
    {
      if( !found_20pmin2 )
      {
        if( this->pvoltage.at(j) >= _20pmin ) //stop after crossing zero
        {
          _20pmin_time2 = xlinearInter( this->ptime.at(j), this->pvoltage.at(j), this->ptime.at(j+1), this->pvoltage.at(j+1), _20pmin );
          found_20pmin2 = true;
        }
      }

      if( !found_10pmin2 )
      {
        if( this->pvoltage.at(j) >= _10pmin ) //stop after crossing zero
        {
          _10pmin_time2 = xlinearInter( this->ptime.at(j), this->pvoltage.at(j), this->ptime.at(j+1), this->pvoltage.at(j+1), _10pmin );
          found_10pmin2 = true;
        }
      }

      if( found_10pmin2 && found_20pmin2 ) break;
    }

    start_time = xlinearInter( _10pmin_time2, _10pmin, _20pmin_time2, _20pmin, 0.0 );


    for( unsigned int j = imin; j< npoints-1; j++ ) // find index of end of pulse
    {
      if( !found_20pmin )
      {
        if( this->pvoltage.at(j) >= _20pmin ) //stop after crossing zero
        {
          _20pmin_time = xlinearInter( this->ptime.at(j), this->pvoltage.at(j), this->ptime.at(j-1), this->pvoltage.at(j-1), _20pmin );
          found_20pmin = true;
        }
      }

      if( !found_10pmin )
      {
        if( this->pvoltage.at(j) >= _10pmin ) //stop after crossing zero
        {
          _10pmin_time = xlinearInter( this->ptime.at(j), this->pvoltage.at(j), this->ptime.at(j-1), this->pvoltage.at(j-1), _10pmin );
          found_10pmin = true;
        }
      }

      if( found_10pmin && found_20pmin ) break;
    }

    end_time = xlinearInter( _10pmin_time, _10pmin, _20pmin_time, _20pmin, 0.0 );


      for( unsigned int j = 0; j < npoints-1; j++ )
      {

        if( this->ptime.at(j) >= start_time )
        {
          istart = j;
          pulse_area = pulse_area + ((this->ptime.at(j)-start_time)) * this->pvoltage.at(j);
          break;
        }

      }

      for(unsigned int j = imin; j < npoints-1; j++)
      {

        if( this->ptime.at(j) >= end_time )
        {
          iend=j-1;
          break;
        }

      }

      //===========Simpson's rule=======
      if( integration_option.compare("Simpson") == 0 )
      {
        std::vector<double> integration_y;
        for ( unsigned int j = istart; j < iend; j++ )
        {
          integration_y.push_back( this->pvoltage.at(j) );
        }
        for( std::size_t i = 0, max = integration_y.size(); i < max; i++ )
        {
          if( i == 0 ) pulse_area = pulse_area + (time_difference/3.0) * (integration_y.at(i));
          else if( i == integration_y.size()-1 ) pulse_area = pulse_area + (time_difference/3.0) * (integration_y.at(i));
          else if( i % 2 == 0 ) pulse_area = pulse_area + 2 * (time_difference/3.0) * (integration_y.at(i));
          else pulse_area = pulse_area + 4 * (time_difference/3.0) * (integration_y.at(i));
        }
      }

      //===========Rectangular=========
      else
      {
        for( unsigned int j = istart; j < iend; j++ )
        {
          pulse_area = pulse_area + (time_difference) * this->pvoltage.at(j);
        }
      }

      return pulse_area; // collected pulse area, assuming voltage is in volts, time is in seconds

    //}else return -1000;

  }else return -1;

}
