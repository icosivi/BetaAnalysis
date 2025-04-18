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
#include <TH1F.h>
#include <TF1.h>
#include <TGraph.h>
#include <TThread.h>
#include <TROOT.h>
#include <TStyle.h>
#include <TImage.h>
#include <TCanvas.h>






ClassImp(Analyzer)

Analyzer::Analyzer(std::vector<float> voltage, std::vector<float> time):TObject(),
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


float Analyzer::Correct_Baseline( int ptN ){

  float mean =0;

  for(std::size_t j = 0, max = ptN; j < max; j++){mean += this->pvoltage.at(j);}

  mean = mean/ptN;

  for(std::size_t j = 0, max = this->pvoltage.size(); j < max; j++){this->pvoltage.at(j) = this->pvoltage.at(j)- mean;}

  return mean;

}


float Analyzer::Find_Noise( const unsigned int inoise){

  float rms = 0.0, mean = 0.0, var = 0.0;

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


std::pair<float, unsigned int> Analyzer::Find_Signal_Maximum(bool confineSearchRegion, float searchRange[2]){

	  float          pmax       = 0.0;
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


std::array<float, 3> Analyzer::Pmax_with_GausFit(const std::pair<float, unsigned int> Pmax, unsigned int maxIndex, int samples_fit){

  std::array<float, 3> result;
  float pmax, tmax, chi2;
  unsigned int pmaxIndex = Pmax.second;
  float time_bin = this->ptime.at(1)-this->ptime.at(0);

  if( pmaxIndex > 5 && pmaxIndex < maxIndex-5 ){

    float time_min = this->ptime.at(pmaxIndex-( (samples_fit-1)/2 ));
    float time_max = this->ptime.at(pmaxIndex+( (samples_fit-1)/2 ));
    TH1D pmax_histo("pmax_histo","pmax_histo",samples_fit,time_min,time_max);

    bool good_fit = true;
    int points_above_zero_counter = 0;

    for(int i=0; i<samples_fit; i++){

      if( Pmax.first*this->pvoltage.at(pmaxIndex-( (samples_fit-1)/2 )+i)>0 ){ 
        
        pmax_histo.Fill( this->ptime.at(pmaxIndex-( (samples_fit-1)/2 )+i) , this->pvoltage.at(pmaxIndex-( (samples_fit-1)/2 )+i) );
        points_above_zero_counter++ ;

      }

     }

     if( points_above_zero_counter < samples_fit-2 ) good_fit = false;

     /*for(int i=0; i<samples_fit; i++){

      if( Pmax.first*this->pvoltage.at(pmaxIndex-( (samples_fit-1)/2 )+i)>0 ) pmax_histo.Fill( this->ptime.at(pmaxIndex-( (samples_fit-1)/2 )+i) , this->pvoltage.at(pmaxIndex-( (samples_fit-1)/2 )+i) );
      else{

        good_fit = false;
        break;

      }

     }*/


   if(good_fit){

    TF1 f("f","gaus",time_min,time_max);
    f.SetParameter(0,Pmax.first);
    f.SetParameter(1,this->ptime.at(Pmax.second));  //pmaxIndex*time_bin
    f.SetParameter(2,samples_fit*time_bin);
    pmax_histo.Fit("f","RN0Q");
    pmax = f.GetParameter(0);
    tmax = f.GetParameter(1);
    chi2 = f.GetChisquare();

   } else {

     pmax = Pmax.first;
     tmax = this->ptime.at(Pmax.second);
     chi2 = -10000;

   }

  } else {

   pmax = Pmax.first;
   tmax = this->ptime.at(Pmax.second);
   chi2 = -10000;

  }


  //if( fabs(pmax-Pmax.first)<0.2*fabs(Pmax.first) ) return std::make_pair( pmax, tmax);
  //else return std::make_pair( Pmax.first, this->ptime.at(Pmax.second) );
  //return std::make_pair( pmax, tmax);

  result = {pmax,tmax,chi2};
  return result;

}


std::array<float, 3> Analyzer::Pmax_for_samples(const std::pair<float, unsigned int> Pmax, unsigned int maxIndex, int samples_fit){

  std::array<float, 3> result;
  float pmax, tmax, sigma;
  unsigned int pmaxIndex = Pmax.second;
  float time_bin = this->ptime.at(1)-this->ptime.at(0);

  if( pmaxIndex > 5 && pmaxIndex < maxIndex-5 ){

    float time_min = this->ptime.at(pmaxIndex-( (samples_fit-1)/2 ));
    float time_max = this->ptime.at(pmaxIndex+( (samples_fit-1)/2 ));
    TH1D pmax_histo("pmax_histo","pmax_histo",samples_fit,time_min,time_max);

    bool good_fit = true;

    for(int i=0; i<samples_fit; i++){

      if(Pmax.first*this->pvoltage.at(pmaxIndex-( (samples_fit-1)/2 )+i)>0) pmax_histo.Fill( this->ptime.at(pmaxIndex-( (samples_fit-1)/2 )+i) , this->pvoltage.at(pmaxIndex-( (samples_fit-1)/2 )+i) );
      else{

        good_fit = false;
        break;

      }

     }

  if(good_fit){

    TF1 f("f","gaus",time_min,time_max);
    f.SetParameter(0,Pmax.first);
    f.SetParameter(1,this->ptime.at(Pmax.second));  //pmaxIndex*time_bin
    f.SetParameter(2,samples_fit*time_bin);
    pmax_histo.Fit("f","RN0Q");
    pmax = f.GetParameter(0);
    tmax = f.GetParameter(1);
    sigma = f.GetParameter(2);

    result = {pmax, tmax, sigma};

  } else {

    result = {-1000., -1000., -1000.};

  }

  } else {

   result = {-1000., -1000., -1000.};

  }


  if( fabs(pmax-Pmax.first)<0.2*fabs(Pmax.first) ) return result;
  else{

    result = {-1000., -1000., -1000.};
    return result;

  } 

}



std::pair<float, unsigned int> Analyzer::Find_Negative_Signal_Maximum( bool confineSearchRegion, float searchRange[2]){

    float pmax = 0.0;
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


std::array<float, 3> Analyzer::Negative_Pmax_with_GausFit(const std::pair<float, unsigned int> NegPmax, unsigned int maxIndex, int samples_fit){

  std::array<float, 3> result;
  float pmax, tmax, chi2;
  unsigned int pmaxIndex = NegPmax.second;
  float time_bin = this->ptime.at(1)-this->ptime.at(0);

  if( pmaxIndex > 5 && pmaxIndex < maxIndex-5 ){

    float time_min = this->ptime.at(pmaxIndex-3);
    float time_max = this->ptime.at(pmaxIndex+3);
    TH1D pmax_histo("pmax_histo","pmax_histo",samples_fit,time_min,time_max);

    bool good_fit = true;

    for(int i=0; i<7; i++){

      if(NegPmax.first*this->pvoltage.at(pmaxIndex-( (samples_fit-1)/2 )+i)>0) pmax_histo.Fill( this->ptime.at(pmaxIndex-( (samples_fit-1)/2 )+i) , -this->pvoltage.at(pmaxIndex-( (samples_fit-1)/2 )+i) );
      else{

        good_fit = false;
        break;

      }

     }

  if(good_fit){

    TF1 f("f","gaus",time_min,time_max);
    f.SetParameter(0,-NegPmax.first);
    f.SetParameter(1,this->ptime.at(NegPmax.second));  //pmaxIndex*time_bin
    f.SetParameter(2,samples_fit*time_bin);
    pmax_histo.Fit("f","RN0Q");
    pmax = -f.GetParameter(0);
    tmax = f.GetParameter(1);
    chi2 = f.GetChisquare();

  } else {

   pmax = NegPmax.first;
   tmax = this->ptime.at(NegPmax.second);
   chi2 = -10000;

  }

  } else {

   pmax = NegPmax.first;
   tmax = this->ptime.at(NegPmax.second);
   chi2 = -10000;

  }

  //return std::make_pair( pmax, tmax);

  result = {pmax,tmax,chi2};
  return result;

}


float Analyzer::Get_Tmax(const std::pair<float, unsigned int> Pmax){

  float tmax = this->ptime.at(Pmax.second);
  return tmax;


}


float Analyzer::Get_Negative_Tmax(const std::pair<float, unsigned int> NegPmax){

  float tmax = this->ptime.at(NegPmax.second);
  return tmax;


}


float Analyzer::DC_Area(float baseline_correction){ //


  float dc_area = 0;
  float time_difference = this->ptime.at(1) - this->ptime.at(0);

  for(int j=5; j<int(this->pvoltage.size()-10); j++ ){ 

    dc_area += (this->pvoltage.at(j))*time_difference ;  //-baseline_correction
    //std::cout<<(this->pvoltage.at(j))*time_difference<<std::endl;

  }

  return dc_area ;

}


float Analyzer::Find_Pulse_Area(const std::pair<float, unsigned int> Pmax){

    float pulse_area = 0.0;
    const float time_difference = this->ptime.at(1) - this->ptime.at(0);

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


float Analyzer::Find_Undershoot_Area(const std::pair<float, unsigned int> Pmax){

    float undershoot_area = 0.0;
    const float time_difference = this->ptime.at(1) - this->ptime.at(0);

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


float Analyzer::Area_NC(const std::pair<float,unsigned int> Pmax, int t_beforeSignal, int t_afterSignal, float rms){

  float area = 0;
  float area_pos = 0;
  float tempinch= 0;
  float tempminch2= 0;
  float bck= 0;
  float TMax = 0;
  float AMax = 0;
  int NFlag = 0;

  //float RMSbck = 2;

  unsigned int imax = Pmax.second;

  std::size_t npoints = this->pvoltage.size();
  if( imax == npoints-1 )imax = imax - 1;//preventing out of range.

  int t_start = imax - t_beforeSignal ;
  if(t_start<=0) t_start = 1;

  int t_stop = imax + t_afterSignal ;
  if(t_stop>=int(npoints)) t_stop = npoints-1;

  for(int j=0;j<t_start;j++)
    {
      bck += this->pvoltage.at(j);
      tempinch++;
    }

  bck *=1/tempinch;
  tempinch=0;

  for(int j=t_start;j<t_stop;j++){
      
      if( this->pvoltage.at(j) > AMax){
        TMax = j;
        AMax = this->pvoltage.at(j);
      }

  }


  for(int j=t_start;j<t_stop;j++){
      //   if there is a signal, do not sum the undershoot
      if (AMax> 3*rms){
        if (j<=TMax) tempminch2 += this->pvoltage.at(j)-bck;
        else  if (j>TMax && this->pvoltage.at(j)-bck > 0 && NFlag ==0){
          tempminch2 += this->pvoltage.at(j)-bck;
        }
        else if (j>TMax && this->pvoltage.at(j)-bck < 0){
          NFlag = 1;
        }
      } else tempminch2 += this->pvoltage.at(j)-bck;
     
      tempinch += this->pvoltage.at(j)-bck;
    }


    area_pos=tempminch2*(this->ptime.at(2)-this->ptime.at(1));;
    area=tempinch*(this->ptime.at(2)-this->ptime.at(1));;
 
  
  return area;

}


float Analyzer::Area_NC_pos(const std::pair<float,unsigned int> Pmax, int t_beforeSignal, int t_afterSignal, float rms){

  float area = 0;
  float area_pos = 0;
  float tempinch= 0;
  float tempminch2= 0;
  float bck= 0;
  float TMax = 0;
  float AMax = 0;
  int NFlag = 0;

  //float RMSbck = 2;

  unsigned int imax = Pmax.second;

  std::size_t npoints = this->pvoltage.size();
  if( imax == npoints-1 )imax = imax - 1;//preventing out of range.

  int t_start = imax - t_beforeSignal ;
  if(t_start<=0) t_start = 1;

  int t_stop = imax + t_afterSignal ;
  if(t_stop>=int(npoints)) t_stop = npoints-1;

  for(int j=0;j<t_start;j++)
    {
      bck += this->pvoltage.at(j);
      tempinch++;
    }

  bck *=1/tempinch;
  tempinch=0;

  for(int j=t_start;j<t_stop;j++){
      
      if( this->pvoltage.at(j) > AMax){
        TMax = j;
        AMax = this->pvoltage.at(j);
      }

  }


  for(int j=t_start;j<t_stop;j++){
      //   if there is a signal, do not sum the undershoot
      if (AMax> 3*rms){
        if (j<=TMax) tempminch2 += this->pvoltage.at(j)-bck;
        else  if (j>TMax && this->pvoltage.at(j)-bck > 0 && NFlag ==0){
          tempminch2 += this->pvoltage.at(j)-bck;
        }
        else if (j>TMax && this->pvoltage.at(j)-bck < 0){
          NFlag = 1;
        }
      } else tempminch2 += this->pvoltage.at(j)-bck;
     
      tempinch += this->pvoltage.at(j)-bck;
    }


    area_pos=tempminch2*(this->ptime.at(2)-this->ptime.at(1));;
    area=tempinch*(this->ptime.at(2)-this->ptime.at(1));;
 
  
  return area_pos;

}


float Analyzer::Pulse_Integration_with_Fixed_Window_Size(const std::pair<float,unsigned int> Pmax, std::string integration_option, 
                                                          float t_beforeSignal, float t_afterSignal){
  
  float pulse_area = 0.0;
  const float time_difference = this->ptime.at(1) - this->ptime.at(0);
  //float tRange[2] = {t_beforeSignal*10e-9, t_afterSignal*10e-9};
  float tRange[2] = {t_beforeSignal, t_afterSignal};

  unsigned int imax = Pmax.second;

  float timeOfMaximum = this->ptime.at(imax);
  std::size_t npoints = this->pvoltage.size();

  if( imax == npoints-1 ) imax = imax - 1;//preventing out of range.

  const float _20pmax = Pmax.first * 0.20;
  const float _10pmax = Pmax.first * 0.10;
  float _20pmax_time = this->ptime.at(0);
  float _10pmax_time = -this->ptime.at(0);
  float start_time = -this->ptime.at(0);
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

  std::vector<float> integration_voltage_vector;
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


float Analyzer::Pulse_Integration_with_Fixed_Window_Size_with_GausFit(const std::pair<float,float> Pmax, unsigned int imax, 
                                                                       std::string integration_option, float t_beforeSignal, float t_afterSignal){
  
  float pulse_area = 0.0;
  const float time_difference = this->ptime.at(1) - this->ptime.at(0);
  
  //float tRange[2] = {t_beforeSignal*10e-9, t_afterSignal*10e-9};
  float tRange[2] = {t_beforeSignal, t_afterSignal};

  //unsigned int imax = Pmax.second;

  float timeOfMaximum = Pmax.second;
  std::size_t npoints = this->pvoltage.size();

  if( imax == npoints-1 ) imax = imax - 1;//preventing out of range.

  const float _20pmax = Pmax.first * 0.20;
  const float _10pmax = Pmax.first * 0.10;
  float _20pmax_time = this->ptime.at(0);
  float _10pmax_time = -this->ptime.at(0);
  float start_time = -this->ptime.at(0);
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

  std::vector<float> integration_voltage_vector;
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


float Analyzer::Pulse_Area_With_Linear_Interpolate_Edge( const std::pair<float,unsigned int> Pmax, std::string integration_option, 
                                                          bool relativeTimeWindow, float StopTime ){

  float pulse_area = 0.0;
  const float time_difference = this->ptime.at(1) - this->ptime.at(0);

  unsigned int imax = Pmax.second;

  float timeOfMaximum = this->ptime.at(imax);
  std::size_t npoints = this->pvoltage.size();

  float stopTime = StopTime*10e-9;

  if( imax == npoints-1 ) imax = imax - 1;//preventing out of range.

  const float _20pmax = Pmax.first * 0.20;
  const float _10pmax = Pmax.first * 0.10;
  float _20pmax_time = 0.0;
  float _10pmax_time = 0.0;
  unsigned int istart = 0;
  unsigned int iend = 0;
  float start_time = 0.0;
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
    std::vector<float> integration_y;
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


float Analyzer::Pulse_Area_With_Linear_Interpolate_Edge_with_GausFit( const std::pair<float,float> Pmax, unsigned int imax, 
                                                                       std::string integration_option, bool relativeTimeWindow, float StopTime ){

  float pulse_area = 0.0;
  const float time_difference = this->ptime.at(1) - this->ptime.at(0);

  //unsigned int imax = Pmax.second;

  float timeOfMaximum = Pmax.second;
  std::size_t npoints = this->pvoltage.size();

  float stopTime = StopTime*10e-9;

  if( imax == npoints-1 ) imax = imax - 1;//preventing out of range.

  const float _20pmax = Pmax.first * 0.20;
  const float _10pmax = Pmax.first * 0.10;
  float _20pmax_time = 0.0;
  float _10pmax_time = 0.0;
  unsigned int istart = 0;
  unsigned int iend = 0;
  float start_time = 0.0;
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
    std::vector<float> integration_y;
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


float Analyzer::Find_Rise_Time(const std::pair<float, unsigned int> Pmax, float bottom , float top){


float rise = 0.0;

  unsigned int itop = this->pvoltage.size()-2, ibottom = 0;

  bool ten = true, ninety = true;

  unsigned int imax = Pmax.second;
  float pmax = Pmax.first;

  float lowerval = pmax * bottom;
  float upperval = pmax * top;

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
  float tb = this->ptime.at(ibottom);
  float pb = this->pvoltage.at(ibottom);
  float tb_1 =  this->ptime.at(ibottom + 1);
  float pb_1 = this->pvoltage.at(ibottom + 1);

  float tt = this->ptime.at(itop);
  float pt = this->pvoltage.at(itop);
  float tt_1 =  this->ptime.at(itop + 1);
  float pt_1 = this->pvoltage.at(itop + 1);

  float tbottom = xlinearInter( tb, pb, tb_1, pb_1, lowerval);
  float ttop    = xlinearInter( tt, pt, tt_1, pt_1, upperval);

  rise  = ttop - tbottom; // rise
  return rise;


}


float Analyzer::Find_Rise_Time_with_GausFit(const std::pair<float, float> Pmax, unsigned int imax, float bottom , float top){


float rise = 0.0;

  //unsigned int itop = this->pvoltage.size()-2, ibottom = 0;
  unsigned int itop = 500, ibottom = 500;

  bool ten = true, ninety = true;

  //unsigned int imax = Pmax.second;
  float pmax = Pmax.first;

  float lowerval = pmax * bottom;
  float upperval = pmax * top;

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
  float tb = this->ptime.at(ibottom);
  float pb = this->pvoltage.at(ibottom);
  float tb_1 =  this->ptime.at(ibottom + 1);
  float pb_1 = this->pvoltage.at(ibottom + 1);

  float tt = this->ptime.at(itop);
  float pt = this->pvoltage.at(itop);
  float tt_1 =  this->ptime.at(itop + 1);
  float pt_1 = this->pvoltage.at(itop + 1);

  float tbottom = xlinearInter( tb, pb, tb_1, pb_1, lowerval);
  float ttop    = xlinearInter( tt, pt, tt_1, pt_1, upperval);

  rise  = ttop - tbottom; // rise
  return rise;


}


float Analyzer::Find_Fall_Time_with_GausFit(const std::pair<float, float> Pmax, unsigned int imax, float bottom , float top){


float rise = 0.0;

  //unsigned int itop = this->pvoltage.size()-2, ibottom = 0;
  unsigned int itop = 500, ibottom = 500;

  bool ten = true, ninety = true;
  std::size_t npoints = this->pvoltage.size()-1;

  //unsigned int imax = Pmax.second;
  float pmax = Pmax.first;

  float lowerval = pmax * bottom;
  float upperval = pmax * top;

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
  float tb = this->ptime.at(ibottom);
  float pb = this->pvoltage.at(ibottom);
  float tb_1 =  this->ptime.at(ibottom - 1);
  float pb_1 = this->pvoltage.at(ibottom - 1);

  float tt = this->ptime.at(itop);
  float pt = this->pvoltage.at(itop);
  float tt_1 =  this->ptime.at(itop - 1);
  float pt_1 = this->pvoltage.at(itop - 1);

  float tbottom = xlinearInter( tb, pb, tb_1, pb_1, lowerval);
  float ttop    = xlinearInter( tt, pt, tt_1, pt_1, upperval);

  rise  = ttop - tbottom; // rise
  return rise;


}


float Analyzer::Find_Dvdt(const int fraction, const int ndif, const std::pair<float,unsigned int> Pmax){

    float time_difference = 0.0;
    float dvdt = 0.0;
    unsigned int ifraction = 0;

    time_difference = this->ptime.at(1) - this->ptime.at(0);

    float pmax = Pmax.first;
    unsigned int imax = Pmax.second;

    for( int j = imax; j>-1; j--)
    {
      if( this->pvoltage.at(j) <= pmax*float(fraction)/100)
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


float Analyzer::Find_Dvdt_with_GausFit(const int fraction, const int ndif, const std::pair<float,float> Pmax, unsigned int imax){

    float time_difference = 0.0;
    float dvdt = 0.0;
    //unsigned int ifraction = 0;
    int ifraction = 0;

    time_difference = this->ptime.at(1) - this->ptime.at(0);

    float pmax = Pmax.first;
    //unsigned int imax = Pmax.second;

    for( int j = imax; j>-1; j--)
    {
      if( this->pvoltage.at(j) <= pmax*float(fraction)/100)
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


float Analyzer::Find_Dvdt2080_with_GausFit(const int ndif, const std::pair<float,float> Pmax, unsigned int imax){

    float time_difference = 0.0;
    float dvdt = 0.0;
    //unsigned int ifraction = 0;
    int ifraction = 0;
    int ifraction2 = 0;

    float _20pmax_time = 0;
    float _80pmax_time = 0;

    time_difference = this->ptime.at(1) - this->ptime.at(0);

    float pmax = Pmax.first;
    //unsigned int imax = Pmax.second;

    for( int j = imax; j>-1; j--)
    {
      if( this->pvoltage.at(j) <= pmax*float(20)/100)
      {
        ifraction = j;

        break;
      }
    }//find index of first point before constant fraction of pulse

    for( int j = imax; j>-1; j--)
    {
      if( this->pvoltage.at(j) <= pmax*float(80)/100)
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


float Analyzer::Rising_Edge_CFD_Time(const float fraction, const std::pair<float,unsigned int> Pmax){

    float pmax = Pmax.first;
    unsigned int imax = Pmax.second;

    float time_fraction = 0.0;
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


float Analyzer::Rising_Edge_CFD_Time_with_GausFit(const float fraction, const std::pair<float, float> Pmax, unsigned int imax){

    float pmax = Pmax.first;
    //unsigned int imax = Pmax.second;

    float time_fraction = 0.0;
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


float Analyzer::Falling_Edge_CFD_Time_with_GausFit(const float fraction, const std::pair<float, float> Pmax, unsigned int imax){

    float pmax = Pmax.first;
    //unsigned int imax = Pmax.second;

    float time_fraction = 0.0;
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



float Analyzer::Find_Time_At_Threshold_with_GausFit(const float thresholdLevel, const std::pair<float,float> Pmax, unsigned int imax){

  float thr = thresholdLevel;

  float timeAtThreshold = 0.0, timeBelowThreshold = 0.0;

  unsigned int timeBelowThreshold_index = 0;

  unsigned int pmax_index = imax;
  float pmax = Pmax.first;
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


float Analyzer::Find_Time_At_Threshold_Falling_Edge_with_GausFit(const float thresholdLevel, const std::pair<float,float> Pmax, unsigned int imax){

  float thr = thresholdLevel;

  float timeAtThreshold = 0.0, timeBelowThreshold = 0.0;

  unsigned int timeBelowThreshold_index = 0;

  unsigned int pmax_index = imax;
  float pmax = Pmax.first;
  std::size_t npoints = this->pvoltage.size();

  bool failure = true;

  if( pmax_index == npoints-1 ) pmax_index = pmax_index - 1;//preventing out of range

    if( pmax < thr ){ return -1000.0;}
    else
    {
      for( int i = pmax_index; i < int(npoints); i++ )
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
float Analyzer::Find_Time_Over_Threshold(const float first_thresholdLevel, const std::pair<float,unsigned int> Pmax, const float second_thresholdLevel){

  float thr1 = first_thresholdLevel;
  float thr2 = second_thresholdLevel;

  float timeAtThreshold1 = 0.0, timeBelowThreshold1 = 0.0;
  float timeAtThreshold2 = 0.0, timeBelowThreshold2 = 0.0;

  unsigned int timeBelowThreshold1_index = 0;
  unsigned int timeBelowThreshold2_index = 0;

  unsigned int pmax_index = Pmax.second;
  float pmax = Pmax.first;
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
float Analyzer::New_Pulse_Area( const std::pair<float,float> Pmax, unsigned int imax, std::string integration_option, float range[2] ){

  if(Pmax.second > range[0] && Pmax.second < range[1]){

    float pulse_area = 0.0;
    const float time_difference = this->ptime.at(1) - this->ptime.at(0);

    //unsigned int imax = Pmax.second;

    float timeOfMaximum = Pmax.second;
    std::size_t npoints = this->pvoltage.size();

    //if( imax == npoints-1 ) imax = imax - 1;//preventing out of range.

    const float _20pmax = Pmax.first * 0.20;
    const float _10pmax = Pmax.first * 0.10;
    float _20pmax_time = 0.0;
    float _10pmax_time = 0.0;
    float _20pmax_time_2 = 0.0;
    float _10pmax_time_2 = 0.0;
    unsigned int istart = 0;
    unsigned int iend = 0;
    float start_time = 0.0;
    float end_time = 0.0;
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
        std::vector<float> integration_y;
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
float Analyzer::New_Undershoot_Area( const std::pair<float,float> Pmax, const std::pair<float,float> Pmin, unsigned int imin, 
                                      std::string integration_option, float range[2]){

  if(Pmax.second > range[0] && Pmax.second < range[1]){

    float pulse_area = 0.0;
    const float time_difference = this->ptime.at(1) - this->ptime.at(0);

    //unsigned int imin = Pmin.second;
    std::size_t npoints = this->pvoltage.size();

    //if( imin == npoints-1 ) imin = imin - 1;//preventing out of range.

    const float _20pmin = Pmin.first * 0.20;
    const float _10pmin = Pmin.first * 0.10;
    float _20pmin_time2 = 0.0;
    float _10pmin_time2 = 0.0;
    float _20pmin_time = 0.0;
    float _10pmin_time = 0.0;
    unsigned int istart = 0;
    unsigned int iend = 0;
    float start_time = 0.0;
    float end_time = 0.0;
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
        std::vector<float> integration_y;
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
