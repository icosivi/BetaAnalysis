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

#include "../include/Analyzer.hpp"
#include "../include/general.hpp"



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


void Analyzer::Correct_Baseline( int ptN ){

  double mean =0;

  for(std::size_t j = 0, max = ptN; j < max; j++){mean += this->pvoltage.at(j);}

  mean = mean/ptN;

  for(std::size_t j = 0, max = this->pvoltage.size(); j < max; j++){this->pvoltage.at(j) = this->pvoltage.at(j)- mean;}

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


// Careful with Gaussian fit option
std::pair<double, unsigned int> Analyzer::Find_Signal_Maximum(bool confineSearchRegion, double searchRange[2], bool gaussian_fit){

	  double          pmax       = 0.0;
    unsigned int    pmaxIndex  = 0;
    //bool   strangePeak = true;
    bool   firstPoint  = true;
    std::size_t npoints = this->pvoltage.size();
    double time_bin = this->ptime.at(1)-this->ptime.at(0);

    if( confineSearchRegion )
    {
      for( std::size_t j = 0; j < npoints; j++)
      {
        if( searchRange[0] <= this->ptime.at(j) && this->ptime.at(j) <= searchRange[1] ) //zoom in to find the Pmax
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

    // Gaussian fit to find Pmax. Only for signals in the center of the time window, with an amplitude well above noise threshold
    if( gaussian_fit == true && pmaxIndex>5 && pmaxIndex<995 ){

      double time_min = time_bin*(pmaxIndex-3);
      double time_max = time_bin*(pmaxIndex+3);
      TH1D pmax_histo("pmax_histo","pmax_histo",7,time_min,time_max);

      bool good_fit = true;

      for(int i=0; i<7; i++){

        if(pmax*this->pvoltage.at(pmaxIndex-3+i)>0) pmax_histo.Fill( time_bin*(pmaxIndex-3+i) , this->pvoltage.at(pmaxIndex-3+i) );
        else{

          good_fit = false;
          break;

        }

       }

    if(good_fit){

      TF1 f("f","gaus",time_min,time_max);
      f.SetParameter(0,pmax);
      f.SetParameter(1,pmaxIndex*time_bin);
      f.SetParameter(2,7*time_bin);
      pmax_histo.Fit("f","RN0Q");
      pmax = f.GetParameter(0);

    }

    }

    return std::make_pair( pmax, pmaxIndex);



}



double Analyzer::Get_Tmax(const std::pair<double, unsigned int> Pmax){


  /*if( Pmax.first == this->pvoltage.at(Pmax.second)) double tmax = this->ptime.at(Pmax.second);
  else{



  }*/

  double tmax = this->ptime.at(Pmax.second);
  return tmax;


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
    else time_fraction = time_fraction + (this->ptime.at(ifraction+1) - this->ptime.at(ifraction))* (pmax*fraction/100.0 - this->pvoltage.at(ifraction)) /(this->pvoltage.at(ifraction+1) - this->pvoltage.at(ifraction));

    return time_fraction;

}


double Analyzer::Find_Time_At_Threshold(const double thresholdLevel, const std::pair<double,unsigned int> Pmax){

  double thr = thresholdLevel/1000;

  double timeAtThreshold = 0.0, timeBelowThreshold = 0.0;

  unsigned int timeBelowThreshold_index = 0;

  unsigned int pmax_index = Pmax.second;
  double pmax = Pmax.first;
  std::size_t npoints = this->pvoltage.size();

  if( pmax_index == npoints-1 ) pmax_index = pmax_index - 1;//preventing out of range

    if( pmax < thr ){ return 9999.0;}
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

      timeAtThreshold = xlinearInter( timeBelowThreshold, this->pvoltage.at(timeBelowThreshold_index), this->pvoltage.at(timeBelowThreshold_index+1), this->pvoltage.at(timeBelowThreshold_index+1), thr );

      return timeAtThreshold;
    }

}





