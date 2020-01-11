/////////////////////////////////////////////
////                                      //
////  Baseline Correction implementation //
////                                    //
/////////////////////////////////////////


//==============================================================================
// Headers

#include "Waveform_Analysis.hpp"

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

/*==============================================================================
Baseline correction using point from 0 to N.

  std::vector<double>& voltageVec := (stl) vector contains the voltage.
  int ptN := point N

  return: None.
=============================================================================*/
void Correct_Baseline(
  std::vector<double>& voltageVec,
  int ptN)
{
  double mean =0;

  for(std::size_t j = 0, max = ptN; j < max; j++){mean += voltageVec.at(j);}

  mean = mean/ptN;

  for(std::size_t j = 0, max = voltageVec.size(); j < max; j++){voltageVec.at(j) = voltageVec.at(j)- mean;}
}

/*==============================================================================
Alternate implementation of baseline correction.(2)

  std::vector<double>& voltageVec := (stl) vector contains the voltage.
  double fractional_pts := using fraction of the points for baseline correction.

  return: None.
==============================================================================*/
void Correct_Baseline2( std::vector<double>& voltageVec, double fractional_pts )
{
  double mean =0;

  int noise_pts = fractional_pts * voltageVec.size();

  for(std::size_t j = 0, max = noise_pts; j < max; j++){mean += voltageVec.at(j);}

  mean = mean/noise_pts;

  for(std::size_t j = 0, max = voltageVec.size(); j < max; j++){ voltageVec.at(j) = voltageVec.at(j)- mean; }
}

/*==============================================================================
Alternate implementation of baseline correction.(2)

  std::vector<double>& voltageVec := (stl) vector contains the voltage.
  double fractional_pts := using fraction of the points for baseline correction.

  return: None.
==============================================================================*/
void Correct_Baseline3(
  std::vector<double>& voltageVec,
  std::vector<double> timeVec,
  double tRange[2]
)
{
  double mean =0;
  int counter =0;
  std::size_t npoints = voltageVec.size();

  for(std::size_t j = 0, max = npoints; j < max; j++)
  {
    if( timeVec.at(j) >= tRange[0] && timeVec.at(j) <= tRange[1] )
    {
      mean += voltageVec.at(j);
      counter++;
    }
  }

  if(counter != 0)
  {
    mean = mean/counter;
    for(std::size_t j = 0; j < npoints; j++){ voltageVec.at(j) = voltageVec.at(j)- mean;}
  }
  else{ Correct_Baseline2(voltageVec, 0.3); }
}

//==============================================================================
//==============================================================================
//find the baseline correction for SSRL data.

void SSRL_Baseline( std::vector<double> &w, std::string workerID )
{
  //TThread::Lock();
  std::string histoName = workerID+"voltage_histogram";
  std::string fitName = workerID+"gausFit";
  delete gROOT->FindObject(histoName.c_str());
  TH1D *voltage_histogram = new TH1D( histoName.c_str(), "", 100, 1, 1);
  for( std::size_t i = 0, max = w.size(); i < max; i++ )
  {
    voltage_histogram->Fill( w.at(i) );
  }
  TThread::Lock();
  TF1 *gausFit = new TF1( fitName.c_str(), "gaus" );
  gausFit->SetParameter( 1, voltage_histogram->GetBinCenter(voltage_histogram->GetMaximumBin() ) );
  voltage_histogram->Fit( gausFit, "Q0" );
  double base_line_correction = gausFit->GetParameter( 1 );
  delete voltage_histogram;
  //TThread::UnLock();
  for( std::size_t j = 0, max = w.size(); j < max; j++){ w.at(j) = w.at(j)- base_line_correction; }
}

/*==============================================================================
find the baseline correction for SSRL data. Alternative implementation.
==============================================================================*/
void SSRL_Baseline( std::vector<double> &w, double &RMS )
{
  TThread::Lock();
  delete gROOT->FindObject("voltage_histogram");
  TH1D *voltage_histogram = new TH1D( "voltage_histogram", "", 100, 1, 1);
  for( std::size_t i = 0, max = w.size(); i < max; i++ )
  {
    voltage_histogram->Fill( w.at(i) );
  }
  double sample_mean = voltage_histogram->GetMean(1);
  double sample_sigma = voltage_histogram->GetStdDev(1);
  double peak = voltage_histogram->GetBinCenter(voltage_histogram->GetMaximumBin() );
  double fwhm = 2.32 * sample_sigma;
  TF1 *gausFit = new TF1( "gausFit", "gaus", peak - fwhm/2.0, peak + fwhm/2.0 );
  gausFit->SetParameter( 1, peak );
  voltage_histogram->Fit( "gausFit", "Q0R" );
  double base_line_correction = gausFit->GetParameter( 1 );
  RMS = gausFit->GetParameter(2);
  delete voltage_histogram;
  TThread::UnLock();
  for( std::size_t j = 0, max = w.size(); j < max; j++){ w.at(j) = w.at(j)- base_line_correction; }
}


/*==============================================================================
SSRL noise correction and find the nosie.
This is dynamic method that looks for empty gap between bunches
in order to reduce the influence from the previous bunches
  param

===============================================================================*/
void SSRL_Dynamic_Noise_And_Baseline(
  std::vector<double> &w, //voltage
  std::vector<double> t, //time
  double assisting_threshold,
  const double time_separation,
  double &RMS
)
{
  std::vector<double> multiple_singal_pmax_v = {};
  multiple_singal_pmax_v.reserve(t.size());
  std::vector<double> multiple_singal_tmax_v = {};
  multiple_singal_tmax_v.reserve(t.size());
  std::vector<int> indexing_v;
  indexing_v.reserve(t.size());
  std::vector<double> noise_pt;
  noise_pt.reserve(t.size());

  std::size_t npoints = w.size();

  //std::cout<<w.at(0)<<"  sssshhhhh\n";
  double pmax = 0.0;
  unsigned int pmax_index = 0;

  bool candidate_signal = false;
  bool noisy_event = true;

  int fail_counter = 0;

  //double ss_pt;
  //double ee_pt;

  for( unsigned int i = 0; i < npoints; i++ )
  {
    if( !candidate_signal )
    {
      if( w.at(i) >= assisting_threshold )
      {
        if( w.at(i) > pmax )
        {
          pmax = w.at(i);
          pmax_index = i;
          candidate_signal = true;
          //if( noisy_event ) noisy_event = false;
          //std::cout<<"pass noisy \n";
        }
      }
    }
    else
    {
      if( w.at(i) > pmax )
      {
        pmax = w.at(i);
        pmax_index = i;
      }
      else if( (w.at(i) < assisting_threshold) && (assisting_threshold - w.at(i)) <= (assisting_threshold) )
      {
        multiple_singal_pmax_v.push_back( pmax );
        multiple_singal_tmax_v.push_back( t.at(pmax_index) );
        indexing_v.push_back( pmax_index );
        pmax = w.at(i);
        pmax_index = i;
        candidate_signal = false;
        //std::cout<<"pass noisy 22 \n";
        if( noisy_event ) noisy_event = false;
      }
      else{}
    }
  }

  if( noisy_event || multiple_singal_tmax_v.size()<2 )
  {
    if( noisy_event )
    {
      //std::cout<<"noisey ";
      //std::cout<<multiple_singal_tmax_v.size() << std::endl;
      for( std::size_t i = 0, max = t.size(); i < max; i++ )
      {
        //std::cout<<"swhat " << max << "\n";
        noise_pt.push_back(w.at(i));
      }
    }
    else if( multiple_singal_tmax_v.size() == 1 )
    {
      //std::cout<<"not noisey \n";
      //unsigned int check = 0;
      double start_time = t.at(0);
      double end_time = start_time+5000.0;
      double max_time = t.at(t.size()-1);
      double increment_time = 5000.0;
      bool _fail_ = true;

      while( true )
      {
        if( ( ( (multiple_singal_tmax_v.at(0)+5000) > start_time ) && ( (multiple_singal_tmax_v.at(0)-5000) < start_time ) ) || ( ( (multiple_singal_tmax_v.at(0)+5000) > end_time ) && ( (multiple_singal_tmax_v.at(0)-5000) < end_time ) ) )
        {
          //std::cout<<"incre\n";
          start_time = start_time + increment_time;
          end_time = end_time + increment_time;
          if( end_time >= max_time ) end_time = max_time;
          if( start_time >= max_time ) break;
        }
        else
        {
          //ss_pt = start_time;
          //ee_pt = end_time;
          _fail_ = false;
          //std::cout<< "dd " << start_time << std::endl;
          //std::cout<< "dd--s " << end_time << std::endl;
          break;
        }
      }
      if( !_fail_ )
      {
        //std::cout<<"not fail "<< start_time<< " end " << end_time << std::endl;
        for( std::size_t i = 0, max = t.size(); i < max; i++ )
        {
          if( t.at(i) >= start_time && t.at(i) <= end_time ) noise_pt.push_back(w.at(i));
        }
      }
      else
      {
        //std::cout<<"fail"<<std::endl;
        for( std::size_t i = 0, max = t.size(); i < max; i++ )
        {
          if( w.at(i) < assisting_threshold ) noise_pt.push_back(w.at(i));
        }
      }
    }
    else{}

    double mean = std::accumulate(std::begin(noise_pt), std::end(noise_pt), 0.0 )/noise_pt.size();
    //std::cout<<"Mean " << mean << " n "<< double(noise_pt.size())<<  std::endl;
    double accum = 0.0;
    for( std::size_t i = 0, max = noise_pt.size(); i < max; i++ )
    {
      accum += (noise_pt.at(i)-mean)*(noise_pt.at(i)-mean);
    }
    double noise = std::sqrt(accum/noise_pt.size());
    RMS = noise;
    //if(RMS > 20.0)std::cout<< "why? " << RMS << " is big " << ss_pt << " " << ee_pt << " "  << std::endl;
    //if(indexing_v.size()!=0)std::cout<<"size " << multiple_singal_tmax_v.at(0) << " " << multiple_singal_pmax_v.at(0) << std::endl;
    for( std::size_t i = 0, max = w.size(); i < max; i++ )
    {
      w.at(i) = w.at(i) - mean;
    }
  }
  else
  {
    bool satisfied = false;
    for( std::size_t i = 0, max = multiple_singal_pmax_v.size()-1; i < max; i++ )
    {
      if( (multiple_singal_tmax_v.at(i+1)-multiple_singal_tmax_v.at(i)) >=  time_separation )
      {
        std::pair<double,int> sec_pmax = std::make_pair( multiple_singal_pmax_v.at(i+1), indexing_v.at(i+1) );
        double rise = Find_Rise_Time( w, t, sec_pmax, 0.1, 0.9 );

        for( std::size_t j=indexing_v.at(i), jmax = indexing_v.at(i+1); j<jmax; j++ )
        {
          //ss_pt = multiple_singal_tmax_v.at(i)+5000.0;
          //ee_pt = multiple_singal_tmax_v.at(i+1)-(rise/0.9);
          if( t.at(j) >= multiple_singal_tmax_v.at(i)+5000.0 && t.at(j) <= multiple_singal_tmax_v.at(i+1)-(rise/0.9+500) )
          {
            noise_pt.push_back( w.at(j) );
          }
        }
        satisfied = true;
        break;
      }
    }

    if( !satisfied )
    {
      unsigned int check = 0;
      double start_time = t.at(0);
      double end_time = start_time+5000.0;
      double max_time = t.at(t.size()-1);
      double increment_time = 5000.0;
      double selection_width = 5000;
      bool _fail_ = false;
      //std::cout<<"not satifsid\n";
      //std::cin.get();
      while(true)
      {
        //if( start_pt <= indexing_v.at(check) && indexing_v.at(check) <= end_pt )
        if( check == multiple_singal_tmax_v.size() ) break;
        if( ( ( (multiple_singal_tmax_v.at(check)+selection_width) > start_time ) && ( (multiple_singal_tmax_v.at(check)-selection_width) < start_time ) ) || ( ( (multiple_singal_tmax_v.at(check)+selection_width) > end_time ) && ( (multiple_singal_tmax_v.at(check)-selection_width) < end_time ) ) )
        {
          start_time = start_time + increment_time;
          end_time = end_time + increment_time;
          //ss_pt = start_time;
          //ee_pt = end_time;
          if( end_time >= max_time ) end_time = max_time;
          if( start_time >= max_time )
          {
            if(fail_counter != 4 )
            {
              fail_counter++;
              increment_time = increment_time - 500;
              start_time = t.at(0);
              end_time = start_time + increment_time;
              selection_width = selection_width - 500;
            }
            else{ _fail_ = true; break;}
          }
        }
        else check++;
      }

      if( _fail_ )
      {
        for( std::size_t i = 0, max = w.size(); i < max; i++ )
        {
          if( w.at(i) < assisting_threshold ) noise_pt.push_back(w.at(i));
        }
      }
      else
      {
        for( std::size_t i = 0, max = t.size(); i < max; i++ )
        {
          if( t.at(i) >= start_time && t.at(i) <= end_time ) noise_pt.push_back(w.at(i));
        }
      }
    }
    else{}


    double mean = std::accumulate(noise_pt.begin(), noise_pt.end(), 0.0 )/noise_pt.size();
    double accum = 0.0;
    //std::cout<<"stop "<<mean<< "  "<< noise_pt.size()<<  " \n";
    for( std::size_t i = 0, max = noise_pt.size(); i < max; i++ )
    {
      accum += (noise_pt.at(i)-mean)*(noise_pt.at(i)-mean);
    }
    double noise = std::sqrt(accum/noise_pt.size());
    RMS = noise;
    /*if(RMS > 0.0)
    {
      std::cout<< "why? " << noise << " is big " << multiple_singal_tmax_v.size() << " ss" << ss_pt << " ee " << ee_pt << " "<< noise_pt.size()<< std::endl;
      for(std::size_t ii = 0, mm = multiple_singal_tmax_v.size(); ii<mm;ii++ )std::cout<<" ssees " << multiple_singal_tmax_v.at(ii) << std::endl;
    }*/

    //std::cout<<"stop22 "<<mean<< "  "<< noise_pt.size()<<  " \n";
    for( std::size_t i = 0, max = w.size(); i < max; i++ )
    {
      w.at(i) = w.at(i) - mean;
    }
  }
}
