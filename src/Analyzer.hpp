#ifndef ANALYZER_H
#define ANALYZER_H

#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <stdio.h>

#include <Riostream.h>
#include <TTreeReader.h>
#include "TTreeReaderValue.h"
#include "TTreeReaderArray.h"
#include "TObject.h"


class Analyzer : public TObject
{

 public:

  Analyzer(std::vector<double> voltage, std::vector<double> time);
  Analyzer();
  Analyzer(const Analyzer &a);
  
  virtual ~Analyzer();

  
  // Baseline Correction
double Correct_Baseline( int ptN );
//void Correct_Baseline2( std::vector<double> &voltageVec, double fractional_pts );
//void Correct_Baseline3( std::vector<double> &voltageVec, std::vector<double> timeVec, double tRange[2] ); //special treatment of ill-signal baseline.

//==========================================================================
// Pulse Maximum (Pmax)
std::pair<double, unsigned int> Find_Signal_Maximum(bool confineSearchRegion, double searchRange[2]);
std::pair<double, double> Pmax_with_GausFit(const std::pair<double, unsigned int> Pmax, unsigned int maxIndex);
std::pair<double, unsigned int> Find_Negative_Signal_Maximum(bool confineSearchRegion,double searchRange[2]);
std::pair<double, double> Negative_Pmax_with_GausFit(const std::pair<double, unsigned int> NegPmax, unsigned int maxIndex);

// Time corresponding to Pmax
double Get_Tmax(const std::pair<double, unsigned int> Pmax);
double Get_Negative_Tmax(const std::pair<double, unsigned int> NegPmax);

//==========================================================================
// Pulse Area
double DC_Area(double baseline_correction);
double Find_Pulse_Area(const std::pair<double,unsigned int> Pmax);
double Find_Undershoot_Area(const std::pair<double, unsigned int> Pmax);

double Pulse_Integration_with_Fixed_Window_Size(const std::pair<double,unsigned int> pmax_holder, std::string integration_option , double t_beforeSignal, double t_afterSignal);

double Pulse_Integration_with_Fixed_Window_Size_with_GausFit(const std::pair<double,double> pmax_holder, unsigned int imax, std::string integration_option , double t_beforeSignal, double t_afterSignal);

double Pulse_Area_With_Linear_Interpolate_Edge(const std::pair<double,unsigned int> Pmax, std::string integration_option, bool relativeTimeWindow, double stopTime);

double Pulse_Area_With_Linear_Interpolate_Edge_with_GausFit(const std::pair<double,double> Pmax, unsigned int imax, std::string integration_option, bool relativeTimeWindow, double stopTime);

double New_Pulse_Area(const std::pair<double,double> Pmax, unsigned int imax, std::string integration_option, double range[2]);

double New_Undershoot_Area(const std::pair<double,double> Pmax, const std::pair<double,double> Pmin, unsigned int imin, std::string integration_option, double range[2]);


//==========================================================================
// CFD
double Rising_Edge_CFD_Time(const double fraction, const std::pair<double,unsigned int> Pmax);
double Rising_Edge_CFD_Time_with_GausFit(const double fraction, const std::pair<double,double> Pmax, unsigned int imax);
double Falling_Edge_CFD_Time_with_GausFit(const double fraction, const std::pair<double,double> Pmax, unsigned int imax);

//==========================================================================
// Rise Time
double Find_Rise_Time(const std::pair<double, unsigned int> Pmax, double bottom = 0.1, double top = 0.9);
double Find_Rise_Time_with_GausFit(const std::pair<double, double> Pmax, unsigned int imax, double bottom = 0.1, double top = 0.9);
double Find_Fall_Time_with_GausFit(const std::pair<double, double> Pmax, unsigned int imax, double bottom = 0.1, double top = 0.9);

//==========================================================================
// Noise
double Find_Noise(const unsigned int inoise);

/*double Find_Noise2(
  std::vector<double> voltageVec,
  const double fractional_pts
);

double Find_Noise_On_Back_Baseline(
  std::vector<double> voltageVec,
  std::vector<double> timeVec,
  double start_time,
  double end_time
);*/

//==========================================================================
// Dvdt

double Find_Dvdt(const int fraction, const int ndif, const std::pair<double,unsigned int> Pmax);
double Find_Dvdt_with_GausFit(const int fraction, const int ndif, const std::pair<double,double> Pmax, unsigned int imax);
double Find_Dvdt2080_with_GausFit(const int ndif, const std::pair<double,double> Pmax, unsigned int imax);

//==========================================================================
// Multiple Peak

/*std::pair <double, unsigned int> Find_Identical_Peak(
  std::vector<double> voltageVec,
  std::vector<double> timeVec,
  unsigned int StartIndex,
  bool limitSearchRegion = false,
  double min_search_range = -1000.0,
  double max_search_range = 1000.0
);

void Get_PmaxTmax_Of_Multiple_Singal(
  const double        assist_threshold,
  std::vector<double> voltageVec,
  std::vector<double> timeVec,
  std::vector<double> &multiple_singal_pmax_v,
  std::vector<double> &multiple_singal_tmax_v,
  std::vector<int>    &indexing_v
);

int Signal_Peak_Counter(
  std::vector<double> &voltageVec,
  std::vector<double> timeVec,
  double assisting_threshold
);*/

//==========================================================================
// time at threshold
double Find_Time_At_Threshold_with_GausFit(const double thresholdLevel, const std::pair<double,double> Pmax, unsigned int imax);
double Find_Time_At_Threshold_Falling_Edge_with_GausFit(const double thresholdLevel, const std::pair<double,double> Pmax, unsigned int imax);
double Find_Time_Over_Threshold(const double thresholdLevel, const std::pair<double,unsigned int> Pmax, const double second_thresholdLevel);

/* void Get_TimeAcrossThreshold(
  const double        thresholdLevel,
  std::vector<double> voltageVec,
  std::vector<double> timeVec,
  std::vector<double> &time_at_threshold_v,
  const unsigned int  expect_count = 6
);*/

  
 private:

 std::vector<double> pvoltage;
 std::vector<double> ptime;
  
  ClassDef(Analyzer,1)
  
};

#endif
