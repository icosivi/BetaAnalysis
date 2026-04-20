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

  Analyzer(std::vector<float> voltage, std::vector<float> time);
  Analyzer();
  Analyzer(const Analyzer &a);
  
  virtual ~Analyzer();

  
  // Baseline Correction
float Correct_Baseline( int ptN );
//void Correct_Baseline2( std::vector<float> &voltageVec, float fractional_pts );
//void Correct_Baseline3( std::vector<float> &voltageVec, std::vector<float> timeVec, float tRange[2] ); //special treatment of ill-signal baseline.

//==========================================================================
// Pulse Maximum (Pmax)
std::pair<float, unsigned int> Find_Signal_Maximum(bool confineSearchRegion, float searchRange[2]);
std::pair<float, float> Pmax_with_GausFit(const std::pair<float, unsigned int> Pmax, unsigned int maxIndex);
std::pair<float, unsigned int> Find_Negative_Signal_Maximum(bool confineSearchRegion,float searchRange[2]);
std::pair<float, float> Negative_Pmax_with_GausFit(const std::pair<float, unsigned int> NegPmax, unsigned int maxIndex);

// Time corresponding to Pmax
float Get_Tmax(const std::pair<float, unsigned int> Pmax);
float Get_Negative_Tmax(const std::pair<float, unsigned int> NegPmax);

//==========================================================================
// Pulse Area
float DC_Area(float baseline_correction);
float Find_Pulse_Area(const std::pair<float,unsigned int> Pmax);
float Find_Undershoot_Area(const std::pair<float, unsigned int> Pmax);

float Pulse_Integration_with_Fixed_Window_Size(const std::pair<float,unsigned int> pmax_holder, std::string integration_option , float t_beforeSignal, float t_afterSignal);

float Pulse_Integration_with_Fixed_Window_Size_with_GausFit(const std::pair<float,float> pmax_holder, unsigned int imax, std::string integration_option , float t_beforeSignal, float t_afterSignal);

float Pulse_Area_With_Linear_Interpolate_Edge(const std::pair<float,unsigned int> Pmax, std::string integration_option, bool relativeTimeWindow, float stopTime);

float Pulse_Area_With_Linear_Interpolate_Edge_with_GausFit(const std::pair<float,float> Pmax, unsigned int imax, std::string integration_option, bool relativeTimeWindow, float stopTime);

float New_Pulse_Area(const std::pair<float,float> Pmax, unsigned int imax, std::string integration_option, float range[2]);

float New_Undershoot_Area(const std::pair<float,float> Pmax, const std::pair<float,float> Pmin, unsigned int imin, std::string integration_option, float range[2]);


//==========================================================================
// CFD
float Rising_Edge_CFD_Time(const float fraction, const std::pair<float,unsigned int> Pmax);
float Rising_Edge_CFD_Time_with_GausFit(const float fraction, const std::pair<float,float> Pmax, unsigned int imax);
float Falling_Edge_CFD_Time_with_GausFit(const float fraction, const std::pair<float,float> Pmax, unsigned int imax);

//==========================================================================
// Rise Time
float Find_Rise_Time(const std::pair<float, unsigned int> Pmax, float bottom = 0.1, float top = 0.9);
float Find_Rise_Time_with_GausFit(const std::pair<float, float> Pmax, unsigned int imax, float bottom = 0.1, float top = 0.9);
float Find_Rise_Time_with_RELU_fit(const std::pair<float, float> Pmax, unsigned int imax, float bottom = 0.1, float top = 0.9);
float Find_Rise_Time_with_LinFit_Rob(const std::pair<float, float> Pmax, unsigned int imax, float bottom = 0.1, float top = 0.9);
float Find_Fall_Time_with_GausFit(const std::pair<float, float> Pmax, unsigned int imax, float bottom = 0.1, float top = 0.9);

//==========================================================================
// Noise
float Find_Noise(const unsigned int inoise);

/*float Find_Noise2(
  std::vector<float> voltageVec,
  const float fractional_pts
);

float Find_Noise_On_Back_Baseline(
  std::vector<float> voltageVec,
  std::vector<float> timeVec,
  float start_time,
  float end_time
);*/

//==========================================================================
// Dvdt

float Find_Dvdt(const int fraction, const int ndif, const std::pair<float,unsigned int> Pmax);
float Find_Dvdt_with_GausFit(const int fraction, const int ndif, const std::pair<float,float> Pmax, unsigned int imax);
float Find_Dvdt2080_with_GausFit(const int ndif, const std::pair<float,float> Pmax, unsigned int imax);

//==========================================================================
// Multiple Peak

/*std::pair <float, unsigned int> Find_Identical_Peak(
  std::vector<float> voltageVec,
  std::vector<float> timeVec,
  unsigned int StartIndex,
  bool limitSearchRegion = false,
  float min_search_range = -1000.0,
  float max_search_range = 1000.0
);

void Get_PmaxTmax_Of_Multiple_Singal(
  const float        assist_threshold,
  std::vector<float> voltageVec,
  std::vector<float> timeVec,
  std::vector<float> &multiple_singal_pmax_v,
  std::vector<float> &multiple_singal_tmax_v,
  std::vector<int>    &indexing_v
);

int Signal_Peak_Counter(
  std::vector<float> &voltageVec,
  std::vector<float> timeVec,
  float assisting_threshold
);*/

//==========================================================================
// time at threshold
float Find_Time_At_Threshold_with_GausFit(const float thresholdLevel, const std::pair<float,float> Pmax, unsigned int imax);
float Find_Time_At_Threshold_Falling_Edge_with_GausFit(const float thresholdLevel, const std::pair<float,float> Pmax, unsigned int imax);
float Find_Time_Over_Threshold(const float thresholdLevel, const std::pair<float,unsigned int> Pmax, const float second_thresholdLevel);

/* void Get_TimeAcrossThreshold(
  const float        thresholdLevel,
  std::vector<float> voltageVec,
  std::vector<float> timeVec,
  std::vector<float> &time_at_threshold_v,
  const unsigned int  expect_count = 6
);*/

//==========================================================================
// getters
inline std::vector<float> const& getVoltages() const { return pvoltage; }
inline std::vector<float> const& getTimes() const { return ptime; }
  
 private:

 std::vector<float> pvoltage;
 std::vector<float> ptime;
  
  ClassDef(Analyzer,1)
  
};

#endif
