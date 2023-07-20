#ifndef BETACONFIG_H
#define BETACONFIG_H

#include <glib-2.0/glib.h>
#include <glib-2.0/glib/gprintf.h>
//#include <glib.h>
//#include <glib/gprintf.h>
#include <gtk/gtk.h>
#include <stdio.h>
#include <stdlib.h>
#include <iostream>

struct BetaConfig
{
  //GKeyFile *gkf;
  //gkf = g_key_file_new();

  const char* header_group          = "HEADER";
  const char* active_channel_group  = "ACTIVE_CHANNEL";
  const char* invert_signal_group   = "INVERT_SIGNAL";
  //const char* directory_group       = "DIRECTORY";

  const char* voltage_scalar_key = "voltage_scalar";
  const char* time_scalar_key    = "time_scalar";

  const char* ch1_key = "ch1";
  const char* ch2_key = "ch2";
  const char* ch3_key = "ch3";
  const char* ch4_key = "ch4";
  /*const char* ch5_key = "ch5";
  const char* ch6_key = "ch6";
  const char* ch7_key = "ch7";
  const char* ch8_key = "ch8";*/

  //const char* directory_key = "location";
  //const char* parallel_key  = "parallel";

  //const char* timestamp_key = "add_scope_timestamp";

  const char* trimTimeWindow_key = "trimTimeWindow";
  const char* trimTimeRangeMin_key = "trimTimeRangeMin";
  const char* trimTimeRangeMax_key = "trimTimeRangeMax";

  const char* limiting_Pmax_search_window_key = "limiting_Pmax_search_window";
  const char* Pmax_search_MinRange_key = "Pmax_search_MinRange";
  const char* Pmax_search_MaxRange_key = "Pmax_search_MaxRange";

  double time_scalar;
  double voltage_scalar;

  //bool ScopeTimestamp;

  bool enable_channel_1;
  bool enable_channel_2;
  bool enable_channel_3;
  bool enable_channel_4;
  /*bool enable_channel_5;
  bool enable_channel_6;
  bool enable_channel_7;
  bool enable_channel_8;*/

  bool invert_channel_1;
  bool invert_channel_2;
  bool invert_channel_3;
  bool invert_channel_4;
  /*bool invert_channel_5;
  bool invert_channel_6;
  bool invert_channel_7;
  bool invert_channel_8;*/

  //std::string directory_location;

  bool trimTimeWindow;
  double trimTimeRangeMin;
  double trimTimeRangeMax;

  bool limiting_Pmax_search_window;
  double Pmax_search_MinRange;
  double Pmax_search_MaxRange;

  void ReadConfig( std::string configFile );
  //void CreateConfig();
};

#endif //BETACONFIG_H
