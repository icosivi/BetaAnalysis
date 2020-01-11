#include "../include/BetaConfig.hpp"

#include <iostream>
#include <string>
#include <sstream>
#include <fstream>
#include <stdlib.h>
#include <stdio.h>
#include <vector>

#include <glib-2.0/glib.h>
#include <glib-2.0/glib/gprintf.h>
//#include <glib.h>
//#include <glib/gprintf.h>
#include <gtk/gtk.h>

void BetaConfig::ReadConfig( std::string configFile )
{
  GKeyFile *gkf;
  gkf = g_key_file_new();
  if( ! g_key_file_load_from_file( gkf, configFile.c_str(), G_KEY_FILE_NONE, NULL) ) std:: cout << "cannot read the configure file. \n";
  else
  {
    this->voltage_scalar   = g_key_file_get_double( gkf, this->header_group, this->voltage_scalar_key, NULL);
    this->time_scalar      = g_key_file_get_double( gkf, this->header_group, this->time_scalar_key, NULL);
    //this->num_core         = g_key_file_get_int64( gkf, this->header_group, this->num_core_key, NULL );
    this->enable_channel_1 = g_key_file_get_boolean( gkf, this->active_channel_group, this->ch1_key, NULL);
    this->enable_channel_2 = g_key_file_get_boolean( gkf, this->active_channel_group, this->ch2_key, NULL);
    this->enable_channel_3 = g_key_file_get_boolean( gkf, this->active_channel_group, this->ch3_key, NULL);
    this->enable_channel_4 = g_key_file_get_boolean( gkf, this->active_channel_group, this->ch4_key, NULL);
    /*this->enable_channel_5 = g_key_file_get_boolean( gkf, this->active_channel_group, this->ch5_key, NULL);
    this->enable_channel_6 = g_key_file_get_boolean( gkf, this->active_channel_group, this->ch6_key, NULL);
    this->enable_channel_7 = g_key_file_get_boolean( gkf, this->active_channel_group, this->ch7_key, NULL);
    this->enable_channel_8 = g_key_file_get_boolean( gkf, this->active_channel_group, this->ch8_key, NULL);*/

    this->invert_channel_1 = g_key_file_get_boolean( gkf, this->invert_signal_group, this->ch1_key, NULL);
    this->invert_channel_2 = g_key_file_get_boolean( gkf, this->invert_signal_group, this->ch2_key, NULL);
    this->invert_channel_3 = g_key_file_get_boolean( gkf, this->invert_signal_group, this->ch3_key, NULL);
    this->invert_channel_4 = g_key_file_get_boolean( gkf, this->invert_signal_group, this->ch4_key, NULL);
    /*this->invert_channel_5 = g_key_file_get_boolean( gkf, this->invert_signal_group, this->ch5_key, NULL);
    this->invert_channel_6 = g_key_file_get_boolean( gkf, this->invert_signal_group, this->ch6_key, NULL);
    this->invert_channel_7 = g_key_file_get_boolean( gkf, this->invert_signal_group, this->ch7_key, NULL);
    this->invert_channel_8 = g_key_file_get_boolean( gkf, this->invert_signal_group, this->ch8_key, NULL);*/

    //this->parallel_switch = g_key_file_get_boolean( gkf, this->parallel_group, this->parallel_key, NULL);

    //this->directory_location = g_key_file_get_string( gkf, this->directory_group, this->directory_key, NULL);

    //this->ScopeTimestamp = g_key_file_get_boolean( gkf, this->header_group, this->timestamp_key, NULL );

    this->trimTimeWindow = g_key_file_get_boolean( gkf, this->header_group, this->trimTimeWindow_key, NULL );
    this->trimTimeRangeMin = g_key_file_get_double( gkf, this->header_group, this->trimTimeRangeMin_key, NULL );
    this->trimTimeRangeMax = g_key_file_get_double( gkf, this->header_group, this->trimTimeRangeMax_key, NULL );

    this->limiting_Pmax_search_window = g_key_file_get_boolean( gkf, this->header_group, this->limiting_Pmax_search_window_key, NULL );
    this->Pmax_search_MinRange = g_key_file_get_double( gkf, this->header_group, this->Pmax_search_MinRange_key, NULL );
    this->Pmax_search_MaxRange = g_key_file_get_double( gkf, this->header_group, this->Pmax_search_MaxRange_key, NULL );
  }
}

/*void BetaConfig::CreateConfig()
{
  std::ofstream beta_config;
  beta_config.open("beta_config_v20180810.ini");

  beta_config << "[" << this->header_group << "]" << std::endl;
  beta_config << this->voltage_scalar_key << "=1000.0" << std::endl;
  beta_config << this->time_scalar_key << "=1.0e12" << std::endl;
  beta_config << this->num_core_key << "=16" << std::endl;
  beta_config << "SLAC_raw_data=false" << std::endl;
  beta_config << "add_scope_timestamp=false" << std::endl;
  beta_config << "trimTimeWindow=false" << std::endl;
  beta_config << "trimTimeRangeMin=-2000" << std::endl;
  beta_config << "trimTimeRangeMax=2000" << std::endl;
  beta_config << this->limiting_Pmax_search_window_key << "=false" << std::endl;
  beta_config << this->Pmax_search_MinRange_key << "=-1000" << std::endl;
  beta_config << this->Pmax_search_MaxRange_key << "=1000" << std::endl;

  beta_config << std::endl;

  beta_config << "[" << this->active_channel_group << "]" << std::endl;
  beta_config << this->ch1_key << "=false" << std::endl;
  beta_config << this->ch2_key << "=true"  << std::endl;
  beta_config << this->ch3_key << "=true"  << std::endl;
  beta_config << this->ch4_key << "=false" << std::endl;
  beta_config << this->ch5_key << "=false"  << std::endl;
  beta_config << this->ch6_key << "=false"  << std::endl;
  beta_config << this->ch7_key << "=false"  << std::endl;
  beta_config << this->ch8_key << "=false"  << std::endl;

  beta_config << std::endl;

  beta_config << "[" << this->invert_signal_group << "]" << std::endl;
  beta_config << this->ch1_key << "=false"  << std::endl;
  beta_config << this->ch2_key << "=true"   << std::endl;
  beta_config << this->ch3_key << "=false"  << std::endl;
  beta_config << this->ch4_key << "=false"  << std::endl;
  beta_config << this->ch5_key << "=false"  << std::endl;
  beta_config << this->ch6_key << "=false"  << std::endl;
  beta_config << this->ch7_key << "=false"  << std::endl;
  beta_config << this->ch8_key << "=false"  << std::endl;

  beta_config << std::endl;

  beta_config << "[" << this->parallel_group << "]" << std::endl;
  beta_config << this->parallel_key << "=true"  << std::endl;

  beta_config << std::endl;
  beta_config << "[" << this->SSRL_group << "]" << std::endl;
  beta_config << this->SSRL_enable_key << "=false"  << std::endl;
  beta_config << this->SSRL_assist_threshold_key << "=15.0"  << std::endl;
  beta_config << this->SSRL_AC_position_key << "=false" << std::endl;
  beta_config << this->SSRL_Dynamic_Noise_Enable_key <<"=false" << std::endl;
  beta_config << this->SSRL_Brute_Forced_Noise_Corr_key <<"=false" << std::endl;

  beta_config << std::endl;

  beta_config << "[" << this->directory_group << "]" << std::endl;
  beta_config << this->directory_key << "=raw/"  << std::endl;

  beta_config.close();
}*/
