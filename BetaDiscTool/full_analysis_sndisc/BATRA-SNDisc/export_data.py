# export_data.py

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import minimize
from scipy.stats import poisson
import ROOT as root
from ROOT import TF1
from scipy.special import gammaln
import math
from math import exp, sqrt, pi
import pandas as pd
import argparse
import glob
import re
import os
import csv
import math
import sys

from firstPass_ProcTools import getBias

def direct_to_table(name_and_df_couples, channel_configs, output_savename, thickness_info):

  number_of_duts = sum(1 for element in channel_configs if element[0] == 1)
  number_of_bias_pts = int(len(name_and_df_couples[0][1]) / number_of_duts)
  thickness_info = list(map(int, thickness_info))
  thickness_col = np.repeat(np.array(thickness_info), number_of_bias_pts)

  pmax_low = []
  pmax_high = []
  nmax_low = []
  tmax_low = []
  tmax_high = []

  pmax_low_mcp = []
  pmax_high_mcp = []
  nmax_low_mcp = []
  tmax_low_mcp = []
  tmax_high_mcp = []

  mcp_channel = False
  for T, _, S in channel_configs:
    if T == 2:
      mcp_channel = True
      pmax_low_mcp = S[0]
      pmax_high_mcp = S[1]

  dfs_to_concat = []
  first_df_found = False
  for index, (var, df) in enumerate(name_and_df_couples):
    if var == "amplitude":
      df_ampl = df[['Channel','Bias','Amplitude MPV']]
      first_df_found = True
      df_ampl.loc[:, 'Amplitude MPV'] = df_ampl['Amplitude MPV'].round(1)
      df_ampl = df_ampl.rename(columns={'Amplitude MPV':'Amplitude / mV'})
      dfs_to_concat.append(df_ampl)
    elif var == "risetime":
      if first_df_found == False:
        first_df_found = True
        df_rt = df[['Channel','Bias','Mean','Sigma']]
      else:
        df_rt = df[['Mean','Sigma']]
      df_rt.loc[:, 'Mean'] = (1000*df_rt['Mean']).round(0)
      df_rt.loc[:, 'Sigma'] = (1000*df_rt['Sigma']).round(0)
      df_rt = df_rt.rename(columns={'Mean':'Rise time / ps','Sigma':'Rise time Unc / ps'})
      dfs_to_concat.append(df_rt)
    elif var == "charge":
      if first_df_found == False:
        first_df_found = True
        df_charge = df[['Channel','Bias','Charge MPV','Landau width','Gaussian sigma','Frac above 1p5 MPV', 'Frac above 1p5 Qmax']]
      else:
        df_charge = df[['Charge MPV','Landau width','Gaussian sigma','Frac above 1p5 MPV', 'Frac above 1p5 Qmax']]
      df_charge.loc[:, 'Charge MPV'] = df_charge['Charge MPV'].round(1)
      df_charge.loc[:, 'Landau width'] = df_charge['Landau width'].round(3)
      df_charge.loc[:, 'Gaussian sigma'] = df_charge['Gaussian sigma'].round(3)
      df_charge.loc[:, 'Frac above 1p5 MPV'] = df_charge['Frac above 1p5 MPV'].round(3)
      df_charge.loc[:, 'Frac above 1p5 Qmax'] = df_charge['Frac above 1p5 Qmax'].round(3)
      df_charge = df_charge.rename(columns={'Charge MPV':'Charge / fC','Landau width':'Landau Cpt Charge','Gaussian sigma':'Gaussian Cpt Charge','Frac above 1p5 MPV':'Frac Charge >1.5xMPV','Frac above 1p5 Qmax':'Frac Charge >1.5xQmax'})
      df_charge['Gain'] = 100*(df_charge['Charge / fC'] / thickness_col).round(2)
      dfs_to_concat.append(df_charge)
    elif var == "rms":
      if first_df_found == False:
        first_df_found = True
        df_rms = df[['Channel','Bias','Mean','Sigma']]
      else:
        df_rms = df[['Mean','Sigma']]
      df_rms.loc[:, 'Sigma'] = df_rms['Sigma'].round(2)
      df_rms = df_rms.rename(columns={'Mean':'RMS Noise / mV', 'Sigma':'RMS Unc / mV'})
      dfs_to_concat.append(df_rms)
    elif var == "timeres":
      if first_df_found == False:
        first_df_found = True
        df_tr = df[['Channel','Bias','Resolution @ 30%','Uncertainty @ 30%']]
      else:
        df_tr = df[['Resolution @ 30%','Uncertainty @ 30%']]
      df_tr = df_tr.rename(columns={'Resolution @ 30%':'TR @ 30% / ps', 'Uncertainty @ 30%':'TR Unc @ 30% / ps'})
      dfs_to_concat.append(df_tr)
  
  dfs_comb = pd.concat(dfs_to_concat, axis=1)
  dfs_comb.loc[:, 'Bias'] = dfs_comb['Bias'].str[:-1].astype(int)
  dfs_comb['Thickness / um'] = thickness_col
  dfs_comb['E field / V/cm'] = 10000*(dfs_comb['Bias'] / dfs_comb['Thickness / um'])
  dfs_comb.loc[:, 'E field / V/cm'] = dfs_comb['E field / V/cm'] // 1
  if ('Rise time / ps' in dfs_comb.columns) and ('Amplitude / mV' in dfs_comb.columns) and ('RMS Noise / mV' in dfs_comb.columns):
    dfs_comb['Jitter / ps'] = dfs_comb['RMS Noise / mV'] / (dfs_comb['Amplitude / mV'] / dfs_comb['Rise time / ps'])
    dfs_comb.loc[:, 'Jitter / ps'] = dfs_comb['Jitter / ps'].round(1)
    unc_cpt_rms = dfs_comb['RMS Unc / mV'] / dfs_comb['RMS Noise / mV']
    unc_cpt_risetime = dfs_comb['Rise time Unc / ps'] / dfs_comb['Rise time / ps']
    unc_cpt_ampl = 0 # idk the unc for a Langaus fit
    dfs_comb['Jitter Unc / ps'] = dfs_comb['Jitter / ps'] * np.sqrt(unc_cpt_rms**2 + unc_cpt_risetime**2 + unc_cpt_ampl**2)
    dfs_comb.loc[:, 'Jitter Unc / ps'] = dfs_comb['Jitter Unc / ps'].round(1)
    if 'TR @ 30% / ps' in dfs_comb.columns:
      dfs_comb['Landau TR Cpt / ps'] = np.sqrt(dfs_comb['TR @ 30% / ps']**2 - dfs_comb['Jitter / ps']**2)
      unc_cpt_jit = dfs_comb['Jitter / ps']*dfs_comb['Jitter Unc / ps']
      unc_cpt_tr = dfs_comb['TR @ 30% / ps']*dfs_comb['TR Unc @ 30% / ps']
      dfs_comb['Landau TR Unc / ps'] = np.sqrt(unc_cpt_jit**2 + unc_cpt_tr**2) / dfs_comb['Landau TR Cpt / ps']
      dfs_comb.loc[:, 'Landau TR Cpt / ps'] = dfs_comb['Landau TR Cpt / ps'].round(1)
      dfs_comb.loc[:, 'Landau TR Unc / ps'] = dfs_comb['Landau TR Unc / ps'].round(1)
      dfs_comb['WF6 Param / ps/um'] = dfs_comb['Landau TR Cpt / ps'] / dfs_comb['Thickness / um']
      dfs_comb['WF6 Param Unc / ps/um'] = dfs_comb['WF6 Param / ps/um'] * dfs_comb['Landau TR Unc / ps'] / dfs_comb['Landau TR Cpt / ps']
      dfs_comb.loc[:, 'WF6 Param / ps/um'] = dfs_comb['WF6 Param / ps/um'].round(2)
      dfs_comb.loc[:, 'WF6 Param Unc / ps/um'] = dfs_comb['WF6 Param Unc / ps/um'].round(2)

  columns = dfs_comb.columns.tolist()
  for col_to_move in ['Thickness / um','E field / V/cm']:
    columns.remove(col_to_move)
  columns[2:2] = ['Thickness / um','E field / V/cm']
  dfs_comb = dfs_comb[columns]

  if mcp_channel == True:
    print(pmax_low_mcp)
    dfs_comb['MCP PMAX low / mV'] = np.tile(pmax_low_mcp, number_of_duts)
    dfs_comb['MCP PMAX high / mV'] = np.tile(pmax_high_mcp, number_of_duts)

  dfs_comb['Bias'] = pd.to_numeric(dfs_comb['Bias'], errors='coerce')
  dfs_comb = dfs_comb.sort_values(by=['Channel','Bias'])

  print(f"[BETA ANALYSIS] : [DATA COLLATOR] Writing data to {output_savename}.csv.")
  dfs_comb.to_csv(output_savename+'.csv', index=False)
