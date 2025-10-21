import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import minimize
from scipy.stats import poisson
import ROOT as root
from ROOT import TF1
from scipy.special import gammaln
import math
import pandas as pd
import argparse
import glob
import re
import os
import csv

var_dict = {"tmax":"t_{max} / ns" , "pmax":"p_max / mV" , "negpmax":"-p_max / mV", "charge":"Q / fC", "area_new":"Area / pWb" , "rms":"RMS / mV"}

def main():
  parser = argparse.ArgumentParser(description='Read CSV fit data against bias and plot charge dist, RMS noise, and time resolution.')
  parser.add_argument('incsv', type=str, help='Input CSV file.')
  #parser.add_argument('--doWaveform', action='store_true', help='Enable the doWaveform option')
  #parser.add_argument('--csvOut', action='store_true', help='Output important data for plots against bias')
  args = parser.parse_args()
  
  incsv = args.incsv
  if os.path.exists(incsv):
    df = pd.read_csv(incsv)
    print(f"Successfully read {incsv}")
  else:
    print(f"File {incsv} doesn't exist")

  #if args.doWaveform: waveform(file_array,tree_array,args.ch-1,total_number_channels)
  #if args.csvOut: makeCSV(file_array,tree_array,total_number_channels)

  ufs = 12
  max_charge = 0.0
  max_rms = 0.0
  max_time_res = 0.0
  fig, axs = plt.subplots(3, 1, figsize=(8,12))
  channels = df['index'].unique()
  for channel in channels:
    dfsub = df.loc[df['index'] == channel]
    bias_str = dfsub['Bias'].to_numpy()
    bias = [int(this_bias_val[:-1]) for this_bias_val in bias_str]
    charge = dfsub['Charge'].to_numpy()
    err_q = dfsub['err_Q'].to_numpy()
    rms = dfsub['RMS'].to_numpy()
    err_rms = dfsub['err_rms'].to_numpy()
    time_res = dfsub['Sigma_value'].to_numpy()

    max_charge = np.max([max_charge, np.max(charge)])
    max_rms = np.max([max_rms, np.max(rms)])
    max_time_res = np.max([max_time_res, 1000*np.max(time_res)])

    axs[0].errorbar(bias, charge, err_q, label=f"Channel {channel}", fmt='-o', markersize=10, markeredgewidth=2, markeredgecolor='k')
    axs[0].set_ylabel('Charge / fC', fontsize=ufs)
    axs[0].set_xlabel('Bias / V', fontsize=ufs)
    axs[0].set_xlim(200,250)
    axs[0].set_ylim(0,max_charge*1.2)
    axs[0].grid(True)
    axs[0].axhline(y=15,color='k',linestyle='--',linewidth='2')
    axs[0].legend(loc="upper right", fontsize=ufs)

    axs[1].errorbar(bias, rms, err_rms, label=f"Channel {channel}", fmt='-o', markersize=10, markeredgewidth=2, markeredgecolor='k')
    axs[1].set_ylabel('RMS noise / mV', fontsize=ufs)
    axs[1].set_xlabel('Bias / V', fontsize=ufs)
    axs[1].set_xlim(200,250)
    axs[1].set_ylim(0,max_rms*1.2)
    axs[1].grid(True)
    axs[1].legend(loc="upper right", fontsize=ufs)

    axs[2].plot(bias, 1000*time_res, label=f"Channel {channel}", linestyle='-', marker='o', markersize=10, markeredgewidth=2, markeredgecolor='k')
    axs[2].set_ylabel('Time resolution / ps', fontsize=ufs)
    axs[2].set_xlabel('Bias / V', fontsize=ufs)
    axs[2].set_xlim(200,250)
    axs[2].set_ylim(0,max_time_res*1.2)
    axs[2].grid(True) 
    axs[2].axhline(y=30,color='k',linestyle='--',linewidth='2')
    axs[2].legend(loc="upper right", fontsize=ufs)

  plt.tight_layout()
  plt.savefig("bias_scan_plots.png",facecolor='w')
  print(f"Saved plot of bias scans to bias_scan_plots.png")

if __name__ == "__main__":
    main()


