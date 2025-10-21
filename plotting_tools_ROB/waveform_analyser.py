import numpy as np
import matplotlib.pyplot as plt
import scipy
from scipy.optimize import minimize
from scipy.stats import poisson, median_abs_deviation
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
import plotly.graph_objects as go

import sys
from datetime import datetime
import matplotlib.pylab as plt
import matplotlib.axes as axes
#from langaus import LanGausFit
from array import array
from landaupy import langauss
from scipy.optimize import curve_fit

var_dict = {"tmax":"t_{max} / 10 ns" , "pmax":"p_max / mV" , "negpmax":"-p_max / mV", "charge":"Q / fC", "area_new":"Area / pWb" , "rms":"RMS / mV"}

def round_to_sig_figs(x, sig):
  if x == 0:
    return 0
  else:
    return round(x, sig - int(math.floor(math.log10(abs(x)))) - 1)

def main():
  parser = argparse.ArgumentParser(description='Read .root files into an array.')
  parser.add_argument('files', metavar='F', type=str, nargs='+',
                      help='List of .root files or wildcard pattern (*.root)')
  args = parser.parse_args()

  files = []
  trees = []

  for pattern in args.files:
    root_files = glob.glob(pattern)
    for root_file in root_files:
      theFile = root.TFile.Open(root_file)
      theTree = theFile.Get("Analysis")
      files.append(theFile)
      trees.append(theTree)

  t_data = [[] for _ in range(len(trees))]
  w_data = [[] for _ in range(len(trees))]

  for j in range(len(trees)):
    tree = trees[j]
    i = 0
    ev_true_count = 0
    for entry in tree:
      i += 1
      if i > 10000:
        continue
      pmax_sig = entry.pmax[0]
      negpmax_sig = entry.negpmax[0]
      pmax_mcp = entry.pmax[1]
      peakfind = entry.cfd[0][1]
      if (pmax_sig > 20) and (pmax_sig < 80) and (negpmax_sig > -30) and (pmax_mcp < 540):
        # W12 15e14/25e14 (pmax_sig > 10) and (pmax_sig < 30) and (negpmax_sig > -30) and (pmax_mcp < 120) and (peakfind > 9) and (peakfind < 14)
        # W13 35e14 (pmax_sig > 55) and (pmax_sig < 80) and (negpmax_sig > -30) and (pmax_mcp < 120) and (peakfind > 9) and (peakfind < 14)
        w_sig = entry.w[0]
        t_sig = entry.t[0]
        w_data[j].extend(w_sig)
        t_data[j].extend(t_sig)
        ev_true_count += 1
        if ev_true_count >= 10: 
          print("Ten values/curves")
          break

  t_data = [np.array(bias) for bias in t_data]
  w_data = [np.array(bias) for bias in w_data]
  print(len(t_data))

  plt.figure(figsize=(10, 6))

  colours = ['black','dodgerblue']
  #colours = ['black','blue']
  alphas = [1.0,0.8]
  edges = ['black','none']
  #labels = ["W5 non irradiated (120 V) G = 8","W13 35e14 (420 V) G = 35"]
  labels = ["W5 non irradiated (120 V) G = 8","W12 15e14 (210 V) G = 4"]

  for i in range(2):
    if i == 0: continue
    time_data = t_data[i]*(10**9)
    print(len(time_data))
    reshaped_time_data = time_data.reshape(10,502)
    reshaped_ampl_data = w_data[i].reshape(10,502)
    #scatter = plt.scatter(time_data,w_data[i],s=10,c=colours[i],marker='o',edgecolor=edges[i],linewidth=0.5,alpha=alphas[i],label=labels[i])
    plt.plot(reshaped_time_data[0],reshaped_ampl_data[0],linestyle='-',color=colours[i],linewidth=1.5,alpha=alphas[i],label=labels[i])
    for j in range(len(reshaped_time_data)-1):
      plt.plot(reshaped_time_data[j+1],reshaped_ampl_data[j+1],linestyle='-',color=colours[i],linewidth=1.5,alpha=alphas[i])

  plt.xlabel('Time [ns]',fontsize=14)
  plt.ylabel('Amplitude [mV]',fontsize=14)
  plt.xticks(fontsize=14)
  plt.yticks(fontsize=14)
  plt.xlim(-3,4)
  #plt.legend(fontsize=14)
  #plt.yscale('log')
  plt.grid(True, linestyle='--', alpha=0.5)
  plt.savefig("./amplitude_analysis.png",dpi=300,facecolor='w')
  plt.show()


if __name__ == "__main__":
    main()
