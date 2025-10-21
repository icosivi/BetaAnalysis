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

def binned_fit_langauss(samples, bins, max_x_val, nan='remove'):
  if nan == 'remove':
    samples = samples[~np.isnan(samples)]

  hist, bin_edges = np.histogram(samples, bins, range=(0,max_x_val), density=True)
  bin_centers = bin_edges[:-1] + np.diff(bin_edges) / 2

  hist = np.insert(hist, 0, sum(samples < bin_edges[0]))
  bin_centers = np.insert(bin_centers, 0, bin_centers[0] - np.diff(bin_edges)[0])
  hist = np.append(hist, sum(samples > bin_edges[-1]))
  bin_centers = np.append(bin_centers, bin_centers[-1] + np.diff(bin_edges)[0])

  hist = hist[1:]
  bin_centers = bin_centers[1:]

  landau_x_mpv_guess = bin_centers[np.argmax(hist)]
  landau_xi_guess = median_abs_deviation(samples) / 5
  gauss_sigma_guess = landau_xi_guess / 10

  popt, pcov = curve_fit(
    lambda x, mpv, xi, sigma: langauss.pdf(x, mpv, xi, sigma),
    xdata=bin_centers,
    ydata=hist,
    p0=[landau_x_mpv_guess, landau_xi_guess, gauss_sigma_guess],
  )
  print("Langaus fit parameters")
  print(popt)
  print("Langaus covariance matrix")
  print(pcov)
  return popt, pcov, hist, bin_centers

def round_to_sig_figs(x, sig):
  if x == 0:
    return 0
  else:
    return round(x, sig - int(math.floor(math.log10(abs(x)))) - 1)

def main():
  parser = argparse.ArgumentParser(description='Read .root files into an array.')
  parser.add_argument('files', metavar='F', type=str, nargs='+',
                      help='List of .root files or wildcard pattern (*.root)')
  parser.add_argument("--ch", type=int, required=True, help='Channel number to which to fit Langauss')
  args = parser.parse_args()

  files = []
  trees = []

  ch_ind = args.ch - 1

  for pattern in args.files:
    root_files = glob.glob(pattern)
    for root_file in root_files:
      theFile = root.TFile.Open(root_file)
      theTree = theFile.Get("Analysis")
      files.append(theFile)
      trees.append(theTree)

  pmax_list = []
  area_list = []

  for entry in theTree:
    pmax_sig = entry.pmax[ch_ind]
    negpmax_sig = entry.negpmax[ch_ind]
    if (pmax_sig < 20) or (pmax_sig > 200) or (negpmax_sig < -25):
      continue
    else:
      area_sig = entry.area_new[ch_ind]
      pmax_list.append(pmax_sig)
      area_list.append(area_sig)

  pmax = np.array(pmax_list)
  area = np.array(area_list)
  area = area/4.7

  max_val = 35
  bins_tot = max_val*4
  area = area[(area>=0) & (area<=max_val)]

  plt.figure(figsize=(10, 6))
  histo, bins, _ = plt.hist(area, bins=bins_tot, range=(0,max_val), color='blue', edgecolor='black', alpha=0.6, density=True)
  #plt.xlabel('Area [pWb]')
  plt.xlabel('Charge [fC]')
  plt.ylabel('Frequency')
  plt.title('Histogram of signal area values')
  plt.yscale('log')
  plt.show()

  bin_centers = bins[:-1] + np.diff(bins) / 2
  popt, pcov, fitted_hist, bin_centers = binned_fit_langauss(area, bins_tot, max_val)

  count_1p0mpv = sum(1 for value in area if value > popt[0])
  count_1p5mpv = sum(1 for value in area if value > 1.5*popt[0])
  print("fEv@1.5:")
  print(count_1p5mpv/count_1p0mpv)

  fig = go.Figure()
  fig.update_layout(
    #xaxis_title='Area [pWb]',
    xaxis_title='Charge [fC]',
    yaxis_title='Probability Density',
    title='Langauss Fit to Histogram of Signal Area Values',
  )

  fig.add_trace(
    go.Histogram(
      x=area,
      name='Histogram of area',
      histnorm='probability density',
      nbinsx=bins_tot,
      opacity=0.6,
      marker=dict(color='blue', line=dict(color='black', width=1)),
    )
  )

  x_axis = np.linspace(0, max_val, 999)
  fig.add_trace(
    go.Scatter(
      x=x_axis,
      y=langauss.pdf(x_axis, *popt),
      name=f'Langauss Fit<br>MPV={popt[0]:.2e}<br>ξ={popt[1]:.2e}<br>σ={popt[2]:.2e}',
      mode='lines',
    )
  )

  fig.show()


if __name__ == "__main__":
    main()
