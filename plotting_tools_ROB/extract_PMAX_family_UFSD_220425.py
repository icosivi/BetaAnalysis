import numpy as np
import matplotlib.pyplot as plt
import scipy
import scipy.optimize as opt
import mpmath
import scipy.special as sp
from scipy.optimize import minimize, curve_fit
from scipy.stats import poisson, median_abs_deviation, norm
import scipy.interpolate as interp
from scipy.interpolate import CubicSpline
from lmfit.models import VoigtModel
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
from sklearn.metrics import mean_absolute_error, mean_squared_error
import matplotlib.gridspec as gridspec

import sys
from datetime import datetime
import matplotlib.pylab as plt
import matplotlib.axes as axes
#from langaus import LanGausFit
from array import array
from landaupy import langauss

var_dict = {"tmax":"t_{max} / 10 ns" , "pmax":"p_max / mV" , "negpmax":"-p_max / mV", "charge":"Q / fC", "area_new":"Area / pWb" , "rms":"RMS / mV"}

def round_to_sig_figs(x, sig):
  if x == 0:
    return 0
  else:
    return round(x, sig - int(math.floor(math.log10(abs(x)))) - 1)

def gaussian(x, A, mu, sigma):
  return A * np.exp(-((x - mu) ** 2) / (2 * sigma ** 2))

def gaussian_fit_binned_data(data_to_bin, fittype):
  if len(data_to_bin) == 0:
    print(f"No data to fit")
    return 0, 0, 0
  counts, bin_edges = np.histogram(data_to_bin, bins=60, range=(-0.4, 0.2))
  bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2
  p0 = [max(counts), np.mean(data_to_bin), np.std(data_to_bin)]
  try:
    params, _ = curve_fit(gaussian, bin_centers, counts, p0=p0)
  except RuntimeError:
    print(f"{fittype} fit failed")
    return 0, 0, 0
  return params

def lorentzian(x, A, x0, gamma):
  return A * gamma**2 / ((x - x0)**2 + gamma**2)

def fit_lorentzian(x, y):
  p0 = [max(y), x[np.argmax(y)], (max(x) - min(x)) / 10]  # Initial guesses
  popt, _ = opt.curve_fit(lorentzian, x, y, p0=p0)
  A, x0, gamma = popt
  y_fit = lorentzian(x, *popt)
  return A, y_fit, lambda x_new: lorentzian(x_new, *popt), compute_errors(y, y_fit)

def fit_voigt(x, y):
  model = VoigtModel()
  params = model.guess(y, x=x)
  result = model.fit(y, params, x=x)
  y_fit = result.best_fit

  def voigt_func(x_new):
    return model.eval(result.params, x=x_new)
    
  return max(y_fit), y_fit, voigt_func, compute_errors(y, y_fit)

def fit_cubic_spline(x, y):
  spline = interp.CubicSpline(x, y)
  x_spline_fine = np.linspace(min(x), max(x), 1000)
  y_fine = spline(x_spline_fine)
  peak = max(y_fine)
  return peak, spline(x), spline, compute_errors(y, spline(x))

def landau_pdf(x, A, x0, sigma):
  v = (x - x0) / sigma
  return A * np.exp(-0.5 * (v + np.exp(-v)))

def fit_landau(x, y):
  p0 = [max(y), x[np.argmax(y)], (max(x) - min(x)) / 10]
  try:
    popt, _ = opt.curve_fit(landau_pdf, x, y, p0=p0)
    A, x0, sigma = popt
    y_fit = landau_pdf(x, *popt)
    return max(y_fit), y_fit, lambda x_new: landau_pdf(x_new, *popt), compute_errors(y, y_fit)
  except RuntimeError:
    print("Fucking fit fail")
    return None, None, None, None

def skewed_gaussian(x, A, mu, sigma, alpha):
  return 2 * A * np.exp(-((x - mu) ** 2) / (2 * sigma ** 2)) * norm.cdf(alpha * (x - mu))

def fit_skewed_gaussian(x, y):
  A0 = max(y)
  mu0 = x[np.argmax(y)]
  sigma0 = (max(x) - min(x)) / 10
  alpha0 = 1  # skew factor

  p0 = [A0, mu0, sigma0, alpha0]

  try:
    popt, _ = opt.curve_fit(skewed_gaussian, x, y, p0=p0)
    A, mu, sigma, alpha = popt
    y_fit = skewed_gaussian(x, *popt)
    return A, y_fit, lambda x_new: skewed_gaussian(x_new, *popt), compute_errors(y, y_fit)
  except RuntimeError:
    return 0, 0, 0, [0, 0]

def compute_errors(y_true, y_fit):
  mae = np.mean(np.abs(y_true - y_fit))
  rmse = np.sqrt(np.mean((y_true - y_fit) ** 2))
  return mae, rmse

def find_CFD_time_with_threshold(x, y, amplitude_value, fraction):
  idx = np.where(y >= amplitude_value*fraction)[0]
  if idx.size > 0:
    x1, y1 = x[idx[0] - 1], y[idx[0] - 1]
    x2, y2 = x[idx[0]], y[idx[0]]
    x_CFD = x1 + ((fraction*amplitude_value - y1) / (y2 - y1)) * (x2 - x1)
    return x_CFD
  else:
    print("Fit failed")
    return fraction*max(y)

def find_CFD_time_with_threshold_spline(x, y, amplitude_value, fraction):
  threshold = fraction * amplitude_value
  spline = CubicSpline(x, y)
  x_fine = np.linspace(x[0], x[-1], 1000)
  y_fine = spline(x_fine)
  idx = np.where((y_fine[:-1] < threshold) & (y_fine[1:] >= threshold))[0]
  if len(idx) == 0:
    print("Fit failed")
    return threshold
  i = idx[0]
  x0, x1 = x_fine[i], x_fine[i+1]
  y0, y1 = y_fine[i], y_fine[i+1]
  x_CFD = x0 + (threshold - y0) * (x1 - x0) / (y1 - y0)
  return x_CFD

def linear_interpolation(x1, y1, x2, y2, y_prime):
  m = (y2 - y1) / (x2 - x1)
  c = y1 - m * x1
  x_prime = (y_prime - c) / m
  return x_prime

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
  t_data_mcp = [[] for _ in range(len(trees))]
  w_data_mcp = [[] for _ in range(len(trees))]

  num_curves = 22798
  ch_sig = 0
  ch_mcp = 0

  for j in range(len(trees)):
    tree = trees[j]
    i = 0
    ev_true_count = 0
    for entry in tree:
      i += 1
      #if i > 10000:
      #  continue
      #entry_index = entry.event
      pmax_sig = entry.Pmax1
      pmax_mcp = entry.Pmax2
      if (pmax_sig > 80) and (pmax_sig < 280) and (pmax_mcp < 300) and (pmax_mcp > 40):
        # W12 15e14/25e14 (pmax_sig > 10) and (pmax_sig < 30) and (negpmax_sig > -30) and (pmax_mcp < 120) and (peakfind > 9) and (peakfind < 14)
        # W13 35e14 (pmax_sig > 55) and (pmax_sig < 80) and (negpmax_sig > -30) and (pmax_mcp < 120) and (peakfind > 9) and (peakfind < 14)
        w_sig = entry.w1
        t_sig = entry.t1
        w_data[j].extend(w_sig)
        t_data[j].extend(t_sig)
        w_mcp = entry.w2
        t_mcp = entry.t2
        w_data_mcp[j].extend(w_mcp)
        t_data_mcp[j].extend(t_mcp)
        ev_true_count += 1
        if ev_true_count >= num_curves: 
          print(f"{num_curves} values/curves")
          break

  t_data = [np.array(bias) for bias in t_data]
  w_data = [np.array(bias) for bias in w_data]
  t_data_mcp = [np.array(bias) for bias in t_data_mcp]
  w_data_mcp = [np.array(bias) for bias in w_data_mcp]

  colours = ['black','dodgerblue']
  #colours = ['black','blue']
  alphas = [1.0,0.8]
  edges = ['black','none']
  #labels = ["W5 new (140 V)","W13 35e14 (420 V) G = 35"]
  #labels = ["W16 new (210 V)","W12 15e14 (210 V) G = 4"]
  labels = ["W17 new (160 V)","W12 15e14 (210 V) G = 4"]

  a_max = []
  a_max_mcp = []
  numptseitherside = 3
  add_noise = False

  amp_fam_1_ind = []
  amp_fam_2_ind = []
  amp_fam_3_ind = []

  for i in range(1):
    #if i == 0: continue
    time_data = t_data[i]*(10**9)
    time_data_mcp = t_data_mcp[i]*(10**9)

    reshaped_time_data = time_data.reshape(num_curves,1002)
    reshaped_ampl_data = w_data[i].reshape(num_curves,1002)
    rtd_mcp = time_data_mcp.reshape(num_curves,1002)
    rad_mcp = w_data_mcp[i].reshape(num_curves,1002)

    if add_noise:
      mean = 0.0
      std_dev = 1.0
      noise = np.random.normal(mean, std_dev, size=1002)
      #noise = np.clip(noise, -1, 1)
      for j in range(num_curves):
        reshaped_ampl_data[j] = reshaped_ampl_data[j] + 0.001*noise
        rad_mcp[j] = rad_mcp[j] + 0.001*noise

    for j in range(num_curves):

      pmax = reshaped_ampl_data[j].max()
      peak_idx = np.argmax(reshaped_ampl_data[j])
      start_idx = max(0, peak_idx - numptseitherside)
      end_idx = min(len(reshaped_ampl_data[j]), peak_idx + numptseitherside + 1)

      x_peak = reshaped_time_data[j][start_idx:end_idx]
      y_peak = reshaped_ampl_data[j][start_idx:end_idx]

      fam_disc = reshaped_ampl_data[j][start_idx+1] - reshaped_ampl_data[j][end_idx-2]
      if fam_disc < -0.002:
        amp_fam_1_ind.append(j)
      elif (fam_disc > -0.002) & (fam_disc < 0.002):
        amp_fam_2_ind.append(j)
      else:
        amp_fam_3_ind.append(j)
      a_max.append(1000*pmax)

  pop_fam1 = len(amp_fam_1_ind)
  pop_fam2 = len(amp_fam_2_ind)
  pop_fam3 = len(amp_fam_3_ind)

  a_max_fam1 = [a_max[i] for i in amp_fam_1_ind]
  a_max_fam2 = [a_max[i] for i in amp_fam_2_ind]
  a_max_fam3 = [a_max[i] for i in amp_fam_3_ind]

  fig, axes = plt.subplots(1, 3, figsize=(20, 6))
  bins = 200
  axes[0].hist(a_max_fam1, bins=bins, range=(80, 280), color='skyblue', edgecolor='black')
  axes[1].hist(a_max_fam2, bins=bins, range=(80, 280), color='skyblue', edgecolor='black')
  axes[2].hist(a_max_fam3, bins=bins, range=(80, 280), color='skyblue', edgecolor='black')
  axes[0].set_title(f'Family 1 : {pop_fam1} Ev')
  axes[1].set_title(f'Family 2 : {pop_fam2} Ev')
  axes[2].set_title(f'Family 3 : {pop_fam3} Ev')
  x_limits = (75, 285)
  y_limits = (0, 150)
  for ax in axes:
    ax.set_xlim(x_limits)
    ax.set_ylim(y_limits)
    ax.grid(True)

  plt.tight_layout()
  plt.savefig("./PMAX_dist_wideselection.png",dpi=300,facecolor='w')
  plt.clf()


if __name__ == "__main__":
  main()
