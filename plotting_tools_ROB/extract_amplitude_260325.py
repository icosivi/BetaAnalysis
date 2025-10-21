import numpy as np
import matplotlib.pyplot as plt
import scipy
import scipy.optimize as opt
import mpmath
import scipy.special as sp
from scipy.optimize import minimize
from scipy.stats import poisson, median_abs_deviation, norm
import scipy.interpolate as interp
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
from scipy.optimize import curve_fit

var_dict = {"tmax":"t_{max} / 10 ns" , "pmax":"p_max / mV" , "negpmax":"-p_max / mV", "charge":"Q / fC", "area_new":"Area / pWb" , "rms":"RMS / mV"}

def round_to_sig_figs(x, sig):
  if x == 0:
    return 0
  else:
    return round(x, sig - int(math.floor(math.log10(abs(x)))) - 1)

def gaussian(x, A, mu, sigma):
  return A * np.exp(-((x - mu) ** 2) / (2 * sigma ** 2))

def gaussian_fit_binned_data(data_to_bin, fittype):
  counts, bin_edges = np.histogram(data_to_bin, bins=150, range=(-2, 1))
  bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2
  p0 = [max(counts), np.mean(data_to_bin), np.std(data_to_bin)]
  figoneoff, axoneoff = plt.subplots(figsize=(8, 6))
  axoneoff.hist(data_to_bin, bins=150, range=(-2, 1), alpha=0.6, color='red', edgecolor='black', label=f'No fit {fittype}', histtype='stepfilled', linewidth=1.5)
  axoneoff.legend()
  figoneoff.savefig(f"{fittype}_nongaussian.png", dpi=300, bbox_inches='tight')
  try:
    params, _ = curve_fit(gaussian, bin_centers, counts, p0=p0)
    print(params)
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
  x_fine = np.linspace(min(x), max(x), 1000)
  y_fine = spline(x_fine)
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
    print("Landau fit failed—perhaps your data isn't 'Landau-y' enough? Try better initial guesses.")
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
    return max(y_fit), y_fit, lambda x_new: skewed_gaussian(x_new, *popt), compute_errors(y, y_fit)
  except RuntimeError:
    return 0, 0, 0, [0, 0]

def compute_errors(y_true, y_fit):
  mae = np.mean(np.abs(y_true - y_fit))
  rmse = np.sqrt(np.mean((y_true - y_fit) ** 2))
  return mae, rmse

def extract_CFD_time(t_vals, a_vals, a_peak, frac):
  A_CFD = frac * a_peak / 1000
  for i in range(len(a_vals) - 1):
    if a_vals[i] < A_CFD and a_vals[i + 1] >= A_CFD:
      t_CFD = t_vals[i] + (A_CFD - a_vals[i]) * (t_vals[i + 1] - t_vals[i]) / (a_vals[i + 1] - a_vals[i])
      return t_CFD

  return None 

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

  cfd10_data = []
  cfd20_data = []
  cfd30_data = []
  cfd30_mcp_data = []
  num_curves = 1421
  ch_sig = 2
  ch_mcp = 3

  for j in range(len(trees)):
    tree = trees[j]
    i = 0
    ev_true_count = 0
    for entry in tree:
      i += 1
      #if i > 10000:
      #  continue
      #entry_index = entry.event
      pmax_sig = entry.pmax[ch_sig]
      negpmax_sig = entry.negpmax[ch_sig]
      pmax_mcp = entry.pmax[ch_mcp]
      cfd10_sig = entry.cfd[ch_sig][0] # 10%
      cfd20_sig = entry.cfd[ch_sig][1] # 20%
      cfd30_mcp = entry.cfd[ch_mcp][2]
      cfd30_sig = entry.cfd[ch_sig][2] # 30%
      if (pmax_sig > 25) and (pmax_mcp < 540) and (pmax_mcp > 20):
        # W12 15e14/25e14 (pmax_sig > 10) and (pmax_sig < 30) and (negpmax_sig > -30) and (pmax_mcp < 120) and (peakfind > 9) and (peakfind < 14)
        # W13 35e14 (pmax_sig > 55) and (pmax_sig < 80) and (negpmax_sig > -30) and (pmax_mcp < 120) and (peakfind > 9) and (peakfind < 14)
        w_sig = entry.w[ch_sig]
        t_sig = entry.t[ch_sig]
        cfd10_data.append(cfd10_sig)
        cfd20_data.append(cfd20_sig)
        cfd30_data.append(cfd30_sig)
        cfd30_mcp_data.append(cfd30_mcp)
        w_data[j].extend(w_sig)
        t_data[j].extend(t_sig)
        w_mcp = entry.w[ch_mcp]
        t_mcp = entry.t[ch_mcp]
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
  labels = ["20μm (CBL) 150 V","W12 15e14 (210 V) G = 4"]
  #labels = ["W17 new (160 V)","W12 15e14 (210 V) G = 4"]

  a_max = []
  a_para = []
  a_gaus = []
  a_lorentz = []
  a_voigt = []
  a_spline = []
  a_landau = []
  a_skewG = []

  para_mae = []
  para_rmse = []
  gaus_mae = []
  gaus_rmse = []
  lorentz_mae = []
  lorentz_rmse = []
  voigt_mae = []
  voigt_rmse = []
  spline_mae = []
  spline_rmse = []
  landau_mae = []
  landau_rmse = []
  skewG_mae = []
  skewG_rmse = []

  a_max_mcp = []
  a_para_mcp = []
  a_gaus_mcp = []
  a_lorentz_mcp = []
  a_voigt_mcp = []
  a_spline_mcp = []
  a_landau_mcp = []
  a_skewG_mcp = []

  make_plots = False
  make_populated_plots = True
  cfd_studies = False
  add_noise = False
  time_res_calc = False
  numptseitherside = 3

  if make_plots:
    #plt.figure(figsize=(10, 6))
    fig, ax = plt.subplots(figsize=(16, 10))

  for i in range(1):
    #if i == 0: continue
    time_data = t_data[i]*(10**9)
    time_data_mcp = t_data_mcp[i]*(10**9)

    reshaped_time_data = time_data.reshape(num_curves,502)
    reshaped_ampl_data = w_data[i].reshape(num_curves,502)
    rtd_mcp = time_data_mcp.reshape(num_curves,502)
    rad_mcp = w_data_mcp[i].reshape(num_curves,502)

    if add_noise:
      mean = 0.0
      std_dev = 1.0
      noise = np.random.normal(mean, std_dev, size=502)
      #noise = np.clip(noise, -1, 1)
      for j in range(num_curves):
        reshaped_ampl_data[j] = reshaped_ampl_data[j] + 0.001*noise
        rad_mcp[j] = rad_mcp[j] + 0.001*noise

    for j in range(num_curves):

      pmax = reshaped_ampl_data[j].max()
      peak_idx = np.argmax(reshaped_ampl_data[j])
      tmax_for_plot = reshaped_time_data[j][peak_idx]
      start_idx = max(0, peak_idx - numptseitherside)
      end_idx = min(len(reshaped_ampl_data[j]), peak_idx + numptseitherside + 1)

      x_peak = reshaped_time_data[j][start_idx:end_idx]
      y_peak = reshaped_ampl_data[j][start_idx:end_idx]

      pmax_mcp = rad_mcp[j].max()
      peak_idx_mcp = np.argmax(rad_mcp[j])
      start_idx_mcp = max(0, peak_idx_mcp - numptseitherside)
      end_idx_mcp = min(len(rad_mcp[j]), peak_idx_mcp + numptseitherside + 1)

      x_peak_mcp = rtd_mcp[j][start_idx_mcp:end_idx_mcp]
      y_peak_mcp = rad_mcp[j][start_idx_mcp:end_idx_mcp]

      if len(x_peak) == 0 or len(x_peak_mcp) == 0:
        continue
      x_fine = np.linspace(min(x_peak), max(x_peak), 1000)

      # parabola
      parabola_coeffs = np.polyfit(x_peak, y_peak, 2)
      a, b, _ = parabola_coeffs
      x_parabola_max = -b / (2 * a)
      y_parabola_max = np.polyval(parabola_coeffs, x_parabola_max)

      y_parabolic_fit = np.polyval(parabola_coeffs, x_peak)
      parabolic_residuals = y_peak - y_parabolic_fit
      parabolic_mae = mean_absolute_error(y_peak, y_parabolic_fit)
      parabolic_rmse = np.sqrt(mean_squared_error(y_peak, y_parabolic_fit))
      parabolic_max_error = np.max(np.abs(parabolic_residuals))

      parabola_coeffs_mcp = np.polyfit(x_peak_mcp, y_peak_mcp, 2)
      a_mcp, b_mcp, _ = parabola_coeffs_mcp
      x_parabola_max_mcp = -b_mcp / (2 * a_mcp)
      y_parabola_max_mcp = np.polyval(parabola_coeffs_mcp, x_parabola_max_mcp)

      idx_below = np.where((y_peak < y_parabola_max) & (x_peak < x_parabola_max))[0][-1]
      t_w_below_Amax_para = x_peak[idx_below]

      # gaussian
      p0 = [np.max(y_peak), x_peak[np.argmax(y_peak)], 1]
      p0_mcp = [np.max(y_peak_mcp), x_peak_mcp[np.argmax(y_peak_mcp)], 1]
      try:
        params, _ = opt.curve_fit(gaussian, x_peak, y_peak, p0=p0)
        A, mu, sigma = params

        y_gaussian_fit = gaussian(x_peak, *params)
        gaussian_residuals = y_peak - y_gaussian_fit
        gaussian_mae = mean_absolute_error(y_peak, y_gaussian_fit)
        gaussian_rmse = np.sqrt(mean_squared_error(y_peak, y_gaussian_fit))
        gaussian_max_error = np.max(np.abs(gaussian_residuals))

        params_mcp, _ = opt.curve_fit(gaussian, x_peak_mcp, y_peak_mcp, p0=p0_mcp)
        _, mu_mcp, _ = params_mcp
      except RuntimeError:
        continue

      # Lorentzian
      try:
        lorentz_peak, lorentz_fit, lorentz_func, lorentz_errors = fit_lorentzian(x_peak, y_peak)
        lorentz_peak_mcp, lorentz_fit_mcp, lorentz_func_mcp, lorentz_errors_mcp = fit_lorentzian(x_peak_mcp, y_peak_mcp)
        y_fine = lorentz_func(x_fine)
        t_Amax_lorentz = x_fine[np.argmin(np.abs(y_fine - lorentz_peak))]
        idx_below = np.where((y_peak < lorentz_peak) & (x_peak < x_fine[np.argmin(np.abs(y_fine - lorentz_peak))]))[0][-1]
        t_w_below_Amax_lorentz = x_peak[idx_below]
      except RuntimeError:
        continue

      # Voigt
      try:
        voigt_peak, voigt_fit, voigt_func, voigt_errors = fit_voigt(x_peak, y_peak)
        voigt_peak_mcp, voigt_fit_mcp, voigt_func_mcp, voigt_errors_mcp = fit_voigt(x_peak_mcp, y_peak_mcp)
      except RuntimeError:
        continue

      # interp_spline
      try:
        spline_peak, spline_fit, spline_func, spline_errors = fit_cubic_spline(x_peak, y_peak)
        spline_peak_mcp, spline_fit_mcp, spline_func_mcp, spline_errors_mcp = fit_cubic_spline(x_peak_mcp, y_peak_mcp)
      except RuntimeError:
        continue

      # landau
      landau_peak, landau_fit, landau_func, landau_errors = fit_landau(x_peak, y_peak)
      landau_peak_mcp, landau_fit_mcp, landau_func_mcp, landau_errors_mcp = fit_landau(x_peak_mcp, y_peak_mcp)
      y_fine = landau_func(x_fine)
      landau_peak = max(y_fine)
      t_Amax_landau = x_fine[np.argmin(np.abs(y_fine - landau_peak))]
      idx_below = np.where((y_peak < landau_peak) & (x_peak < x_fine[np.argmin(np.abs(y_fine - landau_peak))]))[0][-1]
      t_w_below_Amax_landau = x_peak[idx_below]

      # skewed_gaus
      try:
        skewG_peak, skewG_fit, skewG_func, skewG_errors = fit_skewed_gaussian(x_peak, y_peak)
        skewG_peak_mcp, skewG_fit_mcp, skewG_func_mcp, skewG_errors_mcp = fit_skewed_gaussian(x_peak_mcp, y_peak_mcp)
        y_fine = skewG_func(x_fine)
        skewG_peak = max(y_fine)
        t_Amax_skewG = x_fine[np.argmin(np.abs(y_fine - skewG_peak))]
        idx_below = np.where((y_peak < skewG_peak) & (x_peak < x_fine[np.argmin(np.abs(y_fine - skewG_peak))]))[0][-1]
        t_w_below_Amax_skewG = x_peak[idx_below]
      #except RuntimeError:
      except:
        continue

      a_max.append(1000*pmax)
      a_para.append(1000*y_parabola_max)
      a_gaus.append(1000*gaussian(mu, *params))
      a_lorentz.append(1000*lorentz_peak)
      a_voigt.append(1000*voigt_peak)
      a_spline.append(1000*spline_peak)
      a_landau.append(1000*landau_peak)
      a_skewG.append(1000*skewG_peak)

      para_mae.append(parabolic_mae)
      para_rmse.append(parabolic_rmse)
      gaus_mae.append(gaussian_mae)
      gaus_rmse.append(gaussian_rmse)
      lorentz_mae.append(lorentz_errors[0])
      lorentz_rmse.append(lorentz_errors[1])
      voigt_mae.append(voigt_errors[0])
      voigt_rmse.append(voigt_errors[1])
      spline_mae.append(spline_errors[0])
      spline_rmse.append(spline_errors[1])
      landau_mae.append(landau_errors[0])
      landau_rmse.append(landau_errors[1])
      skewG_mae.append(skewG_errors[0])
      skewG_rmse.append(skewG_errors[1])

      a_max_mcp.append(1000*pmax_mcp)
      a_para_mcp.append(1000*y_parabola_max_mcp)
      a_gaus_mcp.append(1000*gaussian(mu_mcp, *params_mcp))
      a_lorentz_mcp.append(1000*lorentz_peak_mcp)
      a_voigt_mcp.append(1000*voigt_peak_mcp)
      a_spline_mcp.append(1000*spline_peak_mcp)
      a_landau_mcp.append(1000*landau_peak_mcp)
      a_skewG_mcp.append(1000*skewG_peak_mcp)

      if make_plots:
        x_linspace = np.linspace(x_peak.min(), x_peak.max(), 400)

        ax.plot(x_linspace, 1000*np.polyval(parabola_coeffs, x_linspace), 'r', label=r"$A_{\rm{para}}$: "+str(round(1000*y_parabola_max, 2))+" mV", linewidth = 3, zorder=2)
        ax.plot(x_linspace, 1000*gaussian(x_linspace, *params), 'g', label=r"$A_{\rm{Gaus}}$: "+str(round(1000*gaussian(mu, *params), 2))+" mV", linewidth = 3, zorder=2)
        ax.plot(x_linspace, 1000*lorentz_func(x_linspace), 'blue', label=r"$A_{\rm{Lorentz}}$: "+str(round(1000*lorentz_peak, 2))+" mV", linewidth = 3, zorder=2)
        ax.plot(x_linspace, 1000*voigt_func(x_linspace), 'orange', label=r"$A_{\rm{Voigt}}$: "+str(round(1000*voigt_peak, 2))+" mV", linewidth = 3, zorder=2)
        ax.plot(x_linspace, 1000*spline_func(x_linspace), 'purple', label=r"$A_{\rm{spline}}$: "+str(round(1000*spline_peak, 2))+" mV", linewidth = 3, zorder=2)
        #ax.plot(x_linspace, 1000*landau_func(x_linspace), 'brown', label=r"$A_{\rm{Landau}}$: "+str(round(1000*landau_peak, 2))+" mV", linewidth = 3, zorder=2)
        ax.plot(x_linspace, 1000*skewG_func(x_linspace), 'cyan', label=r"$A_{\rm{skewed}}$: "+str(round(1000*skewG_peak, 2))+" mV", linewidth = 3, zorder=2)
        ax.axhline(1000*y_parabola_max, color='r', linestyle=':', linewidth = 3, zorder=3)
        ax.axhline(1000*gaussian(mu, *params), color='g', linestyle=':', linewidth = 3, zorder=3)
        ax.axhline(1000*lorentz_peak, color='blue', linestyle=':', linewidth = 3, zorder=3)
        ax.axhline(1000*voigt_peak, color='orange', linestyle=':', linewidth = 3, zorder=3)
        ax.axhline(1000*spline_peak, color='purple', linestyle=':', linewidth = 3, zorder=3)
        #ax.axhline(1000*landau_peak, color='brown', linestyle=':', linewidth = 3, zorder=3)
        ax.axhline(1000*skewG_peak, color='cyan', linestyle=':', linewidth = 3, zorder=3)
        ax.scatter(reshaped_time_data[j],1000*reshaped_ampl_data[j],s=140,c="k",marker='o',edgecolor=edges[i],linewidth=0.5,alpha=1.0, zorder=4)
        ax.scatter(reshaped_time_data[j][peak_idx-3:peak_idx+4],1000*reshaped_ampl_data[j][peak_idx-3:peak_idx+4],s=200,c="white",marker='o',edgecolor=edges[i],linewidth=3,alpha=1.0,
                   label=labels[i]+",\n" + r"$A_{\rm{Sa}}$: "+str(round(1000*pmax, 2))+" mV", zorder=5)
        ax.axvline(tmax_for_plot, color='k', alpha = 0.8, linestyle='--', linewidth = 3)
        #ax.axvline(x_parabola_max, color='k', linestyle=':', label=r"t(A$_{para}$) = "+str(round(1000*x_parabola_max, 2))+" ps", linewidth = 3, zorder=3)
        #ax.axvline(t_w_below_Amax_para, color='k', linestyle='-.', label=r"t(w < A$_{para}$): "+str(round(1000*t_w_below_Amax_para, 2))+" ps", linewidth = 3, zorder=3)
        #ax.axvline(t_w_below_Amax_para+0.1, color='k', linestyle='-.', label=r"t(w < A$_{para}$): ("+str(round(1000*t_w_below_Amax_para, 2))+" + 100) ps", linewidth = 3, zorder=3, alpha=0.6)
        #ax.fill_betweenx(y=np.linspace(-50, 500, 100), x1=t_w_below_Amax_para, x2=x_parabola_max, color='orange', alpha=0.4, hatch='//')
        #ax.axvline(t_Amax_lorentz, color='k', linestyle=':', label=r"t(A$_{Lorentz}$) = "+str(round(1000*t_Amax_lorentz, 2))+" ps", linewidth = 3, zorder=3)
        #ax.axvline(t_w_below_Amax_lorentz, color='k', linestyle='-.', label=r"t(w < A$_{Lorentz}$): "+str(round(1000*t_w_below_Amax_lorentz, 2))+" ps", linewidth = 3, zorder=3)
        #ax.fill_betweenx(y=np.linspace(-50, 500, 100), x1=t_w_below_Amax_lorentz, x2=t_Amax_lorentz, color='orange', alpha=0.4, hatch='//')
        #ax.axvline(t_Amax_landau, color='k', linestyle=':', label=r"t(A$_{Landau}$) = "+str(round(1000*t_Amax_landau, 2))+" ps", linewidth = 3, zorder=3)
        #ax.axvline(t_w_below_Amax_landau, color='k', linestyle='-.', label=r"t(w < A$_{Landau}$): "+str(round(1000*t_w_below_Amax_landau, 2))+" ps", linewidth = 3, zorder=3)
        #ax.fill_betweenx(y=np.linspace(-50, 500, 100), x1=t_w_below_Amax_landau, x2=t_Amax_landau, color='orange', alpha=0.4, hatch='//')
        if num_curves == 1:
          print(f"A_max = {(1000*pmax):.2f} mV")
          print(f"A_para = {(1000*y_parabola_max):.2f} mV")
          print(f"A_Gaus = {(1000*gaussian(mu, *params)):.2f} mV")
          print(f"A_Lorentz = {(1000*lorentz_peak):.2f} mV")
          print(f"A_Voigt = {(1000*voigt_peak):.2f} mV")
          print(f"A_spline = {(1000*spline_peak):.2f} mV")
          print(f"A_Landau = {(1000*landau_peak):.2f} mV")
          print(f"A_skewed = {(1000*skewG_peak):.2f} mV")

          print("\nParabolic Fit Errors:")
          print(f"  MAE  = {parabolic_mae:.4f}")
          print(f"  RMSE = {parabolic_rmse:.4f}")

          print("\nGaussian Fit Errors:")
          print(f"  MAE  = {gaussian_mae:.4f}")
          print(f"  RMSE = {gaussian_rmse:.4f}")

          print("\nLorentz Fit Errors:")
          print(f"  MAE  = {lorentz_errors[0]:.4f}")
          print(f"  RMSE = {lorentz_errors[1]:.4f}")

          print("\nVoigt Fit Errors:")
          print(f"  MAE  = {voigt_errors[0]:.4f}")
          print(f"  RMSE = {voigt_errors[1]:.4f}")

          print("\nInterpolated Spline Errors:")
          print(f"  MAE  = {spline_errors[0]:.4f}")
          print(f"  RMSE = {spline_errors[1]:.4f}")

          print("\nLandau Fit Errors:")
          print(f"  MAE  = {landau_errors[0]:.4f}")
          print(f"  RMSE = {landau_errors[1]:.4f}")

          print("\nSkewed Gaus Fit Errors:")
          print(f"  MAE  = {skewG_errors[0]:.4f}")
          print(f"  RMSE = {skewG_errors[1]:.4f}")

    if make_plots:
      handles, labels = ax.get_legend_handles_labels()
      handles = handles[-1:] + handles[:-1]
      labels = labels[-1:] + labels[:-1]

      ax.set_xlabel('Time [ns]',fontsize=24)
      ax.set_ylabel('Amplitude [mV]',fontsize=24)
      plt.xticks(fontsize=24)
      plt.yticks(fontsize=24)
      ax.tick_params(axis='both', labelsize=24)
      ax.set_xlim(-0.5,0.1)
      ax.set_ylim(20,70)
      ax.legend(handles, labels, fontsize=24)
      #ax.set_yscale('log')
      ax.grid(True, linestyle='--', alpha=0.5)
      plt.tight_layout()
      plt.savefig("./amplitude_analysis.png",dpi=300,facecolor='w')
      #plt.show()
      plt.clf()

  if make_populated_plots:
    fig = plt.figure(figsize=(16, 10))
    gs = gridspec.GridSpec(2, 3, height_ratios=[1, 1], width_ratios=[1, 1, 1])

    ax1 = fig.add_subplot(gs[0, 0])
    ax2 = fig.add_subplot(gs[0, 1])
    ax3 = fig.add_subplot(gs[0, 2])
    ax4 = fig.add_subplot(gs[1, 0])
    ax5 = fig.add_subplot(gs[1, 1])
    ax6 = fig.add_subplot(gs[1, 2])

    ratio_gaus = np.array(a_gaus) / np.array(a_max)
    ratio_para = np.array(a_para) / np.array(a_max)
    ratio_lorentz = np.array(a_lorentz) / np.array(a_max)
    ratio_voigt = np.array(a_voigt) / np.array(a_max)
    ratio_spline = np.array(a_spline) / np.array(a_max)
    ratio_skewG = np.array(a_skewG) / np.array(a_max)

    mae_tot_para = np.mean(para_mae)
    mae_tot_gaus = np.mean(gaus_mae)
    mae_tot_lorentz = np.mean(lorentz_mae)
    mae_tot_voigt = np.mean(voigt_mae)
    mae_tot_spline = np.mean(spline_mae)
    mae_tot_skewG = np.mean(skewG_mae)

    rmse_tot_para = np.mean(para_rmse)
    rmse_tot_gaus = np.mean(gaus_rmse)
    rmse_tot_lorentz = np.mean(lorentz_rmse)
    rmse_tot_voigt = np.mean(voigt_rmse)
    rmse_tot_spline = np.mean(spline_rmse)
    rmse_tot_skewG = np.mean(skewG_rmse)

    #ax1.scatter(a_max, ratio_para, c='r', marker='d', s=20, edgecolors='black', label=r'Parabola / A$_{max}$' + '\nMAE = ' + str(round(mae_tot_para,4)) + " RMSE = " + str(round(rmse_tot_para,4)))
    ax1.scatter(a_max, ratio_para, c='r', marker='D', s=20, edgecolors='black', label=r'Parabola' + '\nRMSE = ' + str(round(rmse_tot_para,4)))
    ax1.axhline(1, color='black', linestyle='dashed', linewidth=1)
    ax1.set_xlabel(r"$A_{\rm{Sa}}$ / mV", fontsize=24)
    ax1.set_ylabel(r"$A_{\rm{para}}$ / $A_{\rm{Sa}}$", fontsize=24)
    ax1.set_ylim(0.9, 1.1)
    ax1.set_xlim(0,600)
    ax1.grid()
    ax1.tick_params(axis="both", labelsize=20)
    ax1.legend(fontsize=20, loc="upper left")

    #ax2.scatter(a_max, ratio_gaus, c='g', marker='d', s=20, edgecolors='black', label=r'Gaussian / A$_{\rm{Sa}}$' + '\nMAE = ' + str(round(mae_tot_para,4)) + " RMSE = " + str(round(rmse_tot_para,4)))
    ax2.scatter(a_max, ratio_gaus, c='g', marker='D', s=20, edgecolors='black', label=r'Gaussian' + '\nRMSE = ' + str(round(rmse_tot_para,4)))
    ax2.axhline(1, color='black', linestyle='dashed', linewidth=1)
    ax2.set_xlabel(r"$A_{\rm{Sa}}$ / mV", fontsize=24)
    ax2.set_ylabel(r"$A_{\rm{Gaus}}$ / $A_{\rm{Sa}}$", fontsize=24)
    ax2.set_ylim(0.9, 1.1)
    ax2.set_xlim(0,600)
    ax2.grid()
    ax2.tick_params(axis="both", labelsize=20)
    ax2.legend(fontsize=20, loc="lower right")

    #ax3.scatter(a_max, ratio_lorentz, c='blue', marker='d', s=20, edgecolors='black', label=r'Lorentz / $A_{\rm{Sa}}$' + '\nMAE = ' + str(round(mae_tot_lorentz,4)) + " RMSE = " + str(round(rmse_tot_lorentz,4)))
    ax3.scatter(a_max, ratio_lorentz, c='blue', marker='D', s=20, edgecolors='black', label=r'Lorentz' + '\nRMSE = ' + str(round(rmse_tot_lorentz,4)))
    ax3.axhline(1, color='black', linestyle='dashed', linewidth=1)
    ax3.set_xlabel(r"$A_{\rm{Sa}}$ / mV", fontsize=24)
    ax3.set_ylabel(r"$A_{\rm{Lorentz}}$ / $A_{\rm{Sa}}$", fontsize=24)
    ax3.set_ylim(0.9, 1.1)
    ax3.set_xlim(0,600)
    ax3.grid()
    ax3.tick_params(axis="both", labelsize=20)
    ax3.legend(fontsize=20, loc="lower right")

    #ax4.scatter(a_max, ratio_voigt, c='orange', marker='d', s=20, edgecolors='black', label=r'Voigt / $A_{\rm{Sa}}$' + '\nMAE = ' + str(round(mae_tot_voigt,4)) + " RMSE = " + str(round(rmse_tot_voigt,4)))
    ax4.scatter(a_max, ratio_voigt, c='orange', marker='D', s=20, edgecolors='black', label=r'Voigt' + '\nRMSE = ' + str(round(rmse_tot_voigt,4)))
    ax4.axhline(1, color='black', linestyle='dashed', linewidth=1)
    ax4.set_xlabel(r"$A_{\rm{Sa}}$ / mV", fontsize=24)
    ax4.set_ylabel(r"$A_{\rm{Voigt}}$ / $A_{\rm{Sa}}$", fontsize=24)
    ax4.set_ylim(0.9, 1.1)
    ax4.set_xlim(0,600)
    ax4.grid()
    ax4.tick_params(axis="both", labelsize=20)
    ax4.legend(fontsize=20, loc="lower right")

    #ax5.scatter(a_max, ratio_spline, c='purple', marker='d', s=20, edgecolors='black', label=r'Spline / $A_{\rm{Sa}}$' + '\nMAE = ' + str(round(mae_tot_spline,4)) + " RMSE = " + str(round(rmse_tot_spline,4)))
    ax5.scatter(a_max, ratio_spline, c='purple', marker='D', s=20, edgecolors='black', label=r'Interpolated spline' + '\nRMSE ≝ ' + str(round(rmse_tot_spline,4)))
    ax5.axhline(1, color='black', linestyle='dashed', linewidth=1)
    ax5.set_xlabel(r"$A_{\rm{Sa}}$ / mV", fontsize=24)
    ax5.set_ylabel(r"$A_{\rm{spline}}$ / $A_{\rm{Sa}}$", fontsize=24)
    ax5.set_ylim(0.9, 1.1)
    ax5.set_xlim(0,600)
    ax5.grid()
    ax5.tick_params(axis="both", labelsize=20)
    ax5.legend(fontsize=20, loc="lower right")

    #ax6.scatter(a_max, ratio_skewG, c='brown', marker='D', s=20, edgecolors='black', label=r'Skewed Gaussian / $A_{\rm{Sa}}$' + '\nMAE = ' + str(round(mae_tot_skewG,4)) + " RMSE = " + str(round(rmse_tot_skewG,4)))
    ax6.scatter(a_max, ratio_skewG, c='brown', marker='D', s=20, edgecolors='black', label=r'Skewed Gaussian' + '\nRMSE = ' + str(round(rmse_tot_skewG,4)))
    ax6.axhline(1, color='black', linestyle='dashed', linewidth=1)
    ax6.set_xlabel(r"$A_{\rm{Sa}}$ / mV", fontsize=24)
    ax6.set_ylabel(r"$A_{\rm{skewed}}$ / $A_{\rm{Sa}}$", fontsize=24)
    ax6.set_ylim(0.9, 1.1)
    ax6.set_xlim(0,600)
    ax6.grid()
    ax6.tick_params(axis="both", labelsize=20)
    ax6.legend(fontsize=20, loc="lower right")

    #fig.suptitle(f"Total {len(a_max)} signal events", fontsize=16, fontweight='bold')
    plt.tight_layout(rect=[0, 0, 1, 0.96])
    plt.savefig("./amplitude_to_sampling_trends.png",dpi=300,facecolor='w')
    plt.clf()

    '''
    fig = plt.figure(figsize=(24, 18))
    gs = gridspec.GridSpec(6, 3, height_ratios=[1, 1, 1, 1, 1, 1], width_ratios=[1, 1, 1])

    ax1 = fig.add_subplot(gs[0, 0])
    ax2 = fig.add_subplot(gs[1, 0])
    ax3 = fig.add_subplot(gs[2, 0])
    ax4 = fig.add_subplot(gs[3, 0])
    ax5 = fig.add_subplot(gs[4, 0])
    ax6 = fig.add_subplot(gs[5, 0])

    ax7 = fig.add_subplot(gs[0, 1])
    ax8 = fig.add_subplot(gs[1, 1])
    ax9 = fig.add_subplot(gs[2, 1])
    ax10 = fig.add_subplot(gs[3, 1])
    ax11 = fig.add_subplot(gs[4, 1])
    ax12 = fig.add_subplot(gs[5, 1])

    ax13 = fig.add_subplot(gs[0, 2])
    ax14 = fig.add_subplot(gs[1, 2])
    ax15 = fig.add_subplot(gs[2, 2])
    ax16 = fig.add_subplot(gs[3, 2])
    ax17 = fig.add_subplot(gs[4, 2])
    ax18 = fig.add_subplot(gs[5, 2])

    ratio_gaus = np.array(a_gaus) / np.array(a_max)
    ratio_para = np.array(a_para) / np.array(a_max)
    ratio_lorentz = np.array(a_lorentz) / np.array(a_max)
    ratio_voigt = np.array(a_voigt) / np.array(a_max)
    ratio_spline = np.array(a_spline) / np.array(a_max)
    ratio_landau = np.array(a_landau) / np.array(a_max)

    mae_tot_para = np.mean(para_mae)
    mae_tot_gaus = np.mean(gaus_mae)
    mae_tot_lorentz = np.mean(lorentz_mae)
    mae_tot_voigt = np.mean(voigt_mae)
    mae_tot_spline = np.mean(spline_mae)
    mae_tot_landau = np.mean(landau_mae)

    rmse_tot_para = np.mean(para_rmse)
    rmse_tot_gaus = np.mean(gaus_rmse)
    rmse_tot_lorentz = np.mean(lorentz_rmse)
    rmse_tot_voigt = np.mean(voigt_rmse)
    rmse_tot_spline = np.mean(spline_rmse)
    rmse_tot_landau = np.mean(landau_rmse)

    ax1.scatter(a_max, ratio_para, c='r', marker='d', s=20, edgecolors='black', label=r'Parabola / A$_{max}$' + '\nMAE = ' + str(round(mae_tot_para,4)) + " RMSE = " + str(round(rmse_tot_para,4)))
    ax1.axhline(1, color='black', linestyle='dashed', linewidth=1)
    ax1.set_xlabel(r"A$_{max}$ / mV", fontsize=14)
    ax1.set_ylabel(r"A$_{para}$ / A$_{max}$", fontsize=14)
    ax1.set_ylim(0.9, 1.1)
    ax1.legend(fontsize=14)

    ax2.scatter(a_max, ratio_gaus, c='g', marker='d', s=20, edgecolors='black', label=r'Gaussian / A$_{max}$' + '\nMAE = ' + str(round(mae_tot_para,4)) + " RMSE = " + str(round(rmse_tot_para,4)))
    ax2.axhline(1, color='black', linestyle='dashed', linewidth=1)
    ax2.set_xlabel(r"A$_{max}$ / mV", fontsize=14)
    ax2.set_ylabel(r"A$_{gaus}$ / A$_{max}$", fontsize=14)
    ax2.set_ylim(0.9, 1.1)
    ax2.legend(fontsize=14)

    ax3.scatter(a_max, ratio_lorentz, c='blue', marker='d', s=20, edgecolors='black', label=r'Lorentz / A$_{max}$' + '\nMAE = ' + str(round(mae_tot_lorentz,4)) + " RMSE = " + str(round(rmse_tot_lorentz,4)))
    ax3.axhline(1, color='black', linestyle='dashed', linewidth=1)
    ax3.set_xlabel(r"A$_{max}$ / mV", fontsize=14)
    ax3.set_ylabel(r"A$_{lorentz}$ / A$_{max}$", fontsize=14)
    ax3.set_ylim(0.9, 1.1)
    ax3.legend(fontsize=14)

    ax4.scatter(a_max, ratio_voigt, c='orange', marker='d', s=20, edgecolors='black', label=r'Voigt / A$_{max}$' + '\nMAE = ' + str(round(mae_tot_voigt,4)) + " RMSE = " + str(round(rmse_tot_voigt,4)))
    ax4.axhline(1, color='black', linestyle='dashed', linewidth=1)
    ax4.set_xlabel(r"A$_{max}$ / mV", fontsize=14)
    ax4.set_ylabel(r"A$_{voigt}$ / A$_{max}$", fontsize=14)
    ax4.set_ylim(0.9, 1.1)
    ax4.legend(fontsize=14)

    ax5.scatter(a_max, ratio_spline, c='purple', marker='d', s=20, edgecolors='black', label=r'Spline / A$_{max}$' + '\nMAE = ' + str(round(mae_tot_spline,4)) + " RMSE = " + str(round(rmse_tot_spline,4)))
    ax5.axhline(1, color='black', linestyle='dashed', linewidth=1)
    ax5.set_xlabel(r"A$_{max}$ / mV", fontsize=14)
    ax5.set_ylabel(r"A$_{spline}$ / A$_{max}$", fontsize=14)
    ax5.set_ylim(0.9, 1.1)
    ax5.legend(fontsize=14)

    ax6.scatter(a_max, ratio_landau, c='brown', marker='D', s=20, edgecolors='black', label=r'Landau / A$_{max}$' + '\nMAE = ' + str(round(mae_tot_landau,4)) + " RMSE = " + str(round(rmse_tot_landau,4)))
    ax6.axhline(1, color='black', linestyle='dashed', linewidth=1)
    ax6.set_xlabel(r"A$_{max}$ / mV", fontsize=14)
    ax6.set_ylabel(r"A$_{landau}$ / A$_{max}$", fontsize=14)
    ax6.set_ylim(0.9, 1.1)
    ax6.legend(fontsize=14)

    ax7.hist(para_mae, bins=50, density=False, alpha=0.9, color='r', edgecolor='black', label='Parabolic', histtype='stepfilled', linewidth=1.5)
    ax7.set_xlabel(r"Mean Absolute Error",fontsize=14)
    ax7.set_ylabel(r"Counts",fontsize=14)
    ax7.set_xlim(0,0.02)
    ax7.set_yscale("log")
    ax7.legend(fontsize=14)

    ax8.hist(gaus_mae, bins=50, density=False, alpha=0.9, color='g', edgecolor='black', label='Gaussian', histtype='stepfilled', linewidth=1.5)
    ax8.set_xlabel(r"Mean Absolute Error",fontsize=14)
    ax8.set_ylabel(r"Counts",fontsize=14)
    ax8.set_xlim(0,0.02)
    ax8.set_yscale("log")
    ax8.legend(fontsize=14)

    ax9.hist(lorentz_mae, bins=50, density=False, alpha=0.9, color='blue', edgecolor='black', label='Lorentz', histtype='stepfilled', linewidth=1.5)
    ax9.set_xlabel(r"Mean Absolute Error",fontsize=14)
    ax9.set_ylabel(r"Counts",fontsize=14)
    ax9.set_xlim(0,0.02)
    ax9.set_yscale("log")
    ax9.legend(fontsize=14)

    ax10.hist(voigt_mae, bins=50, density=False, alpha=0.9, color='orange', edgecolor='black', label='Voigt', histtype='stepfilled', linewidth=1.5)
    ax10.set_xlabel(r"Mean Absolute Error",fontsize=14)
    ax10.set_ylabel(r"Counts",fontsize=14)
    ax10.set_xlim(0,0.02)
    ax10.set_yscale("log")
    ax10.legend(fontsize=14)

    ax11.hist(spline_mae, bins=50, density=False, alpha=0.9, color='purple', edgecolor='black', label='Spline', histtype='stepfilled', linewidth=1.5)
    ax11.set_xlabel(r"Mean Absolute Error",fontsize=14)
    ax11.set_ylabel(r"Counts",fontsize=14)
    ax11.set_xlim(0,0.02)
    ax11.set_yscale("log")
    ax11.legend(fontsize=14)

    ax12.hist(landau_mae, bins=50, density=False, alpha=0.9, color='brown', edgecolor='black', label='Landau', histtype='stepfilled', linewidth=1.5)
    ax12.set_xlabel(r"Mean Absolute Error",fontsize=14)
    ax12.set_ylabel(r"Counts",fontsize=14)
    ax12.set_xlim(0,0.02)
    ax12.set_yscale("log")
    ax12.legend(fontsize=14)
    
    ax13.hist(para_rmse, bins=50, density=False, alpha=0.9, color='r', edgecolor='black', label='Parabolic', histtype='stepfilled', linewidth=1.5)
    ax13.set_xlabel(r"RMS Error",fontsize=14)
    ax13.set_ylabel(r"Counts",fontsize=14)
    ax13.set_xlim(0,0.02)
    ax13.set_yscale("log")
    ax13.legend(fontsize=14)

    ax14.hist(gaus_rmse, bins=50, density=False, alpha=0.9, color='g', edgecolor='black', label='Gaussian', histtype='stepfilled', linewidth=1.5)
    ax14.set_xlabel(r"RMS Error",fontsize=14)
    ax14.set_ylabel(r"Counts",fontsize=14)
    ax14.set_xlim(0,0.02)
    ax14.set_yscale("log")
    ax14.legend(fontsize=14)

    ax15.hist(lorentz_rmse, bins=50, density=False, alpha=0.9, color='blue', edgecolor='black', label='Lorentz', histtype='stepfilled', linewidth=1.5)
    ax15.set_xlabel(r"RMS Error",fontsize=14)
    ax15.set_ylabel(r"Counts",fontsize=14)
    ax15.set_xlim(0,0.02)
    ax15.set_yscale("log")
    ax15.legend(fontsize=14)

    ax16.hist(voigt_rmse, bins=50, density=False, alpha=0.9, color='orange', edgecolor='black', label='Voigt', histtype='stepfilled', linewidth=1.5)
    ax16.set_xlabel(r"RMS Error",fontsize=14)
    ax16.set_ylabel(r"Counts",fontsize=14)
    ax16.set_xlim(0,0.02)
    ax16.set_yscale("log")
    ax16.legend(fontsize=14)

    ax17.hist(spline_rmse, bins=50, density=False, alpha=0.9, color='purple', edgecolor='black', label='Spline', histtype='stepfilled', linewidth=1.5)
    ax17.set_xlabel(r"RMS Error",fontsize=14)
    ax17.set_ylabel(r"Counts",fontsize=14)
    ax17.set_xlim(0,0.02)
    ax17.set_yscale("log")
    ax17.legend(fontsize=14)

    ax18.hist(landau_rmse, bins=50, density=False, alpha=0.9, color='brown', edgecolor='black', label='Landau', histtype='stepfilled', linewidth=1.5)
    ax18.set_xlabel(r"RMS Error",fontsize=14)
    ax18.set_ylabel(r"Counts",fontsize=14)
    ax18.set_xlim(0,0.02)
    ax18.set_yscale("log")
    ax18.legend(fontsize=14)

    fig.suptitle(f"Total {len(a_max)} signal events", fontsize=16, fontweight='bold')
    plt.tight_layout(rect=[0, 0, 1, 0.96])
    plt.savefig("./amplitude_big_data_analysis.png",dpi=300,facecolor='w')
    plt.clf()
    '''
  if cfd_studies:

    cfd30_para = []
    cfd30_gaus = []
    cfd30_lorentz = []
    cfd30_voigt = []
    cfd30_spline = []
    cfd30_landau = []

    cfd30_para_mcp = []
    cfd30_gaus_mcp = []
    cfd30_lorentz_mcp = []
    cfd30_voigt_mcp = []
    cfd30_spline_mcp = []
    cfd30_landau_mcp = []

    for j in range(num_curves):
      idx_para = np.where(np.array(reshaped_ampl_data[j]) > 0.0002*a_para[j])[0][0]
      print(len(reshaped_ampl_data[j]))
      print(len(a_gaus[j]))
      idx_gaus = np.where(np.array(reshaped_ampl_data[j]) > 0.0002*a_gaus[j])[0][0]
      idx_lorentz = np.where(np.array(reshaped_ampl_data[j]) > 0.0002*a_lorentz[j])[0][0]
      idx_voigt = np.where(np.array(reshaped_ampl_data[j]) > 0.0002*a_voigt[j])[0][0]
      idx_spline = np.where(np.array(reshaped_ampl_data[j]) > 0.0002*a_spline[j])[0][0]
      idx_landau = np.where(np.array(reshaped_ampl_data[j]) > 0.0002*a_landau[j])[0][0]
      time_array_event = np.array(reshaped_time_data[j])
      mean_time_para = np.mean(time_array_event[[idx_para-1, idx_para]])
      mean_time_gaus = np.mean(time_array_event[[idx_gaus-1, idx_gaus]])
      mean_time_lorentz = np.mean(time_array_event[[idx_lorentz-1, idx_lorentz]])
      mean_time_voigt = np.mean(time_array_event[[idx_voigt-1, idx_voigt]])
      mean_time_spline = np.mean(time_array_event[[idx_spline-1, idx_spline]])
      mean_time_landau = np.mean(time_array_event[[idx_landau-1, idx_landau]])
      cfd30_para.append(mean_time_para)
      cfd30_gaus.append(mean_time_gaus)
      cfd30_lorentz.append(mean_time_lorentz)
      cfd30_voigt.append(mean_time_voigt)
      cfd30_spline.append(mean_time_spline)
      cfd30_landau.append(mean_time_landau)

      if time_res_calc:
        idx_para_mcp = np.where(np.array(rad_mcp[j]) > 0.0002*a_para_mcp[j])[0][0]
        idx_gaus_mcp = np.where(np.array(rad_mcp[j]) > 0.0002*a_gaus_mcp[j])[0][0]
        idx_lorentz_mcp = np.where(np.array(rad_mcp[j]) > 0.0002*a_lorentz_mcp[j])[0][0]
        idx_voigt_mcp = np.where(np.array(rad_mcp[j]) > 0.0002*a_voigt_mcp[j])[0][0]
        idx_spline_mcp = np.where(np.array(rad_mcp[j]) > 0.0002*a_spline_mcp[j])[0][0]
        idx_landau_mcp = np.where(np.array(rad_mcp[j]) > 0.0002*a_landau_mcp[j])[0][0]
        time_array_event_mcp = np.array(rtd_mcp[j])
        mean_time_para_mcp = 0.5*(time_array_event_mcp[idx_para_mcp-1] + time_array_event_mcp[idx_para_mcp]) #np.mean(time_array_event_mcp[[idx_para_mcp-1, idx_para_mcp]])
        mean_time_gaus_mcp = np.mean(time_array_event_mcp[[idx_gaus_mcp-1, idx_gaus_mcp]])
        mean_time_lorentz_mcp = np.mean(time_array_event_mcp[[idx_lorentz_mcp-1, idx_lorentz_mcp]])
        mean_time_voigt_mcp = np.mean(time_array_event_mcp[[idx_voigt_mcp-1, idx_voigt_mcp]])
        mean_time_spline_mcp = np.mean(time_array_event_mcp[[idx_spline_mcp-1, idx_spline_mcp]])
        mean_time_landau_mcp = np.mean(time_array_event_mcp[[idx_landau_mcp-1, idx_landau_mcp]])
        cfd30_para_mcp.append(mean_time_para_mcp)
        cfd30_gaus_mcp.append(mean_time_gaus_mcp)
        cfd30_lorentz_mcp.append(mean_time_lorentz_mcp)
        cfd30_voigt_mcp.append(mean_time_voigt_mcp)
        cfd30_spline_mcp.append(mean_time_spline_mcp)
        cfd30_landau_mcp.append(mean_time_landau_mcp)

    label_cfd = "CFD@30%"

    if time_res_calc:
      
      fig2oneoff, ax2oneoff = plt.subplots(figsize=(8, 6))
      ax2oneoff.hist([cfd30_para,cfd30_para_mcp], bins=150, range=(-2, 1), alpha=0.6, color=['blue','purple'], edgecolor='black', label=['cfd30_para','cfd30_para_mcp'], histtype='stepfilled', linewidth=1.5)
      dut_fit_para = gaussian_fit_binned_data(cfd30_para, "PARA_DUT")
      mcp_fit_para = gaussian_fit_binned_data(cfd30_para_mcp, "PARA_MCP")
      #ax2oneoff.hist([cfd30_data,cfd30_mcp_data], bins=150, range=(-2, 1), color=['red','orange'], edgecolor='black', label=['cfd30_data','cfd30_mcp_data'], histtype='stepfilled', linewidth=1.5, alpha=0.2)
      #ax2oneoff.legend()
      #fig2oneoff.savefig(f"datavpara.png", dpi=300, bbox_inches='tight')

      tr_data = np.array(cfd30_data) - np.array(cfd30_mcp_data)
      tr_para = np.array(cfd30_para) - np.array(cfd30_para_mcp)
      my_proper_dist = np.subtract(cfd30_para, cfd30_para_mcp)
      #ax2oneoff.hist(my_proper_dist, bins=150, range=(-2, 1), alpha=0.6, color='yellow', edgecolor='black', label='Time resolution', histtype='stepfilled', linewidth=1.5)
      ax2oneoff.legend()
      fig2oneoff.savefig(f"datavpara.png", dpi=300, bbox_inches='tight')
      '''
      for i in range(len(cfd30_para)):
        print("\n")
        print(str(cfd30_para[i]))
        print(str(cfd30_para_mcp[i]) + " -")
        print("---------------------")
        print(str(cfd30_para[i] - cfd30_para_mcp[i]))
      '''
      tr_gaus = np.array(cfd30_gaus) - np.array(cfd30_gaus_mcp)
      tr_lorentz = np.array(cfd30_lorentz) - np.array(cfd30_lorentz_mcp)
      tr_voigt = np.array(cfd30_voigt) - np.array(cfd30_voigt_mcp)
      tr_spline = np.array(cfd30_spline) - np.array(cfd30_spline_mcp)
      tr_landau = np.array(cfd30_landau) - np.array(cfd30_landau_mcp)

      label_cfd = r"$\sigma_{t}^{30%}$"

      data_tr_params = gaussian_fit_binned_data(tr_data, "Data")
      para_tr_params = gaussian_fit_binned_data(tr_para, "Parabolic")
      gaus_tr_params = gaussian_fit_binned_data(tr_gaus, "Gaussian")
      lorentz_tr_params = gaussian_fit_binned_data(tr_lorentz, "Lorentz")
      voigt_tr_params = gaussian_fit_binned_data(tr_voigt, "Voigt")
      spline_tr_params = gaussian_fit_binned_data(tr_spline, "Interpolated spline")
      landau_tr_params = gaussian_fit_binned_data(tr_landau, "Landau")

    fig, axes = plt.subplots(2, 3, figsize=(18, 12))
    rms_diff_para = np.sqrt(np.mean((np.array(cfd30_data) - np.array(cfd30_para)) ** 2)).round(3)
    axes[0,0].hist(cfd30_data, bins=400,range=(-0.8,-0.4),color='gray',edgecolor='black',label=label_cfd + r"(A$_{max}$)")
    axes[0,0].hist(cfd30_para, bins=400,range=(-0.8,-0.4),color='r',edgecolor='black',alpha=0.4,label=label_cfd + r"(A$_{para}$)" + "\n" + r"$\Delta_{RMS}$ = " + str(rms_diff_para))

    rms_diff_gaus = np.sqrt(np.mean((np.array(cfd30_data) - np.array(cfd30_gaus)) ** 2)).round(3)
    axes[0,1].hist(cfd30_data, bins=400,range=(-0.8,-0.4),color='gray',edgecolor='black',label=label_cfd + r"(A$_{max}$)")
    axes[0,1].hist(cfd30_gaus, bins=400,range=(-0.8,-0.4),color='g',edgecolor='black',alpha=0.4,label=label_cfd + r"(A$_{Gaus}$)" + "\n" + r"$\Delta_{RMS}$ = " + str(rms_diff_gaus))

    rms_diff_lorentz = np.sqrt(np.mean((np.array(cfd30_data) - np.array(cfd30_lorentz)) ** 2)).round(3)
    axes[1,0].hist(cfd30_data, bins=400,range=(-0.8,-0.4),color='gray',edgecolor='black',label=label_cfd + r"(A$_{max}$)")
    axes[1,0].hist(cfd30_lorentz, bins=400,range=(-0.8,-0.4),color='blue',edgecolor='black',alpha=0.4,label=label_cfd + r"(A$_{Lorentz}$)" + "\n" + r"$\Delta_{RMS}$ = " + str(rms_diff_lorentz))

    rms_diff_voigt = np.sqrt(np.mean((np.array(cfd30_data) - np.array(cfd30_voigt)) ** 2)).round(3)
    axes[1,1].hist(cfd30_data, bins=400,range=(-0.8,-0.4),color='gray',edgecolor='black',label=label_cfd + r"(A$_{max}$)")
    axes[1,1].hist(cfd30_voigt, bins=400,range=(-0.8,-0.4),color='orange',edgecolor='black',alpha=0.4,label=label_cfd + r"(A$_{Voigt}$)" + "\n" + r"$\Delta_{RMS}$ = " + str(rms_diff_voigt))

    rms_diff_spline = np.sqrt(np.mean((np.array(cfd30_data) - np.array(cfd30_spline)) ** 2)).round(3)
    axes[0,2].hist(cfd30_data, bins=400,range=(-0.8,-0.4),color='gray',edgecolor='black',label=label_cfd + r"(A$_{max}$)")
    axes[0,2].hist(cfd30_spline, bins=400,range=(-0.8,-0.4),color='purple',edgecolor='black',alpha=0.4,label=label_cfd + r"(A$_{spline}$)" + "\n" + r"$\Delta_{RMS}$ = " + str(rms_diff_spline))

    rms_diff_landau = np.sqrt(np.mean((np.array(cfd30_data) - np.array(cfd30_landau)) ** 2)).round(3)
    axes[1,2].hist(cfd30_data, bins=400,range=(-0.8,-0.4),color='gray',edgecolor='black',label=label_cfd + r"(A$_{max}$)")
    axes[1,2].hist(cfd30_landau, bins=400,range=(-0.8,-0.4),color='cyan',edgecolor='black',alpha=0.4,label=label_cfd + r"(A$_{Landau}$)" + "\n" + r"$\Delta_{RMS}$ = " + str(rms_diff_landau))

    if time_res_calc:
      x_tr_fit = np.linspace(-0.8, -0.4, 400)
      data_tr_fit = gaussian(x_tr_fit, *data_tr_params)
      para_tr_fit = gaussian(x_tr_fit, *para_tr_params)
      gaus_tr_fit = gaussian(x_tr_fit, *gaus_tr_params)
      lorentz_tr_fit = gaussian(x_tr_fit, *lorentz_tr_params)
      voigt_tr_fit = gaussian(x_tr_fit, *voigt_tr_params)
      spline_tr_fit = gaussian(x_tr_fit, *spline_tr_params)
      landau_tr_fit = gaussian(x_tr_fit, *landau_tr_params)

      mcp_tr_est = 5
      data_tr_val = np.sqrt((1000*data_tr_params[2])**2 - mcp_tr_est**2)
      para_tr_val = np.sqrt((1000*para_tr_params[2])**2 - mcp_tr_est**2)
      gaus_tr_val = np.sqrt((1000*gaus_tr_params[2])**2 - mcp_tr_est**2)
      lorentz_tr_val = np.sqrt((1000*lorentz_tr_params[2])**2 - mcp_tr_est**2)
      voigt_tr_val = np.sqrt((1000*voigt_tr_params[2])**2 - mcp_tr_est**2)
      spline_tr_val = np.sqrt((1000*spline_tr_params[2])**2 - mcp_tr_est**2)
      landau_tr_val = np.sqrt((1000*landau_tr_params[2])**2 - mcp_tr_est**2)

      for i in range(2):
        for j in range(3):
          axes[i,j].plot(x_fit, data_tr_fit, 'k--', linewidth=2, label=r"$\sigma_{tr}$ = " + str(round(data_tr_val, 1)) + " ps")

      axes[0,0].plot(x_fit, para_tr_fit, 'r', linewidth=2, label=r"$\sigma_{tr}^{para}$ = " + str(round(data_tr_val, 1)) + " ps")
      axes[0,1].plot(x_fit, gaus_tr_fit, 'g', linewidth=2, label=r"$\sigma_{tr}^{Gaus}$ = " + str(round(data_tr_val, 1)) + " ps")
      axes[1,0].plot(x_fit, lorentz_tr_fit, 'blue', linewidth=2, label=r"$\sigma_{tr}^{Lorentz}$ = " + str(round(data_tr_val, 1)) + " ps")
      axes[1,1].plot(x_fit, voigt_tr_fit, 'orange', linewidth=2, label=r"$\sigma_{tr}^{Voigt}$ = " + str(round(data_tr_val, 1)) + " ps")
      axes[0,2].plot(x_fit, spline_tr_fit, 'purple', linewidth=2, label=r"$\sigma_{tr}^{spline}$ = " + str(round(data_tr_val, 1)) + " ps")
      axes[1,2].plot(x_fit, landau_tr_fit, 'cyan', linewidth=2, label=r"$\sigma_{tr}^{Landau}$ = " + str(round(data_tr_val, 1)) + " ps")

    for i in range(2):
      for j in range(3):
        axes[i,j].set_xlabel(label_cfd + r"/ ns",fontsize=14)
        axes[i,j].set_ylabel(r"Events",fontsize=14)
        axes[i,j].set_xlim(-2.5, 0.5)
        axes[i,j].legend(fontsize=14)
        axes[i,j].grid(True, axis='both', linestyle='--', alpha=0.5)

    fig.suptitle(f"Total {len(a_max)} signal events", fontsize=16, fontweight='bold')
    plt.tight_layout(rect=[0, 0, 1, 0.96])
    if time_res_calc:
      plt.savefig("./timeres_30_260325.png",dpi=300,facecolor='w')
    else:
      plt.savefig("./cfd_30_260325.png",dpi=300,facecolor='w')
    plt.clf()




if __name__ == "__main__":
  main()
