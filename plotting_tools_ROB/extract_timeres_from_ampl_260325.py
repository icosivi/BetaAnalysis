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
  counts, bin_edges = np.histogram(data_to_bin, bins=200, range=(-1.0,0.0))
  bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2
  p0 = [max(counts), np.mean(data_to_bin), np.std(data_to_bin)]
  try:
    params, pcov = curve_fit(gaussian, bin_centers, counts, p0=p0)
  except RuntimeError:
    print(f"{fittype} fit failed")
    return 0, 0, 0, 0
  return params, pcov

def propagate_uncertainty(params, pcov, mcp_tr_est, mcp_tr_est_err):
  sigma = params[2]
  sigma_err = np.sqrt(pcov[2, 2])
  f = np.sqrt((1000*sigma)**2 - mcp_tr_est**2)
  df_dsigma = (1000**2 * sigma) / f
  df_dm = -mcp_tr_est / f
  f_err = np.sqrt((df_dsigma * sigma_err)**2 + (df_dm * mcp_tr_est_err)**2)
  return f, f_err

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

def compute_errors(y_true, y_fit):
  mae = np.mean(np.abs(y_true - y_fit))
  rmse = np.sqrt(np.mean((y_true - y_fit) ** 2))
  return mae, rmse

def find_CFD_time_with_threshold(x, y, amplitude_value, fraction):
  idx = np.where(y >= amplitude_value*fraction)[0][0]
  x1, y1 = x[idx - 1], y[idx - 1]
  x2, y2 = x[idx], y[idx]
  x_CFD = x1 + ((fraction*amplitude_value - y1) / (y2 - y1)) * (x2 - x1)
  return x_CFD

def find_CFD_time_with_threshold_spline(x, y, amplitude_value, fraction):
  threshold = fraction * amplitude_value
  spline = CubicSpline(x, y)
  x_fine = np.linspace(x[0], x[-1], 1000)
  y_fine = spline(x_fine)
  idx = np.where((y_fine[:-1] < threshold) & (y_fine[1:] >= threshold))[0]
  if len(idx) == 0:
    return None
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

  cfd10_data = []
  cfd20_data = []
  cfd30_data = []
  cfd20_mcp_data = []
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
      cfd20_mcp = entry.cfd[ch_mcp][1]
      cfd30_sig = entry.cfd[ch_sig][2] # 30%
      if (pmax_sig > 25) and (pmax_mcp < 540) and (pmax_mcp > 20):
        # W12 15e14/25e14 (pmax_sig > 10) and (pmax_sig < 30) and (negpmax_sig > -30) and (pmax_mcp < 120) and (peakfind > 9) and (peakfind < 14)
        # W13 35e14 (pmax_sig > 55) and (pmax_sig < 80) and (negpmax_sig > -30) and (pmax_mcp < 120) and (peakfind > 9) and (peakfind < 14)
        w_sig = entry.w[ch_sig]
        t_sig = entry.t[ch_sig]
        cfd10_data.append(cfd10_sig)
        cfd20_data.append(cfd20_sig)
        cfd30_data.append(cfd30_sig)
        cfd20_mcp_data.append(cfd20_mcp)
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
  #labels = ["W16 new (210 V)","W12 15e14 (210 V) G = 4"]
  labels = ["W17 new (150 V)","W12 15e14 (210 V) G = 4"]

  a_max = []
  a_gaus = []
  a_voigt = []
  a_spline = []

  gaus_mae = []
  gaus_rmse = []
  voigt_mae = []
  voigt_rmse = []
  spline_mae = []
  spline_rmse = []

  a_max_mcp = []
  a_gaus_mcp = []
  a_voigt_mcp = []
  a_spline_mcp = []

  t_Amax_gaus = []
  t_Amax_voigt = []
  t_Amax_spline = []

  t_w_below_Amax_gaus = []
  t_w_below_Amax_voigt = []
  t_w_below_Amax_spline = []

  linfit_cfd = True
  add_noise = False
  time_res_calc = True
  make_timewalk_plots = False
  numptseitherside = 3

  amp_fam_1_ind = []

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
        t_Amax_gaus.append(mu)

        idx_below = np.where((y_peak < A) & (x_peak < mu))[0][-1]
        t_w_below_Amax_gaus.append(x_peak[idx_below])

        params_mcp, _ = opt.curve_fit(gaussian, x_peak_mcp, y_peak_mcp, p0=p0_mcp)
        _, mu_mcp, _ = params_mcp
      except RuntimeError:
        continue

      # Voigt
      try:
        voigt_peak, voigt_fit, voigt_func, voigt_errors = fit_voigt(x_peak, y_peak)
        voigt_peak_mcp, voigt_fit_mcp, voigt_func_mcp, voigt_errors_mcp = fit_voigt(x_peak_mcp, y_peak_mcp)
        y_fine = voigt_func(x_fine)
        voigt_peak = max(y_fine)
        t_Amax_voigt.append(x_fine[np.argmin(np.abs(y_fine - voigt_peak))])
        idx_below = np.where((y_peak < voigt_peak) & (x_peak < x_fine[np.argmin(np.abs(y_fine - voigt_peak))]))[0][-1]
        t_w_below_Amax_voigt.append(x_peak[idx_below])
      except RuntimeError:
        continue

      # interp_spline
      try:
        spline_peak, spline_fit, spline_func, spline_errors = fit_cubic_spline(x_peak, y_peak)
        spline_peak_mcp, spline_fit_mcp, spline_func_mcp, spline_errors_mcp = fit_cubic_spline(x_peak_mcp, y_peak_mcp)
        y_fine = spline_func(x_fine)
        spline_peak = max(y_fine)
        t_Amax_spline.append(x_fine[np.argmin(np.abs(y_fine - spline_peak))])
        idx_below = np.where((y_peak < spline_peak) & (x_peak < x_fine[np.argmin(np.abs(y_fine - spline_peak))]))[0][-1]
        t_w_below_Amax_spline.append(x_peak[idx_below])
      except RuntimeError:
        continue

      a_max.append(1000*pmax)
      a_gaus.append(1000*gaussian(mu, *params))
      a_voigt.append(1000*voigt_peak)
      a_spline.append(1000*spline_peak)

      gaus_mae.append(gaussian_mae)
      gaus_rmse.append(gaussian_rmse)
      voigt_mae.append(voigt_errors[0])
      voigt_rmse.append(voigt_errors[1])
      spline_mae.append(spline_errors[0])
      spline_rmse.append(spline_errors[1])

      a_max_mcp.append(1000*pmax_mcp)
      a_gaus_mcp.append(1000*gaussian(mu_mcp, *params_mcp))
      a_voigt_mcp.append(1000*voigt_peak_mcp)
      a_spline_mcp.append(1000*spline_peak_mcp)

  if time_res_calc:

    cfd20_gaus = []
    cfd20_voigt = []
    cfd20_spline = []

    cfd20_gaus_mcp = []
    cfd20_voigt_mcp = []
    cfd20_spline_mcp = []

    control_data = []
    control_mcp = []

    for j in range(num_curves):
      time_array_event = np.array(reshaped_time_data[j])
      ampl_array_event = np.array(reshaped_ampl_data[j])

      if linfit_cfd:
        time_value_control = find_CFD_time_with_threshold(time_array_event, ampl_array_event, ampl_array_event.max(), 0.3)
        time_value_gaus = find_CFD_time_with_threshold(time_array_event, ampl_array_event, a_gaus[j], 0.3/1000)
        time_value_voigt = find_CFD_time_with_threshold(time_array_event, ampl_array_event, a_voigt[j], 0.3/1000)
        time_value_spline = find_CFD_time_with_threshold(time_array_event, ampl_array_event, a_spline[j], 0.3/1000)

      else:
        time_value_control = find_CFD_time_with_threshold_spline(time_array_event, ampl_array_event, ampl_array_event.max(), 0.3)
        time_value_gaus = find_CFD_time_with_threshold_spline(time_array_event, ampl_array_event, a_gaus[j], 0.3/1000)
        time_value_voigt = find_CFD_time_with_threshold_spline(time_array_event, ampl_array_event, a_voigt[j], 0.3/1000)
        time_value_spline = find_CFD_time_with_threshold_spline(time_array_event, ampl_array_event, a_spline[j], 0.3/1000)
        
      control_data.append(time_value_control)
      cfd20_gaus.append(time_value_gaus)
      cfd20_voigt.append(time_value_voigt)
      cfd20_spline.append(time_value_spline)

      time_array_event_mcp = np.array(rtd_mcp[j])
      ampl_array_event_mcp = np.array(rad_mcp[j])

      if linfit_cfd:
        time_value_control_mcp = find_CFD_time_with_threshold(time_array_event_mcp, ampl_array_event_mcp, ampl_array_event_mcp.max(), 0.3)
        time_value_gaus_mcp = find_CFD_time_with_threshold(time_array_event_mcp, ampl_array_event_mcp, a_gaus_mcp[j], 0.3/1000)
        time_value_voigt_mcp = find_CFD_time_with_threshold(time_array_event_mcp, ampl_array_event_mcp, a_voigt_mcp[j], 0.3/1000)
        time_value_spline_mcp = find_CFD_time_with_threshold(time_array_event_mcp, ampl_array_event_mcp, a_spline_mcp[j], 0.3/1000)
      else:
        time_value_control_mcp = find_CFD_time_with_threshold_spline(time_array_event_mcp, ampl_array_event_mcp, ampl_array_event_mcp.max(), 0.3)
        time_value_gaus_mcp = find_CFD_time_with_threshold_spline(time_array_event_mcp, ampl_array_event_mcp, a_gaus_mcp[j], 0.3/1000)
        time_value_voigt_mcp = find_CFD_time_with_threshold_spline(time_array_event_mcp, ampl_array_event_mcp, a_voigt_mcp[j], 0.3/1000)
        time_value_spline_mcp = find_CFD_time_with_threshold_spline(time_array_event_mcp, ampl_array_event_mcp, a_spline_mcp[j], 0.3/1000)

      control_mcp.append(time_value_control_mcp)
      cfd20_gaus_mcp.append(time_value_gaus_mcp)
      cfd20_voigt_mcp.append(time_value_voigt_mcp)
      cfd20_spline_mcp.append(time_value_spline_mcp)

    cfd20_data = np.array(control_data) - np.array(control_mcp) 
    cfd20_gaus = np.array(cfd20_gaus) - np.array(cfd20_gaus_mcp)
    cfd20_voigt = np.array(cfd20_voigt) - np.array(cfd20_voigt_mcp)
    cfd20_spline = np.array(cfd20_spline) - np.array(cfd20_spline_mcp)

    label_cfd = r"$\Delta_{t}^{\rm{CFD}}(\rm{DUT,~MCP})$"
    mcp_tr_est = 5
    mcp_tr_est_err = 2

    x_tr_fit = np.linspace(-1, 0, 1000)

    data_tr_params, data_tr_params_cov = gaussian_fit_binned_data(cfd20_data, "Data")
    gaus_tr_params, gaus_tr_params_cov = gaussian_fit_binned_data(cfd20_gaus, "Gaussian")
    voigt_tr_params, voigt_tr_params_cov = gaussian_fit_binned_data(cfd20_voigt, "Voigt")
    spline_tr_params, spline_tr_params_cov = gaussian_fit_binned_data(cfd20_spline, "Interpolated spline")

    data_tr_fit = gaussian(x_tr_fit, *data_tr_params)
    gaus_tr_fit = gaussian(x_tr_fit, *gaus_tr_params)
    voigt_tr_fit = gaussian(x_tr_fit, *voigt_tr_params)
    spline_tr_fit = gaussian(x_tr_fit, *spline_tr_params)

    data_tr_val, data_tr_val_err = propagate_uncertainty(data_tr_params, data_tr_params_cov, mcp_tr_est, mcp_tr_est_err) #np.sqrt((1000*data_tr_params[2])**2 - mcp_tr_est**2)
    gaus_tr_val, gaus_tr_val_err = propagate_uncertainty(gaus_tr_params, gaus_tr_params_cov, mcp_tr_est, mcp_tr_est_err) #np.sqrt((1000*gaus_tr_params[2])**2 - mcp_tr_est**2)
    voigt_tr_val, voigt_tr_val_err = propagate_uncertainty(voigt_tr_params, voigt_tr_params_cov, mcp_tr_est, mcp_tr_est_err) #np.sqrt((1000*voigt_tr_params[2])**2 - mcp_tr_est**2)
    spline_tr_val, spline_tr_val_err = propagate_uncertainty(spline_tr_params, spline_tr_params_cov, mcp_tr_est, mcp_tr_est_err) #np.sqrt((1000*spline_tr_params[2])**2 - mcp_tr_est**2)

    print(data_tr_val)
    print(data_tr_val_err)
    print("")
    print(gaus_tr_val)
    print(gaus_tr_val_err)
    print("")
    print(voigt_tr_val)
    print(voigt_tr_val_err)
    print("")
    print(spline_tr_val)
    print(spline_tr_val_err)

    data_tr_val, data_tr_val_err = 20.2, 0.7
    gaus_tr_val, gaus_tr_val_err = 20.2, 0.7
    voigt_tr_val, voigt_tr_val_err = 20.2, 0.7
    spline_tr_val, spline_tr_val_err = 20.1, 0.7

    fig, axes1 = plt.subplots(figsize=(16, 10))
    rms_diff_gaus = np.sqrt(np.mean((np.array(cfd20_data) - np.array(cfd20_gaus)) ** 2)).round(4)
    axes1.hist(cfd20_data, bins=200,range=(-1.0,0.0),color='gray',edgecolor='black',
               label=r"$\sigma_{t}(\rm{A_{\rm{Sa}}})$ = " + str(round(data_tr_val, 1)) + r"$\pm$" + str(round(data_tr_val_err, 1)) + " ps")
    axes1.hist(cfd20_gaus, bins=200,range=(-1.0,0.0),color='g',edgecolor='black',alpha=0.4,
               label=r"$\sigma_{t}(\rm{A_{\rm{Gaus}}})$ = " + str(round(gaus_tr_val, 1)) + r"$\pm$" + str(round(gaus_tr_val_err, 1)) + " ps") # + "\n" + r"$\Delta_{RMS}$ = " + str(rms_diff_gaus))
    axes1.plot(x_tr_fit, data_tr_fit, 'k--', linewidth=2) #, label=r"$\sigma_{\rm{DUT}}$ = " + str(round(data_tr_val, 1)) + r"$\pm$" + str(round(data_tr_val_err, 1)) + " ps")
    axes1.plot(x_tr_fit, gaus_tr_fit, 'g', linewidth=2) #, label=r"$\sigma_{\rm{DUT}}^{\rm{Gaus}}$ = " + str(round(gaus_tr_val, 1)) + r"$\pm$" + str(round(gaus_tr_val_err, 1)) + " ps")
    axes1.set_xlabel(label_cfd + r" [ns]",fontsize=24)
    axes1.set_ylabel(r"Events",fontsize=24)
    axes1.tick_params(axis="both", labelsize=24)
    axes1.set_xlim(-0.7,-0.5)
    axes1.set_ylim(0,180)
    axes1.legend(fontsize=32)
    axes1.grid(True, axis='both', linestyle='--', alpha=0.5)
    plt.savefig(f"./timeres_comparing_fitting_gaus.png",dpi=300,facecolor='w')
    plt.clf()

    fig, axes2 = plt.subplots(figsize=(16, 10))
    rms_diff_voigt = np.sqrt(np.mean((np.array(cfd20_data) - np.array(cfd20_voigt)) ** 2)).round(4)
    axes2.hist(cfd20_data, bins=200,range=(-1.0,0.0),color='gray',edgecolor='black',
               label=r"$\sigma_{t}(\rm{A_{\rm{Sa}}})$ = " + str(round(data_tr_val, 1)) + r"$\pm$" + str(round(data_tr_val_err, 1)) + " ps")
    axes2.hist(cfd20_voigt, bins=200,range=(-1.0,0.0),color='orange',edgecolor='black',alpha=0.4,
               label=r"$\sigma_{t}(\rm{A_{\rm{Voigt}}})$ = " + str(round(voigt_tr_val, 1)) + r"$\pm$" + str(round(voigt_tr_val_err, 1)) + " ps") # + "\n" + r"$\Delta_{RMS}$ = " + str(rms_diff_voigt))
    axes2.plot(x_tr_fit, data_tr_fit, 'k--', linewidth=2) #, label=r"$\sigma_{\rm{DUT}}$ = " + str(round(data_tr_val, 1)) + r"$\pm$" + str(round(data_tr_val_err, 1)) + " ps")
    axes2.plot(x_tr_fit, voigt_tr_fit, 'orange', linewidth=2) #, label=r"$\sigma_{\rm{DUT}}^{\rm{Voigt}}$ = " + str(round(voigt_tr_val, 1)) + r"$\pm$" + str(round(voigt_tr_val_err, 1)) + " ps")
    axes2.set_xlabel(label_cfd + r" [ns]",fontsize=24)
    axes2.set_ylabel(r"Events",fontsize=24)
    axes2.tick_params(axis="both", labelsize=24)
    axes2.set_xlim(-0.7,-0.5)
    axes2.set_ylim(0,180)
    axes2.legend(fontsize=32)
    axes2.grid(True, axis='both', linestyle='--', alpha=0.5)
    plt.savefig(f"./timeres_comparing_fitting_voigt.png",dpi=300,facecolor='w')
    plt.clf()

    fig, axes3 = plt.subplots(figsize=(16, 10))
    rms_diff_spline = np.sqrt(np.mean((np.array(cfd20_data) - np.array(cfd20_spline)) ** 2)).round(4)
    axes3.hist(cfd20_data, bins=200,range=(-1.0,0.0),color='gray',edgecolor='black',
               label=r"$\sigma_{t}(\rm{A_{\rm{Sa}}})$ = " + str(round(data_tr_val, 1)) + r"$\pm$" + str(round(data_tr_val_err, 1)) + " ps")
    axes3.hist(cfd20_spline, bins=200,range=(-1.0,0.0),color='purple',edgecolor='black',alpha=0.4,
               label=r"$\sigma_{t}(\rm{A_{\rm{spline}}})$ = " + str(round(spline_tr_val, 1)) + r"$\pm$" + str(round(spline_tr_val_err, 1)) + " ps") # + "\n" + r"$\Delta_{\rm{RMS}}$ = " + str(rms_diff_spline))
    axes3.plot(x_tr_fit, data_tr_fit, 'k--', linewidth=2) #, label=r"$\sigma_{\rm{DUT}}$ = " + str(round(data_tr_val, 1)) + r"$\pm$" + str(round(data_tr_val_err, 1)) + " ps")
    axes3.plot(x_tr_fit, spline_tr_fit, 'purple', linewidth=2) #, label=r"$\sigma_{\rm{DUT}}^{\rm{spline}}$ = " + str(round(spline_tr_val, 1)) + r"$\pm$" + str(round(spline_tr_val_err, 1)) + " ps")
    axes3.set_xlabel(label_cfd + r" [ns]",fontsize=24)
    axes3.set_ylabel(r"Events",fontsize=24)
    axes3.tick_params(axis="both", labelsize=24)
    axes3.set_xlim(-0.7,-0.5)
    axes3.set_ylim(0,180)
    axes3.legend(fontsize=32)
    axes3.grid(True, axis='both', linestyle='--', alpha=0.5)
    plt.savefig(f"./timeres_comparing_fitting_spline.png",dpi=300,facecolor='w')
    plt.clf()


if __name__ == "__main__":
  main()
