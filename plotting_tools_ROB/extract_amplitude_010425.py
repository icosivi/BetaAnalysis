import numpy as np
import matplotlib.pyplot as plt
import scipy
import scipy.optimize as opt
import mpmath
import scipy.special as sp
from scipy.optimize import minimize
from scipy.stats import poisson, median_abs_deviation
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
  counts, bin_edges = np.histogram(data_to_bin, bins=100, range=(-1, 0))
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
    print("Landau fit failed—perhaps your data isn't 'Landau-y' enough? Try better initial guesses.")
    return None, None, None, None

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
  num_curves = 1000
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
      if (pmax_sig > 25) and (pmax_mcp < 540) and (pmax_mcp > 40):
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
  labels = ["W17 new (160 V)","W12 15e14 (210 V) G = 4"]

  a_max = []
  a_para = []
  a_gaus = []
  a_lorentz = []
  a_voigt = []
  a_spline = []
  a_landau = []
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

  a_max_mcp = []
  a_para_mcp = []
  a_gaus_mcp = []
  a_lorentz_mcp = []
  a_voigt_mcp = []
  a_spline_mcp = []
  a_landau_mcp = []

  t_Amax_para = []
  t_Amax_gaus = []
  t_Amax_lorentz = []
  t_Amax_voigt = []
  t_Amax_spline = []
  t_Amax_landau = []

  t_w_below_Amax_para = []
  t_w_below_Amax_gaus = []
  t_w_below_Amax_lorentz = []
  t_w_below_Amax_voigt = []
  t_w_below_Amax_spline = []
  t_w_below_Amax_landau = []

  make_plots = False
  make_populated_plots = False
  cfd_studies = True
  add_noise = False
  time_res_calc = True
  make_timewalk_plots = False
  numptseitherside = 3

  if make_plots:
    plt.figure(figsize=(10, 6))

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

      # parabola
      parabola_coeffs = np.polyfit(x_peak, y_peak, 2)
      a, b, _ = parabola_coeffs
      x_parabola_max = -b / (2 * a)
      y_parabola_max = np.polyval(parabola_coeffs, x_parabola_max)
      res_effect_y = np.polyval(parabola_coeffs, x_fine)
      t_Amax_para.append(x_fine[np.argmin(np.abs(res_effect_y - max(res_effect_y)))])

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
      t_w_below_Amax_para.append(x_peak[idx_below])

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

      # Lorentzian
      try:
        lorentz_peak, lorentz_fit, lorentz_func, lorentz_errors = fit_lorentzian(x_peak, y_peak)
        lorentz_peak_mcp, lorentz_fit_mcp, lorentz_func_mcp, lorentz_errors_mcp = fit_lorentzian(x_peak_mcp, y_peak_mcp)
        y_fine = lorentz_func(x_fine)
        lorentz_peak = max(y_fine)
        t_Amax_lorentz.append(x_fine[np.argmin(np.abs(y_fine - lorentz_peak))])
        idx_below = np.where((y_peak < lorentz_peak) & (x_peak < x_fine[np.argmin(np.abs(y_fine - lorentz_peak))]))[0][-1]
        t_w_below_Amax_lorentz.append(x_peak[idx_below])
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

      # landau
      landau_peak, landau_fit, landau_func, landau_errors = fit_landau(x_peak, y_peak)
      landau_peak_mcp, landau_fit_mcp, landau_func_mcp, landau_errors_mcp = fit_landau(x_peak_mcp, y_peak_mcp)
      y_fine = landau_func(x_fine)
      landau_peak = max(y_fine)
      t_Amax_landau.append(x_fine[np.argmin(np.abs(y_fine - landau_peak))])
      idx_below = np.where((y_peak < landau_peak) & (x_peak < x_fine[np.argmin(np.abs(y_fine - landau_peak))]))[0][-1]
      t_w_below_Amax_landau.append(x_peak[idx_below])

      a_max.append(1000*pmax)
      a_para.append(1000*y_parabola_max)
      a_gaus.append(1000*gaussian(mu, *params))
      a_lorentz.append(1000*lorentz_peak)
      a_voigt.append(1000*voigt_peak)
      a_spline.append(1000*spline_peak)
      a_landau.append(1000*landau_peak)

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

      a_max_mcp.append(1000*pmax_mcp)
      a_para_mcp.append(1000*y_parabola_max_mcp)
      a_gaus_mcp.append(1000*gaussian(mu_mcp, *params_mcp))
      a_lorentz_mcp.append(1000*lorentz_peak_mcp)
      a_voigt_mcp.append(1000*voigt_peak_mcp)
      a_spline_mcp.append(1000*spline_peak_mcp)
      a_landau_mcp.append(1000*landau_peak_mcp)

      if make_plots:
        x_linspace = np.linspace(x_peak.min(), x_peak.max(), 400)

        plt.scatter(reshaped_time_data[j],reshaped_ampl_data[j],s=20,c=colours[i],marker='o',edgecolor=edges[i],linewidth=0.5,alpha=alphas[i],label=labels[i]+" "+str(round(1000*pmax, 2))+" mV", zorder=1)
        plt.plot(x_linspace, np.polyval(parabola_coeffs, x_linspace), 'r', label="Parabolic Fit", linewidth = 2, zorder=2)
        plt.plot(x_linspace, gaussian(x_linspace, *params), 'g', label="Gaussian Fit", linewidth = 2, zorder=2)
        plt.plot(x_linspace, lorentz_func(x_linspace), 'blue', label="Lorentz Fit", linewidth = 2, zorder=2)
        plt.plot(x_linspace, voigt_func(x_linspace), 'purple', label="Voigt Fit", linewidth = 2, zorder=2)
        plt.plot(x_linspace, spline_func(x_linspace), 'orange', label="Interpolated Spline", linewidth = 2, zorder=2)
        plt.plot(x_linspace, landau_func(x_linspace), 'cyan', label="Landau Fit", linewidth = 2, zorder=2)
        plt.axhline(y_parabola_max, color='r', linestyle=':', label=r"A$_{para}$: "+str(round(1000*y_parabola_max, 2))+" mV", linewidth = 2, zorder=3)
        plt.axhline(gaussian(mu, *params), color='g', linestyle=':', label=r"A$_{Gaus}$: "+str(round(1000*gaussian(mu, *params), 2))+" mV", linewidth = 2, zorder=3)
        plt.axhline(lorentz_peak, color='blue', linestyle=':', label=r"A$_{Lorentz}$: "+str(round(1000*lorentz_peak, 2))+" mV", linewidth = 2, zorder=3)
        plt.axhline(voigt_peak, color='purple', linestyle=':', label=r"A$_{Voigt}$: "+str(round(1000*voigt_peak, 2))+" mV", linewidth = 2, zorder=3)
        plt.axhline(spline_peak, color='orange', linestyle=':', label=r"A$_{spline}$: "+str(round(1000*spline_peak, 2))+" mV", linewidth = 2, zorder=3)
        plt.axhline(landau_peak, color='cyan', linestyle=':', label=r"A$_{Landau}$: "+str(round(1000*landau_peak, 2))+" mV", linewidth = 2, zorder=3)
        plt.scatter(reshaped_time_data[j],reshaped_ampl_data[j],s=30,c=colours[i],marker='o',edgecolor=edges[i],linewidth=0.5,alpha=alphas[i], zorder=4)
        #plt.plot(reshaped_time_data[0],reshaped_ampl_data[0],linestyle='-',color=colours[i],linewidth=1.5,alpha=alphas[i],label=labels[i])
        #for j in range(len(reshaped_time_data)-1):
        #  plt.plot(reshaped_time_data[j+1],reshaped_ampl_data[j+1],linestyle='-',color=colours[i],linewidth=1.5,alpha=alphas[i])
        if num_curves == 1:
          print(f"A_max = {(1000*pmax):.2f} mV")
          print(f"A_para = {(1000*y_parabola_max):.2f} mV")
          print(f"A_Gaus = {(1000*gaussian(mu, *params)):.2f} mV")
          print(f"A_Lorentz = {(1000*lorentz_peak):.2f} mV")
          print(f"A_Voigt = {(1000*voigt_peak):.2f} mV")
          print(f"A_spline = {(1000*spline_peak):.2f} mV")
          print(f"A_Landau = {(1000*landau_peak):.2f} mV")

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

    if make_plots:
      plt.xlabel('Time [ns]',fontsize=14)
      plt.ylabel('Amplitude [mV]',fontsize=14)
      plt.xticks(fontsize=14)
      plt.yticks(fontsize=14)
      plt.xlim(-2,1)
      plt.legend(fontsize=12)
      #plt.yscale('log')
      plt.grid(True, linestyle='--', alpha=0.5)
      plt.tight_layout()
      plt.savefig("./amplitude_analysis.png",dpi=300,facecolor='w')
      #plt.show()
      plt.clf()

  if make_populated_plots:
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

  if cfd_studies:

    cfd20_para = []
    cfd20_gaus = []
    cfd20_lorentz = []
    cfd20_voigt = []
    cfd20_spline = []
    cfd20_landau = []

    cfd20_para_mcp = []
    cfd20_gaus_mcp = []
    cfd20_lorentz_mcp = []
    cfd20_voigt_mcp = []
    cfd20_spline_mcp = []
    cfd20_landau_mcp = []

    control_data = []
    control_mcp = []

    for j in range(num_curves):
      time_array_event = np.array(reshaped_time_data[j])
      ampl_array_event = np.array(reshaped_ampl_data[j])

      time_value_control = find_CFD_time_with_threshold(time_array_event, ampl_array_event, ampl_array_event.max(), 0.2)
      time_value_para = find_CFD_time_with_threshold(time_array_event, ampl_array_event, a_para[j], 0.2/1000)
      time_value_gaus = find_CFD_time_with_threshold(time_array_event, ampl_array_event, a_gaus[j], 0.2/1000)
      time_value_lorentz = find_CFD_time_with_threshold(time_array_event, ampl_array_event, a_lorentz[j], 0.2/1000)
      time_value_voigt = find_CFD_time_with_threshold(time_array_event, ampl_array_event, a_voigt[j], 0.2/1000)
      time_value_spline = find_CFD_time_with_threshold(time_array_event, ampl_array_event, a_spline[j], 0.2/1000)
      time_value_landau = find_CFD_time_with_threshold(time_array_event, ampl_array_event, a_landau[j], 0.2/1000)
      control_data.append(time_value_control)
      cfd20_para.append(time_value_para)
      cfd20_gaus.append(time_value_gaus)
      cfd20_lorentz.append(time_value_lorentz)
      cfd20_voigt.append(time_value_voigt)
      cfd20_spline.append(time_value_spline)
      cfd20_landau.append(time_value_landau)

      if time_res_calc:
        time_array_event_mcp = np.array(rtd_mcp[j])
        ampl_array_event_mcp = np.array(rad_mcp[j])

        time_value_control_mcp = find_CFD_time_with_threshold(time_array_event_mcp, ampl_array_event_mcp, ampl_array_event_mcp.max(), 0.2)
        time_value_para_mcp = find_CFD_time_with_threshold(time_array_event_mcp, ampl_array_event_mcp, a_para_mcp[j], 0.2/1000)
        time_value_gaus_mcp = find_CFD_time_with_threshold(time_array_event_mcp, ampl_array_event_mcp, a_gaus_mcp[j], 0.2/1000)
        time_value_lorentz_mcp = find_CFD_time_with_threshold(time_array_event_mcp, ampl_array_event_mcp, a_lorentz_mcp[j], 0.2/1000)
        time_value_voigt_mcp = find_CFD_time_with_threshold(time_array_event_mcp, ampl_array_event_mcp, a_voigt_mcp[j], 0.2/1000)
        time_value_spline_mcp = find_CFD_time_with_threshold(time_array_event_mcp, ampl_array_event_mcp, a_spline_mcp[j], 0.2/1000)
        time_value_landau_mcp = find_CFD_time_with_threshold(time_array_event_mcp, ampl_array_event_mcp, a_landau_mcp[j], 0.2/1000)
        control_mcp.append(time_value_control_mcp)
        cfd20_para_mcp.append(time_value_para_mcp)
        cfd20_gaus_mcp.append(time_value_gaus_mcp)
        cfd20_lorentz_mcp.append(time_value_lorentz_mcp)
        cfd20_voigt_mcp.append(time_value_voigt_mcp)
        cfd20_spline_mcp.append(time_value_spline_mcp)
        cfd20_landau_mcp.append(time_value_landau_mcp)

    label_cfd = "CFD@20%"

    if time_res_calc:
      
      #cfd20_data = np.array(cfd20_data) - np.array(cfd20_mcp_data)
      cfd20_data = np.array(control_data) - np.array(control_mcp) 
      cfd20_para = np.array(cfd20_para) - np.array(cfd20_para_mcp)
      cfd20_gaus = np.array(cfd20_gaus) - np.array(cfd20_gaus_mcp)
      cfd20_lorentz = np.array(cfd20_lorentz) - np.array(cfd20_lorentz_mcp)
      cfd20_voigt = np.array(cfd20_voigt) - np.array(cfd20_voigt_mcp)
      cfd20_spline = np.array(cfd20_spline) - np.array(cfd20_spline_mcp)
      cfd20_landau = np.array(cfd20_landau) - np.array(cfd20_landau_mcp)

      label_cfd = r"$\sigma_{t}^{20\%}$"

      data_tr_params = gaussian_fit_binned_data(cfd20_data, "Data")
      para_tr_params = gaussian_fit_binned_data(cfd20_para, "Parabolic")
      gaus_tr_params = gaussian_fit_binned_data(cfd20_gaus, "Gaussian")
      lorentz_tr_params = gaussian_fit_binned_data(cfd20_lorentz, "Lorentz")
      voigt_tr_params = gaussian_fit_binned_data(cfd20_voigt, "Voigt")
      spline_tr_params = gaussian_fit_binned_data(cfd20_spline, "Interpolated spline")
      landau_tr_params = gaussian_fit_binned_data(cfd20_landau, "Landau")

    fig, axes = plt.subplots(2, 3, figsize=(18, 12))
    rms_diff_para = np.sqrt(np.mean((np.array(cfd20_data) - np.array(cfd20_para)) ** 2)).round(3)
    axes[0,0].hist(cfd20_data, bins=100,range=(-1, 0),color='gray',edgecolor='black',label=label_cfd + r"(A$_{max}$)")
    axes[0,0].hist(cfd20_para, bins=100,range=(-1, 0),color='r',edgecolor='black',alpha=0.4,label=label_cfd + r"(A$_{para}$)" + "\n" + r"$\Delta_{RMS}$ = " + str(rms_diff_para))

    rms_diff_gaus = np.sqrt(np.mean((np.array(cfd20_data) - np.array(cfd20_gaus)) ** 2)).round(3)
    axes[0,1].hist(cfd20_data, bins=100,range=(-1, 0),color='gray',edgecolor='black',label=label_cfd + r"(A$_{max}$)")
    axes[0,1].hist(cfd20_gaus, bins=100,range=(-1, 0),color='g',edgecolor='black',alpha=0.4,label=label_cfd + r"(A$_{Gaus}$)" + "\n" + r"$\Delta_{RMS}$ = " + str(rms_diff_gaus))

    rms_diff_lorentz = np.sqrt(np.mean((np.array(cfd20_data) - np.array(cfd20_lorentz)) ** 2)).round(3)
    axes[1,0].hist(cfd20_data, bins=100,range=(-1, 0),color='gray',edgecolor='black',label=label_cfd + r"(A$_{max}$)")
    axes[1,0].hist(cfd20_lorentz, bins=100,range=(-1, 0),color='blue',edgecolor='black',alpha=0.4,label=label_cfd + r"(A$_{Lorentz}$)" + "\n" + r"$\Delta_{RMS}$ = " + str(rms_diff_lorentz))

    rms_diff_voigt = np.sqrt(np.mean((np.array(cfd20_data) - np.array(cfd20_voigt)) ** 2)).round(3)
    axes[1,1].hist(cfd20_data, bins=100,range=(-1, 0),color='gray',edgecolor='black',label=label_cfd + r"(A$_{max}$)")
    axes[1,1].hist(cfd20_voigt, bins=100,range=(-1, 0),color='orange',edgecolor='black',alpha=0.4,label=label_cfd + r"(A$_{Voigt}$)" + "\n" + r"$\Delta_{RMS}$ = " + str(rms_diff_voigt))

    rms_diff_spline = np.sqrt(np.mean((np.array(cfd20_data) - np.array(cfd20_spline)) ** 2)).round(3)
    axes[0,2].hist(cfd20_data, bins=100,range=(-1, 0),color='gray',edgecolor='black',label=label_cfd + r"(A$_{max}$)")
    axes[0,2].hist(cfd20_spline, bins=100,range=(-1, 0),color='purple',edgecolor='black',alpha=0.4,label=label_cfd + r"(A$_{spline}$)" + "\n" + r"$\Delta_{RMS}$ = " + str(rms_diff_spline))

    rms_diff_landau = np.sqrt(np.mean((np.array(cfd20_data) - np.array(cfd20_landau)) ** 2)).round(3)
    axes[1,2].hist(cfd20_data, bins=100,range=(-1, 0),color='gray',edgecolor='black',label=label_cfd + r"(A$_{max}$)")
    axes[1,2].hist(cfd20_landau, bins=100,range=(-1, 0),color='brown',edgecolor='black',alpha=0.4,label=label_cfd + r"(A$_{Landau}$)" + "\n" + r"$\Delta_{RMS}$ = " + str(rms_diff_landau))

    if time_res_calc:
      x_tr_fit = np.linspace(-1, 0, 1000)
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
          axes[i,j].plot(x_tr_fit, data_tr_fit, 'k--', linewidth=2, label=r"$\sigma_{tr}$ = " + str(round(data_tr_val, 1)) + " ps")

      axes[0,0].plot(x_tr_fit, para_tr_fit, 'r', linewidth=2, label=r"$\sigma_{tr}^{para}$ = " + str(round(para_tr_val, 1)) + " ps")
      axes[0,1].plot(x_tr_fit, gaus_tr_fit, 'g', linewidth=2, label=r"$\sigma_{tr}^{Gaus}$ = " + str(round(gaus_tr_val, 1)) + " ps")
      axes[1,0].plot(x_tr_fit, lorentz_tr_fit, 'blue', linewidth=2, label=r"$\sigma_{tr}^{Lorentz}$ = " + str(round(lorentz_tr_val, 1)) + " ps")
      axes[1,1].plot(x_tr_fit, voigt_tr_fit, 'orange', linewidth=2, label=r"$\sigma_{tr}^{Voigt}$ = " + str(round(voigt_tr_val, 1)) + " ps")
      axes[0,2].plot(x_tr_fit, spline_tr_fit, 'purple', linewidth=2, label=r"$\sigma_{tr}^{spline}$ = " + str(round(spline_tr_val, 1)) + " ps")
      axes[1,2].plot(x_tr_fit, landau_tr_fit, 'brown', linewidth=2, label=r"$\sigma_{tr}^{Landau}$ = " + str(round(landau_tr_val, 1)) + " ps")

    for i in range(2):
      for j in range(3):
        axes[i,j].set_xlabel(label_cfd + r"/ ns",fontsize=14)
        axes[i,j].set_ylabel(r"Events",fontsize=14)
        axes[i,j].set_xlim(-0.8,-0.4)
        axes[i,j].legend(fontsize=14)
        axes[i,j].grid(True, axis='both', linestyle='--', alpha=0.5)

    fig.suptitle(f"Total {len(a_max)} signal events", fontsize=16, fontweight='bold')
    plt.tight_layout(rect=[0, 0, 1, 0.96])
    if time_res_calc:
      print("NEW DATA")
      plt.savefig("./timeres_20_260325.png",dpi=300,facecolor='w')
    else:
      plt.savefig("./cfd_20_260325.png",dpi=300,facecolor='w')
    plt.clf()

  if make_timewalk_plots:
    deltaT_para = np.array(t_Amax_para) - np.array(t_w_below_Amax_para)
    deltaT_gaus = np.array(t_Amax_gaus) - np.array(t_w_below_Amax_gaus)
    deltaT_lorentz = np.array(t_Amax_lorentz) - np.array(t_w_below_Amax_lorentz)
    deltaT_voigt = np.array(t_Amax_voigt) - np.array(t_w_below_Amax_voigt)
    deltaT_spline = np.array(t_Amax_spline) - np.array(t_w_below_Amax_spline)
    deltaT_landau = np.array(t_Amax_landau) - np.array(t_w_below_Amax_landau)

    deltaT_para = np.where(deltaT_para > 0.1, deltaT_para - 0.1, deltaT_para)
    deltaT_gaus = np.where(deltaT_gaus > 0.1, deltaT_gaus - 0.1, deltaT_gaus)
    deltaT_lorentz = np.where(deltaT_lorentz > 0.1, deltaT_lorentz - 0.1, deltaT_lorentz)
    deltaT_voigt = np.where(deltaT_voigt > 0.1, deltaT_voigt - 0.1, deltaT_voigt)
    deltaT_spline = np.where(deltaT_spline > 0.1, deltaT_spline - 0.1, deltaT_spline)
    deltaT_landau = np.where(deltaT_landau > 0.1, deltaT_landau - 0.1, deltaT_landau)

    label_timewalk = r"$\Delta$t(A$_{fit}$,w$<$A$_{fit}$)"
    num_bins_timewalk = 10

    fig, axes = plt.subplots(2, 3, figsize=(18, 12))
    counts_para, bins_para, _ = axes[0,0].hist(deltaT_para, bins=num_bins_timewalk,range=(0.0, 0.1),color='r',edgecolor='black',label=r"$\Delta$t$_{para}$")
    counts_gaus, bins_gaus, _ = axes[0,1].hist(deltaT_gaus, bins=num_bins_timewalk,range=(0.0, 0.1),color='g',edgecolor='black',label=r"$\Delta$t$_{Gaus}$")
    counts_lorentz, bins_lorentz, _ = axes[1,0].hist(deltaT_lorentz, bins=num_bins_timewalk,range=(0.0, 0.1),color='blue',edgecolor='black',label=r"$\Delta$t$_{Lorentz}$")
    counts_voigt, bins_voigt, _ = axes[1,1].hist(deltaT_voigt, bins=num_bins_timewalk,range=(0.0, 0.1),color='orange',edgecolor='black',label=r"$\Delta$t$_{Voigt}$")
    counts_spline, bins_spline, _ = axes[0,2].hist(deltaT_spline, bins=num_bins_timewalk,range=(0.0, 0.1),color='purple',edgecolor='black',label=r"$\Delta$t$_{spline}$")
    counts_landau, bins_landau, _ = axes[1,2].hist(deltaT_landau, bins=num_bins_timewalk,range=(0.0, 0.1),color='brown',edgecolor='black',label=r"$\Delta$t$_{Landau}$")
    
    if time_res_calc:
      axes_RHS = np.empty((2, 3), dtype=object)

    for i in range(2):
      for j in range(3):
        axes[i,j].set_xlabel(label_timewalk + r" / ns",fontsize=14)
        axes[i,j].set_ylabel(r"Events",fontsize=14)
        axes[i,j].set_xlim(0.0,0.1)
        axes[i,j].legend(loc='upper left',fontsize=14)
        axes[i,j].grid(True, axis='both', linestyle='--', alpha=0.5)
        if time_res_calc:
          axes_RHS[i,j] = axes[i,j].twinx()
          axes_RHS[i,j].set_ylabel(r"$\sigma_{tr}$ of events in given bin / ps",fontsize=14)

    if time_res_calc:

      bin_indices_para = np.digitize(deltaT_para, bins_para, right=False) - 1
      bin_indices_gaus = np.digitize(deltaT_gaus, bins_gaus, right=False) - 1
      bin_indices_lorentz = np.digitize(deltaT_lorentz, bins_lorentz, right=False) - 1
      bin_indices_voigt = np.digitize(deltaT_voigt, bins_voigt, right=False) - 1
      bin_indices_spline = np.digitize(deltaT_spline, bins_spline, right=False) - 1
      bin_indices_landau = np.digitize(deltaT_landau, bins_landau, right=False) - 1

      bin_indices_para[bin_indices_para == num_bins_timewalk] = (num_bins_timewalk-1)
      bin_indices_gaus[bin_indices_gaus == num_bins_timewalk] = (num_bins_timewalk-1)
      bin_indices_lorentz[bin_indices_lorentz == num_bins_timewalk] = (num_bins_timewalk-1)
      bin_indices_voigt[bin_indices_voigt == num_bins_timewalk] = (num_bins_timewalk-1)
      bin_indices_spline[bin_indices_spline == num_bins_timewalk] = (num_bins_timewalk-1)
      bin_indices_landau[bin_indices_landau == num_bins_timewalk] = (num_bins_timewalk-1)

      binned_para = [np.where(bin_indices_para == i)[0] for i in range(num_bins_timewalk)]
      binned_gaus = [np.where(bin_indices_gaus == i)[0] for i in range(num_bins_timewalk)]
      binned_lorentz = [np.where(bin_indices_lorentz == i)[0] for i in range(num_bins_timewalk)]
      binned_voigt = [np.where(bin_indices_voigt == i)[0] for i in range(num_bins_timewalk)]
      binned_spline = [np.where(bin_indices_spline == i)[0] for i in range(num_bins_timewalk)]
      binned_landau = [np.where(bin_indices_landau == i)[0] for i in range(num_bins_timewalk)]

      selected_tr_events_para = [cfd20_para[binned_para[i]] for i in range(num_bins_timewalk)]
      selected_tr_events_gaus = [cfd20_gaus[binned_gaus[i]] for i in range(num_bins_timewalk)]
      selected_tr_events_lorentz = [cfd20_lorentz[binned_lorentz[i]] for i in range(num_bins_timewalk)]
      selected_tr_events_voigt = [cfd20_voigt[binned_voigt[i]] for i in range(num_bins_timewalk)]
      selected_tr_events_spline = [cfd20_spline[binned_spline[i]] for i in range(num_bins_timewalk)]
      selected_tr_events_landau = [cfd20_landau[binned_landau[i]] for i in range(num_bins_timewalk)]
      
      selected_para_tr_params = [gaussian_fit_binned_data(selected_tr_events_para[i], "Parabolic") for i in range(num_bins_timewalk)]
      selected_gaus_tr_params = [gaussian_fit_binned_data(selected_tr_events_gaus[i], "Gaussian") for i in range(num_bins_timewalk)]
      selected_lorentz_tr_params = [gaussian_fit_binned_data(selected_tr_events_lorentz[i], "Lorentz") for i in range(num_bins_timewalk)]
      selected_voigt_tr_params = [gaussian_fit_binned_data(selected_tr_events_voigt[i], "Voigt") for i in range(num_bins_timewalk)]
      selected_spline_tr_params = [gaussian_fit_binned_data(selected_tr_events_spline[i], "Interpolated spline") for i in range(num_bins_timewalk)]
      selected_landau_tr_params = [gaussian_fit_binned_data(selected_tr_events_landau[i], "Landau") for i in range(num_bins_timewalk)]

      bin_centres_para = (bins_para[:-1] + bins_para[1:]) / 2
      bin_centres_gaus = (bins_gaus[:-1] + bins_gaus[1:]) / 2
      bin_centres_lorentz = (bins_lorentz[:-1] + bins_lorentz[1:]) / 2
      bin_centres_voigt = (bins_voigt[:-1] + bins_voigt[1:]) / 2
      bin_centres_spline = (bins_spline[:-1] + bins_spline[1:]) / 2
      bin_centres_landau = (bins_landau[:-1] + bins_landau[1:]) / 2

      mcp_tr_est = 5
      per_bin_tr_val_para = np.sqrt(np.where((1000*np.array(selected_para_tr_params)[:,2])**2 - mcp_tr_est**2 < 0, 0, (1000*np.array(selected_para_tr_params)[:,2])**2 - mcp_tr_est**2))
      per_bin_tr_val_para = np.sqrt(np.where((1000*np.array(selected_para_tr_params)[:,2])**2 - mcp_tr_est**2 < 0, 0, (1000*np.array(selected_para_tr_params)[:,2])**2 - mcp_tr_est**2))
      per_bin_tr_val_gaus = np.sqrt(np.where((1000*np.array(selected_gaus_tr_params)[:,2])**2 - mcp_tr_est**2 < 0, 0, (1000*np.array(selected_gaus_tr_params)[:,2])**2 - mcp_tr_est**2))
      per_bin_tr_val_lorentz = np.sqrt(np.where((1000*np.array(selected_lorentz_tr_params)[:,2])**2 - mcp_tr_est**2 < 0, 0, (1000*np.array(selected_lorentz_tr_params)[:,2])**2 - mcp_tr_est**2))
      per_bin_tr_val_voigt = np.sqrt(np.where((1000*np.array(selected_voigt_tr_params)[:,2])**2 - mcp_tr_est**2 < 0, 0, (1000*np.array(selected_voigt_tr_params)[:,2])**2 - mcp_tr_est**2))
      per_bin_tr_val_spline = np.sqrt(np.where((1000*np.array(selected_spline_tr_params)[:,2])**2 - mcp_tr_est**2 < 0, 0, (1000*np.array(selected_spline_tr_params)[:,2])**2 - mcp_tr_est**2))
      per_bin_tr_val_landau = np.sqrt(np.where((1000*np.array(selected_landau_tr_params)[:,2])**2 - mcp_tr_est**2 < 0, 0, (1000*np.array(selected_landau_tr_params)[:,2])**2 - mcp_tr_est**2))

      axes_RHS[0,0].plot(bin_centres_para, per_bin_tr_val_para, 'o-', color='k', label='$\sigma_{tr}^{para}$', markersize=10, linewidth=2, alpha=0.5)
      axes_RHS[0,1].plot(bin_centres_gaus, per_bin_tr_val_gaus, 'o-', color='k', label='$\sigma_{tr}^{Gaus}$', markersize=10, linewidth=2, alpha=0.5)
      axes_RHS[1,0].plot(bin_centres_lorentz, per_bin_tr_val_lorentz, 'o-', color='k', label='$\sigma_{tr}^{Lorentz}$', markersize=10, linewidth=2, alpha=0.5)
      axes_RHS[1,1].plot(bin_centres_voigt, per_bin_tr_val_voigt, 'o-', color='k', label='$\sigma_{tr}^{Voigt}$', markersize=10, linewidth=2, alpha=0.5)
      axes_RHS[0,2].plot(bin_centres_spline, per_bin_tr_val_spline, 'o-', color='k', label='$\sigma_{tr}^{spline}$', markersize=10, linewidth=2, alpha=0.5)
      axes_RHS[1,2].plot(bin_centres_landau, per_bin_tr_val_landau, 'o-', color='k', label='$\sigma_{tr}^{Landau}$', markersize=10, linewidth=2, alpha=0.5)
      for i in range(2):
        for j in range(3):
          axes_RHS[i,j].legend(loc='upper right',fontsize=14)

    fig.suptitle(f"Total {len(a_max)} signal events", fontsize=16, fontweight='bold')
    plt.tight_layout(rect=[0, 0, 1, 0.96])
    plt.savefig("./timewalk_010425.png",dpi=300,facecolor='w')
    plt.clf()




if __name__ == "__main__":
  main()
