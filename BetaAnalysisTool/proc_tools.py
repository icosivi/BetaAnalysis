# proc_tools.py

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import minimize, curve_fit
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

'''
from proc_tools import get_fit_results, hist_tree_file_basics, plot_fit_curves, getBias
'''

root.gErrorIgnoreLevel = root.kWarning
import warnings
warnings.filterwarnings("ignore", category=RuntimeWarning)

def round_to_sig_figs(x, sig):
  if x == 0:
    return 0
  else:
    return round(x, sig - int(math.floor(math.log10(abs(x)))) - 1)

def get_fit_results(arr_of_fits, arr_of_biases, arr_of_nevents, channel_of_dut, decomp_sigma=False,simplified=False,add_threshold_frac=None):
  arr_of_ch = []
  arr_of_biases_fitted = []
  arr_of_mean = []
  arr_of_sigma = []
  arr_of_ampl = []
  arr_langaus_sig = []
  arr_of_red_chi2 = []
  arr_of_nevents_fitted = []
  for channel_i, fit_func in enumerate(arr_of_fits):
    if fit_func.GetNpar() > 3: # NEED TO CHANGE THIS TO TYPE OF TIME RES FIT
      langaus_sig = fit_func.GetParameter(3)
    else:
      langaus_sig = 0

    #ndf = fit_func.GetNDF()
    #npar = fit_func.GetNpar()
    #nevents = ndf + npar
    mean = fit_func.GetParameter(1)  # Mean of the gauss distribution
    sigma = fit_func.GetParameter(2) # Sigma of the gauss distribution
    amplitude = fit_func.GetParameter(0)  # Amplitude of the gauss distribution
    chi2 = fit_func.GetChisquare()  # Chi-squared value of the fit
    ndf = fit_func.GetNDF()  # Number of degrees of freedom
    arr_of_nevents_fitted.append(arr_of_nevents[channel_of_dut[channel_i]])
    arr_of_mean.append(round_to_sig_figs(mean,4))
    arr_of_sigma.append(round_to_sig_figs(sigma,4))
    arr_of_ampl.append(round_to_sig_figs(amplitude,3))
    arr_langaus_sig.append(round_to_sig_figs(langaus_sig,4))
    arr_of_ch.append("Ch" + str(channel_of_dut[channel_i]))
    arr_of_biases_fitted.append(arr_of_biases[channel_i])
    if ndf == 0:
      ndf = 1
    arr_of_red_chi2.append(round_to_sig_figs((chi2/ndf),3))

  df_of_results = pd.DataFrame({
    "Channel": arr_of_ch,
    "Bias": arr_of_biases_fitted,
    "Mean": arr_of_mean,
    "Sigma": arr_of_sigma,
    "Amplitude": arr_of_ampl,
    "RChi2": arr_of_red_chi2,
    "NEvents": arr_of_nevents_fitted
  })

  return df_of_results

def getBias(filename, chnum):
  pattern_ch = f"Ch{chnum}"
  ch_match = re.search(pattern_ch, filename)
  if ch_match:
    start_index = ch_match.end()
    substring_to_search = filename[start_index:]
        
    pattern_bias = r"-(\d{2,4}V)"
    bias_match = re.search(pattern_bias, substring_to_search)
    if bias_match:
      return bias_match.group(1)
    else:
      print("[GetBias] : BIAS NOT FOUND")
      return None
  else:
    return None

def landau_tr_quad_fit(Q_df, tr_df, tre_df):
  Q = Q_df.to_numpy()
  tr = tr_df.to_numpy()
  tre = tre_df.to_numpy()
  def model(Q, A, B):
    return A**2 + (B/Q)**2
  tr2 = tr**2
  tre2 = 2 * tr * tre
  popt, pcov = curve_fit(model, Q, tr2, sigma=tre2, absolute_sigma=True)
  A_fit, B_fit = popt
  A_unc, B_unc = np.sqrt(np.diag(pcov))
  return A_fit, A_unc

def hist_tree_file_basics(tree,file,var,index,nBins,xLower,xUpper,biasVal,cut_cond,ch):
  var_dict = {"tmax":"t_{max} / 10 ns" , "pmax":"p_max / mV" , "negpmax":"-p_max / mV", "risetime":"Rise time / ns", "area_new":"Area / pWb" , "rms":"RMS / mV", "dvdt":"dV/dt / mV/ns", "dvdt_2080": "dV/dt[20%:80%] / mV/ns"}
  thisHist = root.TH1F("CH "+str(ch)+" "+biasVal, var+";"+var_dict[var]+";Events", nBins, xLower, xUpper)
  if (var == "pmax") or (var == "negpmax") or (var == "tmax"):
    tree.Draw(var+"["+str(ch)+"]>>CH "+str(ch)+" "+biasVal,"event>-1")
  else:
    tree.Draw(var+"["+str(ch)+"]>>CH "+str(ch)+" "+biasVal,cut_cond)
  thisHist.SetLineWidth(2)
  thisHist.SetLineColor(index+1)
  return thisHist

def plot_fit_curves(xLower,xUpper,fit_type,hist_to_fit,index,biasVal):
  thisFit = TF1(fit_type+"_hist"+biasVal, fit_type, xLower, xUpper)
  hist_to_fit.Fit(thisFit, "Q")
  thisFit.SetLineWidth(3)
  thisFit.SetLineColor(index+1)
  #thisFit.SetLineStyle(2)
  return thisFit
