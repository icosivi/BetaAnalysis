import numpy as np
import matplotlib.pyplot as plt
import scipy
from scipy.optimize import minimize
from scipy.stats import poisson, median_abs_deviation
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
from array import array
from scipy.optimize import curve_fit

from firstPass_ProcTools import getBias
from secondPass_ProcTools import binned_fit_gaussian, gaussian

def round_to_sig_figs(x, sig):
  if x == 0:
    return 0
  else:
    return round(x, sig - int(math.floor(math.log10(abs(x)))) - 1)

def plot_gaussian(var, file, file_index, tree, channel_array, nBins, SNDisc_signal_events_bool, savename):

  arr_of_ch = []
  arr_of_biases = []
  arr_of_mean = []
  arr_of_sigma = []
  arr_of_sse = []
  arr_of_rchi2 = []

  SNDisc_Index = 0
  var_dict = {"risetime":"Rise time / ns", "rms":"RMS / mV"}
  for ch_ind, ch_val in enumerate(channel_array):
    data_filtered_signal = []
    sensorType, AtQfactor, mcp_selection = ch_val

    if sensorType == 1:
      bias_of_channel = getBias(str(file), ch_ind)
      for i, entry in enumerate(tree):
        if SNDisc_signal_events_bool[SNDisc_Index][i]:
          if (var == "risetime") & (SNDisc_signal_events_bool[SNDisc_Index][i]):
            sig_event = entry.risetime[ch_ind]
            data_filtered_signal.append(sig_event)
          if (var == "rms") & SNDisc_signal_events_bool[SNDisc_Index][i]:
            sig_event = entry.rms[ch_ind]
            data_filtered_signal.append(sig_event)
      SNDisc_Index += 1
    else:
      continue

    plt.figure(figsize=(10, 6))
    data_var = np.array(data_filtered_signal)

    popt, sse, rchi2 = binned_fit_gaussian(data_var, nBins)

    arr_of_ch.append("Ch"+str(ch_ind))
    arr_of_biases.append(bias_of_channel)
    arr_of_mean.append(popt[1])
    arr_of_sigma.append(popt[2])
    arr_of_sse.append(sse)
    arr_of_rchi2.append(rchi2)

    fig = go.Figure()
    fig.update_layout(
      xaxis_title=var_dict[var],
      yaxis_title='Probability Density',
      title='Gaussian Fit',
    )

    fig.add_trace(
      go.Histogram(
        x=data_var,
        xbins=dict(start=0, end=2*np.mean(data_var), size=2*np.mean(data_var)/nBins),
        name='Histogram',
        histnorm='probability density',
        #nbinsx=nBins,
        opacity=0.6,
        marker=dict(color='green', line=dict(color='black', width=1)),
      )
    )

    x_axis = np.linspace(0, 2*np.mean(data_var), 999)
    fig.add_trace(
      go.Scatter(
        x=x_axis,
        y=gaussian(x_axis, *popt[:-1]),
        name=f'Gaussian Fit<br>μ={popt[1]:.2e}<br>σ={popt[2]:.2e}',
        mode='lines',
        line=dict(width=3,dash='dash',color='black'),
      )
    )

    if not os.path.exists(var):
      os.makedirs(var)
    fig.write_image(var+"/"+savename+"_Ch"+str(ch_ind)+".png")
    print("[BETA ANALYSIS]: [GAUSSIAN PLOTTER] Saved file "+var+"/"+savename+"_Ch"+str(ch_ind)+".png")

  df_of_results = pd.DataFrame({
    "Channel": arr_of_ch,
    "Bias": arr_of_biases,
    "Mean": arr_of_mean,
    "Sigma": arr_of_sigma,
    "SSE score": arr_of_sse,
    "Red. Chi2": arr_of_rchi2,
  })
  return df_of_results
