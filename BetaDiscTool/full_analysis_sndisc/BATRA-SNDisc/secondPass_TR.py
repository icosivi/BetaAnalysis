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
from secondPass_ProcTools import binned_fit_gaussian, gaussian, get_fit_results_TR

def round_to_sig_figs(x, sig):
  if x == 0:
    return 0
  else:
    return round(x, sig - int(math.floor(math.log10(abs(x)))) - 1)

def plot_TR(var, file, file_index, tree, channel_array, nBins, SNDisc_signal_events_bool, savename, mcp_specs):

  nBins = 100

  arr_of_ch = []
  arr_of_biases = []
  arr_of_mean = []
  arr_of_sigma = []
  arr_of_sigma_uncs = []
  arr_of_sse = []
  arr_of_rchi2 = []

  SNDisc_Index = 0
  data_filtered_signal_dut_set = []
  data_filtered_signal_mcp = []
  pmax_unselected_mcp = []
  ch_literal_index = []
  mcp_exists = False
  var_dict = {"timeres":"tn-tn+1 / ns"}
  for ch_ind, ch_val in enumerate(channel_array):
    data_filtered_signal_dut = []
    sensorType, AtQfactor, mcp_selection = ch_val
    if sensorType == 1:
      bias_of_channel = getBias(str(file), ch_ind)
      arr_of_biases.append(bias_of_channel)
      ch_literal_index.append(ch_ind)
      for i, entry in enumerate(tree):
        if SNDisc_signal_events_bool[SNDisc_Index][i]:
          sig_event = entry.cfd[ch_ind][2]
          data_filtered_signal_dut.append(sig_event)
      SNDisc_Index += 1
      data_filtered_signal_dut_set.append(data_filtered_signal_dut)
    elif sensorType == 2:
      mcp_exists = True
      mcp_pmax_selection = ch_val[2]
      for i, entry in enumerate(tree):
        pmax_unselected_mcp_event = entry.pmax[ch_ind]
        sig_event = entry.cfd[ch_ind][2]
        data_filtered_signal_mcp.append(sig_event)
        pmax_unselected_mcp.append(pmax_unselected_mcp_event)
    else:
      continue

  if mcp_exists:
    # need to do index selection for each DUT which will depend on data_filtered_signal
    cfd_mcp_indexed_per_dut = []
    pmax_mcp_indexed_per_dut = []
    for dut_val in SNDisc_signal_events_bool:
      filtered_mcp_cfd_for_this_dut = [val for val, keep in zip(data_filtered_signal_mcp, dut_val) if keep]
      filtered_mcp_pmax_for_this_dut = [val for val, keep in zip(pmax_unselected_mcp, dut_val) if keep]
      cfd_mcp_indexed_per_dut.append(filtered_mcp_cfd_for_this_dut)
      pmax_mcp_indexed_per_dut.append(filtered_mcp_pmax_for_this_dut)

    # selection according to PMAX of MCP, applied to MCP and DUT together
    selected_indices_per_dut = []
    for dut_values in pmax_mcp_indexed_per_dut:
      indices = [i for i, val in enumerate(dut_values) if mcp_pmax_selection[0] <= val <= mcp_pmax_selection[1]]
      selected_indices_per_dut.append(indices)

    def filter_list_by_indices(data_list, selected_indices_per_dut):
      filtered = []
      for dut_values, indices in zip(data_list, selected_indices_per_dut):
        filtered.append([dut_values[i] for i in indices])
      return filtered

    final_filter_dut = filter_list_by_indices(data_filtered_signal_dut_set, selected_indices_per_dut)
    final_filter_mcp = filter_list_by_indices(cfd_mcp_indexed_per_dut, selected_indices_per_dut)

    data_dut_low = []
    data_mcp_low = []
    for fi, fv in enumerate(final_filter_dut):
      fv = np.array(fv)
      fv_new = fv[fv < 0]
      filtered_ev_mcp = np.array(final_filter_mcp[fi])
      filtered_ev_mcp_new = filtered_ev_mcp[fv < 0]
      data_dut_low.append(fv_new)
      data_mcp_low.append(filtered_ev_mcp_new)

  final_filter_dut = data_dut_low
  final_filter_mcp = data_mcp_low

  for data_dut_ind, data_dut_val in enumerate(final_filter_dut):

    data_var = np.array(data_dut_val) - np.array(final_filter_mcp[data_dut_ind])

    plt.figure(figsize=(10, 6))

    popt, sse, rchi2 = binned_fit_gaussian(data_var, nBins)

    arr_of_ch.append("Ch"+str(ch_literal_index[data_dut_ind]))
    arr_of_mean.append(popt[1])
    arr_of_sigma.append(popt[2])
    arr_of_sigma_uncs.append(popt[3])
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
        xbins=dict(start=max(-3,min(data_var)), end=max(data_var), size=(max(data_var)-min(data_var))/nBins),
        name='Histogram',
        histnorm='probability density',
        #nbinsx=nBins,
        opacity=0.6,
        marker=dict(color='green', line=dict(color='black', width=1)),
      )
    )

    x_axis = np.linspace(max(-3,min(data_var)), max(data_var), 999)
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
    fig.write_image(var+"/"+savename+"_Ch"+str(ch_literal_index[data_dut_ind])+".png")
    print("[BETA ANALYSIS]: [GAUSSIAN PLOTTER] Saved file "+var+"/"+savename+"_Ch"+str(ch_literal_index[data_dut_ind])+".png")

  df_of_results = pd.DataFrame({
    "Channel": arr_of_ch,
    "Bias": arr_of_biases,
    "Mean": arr_of_mean,
    "Sigma": arr_of_sigma,
    "Uncertainty": arr_of_sigma_uncs,
    "SSE score": arr_of_sse,
    "Red. Chi2": arr_of_rchi2,
  })

  df_of_results =  get_fit_results_TR(df_of_results, mcp_specs)

  return df_of_results
