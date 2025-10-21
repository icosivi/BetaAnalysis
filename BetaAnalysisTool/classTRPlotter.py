# classPlotter.py

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import minimize
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

from proc_tools_TR import get_fit_results_TR, hist_tree_file_timeres, plot_fit_curves
from proc_tools import getBias

class plotTRVar:
  def __init__(self, var, nBins, xLower, xUpper, log_scale, save_name):

    self.var = var
    self.nBins = nBins
    self.xLower = xLower
    self.xUpper = xUpper
    self.log_scale = log_scale
    self.save_name = save_name

  def run(self, file, file_index, tree, channel_array, mcp_tr):

    arr_of_hists = []
    arr_of_biases = []
    arr_of_nevents = []

    root.gErrorIgnoreLevel = root.kWarning

    result = []
    for i, (_, _, (A, B, C, D, E)) in enumerate(channel_array):
      if (A == 0.0) or (A == []):
        condition = f"pmax[{i}] > 0.0 && pmax[{i}] < {B} && negpmax[{i}] > {C} && tmax[{i}] > {D} && tmax[{i}] < {E}"
      else:
        condition = f"pmax[{i}] > {A[file_index]} && pmax[{i}] < {B} && negpmax[{i}] > {C} && tmax[{i}] > {D} && tmax[{i}] < {E}"
      result.append(condition)

    duts_to_analyse = []
    channel_of_dut = []
    mcp_exists = False

    for j in range(len(channel_array)):
      if (channel_array[j][0] == 1):
        bias = getBias(str(file), j)
        arr_of_biases.append(bias)
        duts_to_analyse.append([["cfd["+str(j)+"][2]-cfd[","cfd["+str(j)+"][0]-cfd[","cfd["+str(j)+"][4]-cfd["], result[j], j])
        channel_of_dut.append(j)
      elif (channel_array[j][0] == 2):
        mcp_channel = [[str(j)+"][2]",str(j)+"][0]",str(j)+"][4]"], result[j]]
        mcp_exists = True

    duts_vars_cuts = []
    if (mcp_exists == True) & (mcp_tr != (0, 0)):
      for dut in duts_to_analyse:
        vars_cuts_defined = ([x + y for x, y in zip(dut[0], mcp_channel[0])], dut[1] + " && " + mcp_channel[1], dut[2])
        duts_vars_cuts.append(vars_cuts_defined)
    if (mcp_exists == False) | (mcp_tr == (0, 0)):
      for dut in duts_to_analyse:
        vars_cuts_defined = ([x + y for x, y in zip(dut[0], mcp_channel[0])], dut[1] + " && " + mcp_channel[1], dut[2])
        duts_vars_cuts.append(vars_cuts_defined)
      adapted_dut_channel = [second_plane[4:9] for second_plane in duts_to_analyse[1][0]]
      dual_plane_var = [x + y for x, y in zip(duts_to_analyse[0][0], adapted_dut_channel)]
      dual_plane_duts = [dual_plane_var, duts_to_analyse[0][1]+" && "+duts_to_analyse[1][1], int(mcp_channel[0][1][0])]
      duts_vars_cuts.append(dual_plane_duts)
      arr_of_biases.append("MCP")
      channel_of_dut.append(int(mcp_channel[0][1][0]))

    hists_to_plot = []
    for j in range(len(duts_vars_cuts)):
      hist_down_up_dev = []
      for dut_var_ind, dut_var in enumerate(duts_vars_cuts[j][0]):
        thisHist = hist_tree_file_timeres(tree, file, dut_var, duts_vars_cuts[j][2], self.nBins, self.xLower, self.xUpper, arr_of_biases[j], duts_vars_cuts[j][1])
        num_ev = thisHist.GetEntries()
        if dut_var_ind == 0:
          hists_to_plot.append(thisHist)
        hist_down_up_dev.append(thisHist)
      arr_of_hists.append(hist_down_up_dev)
      arr_of_biases.append(bias)
      arr_of_nevents.append(num_ev)

    c1 = root.TCanvas("c1", f"Distribution time resolution", 800, 600)
    if self.log_scale:
      c1.SetLogy()

    max_y = max(hist.GetMaximum() for hist in hists_to_plot) * 1.05
    hists_to_plot[0].GetYaxis().SetRangeUser(1 if self.log_scale else 0, max_y)
    hists_to_plot[0].SetTitle(f"Distribution time resolution")
    hists_to_plot[0].Draw()
    if len(hists_to_plot) > 1:
      for hist_to_draw in hists_to_plot[1:]:
        hist_to_draw.Draw("SAME")

    for i, thisHist in enumerate(hists_to_plot):
      thisFit = plot_fit_curves(self.xLower, self.xUpper, "gaus", hists_to_plot[i], channel_of_dut[i], arr_of_biases[i])
      thisFit.Draw("SAME")

    legend = root.TLegend(0.7, 0.7, 0.9, 0.9)
    for i in range(len(hists_to_plot)):
      legend.AddEntry(hists_to_plot[i], arr_of_biases[i] + " CH " + str(channel_of_dut[i]), "l")
    legend.Draw()

    if not os.path.exists("timeres"):
      os.makedirs("timeres")
    c1.SaveAs("timeres/"+self.save_name)
    print(f"[BETA ANALYSIS]: [TIME RESOLUTION] Saved time resolution as timeres/"+self.save_name)

    arr_of_fits = []
    for i, nom_up_down_hists in enumerate(arr_of_hists):
      fit_down_up_dev = []
      for j, toa_thresh_hist in enumerate(nom_up_down_hists):
        thisFit = plot_fit_curves(self.xLower, self.xUpper, "gaus", toa_thresh_hist, channel_of_dut[i], arr_of_biases[i])
        fit_down_up_dev.append(thisFit)
      arr_of_fits.append(fit_down_up_dev)

    fit_results = get_fit_results_TR(arr_of_fits, arr_of_biases, arr_of_nevents, channel_of_dut, mcp_tr)
    return fit_results
