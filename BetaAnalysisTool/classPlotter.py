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

from proc_tools import get_fit_results, hist_tree_file_basics, plot_fit_curves, getBias

class plotVar:
  def __init__(self, var, nBins, xLower, xUpper, log_scale, save_name, fit=None):

    self.var = var
    self.nBins = nBins
    self.xLower = xLower
    self.xUpper = xUpper
    self.log_scale = log_scale
    self.save_name = save_name
    self.fit = fit

  def run(self, file, file_index, tree, channel_array):

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

    channel_of_dut = []
    for j in range(len(channel_array)):
      bias = getBias(str(file), j)
      if (channel_array[j][0] == 1) or (channel_array[j][0] == 2):
        if (channel_array[j][0] == 1): channel_of_dut.append(j)
        thisHist = hist_tree_file_basics(tree, file, self.var, j, self.nBins, self.xLower, self.xUpper, bias, result[j], j)
        num_ev = thisHist.GetEntries()
      else:
        thisHist = None
        num_ev = 0
      arr_of_hists.append(thisHist)
      arr_of_biases.append(bias)
      arr_of_nevents.append(num_ev)

    c1 = root.TCanvas("c1", f"Distribution {self.var}", 800, 600)
    if self.log_scale:
      c1.SetLogy()

    valid_hists = [hist for hist in arr_of_hists if hist is not None]
    arr_of_biases = [bias for bias in arr_of_biases if bias is not None]

    max_y = max(hist.GetMaximum() for hist in valid_hists) * 1.05
    valid_hists[0].GetYaxis().SetRangeUser(1 if self.log_scale else 0, max_y)
    valid_hists[0].SetTitle(f"Distribution {self.var}")
    valid_hists[0].Draw()
    if len(valid_hists) > 1:
      for hist_to_draw in valid_hists[1:]:
        hist_to_draw.Draw("SAME")

    if (self.fit) is not None:
      arr_of_fits = []
      for i, thisHist in enumerate(arr_of_hists):
        if (channel_array[i][0] == 1):
          thisFit = plot_fit_curves(self.xLower, self.xUpper, self.fit, arr_of_hists[i], i, arr_of_biases[i])
          arr_of_fits.append(thisFit)
          thisFit.Draw("SAME")
        else:
          arr_of_fits.append(None)

    legend = root.TLegend(0.7, 0.7, 0.9, 0.9)
    for i in range(len(valid_hists)):
      legend.AddEntry(valid_hists[i], arr_of_biases[i] + " CH " + str(i+1), "l")

    legend.Draw()
    if not os.path.exists(self.var):
      os.makedirs(self.var)
    c1.SaveAs(self.var+"/"+self.save_name)
    print(f"[BETA ANALYSIS]: [PLOTTER] Saved {self.var} as {self.var}/"+self.save_name)

    if (self.fit) is not None:
      arr_of_fits = [fit for fit in arr_of_fits if fit is not None]
      fit_results = get_fit_results(arr_of_fits, arr_of_biases, arr_of_nevents, channel_of_dut)
      return fit_results
