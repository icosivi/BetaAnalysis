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

from firstPass_ProcTools import get_fit_results, hist_tree_file_basics, plot_fit_curves, getBias

class plotVar:
  def __init__(self, var, nBins, xUpper, log_scale, save_name):

    self.var = var
    self.nBins = nBins
    self.xUpper = xUpper
    self.log_scale = log_scale
    self.save_name = save_name

  def run(self, file, file_index, tree, channel_array):

    arr_of_hists = []
    arr_of_biases = []

    root.gErrorIgnoreLevel = root.kWarning

    channel_of_dut = []
    for ch_ind, ch_val in enumerate(channel_array):
      bias = getBias(str(file), ch_ind)
      sensorType, AtQfactor, ansatz_pmax = ch_val
      if sensorType == 1:
        channel_of_dut.append(ch_ind)
        condition = f"pmax[{ch_ind}] > 0.0 && pmax[{ch_ind}] < {self.xUpper}"
        thisHist = hist_tree_file_basics(tree, file, self.var, ch_ind, self.nBins, self.xUpper, bias, condition, ch_ind)
      elif sensorType == 2:
        condition = f"pmax[{ch_ind}] > 0.0"
        thisHist = hist_tree_file_basics(tree, file, self.var, ch_ind, self.nBins, self.xUpper, bias, condition, ch_ind)
      else:
        thisHist = None
      arr_of_hists.append(thisHist)
      arr_of_biases.append(bias)

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

    legend = root.TLegend(0.7, 0.7, 0.9, 0.9)
    for i in range(len(valid_hists)):
      legend.AddEntry(valid_hists[i], arr_of_biases[i] + " CH " + str(i+1), "l")

    legend.Draw()
    if not os.path.exists(self.var):
      os.makedirs(self.var)
    c1.SaveAs(self.var+"/"+self.save_name)
    print(f"[BETA ANALYSIS]: [PLOTTER] Saved {self.var} as {self.var}/"+self.save_name)
