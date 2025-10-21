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

from firstPass_ClassPlotter import plotVar
from cardReader import read_text_card
from secondPass_Langaus import plot_langaus
from secondPass_Gaussian import plot_gaussian
from secondPass_TR import plot_TR
from export_data import direct_to_table

from SNDisc import SNDisc_extract_signal

def main():
  parser = argparse.ArgumentParser(description='Read a text card containing information and location of ROOT analysis files and plot distributions of corresponding variables.')
  parser.add_argument('config', type=str, help='Path to the configuration text card (e.g., config.txt)')
  args = parser.parse_args()

  print(f"[BETA ANALYSIS] : [CARD READER] Reading text card {args.config}.")
  config, thicknesses, mcp_specs = read_text_card(args.config)
  thicknesses = [thickness for thickness in thicknesses if thickness != "nDUT"]
  output_name = os.path.splitext(os.path.basename(args.config))[0]

  file_list = config.get('files', [])
  channels = config.get('channels', [[0, 1]] * 8)

  print(f"[BETA ANALYSIS] : [FILE READER] Reading files {file_list}.")

  while config['channels'] and config['channels'][-1][0] == 0:
    config['channels'].pop()

  for ch in config['channels']:
    if ch[0] == 0: ch[1] = 0

  channel_mapping = {
    1: "DUT",
    2: "MCP",
    3: "Reference sensor",
    0: "Unknown or unused channel"
  }

  board_mapping = {
    4.7: "on SC board",
    5: "on Mignone board",
    1: "unmounted",
    0: ""
  }

  root.gROOT.SetBatch(True)

  file_array = []
  tree_array = []
  output_name_array = []

  for pattern in file_list:
    root_files = glob.glob(pattern)
    output_name_const = output_name
    for root_file in root_files:
      try:
        theFile = root.TFile(root_file)
        file_array.append(theFile)
        tree_array.append(theFile.Get("Analysis"))
        output_name_w_bias = ""
        for ch_ind, ch_val in enumerate(config['channels']):
          if ch_val[0] == 1:
            bias_after_channel = re.search(rf"Ch{ch_ind}-(\d+)V_", pattern)
            output_name_w_bias = output_name_w_bias + f"_Ch{ch_ind}-" + bias_after_channel.group(1) +"V"
        output_name_const = "hist_" + output_name_const + output_name_w_bias
        output_name_array.append(output_name_const)
      except Exception as e:
        print(f"Error reading {root_file}: {e}")

  if len(file_array) == 0:
    print(f"[BETA ANALYSIS] : [FILE READER] No files found.")
    sys.exit(0)
  else:
    print(f"[BETA ANALYSIS] : [FILE READER] Total {len(file_array)} input ROOT files read.")

  if config['pass_criteria'][0] == True:
    sentence = "First pass to plot unabridged PMAX distributions to calculate ansatz PMAX cuts between signal and background"
  else:
    sentence = "Second pass to plot AMPLITUDE, RISE TIME, CHARGE, RMS, and TIME RESOLUTION distributions, with specified selections"

  print(f"[BETA ANALYSIS] : [CARD READER] " + sentence + " for the following setup:")
  for i, ch in enumerate(config['channels']):
    print(f"        CH {i} : {channel_mapping.get(ch[0])} {board_mapping.get(ch[1])}")

  data_out = []
  if config['pass_criteria'][0] == True:
    print(f"[BETA ANALYSIS]: [PLOTTER] Plotting PMAX distribution (note that on the first pass no selections are applied to the phase space)")
    for file_ind, file_real in enumerate(file_array):
      plot_pmax = plotVar("pmax", config['pass_criteria'][1], config['pass_criteria'][1], True, output_name_array[file_ind]+"_pmax.png")
      plot_pmax.run(file_real, file_ind, tree_array[file_ind], config['channels'])
  if config['pass_criteria'][0] == False:
    print(f"[BETA ANALYSIS] : [SIGNAL-NOISE DISCRIMINATOR] Running SNDisc NN to extract signal events from the ROOT file")
    signal_event_array = []
    amplitude_dfs = []
    if not os.path.exists("NN_performance"):
      os.makedirs("NN_performance")
    for file_ind, file_real in enumerate(file_array):
      signal_events, df_data = SNDisc_extract_signal(file_real, file_ind, tree_array[file_ind], config['channels'], config['pass_criteria'][1], "NN_performance/"+output_name_array[file_ind])
      signal_event_array.append(signal_events)
      amplitude_dfs.append(df_data)
    amplitude_data = pd.concat(amplitude_dfs, ignore_index=True)
    print(amplitude_data.sort_values(by=['Channel','Bias']))
    data_out.append(('amplitude', amplitude_data.sort_values(by=['Channel','Bias'])))

    # SIGNAL_EVENTS IS INFORMATION CONTAINING TRUE FALSE VALUES FOR THE INDICES OF THE EVENTS TO KEEP, SO NEED TO DO A BOOLEAN MATCH WITH ROOT DATA
    #charge_events = get_root_data(file_real, file_ind, tree_array[file_ind], signal_event_array, charge)
    #timeres_events = get_root_data(file_real, file_ind, tree_array[file_ind], signal_event_array, timeres)
    print(f"[BETA ANALYSIS] : [SIGNAL-NOISE DISCRIMINATOR] Signal events extracted and risetime, charge, RMS, and time resolution computation underway")
    risetime_dfs = []
    for file_ind, file_real in enumerate(file_array):
      df_data = plot_gaussian("risetime", file_real, file_ind, tree_array[file_ind], config['channels'], int(0.5*config['pass_criteria'][1]), signal_event_array[file_ind], output_name_array[file_ind]+"_risetime")
      risetime_dfs.append(df_data)
    risetime_data = pd.concat(risetime_dfs, ignore_index=True)
    print(risetime_data.sort_values(by=['Channel','Bias']))
    data_out.append(('risetime', risetime_data.sort_values(by=['Channel','Bias'])))
    charge_dfs = []
    for file_ind, file_real in enumerate(file_array):
      df_data = plot_langaus('charge', file_real, file_ind, tree_array[file_ind], config['channels'], int(0.5*config['pass_criteria'][1]), 0.2*config['pass_criteria'][1], signal_event_array[file_ind], output_name_array[file_ind]+"_charge")
      charge_dfs.append(df_data)
    charge_data = pd.concat(charge_dfs, ignore_index=True)
    print(charge_data.sort_values(by=['Channel','Bias']))
    data_out.append(('charge', charge_data.sort_values(by=['Channel','Bias'])))
    rms_dfs = []
    for file_ind, file_real in enumerate(file_array):
      df_data = plot_gaussian('rms', file_real, file_ind, tree_array[file_ind], config['channels'], int(0.3*config['pass_criteria'][1]), signal_event_array[file_ind], output_name_array[file_ind]+"_rms")
      rms_dfs.append(df_data)
    rms_data = pd.concat(rms_dfs, ignore_index=True)
    print(rms_data.sort_values(by=['Channel','Bias']))
    data_out.append(('rms', rms_data.sort_values(by=['Channel','Bias'])))
    print(f"[BETA ANALYSIS]: [TIME RESOLUTION] Performing Gaussian fit to DUT-MCP channels")
    time_res_dfs = []
    for file_ind, file_real in enumerate(file_array):
      df_data = plot_TR('timeres', file_real, file_ind, tree_array[file_ind], config['channels'], int(0.2*config['pass_criteria'][1]), signal_event_array[file_ind], output_name_array[file_ind]+"_TR", mcp_specs)
      time_res_dfs.append(df_data)
    time_res_data = pd.concat(time_res_dfs, ignore_index=True)
    print(time_res_data.sort_values(by=['Channel','Bias']))
    data_out.append(('timeres', time_res_data.sort_values(by=['Channel','Bias'])))

  if len(data_out) > 1: direct_to_table(data_out, config['channels'], output_name, thicknesses)

if __name__ == "__main__":
    main()
