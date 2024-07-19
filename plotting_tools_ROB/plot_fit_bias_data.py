import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import minimize
from scipy.stats import poisson
import ROOT as root
from ROOT import TF1
from scipy.special import gammaln
import math
import pandas as pd
import argparse
import glob
import re
import os
import csv

var_dict = {"tmax":"t_{max} / ns" , "pmax":"p_max / mV" , "negpmax":"-p_max / mV", "charge":"Q / fC", "area_new":"Area / pWb" , "rms":"RMS / mV"}

def main():
  parser = argparse.ArgumentParser(description='Read CSV fit data against bias and plot charge dist, RMS noise, and time resolution.')
  parser.add_argument('incsv', metavar='F', type=str, nargs='+',
                      help='Input CSV file.')
  #parser.add_argument('--doWaveform', action='store_true', help='Enable the doWaveform option')
  #parser.add_argument('--csvOut', action='store_true', help='Output important data for plots against bias')
  args = parser.parse_args()

  if os.path.exists(incsv):
    existing_df = pd.read_csv(incsv)
    print(f"Successfully read {incsv}.")
  else:
    print(f"File {incsv} doesn't exist.")

  #if args.doWaveform: waveform(file_array,tree_array,args.ch-1,total_number_channels)
  #if args.csvOut: makeCSV(file_array,tree_array,total_number_channels)

if __name__ == "__main__":
    main()
