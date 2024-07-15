import os
import matplotlib.pyplot as plt
import pandas as pd
import lecroyparser 
import numpy as np
import ROOT
import cppyy
from array import array
import re
import sys
import argparse

#path = "/media/gp19133/EXTERNAL_USB/fastdaqtest/"
#path = "./binary_files/"

def main():
  df = pd.DataFrame()

  parser = argparse.ArgumentParser(description='[L0 pre-proc]: merge ROOT files across channel into a single ROOT file.')
  parser.add_argument('-i', '--input', required=True, help='Path to input directory (do not need to state the files, just the path)')
  parser.add_argument('-o', '--output', required=True, help='Path of output directory (will create it if it doesn\'t exist)')
  parser.add_argument('-ch', '--channel', type=int, required=True, help='Channel number for the event files for that channel to be merged from L0')

  args = parser.parse_args()

  CVAL = "C"+str(args.channel)

  index = 0
  pattern = r'(\d{5})'

  total = int(len(os.listdir(args.input))/2)
  total_minus_false_files = total - 3

  for file in os.listdir(args.input):
    if CVAL+"--Trace--" in file:
      #C1--Trace--00000.trc
      match = re.search(pattern, file)
      if match:
        event_number = match.group(1).zfill(5)
      else:
        print("[L0 pre-proc]: No 5-digit number found in the filename")
      if event_number == "00000" or event_number == "00001" or event_number == "00015":
        continue

      if not os.path.exists(args.output):
        print(f"[L0 pre-proc]: Output directory '{args.output}' does not exist. Creating it.")
        os.makedirs(args.output)

      lecroy_data = lecroyparser.ScopeData(args.input+"/"+file)
      #print(dir(lecroy_data)) # get attributes of the ScopeData class
      xy_data = pd.DataFrame({'event': 0, 't': lecroy_data.x, 'w': lecroy_data.y})

      df = pd.concat([df, xy_data], ignore_index=True)

      index += 1
      print(f"[L0 pre-proc]: {index} files processed of {total}")

      df['event'] = np.repeat(np.arange(1000),10002)
      grouped_df = df.groupby('event').agg({
        't': lambda x: x.tolist(),
        'w': lambda x: x.tolist()
      }).reset_index()

      data = df.groupby('event').agg({'t': list, 'w': list}).reset_index()
      data['event'] = data['event']+int(event_number)*1000
      print(file)
      print(data)

      for idx, row in data.iterrows():
        t_vector = cppyy.gbl.std.vector['double']()
        w_vector = cppyy.gbl.std.vector['double']()
    
        for value in row['t']:
          t_vector.push_back(value)
        for value in row['w']:
          w_vector.push_back(value)

        data.at[idx, 't'] = t_vector
        data.at[idx, 'w'] = w_vector
    
      root_file = ROOT.TFile(f"{args.output}/SPS_TestBeam-2024_{CVAL}_{event_number}.root", "RECREATE")
      tree = ROOT.TTree("wfm", "Skimmed tree containing waveform information")

      event = array('i', [0])
      t_vector = ROOT.std.vector('double')()
      w_vector = ROOT.std.vector('double')()

      tree.Branch("event", event, "event/I")
      tree.Branch("t", t_vector)
      tree.Branch("w", w_vector)

      for idx, row in data.iterrows():
        event[0] = int(row['event'])
        t_vector.clear()
        w_vector.clear()
        for value in row['t']:
          t_vector.push_back(value)
        for value in row['w']:
          w_vector.push_back(value)
        tree.Fill()

      tree.Write()
      root_file.Close()

if __name__ == "__main__":
  main()
