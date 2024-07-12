import os
import matplotlib.pyplot as plt
import pandas as pd
import lecroyparser 
import numpy as np
import ROOT
import cppyy
from array import array
import re
import argparse

#path = "/media/gp19133/Transcend/TB_SPS_ETL_sensors_June2024/220V/"
#path = "./binary_files/"

def main():
  df = pd.DataFrame()

  parser = argparse.ArgumentParser(description='[L0 pre-proc]: merge ROOT files across a channel into a single ROOT file for said channel.')
  parser.add_argument('-i', '--input', required=True, help='Path to output directory (do not need to state the files, just the path)')
  parser.add_argument('-o', '--output', required=True, help='Path of output directory (will create it if it doesn\'t exist)')
  parser.add_argument('-ch', '--channel', type=int, required=True, help='Channel number for the event files for that channel to be merged from L0')

  args = parser.parse_args()

  CVAL = "C"+str(args.channel)
  pattern = r'(\d{5})'

  index = 0
  total = int(len(os.listdir(args.input))/4)

  for file in os.listdir(args.input):
    if CVAL+"--Trace--" in file:
      #C1--Trace--00000.trc
      match = re.search(pattern, file)
      if match:
        event_number = match.group(1).zfill(5)
      else:
        print("[L0 pre-proc]: No 5-digit number found in the filename")

      if not os.path.exists(args.output):
        print(f"[L0 pre-proc]: Output directory '{args.output}' does not exist. Creating it.")
        os.makedirs(args.output)

      lecroy_data = lecroyparser.ScopeData(args.input+"/"+file)
      #print(dir(lecroy_data)) # get attributes of the ScopeData class
      if int(event_number) < 15:
        xy_data = pd.DataFrame({'event': event_number, 't': np.linspace(-1.00720939e-08, 1.00279064e-08, 402), 'w': np.zeros(402)})
      else:
        xy_data = pd.DataFrame({'event': event_number, 't': lecroy_data.x, 'w': lecroy_data.y})
      df = pd.concat([df, xy_data], ignore_index=True)
      index += 1
      if index % 1000 == 0:
        print(f"[L0 pre-proc]: {index} events processed of {total}")

      if index < 10:
        continue

      if index % 10000 == 0:
        data = df.groupby('event').agg({'t': list, 'w': list}).reset_index()
        for idx, row in data.iterrows():
          t_vector = cppyy.gbl.std.vector['double']()
          w_vector = cppyy.gbl.std.vector['double']()
    
          for value in row['t']:
            t_vector.push_back(value)
          for value in row['w']:
            w_vector.push_back(value)

          data.at[idx, 't'] = t_vector
          data.at[idx, 'w'] = w_vector
    
        root_file = ROOT.TFile(f"{args.output}/SPS_TestBeam-2024_{CVAL}_{index}.root", "RECREATE")
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
