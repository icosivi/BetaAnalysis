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

df = pd.DataFrame()
path = "/media/gp19133/EXTERNAL_USB/fastdaqtest/"
#path = "./binary_files/"
pattern = r'(\d{5})'

CVAL = "C2"

index = 0
total = int(len(os.listdir(path))/2)
total_minus_false_files = total - 3

for file in os.listdir(path):
    if CVAL+"--Trace--" in file:
        #C1--Trace--00000.trc
        match = re.search(pattern, file)
        if match:
            event_number = match.group(1).zfill(5)
        else:
            print("No 5-digit number found in the filename")
        if event_number == "00000" or event_number == "00001" or event_number == "00015":
            continue

        lecroy_data = lecroyparser.ScopeData(path+file)
        #print(dir(lecroy_data)) # get attributes of the ScopeData class
        xy_data = pd.DataFrame({'event': 0, 't': lecroy_data.x, 'w': lecroy_data.y})

        df = pd.concat([df, xy_data], ignore_index=True)
        #time_diff = df['t'].diff()
        #myXPoints = df['t'].iloc[[10002,20004,30006,40008,50010,60012,70014,80016,90018,100020]]

        index += 1
        print(f"{index} files processed of {total}")
        #plt.figure(figsize=(25, 6))
        #plt.plot(df['t'].iloc[0:100020],df['w'].iloc[0:100020])
        #for point in myXPoints:
        #    plt.axvline(x = point, linewidth=2, color='k')
        #plt.show()

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
    
        root_file = ROOT.TFile(f"SPS_TestBeam_ROOT_files_fastdaqtest_eventfix/SPS_TestBeam-2024_{CVAL}_{event_number}.root", "RECREATE")
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
        df = pd.DataFrame()
