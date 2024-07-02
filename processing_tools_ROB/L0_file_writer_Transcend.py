import os
import matplotlib.pyplot as plt
import pandas as pd
import lecroyparser 
import numpy as np
import ROOT
import cppyy
from array import array
import re

df = pd.DataFrame()
path = "/media/gp19133/Transcend/TB_SPS_ETL_sensors_June2024/230V/"
#path = "./binary_files/"
pattern = r'(\d{5})'

CVAL = "C4"

index = 0
total = int(len(os.listdir(path))/4)

for file in os.listdir(path):
    if CVAL+"--Trace--" in file:
        #C1--Trace--00000.trc
        match = re.search(pattern, file)
        if match:
            event_number = match.group(1).zfill(5)
        else:
            print("No 5-digit number found in the filename")
        lecroy_data = lecroyparser.ScopeData(path+file)
        #print(dir(lecroy_data)) # get attributes of the ScopeData class
        xy_data = pd.DataFrame({'event': event_number, 't': lecroy_data.x, 'w': lecroy_data.y})
        df = pd.concat([df, xy_data], ignore_index=True)
        index += 1
        if index % 1000 == 0:
            print(f"{index} events processed of {total}")
        #plt.plot(lecroy_data.x,lecroy_data.y)
        #plt.show()

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
    
            root_file = ROOT.TFile(f"SPS_TestBeam_ROOT_files_Transcend_230V/SPS_TestBeam-2024_{CVAL}_{index}.root", "RECREATE")
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
