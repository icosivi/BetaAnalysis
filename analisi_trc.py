import os
import matplotlib.pyplot as plt
import pandas as pd
import lecroyparser 
import numpy as np
import ROOT


#data = lecroyparser.ScopeData(path,parseAll=True)
#print(type(data.x))
#print(data.x)
#plt.plot(data.x,data.y)
#plt.show()


df = pd.DataFrame()
path = "/Volumes/Elements/TB_SPS_ETL_sensors_June2024/220V/"

for event, file in enumerate(os.listdir(path)):
    if event>10:
        event_long = str(int(event)).zfill(5)
        #file = "C1--Trace--"+event_long+".trc" 
        lecroy = lecroyparser.ScopeData(path+file)
        xy_data = pd.DataFrame({'event': event_long, 'w': np.lecroy.x, 't': lecroy.y})
        #df = df.append(xy_data, ignore_index=True)
        df = pd.concat([df, xy_data], ignore_index=True)
    
data = {key: df[key].values for key in ['event', 'w', 't']}
rdf = ROOT.RDF.MakeNumpyDataFrame(data)
rdf.Snapshot('wfm', 'pippo.root')