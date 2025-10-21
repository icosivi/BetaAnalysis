import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

def perform_diagnostics(df_data, ch_index, savename):
  colours = ["r","orange","yellow","lime","green","blue","purple","magenta"]
  markers = ["o","v","s","^","D","p","d","h"]

  fig, axes = plt.subplots(1, 5, figsize=(25, 5))
  for i, y_column in enumerate(['Signal event count', 'Amplitude MPV', 'SSE score', 'Red. Chi2', 'Mod. Chi2']):
    ax = axes[i]
    for j, label in enumerate(df_data['Number of epochs'].unique()):
      subset = df_data[df_data['Number of epochs'] == label]
      ax.plot(subset['Probability threshold'], subset[y_column], label=f'{label} ({y_column})', color=colours[j], marker=markers[j], markersize=8, markeredgecolor='black', markeredgewidth=1)
      ax.set_xlabel('Probability threshold')
      ax.set_ylabel(y_column)
      if i >= 2:
        ax.set_yscale("log")
      ax.legend()

  plt.tight_layout()
  plt.savefig(savename+"_Ch"+str(ch_index)+"_pdperformanceplots.png",facecolor='w')
  df_data.to_csv(savename + "diagnostics_SNDisc_Ch"+str(ch_index)+".csv", index=False)
  return
