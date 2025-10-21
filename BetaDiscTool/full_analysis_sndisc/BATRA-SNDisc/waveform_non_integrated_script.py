# waveform.py

from proc_tools import get_fit_results, hist_tree_file_basics, plot_fit_curves, getBias

def plot_waveform(files,trees,ch,total_number_channels):
  print("Waveform analysis")
  t_data = [[] for _ in range(total_number_channels)]
  w_data = [[] for _ in range(total_number_channels)]
  event_limit = 24500
  event_start = 4500
  for i in range(len(trees)):
    tree = trees[i]
    if ch != -1:
      event_cycle = 0
      for entry_index in range(event_start, event_limit+event_start):
        tree.GetEntry(entry_index)
        if (event_cycle > event_limit): continue
        event_cycle += 1
        t = tree.t
        w = tree.w
        pmax = tree.pmax
        if 0 < pmax[ch] < 200:
          t_data[0].extend(t[ch])
          w_data[0].extend(w[ch])
    else:
      for ch_it in range(total_number_channels):
        print("Channel number:")
        print(ch_it)
        event_cycle = 0
        for entry_index in range(event_start, event_limit+event_start):
          tree.GetEntry(entry_index)
          if (event_cycle > event_limit): continue
          event_cycle += 1
          t = tree.t
          w = tree.w
          pmax = tree.pmax
          if 0 < pmax[ch_it] < 200:
            t_data[ch_it].extend(t[ch_it])
            w_data[ch_it].extend(w[ch_it])
    
  t_data = [np.array(channel) for channel in t_data]
  w_data = [np.array(channel) for channel in w_data]

  fig, ax = plt.subplots(figsize=(10, 6))

  colors = plt.cm.viridis(np.linspace(0, 1, total_number_channels))

  for i in range(3):
    ax.scatter(t_data[i], w_data[i], label=f'Channel {i+1}', color=colors[i], s=1)

  ax.set_xlim(-0.5e-8,1e-8)
  ax.set_ylim(min(min(inner_array) for inner_array in w_data),max(max(inner_array) for inner_array in w_data))
  ax.set_xlabel('t / 10 ns', fontsize=14)
  ax.set_ylabel('Amplitude / V', fontsize=14)
  ax.set_title(f'Waveform for events between Ev No.{event_start} and Ev No.{event_limit}', fontsize=16)
  ax.legend(loc = "upper right")
  ax.grid(True)
  plt.tight_layout()
  plt.savefig("waveform_comparison.png",facecolor='w')
