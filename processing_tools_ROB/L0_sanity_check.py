import ROOT
import numpy as np

file1 = ROOT.TFile("SPS_TestBeam_ROOT_files_Transcend_220V/SPS_TestBeam-2024_C1_80000.root")
file2 = ROOT.TFile("SPS_TestBeam_ROOT_files_Transcend_220V/SPS_TestBeam-2024_C1_60000.root")
df = ROOT.RDataFrame("wfm", file1)
df_comp = ROOT.RDataFrame("wfm", file2)

num_bins = 5000
w_min = df.Min("w").GetValue()
w_max = df.Max("w").GetValue()

hist_model = ROOT.RDF.TH1DModel("w_hist", "w Distribution;w;counts", num_bins, w_min, w_max)
hist_w = df.Histo1D(hist_model, "w")

canvas_hist = ROOT.TCanvas("canvas_hist", "w Distribution", 800, 600)
hist_w.Draw()
canvas_hist.SaveAs("w_distribution.png")

df_numpy = df.AsNumpy(columns=["event", "t", "w"])

t_values = [list(vec) for vec in df_numpy["t"]]
w_values = [list(vec) for vec in df_numpy["w"]]

# Add converted lists back to the numpy dictionary
df_numpy["t"] = t_values
df_numpy["w"] = w_values

# Print the values
print("Event values:", df_numpy["event"])
print("'t' values:", df_numpy["t"][920])
print("'w' values:", df_numpy["w"][920])
print("Total number of events:", len(df_numpy["event"]))

print(len(df_numpy["event"]))
example_event_index = 920
t_values = df_numpy["t"][example_event_index]
w_values = df_numpy["w"][example_event_index]
t_array = np.array(t_values)
w_array = np.array(w_values)

graph = ROOT.TGraph(len(t_array), t_array, w_array)
canvas_graph = ROOT.TCanvas("canvas_graph", "t vs w", 800, 600)
graph.SetTitle("Example Event: t vs w; t / ns; w / pWb")
graph.SetMarkerStyle(20)
graph.SetMarkerSize(1)
graph.Draw("APL")

canvas_graph.SaveAs("t_vs_w_example_event.png")
canvas_graph.Draw()

# Check timing of waveforms that they overlap and don't add consecutively
canvas_tw = ROOT.TCanvas("canvas_graph", "t vs w", 800, 600)
t_1000ev_1 = []
w_1000ev_1 = []
t_1000ev_2 = []
w_1000ev_2 = []

df_numpy2 = df_comp.AsNumpy(columns=["event", "t", "w"])
t_values = [list(vec) for vec in df_numpy2["t"]]
w_values = [list(vec) for vec in df_numpy2["w"]]
df_numpy2["t"] = t_values
df_numpy2["w"] = w_values

for event_index in range(len(df_numpy["t"])):
  t_ev = df_numpy["t"][event_index]
  w_ev = df_numpy["w"][event_index]
  t_arr = np.array(t_ev)
  w_arr = np.array(w_ev)
  if event_index == 1:
    print(len(t_arr))
  t_1000ev_1.append(t_arr.tolist())
  w_1000ev_1.append(w_arr.tolist())

for event_index in range(len(df_numpy2["t"])):
  t_ev = df_numpy2["t"][event_index]
  w_ev = df_numpy2["w"][event_index]
  t_arr = np.array(t_ev)
  w_arr = np.array(w_ev)
  if event_index == 1:
    print(len(t_arr))
  t_1000ev_2.append(t_arr.tolist())
  w_1000ev_2.append(w_arr.tolist())


print(len(t_1000ev_1))
print(len(t_1000ev_2))
t_1000ev_1 = np.array([item for sublist in t_1000ev_1 for item in sublist])
w_1000ev_1 = np.array([item for sublist in w_1000ev_1 for item in sublist])
t_1000ev_2 = np.array([item for sublist in t_1000ev_2 for item in sublist])
w_1000ev_2 = np.array([item for sublist in w_1000ev_2 for item in sublist])
print(len(t_1000ev_1))
graph_tw = ROOT.TGraph(len(t_1000ev_1), t_1000ev_1, w_1000ev_1)
graph_tw.SetTitle("Distribution of 1000 events; t vs w; t / ns; w / pWb")
graph_tw.SetMarkerStyle(20)
graph_tw.SetMarkerSize(1)
graph_tw.SetLineColor(ROOT.kRed)
graph_tw.SetMarkerColor(ROOT.kRed)
graph_tw.Draw("APL")

graph_tw_2 = ROOT.TGraph(len(t_1000ev_2), t_1000ev_2, w_1000ev_2)
graph_tw_2.SetMarkerStyle(20)
graph_tw_2.SetMarkerSize(1)
graph_tw_2.SetLineColor(ROOT.kBlue)
graph_tw_2.SetMarkerColor(ROOT.kBlue)
graph_tw_2.Draw("PL SAME")

legend = ROOT.TLegend(0.7, 0.7, 0.9, 0.9)
legend.AddEntry(graph_tw, "1st set Ev", "lp")
legend.AddEntry(graph_tw_2, "2nd set Ev", "lp")
legend.Draw()

canvas_tw.SaveAs("t_vs_w_1000ev_l0.png")
canvas_tw.Draw()
