import ROOT
import numpy as np
import sys

#file = ROOT.TFile("SPS_TestBeam_ROOT_files_Transcend_220V/SPS_TestBeam-2024_C1_80000.root")
#file = ROOT.TFile("SPS_TestBeam_ROOT_files_Transcend_220V_comb_events/SPS_TestBeam-2024_C3.root")
file = ROOT.TFile("SPS_TestBeam_ROOT_files_fastdaqtest_comb_events/SPS_TestBeam-2024_C1.root")
#file = ROOT.TFile("SPS_TestBeam_ROOT_files_Transcend_prontoFede/SPS_TestBeam-2024_220V.root")
df = ROOT.RDataFrame("wfm", file)

num_bins = 5000
w_min = df.Min("w").GetValue()
w_max = df.Max("w").GetValue()
#w_min = df.Min("w3").GetValue()
#w_max = df.Max("w3").GetValue()

hist_model = ROOT.RDF.TH1DModel("w_hist", "w Distribution;w;counts", num_bins, w_min, w_max)
hist_w = df.Histo1D(hist_model, "w")
#hist_w = df.Histo1D(hist_model, "w3")

canvas_hist = ROOT.TCanvas("canvas_hist", "w Distribution", 800, 600)
hist_w.Draw()
canvas_hist.SaveAs("w_distribution.png")

df_numpy = df.AsNumpy(columns=["event", "t", "w"])
print("Fine here")

t_values = []
w_values = []
for i in range(200):
  t_value = df_numpy["t"][i]
  w_value = df_numpy["w"][i]
  t_values.append(t_value)
  w_values.append(w_value)
#t_values = [list(vec) for vec in df_numpy["t"]]
#w_values = [list(vec) for vec in df_numpy["w"]]

#print("Event values:", df_numpy["event"])
#print("'t' values:", df_numpy["t"][116])
#print("'w' values:", df_numpy["w"][116])
print("Total number of events:", len(df_numpy["event"]))

example_event_index = 16
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
sys.exit()
# Check timing of waveforms that they overlap and don't add consecutively
canvas_tw = ROOT.TCanvas("canvas_graph", "t vs w", 800, 600)
t_1000ev_1 = []
w_1000ev_1 = []
t_1000ev_2 = []
w_1000ev_2 = []
t_1000ev_3 = []
w_1000ev_3 = []

for event_index in range(len(df_numpy["t"])):
  t_ev = df_numpy["t"][event_index]
  w_ev = df_numpy["w"][event_index]
  t_arr = np.array(t_ev)
  w_arr = np.array(w_ev)
  if len(t_arr) != 402:
    print(f"The event {event_index} has an unusual number of data points of {str(len(t_arr))}")
  if (event_index >= 0) and event_index < 10000:
    if event_index == 6:
      print(len(t_arr))
    t_1000ev_1.append(t_arr.tolist())
    w_1000ev_1.append(w_arr.tolist())
  elif (event_index > 10000) and (event_index < 20000):
    if event_index == 10001:
      print(len(t_arr))
    t_1000ev_2.append(t_arr.tolist())
    w_1000ev_2.append(w_arr.tolist())
  elif (event_index > 20000) and (event_index < 30000):
    if event_index == 20001:
      print(len(t_arr))
    t_1000ev_3.append(t_arr.tolist())
    w_1000ev_3.append(w_arr.tolist())


t_1000ev_1 = np.array([item for sublist in t_1000ev_1 for item in sublist])
w_1000ev_1 = np.array([item for sublist in w_1000ev_1 for item in sublist])
t_1000ev_2 = np.array([item for sublist in t_1000ev_2 for item in sublist])
w_1000ev_2 = np.array([item for sublist in w_1000ev_2 for item in sublist])
t_1000ev_3 = np.array([item for sublist in t_1000ev_3 for item in sublist])
w_1000ev_3 = np.array([item for sublist in w_1000ev_3 for item in sublist])


graph_tw_1 = ROOT.TGraph(len(t_1000ev_2), t_1000ev_1[:402], w_1000ev_1[:402])
graph_tw_1.SetTitle("Distribution of 1000 events; t vs w; t / ns; w / pWb")
graph_tw_1.SetMarkerStyle(20)
graph_tw_1.SetMarkerSize(1)
graph_tw_1.SetLineColor(ROOT.kRed)
graph_tw_1.SetMarkerColor(ROOT.kRed)
graph_tw_1.Draw("APL")

graph_tw_2 = ROOT.TGraph(len(t_1000ev_2), t_1000ev_2, w_1000ev_2)
graph_tw_2.SetMarkerStyle(20)
graph_tw_2.SetMarkerSize(1)
graph_tw_2.SetLineColor(ROOT.kBlue)
graph_tw_2.SetMarkerColor(ROOT.kBlue)
graph_tw_2.Draw("PL SAME")

graph_tw_3 = ROOT.TGraph(len(t_1000ev_3), t_1000ev_3, w_1000ev_3)
graph_tw_3.SetMarkerStyle(20)
graph_tw_3.SetMarkerSize(1)
graph_tw_3.SetLineColor(ROOT.kGreen)
graph_tw_3.SetMarkerColor(ROOT.kGreen)
graph_tw_3.Draw("PL SAME")

legend = ROOT.TLegend(0.7, 0.7, 0.9, 0.9)
legend.AddEntry(graph_tw_1, "First 10000 Ev", "lp")
legend.AddEntry(graph_tw_2, "Second 10000 Ev", "lp")
legend.AddEntry(graph_tw_3, "Third 10000 Ev", "lp")
legend.Draw()

canvas_tw.SaveAs("t_vs_w_1000ev.png")
canvas_tw.Draw()
