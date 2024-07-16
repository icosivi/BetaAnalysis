import ROOT
import numpy as np

file = ROOT.TFile("SPS_TestBeam_ROOT_files_Transcend_prontoFede/SPS_TestBeam-2024_220V.root")
df = ROOT.RDataFrame("wfm", file)

num_bins = 5000
w_min = df.Min("w1").GetValue()
w_max = df.Max("w1").GetValue()

hist_model = ROOT.RDF.TH1DModel("w1_hist", "w1 Distribution;w1;counts", num_bins, w_min, w_max)
hist_w = df.Histo1D(hist_model, "w1")

canvas_hist = ROOT.TCanvas("canvas_hist", "w1 Distribution", 800, 600)
hist_w.Draw()
canvas_hist.SaveAs("w1_distribution.png")

df_numpy = df.AsNumpy(columns=["event", "t1", "w1"])

t_values = [list(vec) for vec in df_numpy["t1"]]
w_values = [list(vec) for vec in df_numpy["w1"]]

df_numpy["t1"] = t_values
df_numpy["w1"] = w_values

print("Event values:", df_numpy["event"])
print("'t1' values:", df_numpy["t1"][520])
print("'w1' values:", df_numpy["w1"][520])
print("Total number of events:", len(df_numpy["event"]))

print(len(df_numpy["event"]))

example_event_index = 520
t_values = df_numpy["t1"][example_event_index]
w_values = df_numpy["w1"][example_event_index]
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

for event_index in range(len(df_numpy["t1"])):
  t_ev = df_numpy["t1"][event_index]
  w_ev = df_numpy["w1"][event_index]
  t_arr = np.array(t_ev)
  w_arr = np.array(w_ev)
  if len(t_arr) != 402:
    print(f"The event {event_index} has an unusual number of data points of {str(len(t_arr))}")
