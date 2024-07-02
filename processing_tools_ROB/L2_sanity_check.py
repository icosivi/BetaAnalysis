import ROOT
import numpy as np

file = ROOT.TFile("SPS_TestBeam_ROOT_files_Transcend_prontoFede/SPS_TestBeam-2024_220V.root")
df = ROOT.RDataFrame("wfm", file)

num_bins = 5000
w_min = df.Min("w2").GetValue()
w_max = df.Max("w2").GetValue()

hist_model = ROOT.RDF.TH1DModel("w2_hist", "w2 Distribution;w2;counts", num_bins, w_min, w_max)
hist_w = df.Histo1D(hist_model, "w2")

canvas_hist = ROOT.TCanvas("canvas_hist", "w2 Distribution", 800, 600)
hist_w.Draw()
canvas_hist.SaveAs("w2_distribution.png")

df_numpy = df.AsNumpy(columns=["event", "t2", "w2"])

t_values = [list(vec) for vec in df_numpy["t2"]]
w_values = [list(vec) for vec in df_numpy["w2"]]

df_numpy["t2"] = t_values
df_numpy["w2"] = w_values

print("Event values:", df_numpy["event"])
print("'t2' values:", df_numpy["t2"][520])
print("'w2' values:", df_numpy["w2"][520])
print("Total number of events:", len(df_numpy["event"]))

print(len(df_numpy["event"]))

example_event_index = 520
t_values = df_numpy["t2"][example_event_index]
w_values = df_numpy["w2"][example_event_index]
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
