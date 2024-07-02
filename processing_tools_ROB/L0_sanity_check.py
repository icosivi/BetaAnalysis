import ROOT
import numpy as np

file = ROOT.TFile("SPS_TestBeam_ROOT_files_Transcend_220V/SPS_TestBeam-2024_C1_80000.root")
df = ROOT.RDataFrame("wfm", file)

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
