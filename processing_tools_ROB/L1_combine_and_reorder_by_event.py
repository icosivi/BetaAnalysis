import ROOT
from ROOT import TFile, TTree, std
from array import array  # Import the array module

CVAL = "C1"
path = "SPS_TestBeam_ROOT_files_Transcend_230V/"
#input_files = ["SPS_TestBeam-2024_"+CVAL+"_100000.root", "SPS_TestBeam-2024_"+CVAL+"_80000.root", "SPS_TestBeam-2024_"+CVAL+"_10000.root",
#               "SPS_TestBeam-2024_"+CVAL+"_90000.root", "SPS_TestBeam-2024_"+CVAL+"_20000.root", "SPS_TestBeam-2024_"+CVAL+"_30000.root",
#               "SPS_TestBeam-2024_"+CVAL+"_40000.root", "SPS_TestBeam-2024_"+CVAL+"_50000.root", "SPS_TestBeam-2024_"+CVAL+"_60000.root",
#               "SPS_TestBeam-2024_"+CVAL+"_70000.root"]

input_files = ["SPS_TestBeam-2024_"+CVAL+"_10000.root","SPS_TestBeam-2024_"+CVAL+"_20000.root", "SPS_TestBeam-2024_"+CVAL+"_30000.root",
               "SPS_TestBeam-2024_"+CVAL+"_40000.root", "SPS_TestBeam-2024_"+CVAL+"_50000.root", "SPS_TestBeam-2024_"+CVAL+"_60000.root",
               "SPS_TestBeam-2024_"+CVAL+"_70000.root"]

output_file = "SPS_TestBeam_ROOT_files_Transcend_230V_comb_events/SPS_TestBeam-2024_"+CVAL+".root"

entries = []

# Function to collect entries from a tree
def collect_entries(tree):
    for entry in tree:
        event = entry.event
        w = std.vector('double')(entry.w)
        t = std.vector('double')(entry.t)
        entries.append((event, w, t))

# Collect entries from all input files
for file_name in input_files:
    file = TFile.Open(path + file_name)
    tree = file.Get("wfm")
    collect_entries(tree)
    file.Close()

# Sort entries by 'event'
entries.sort(key=lambda x: x[0])

# Create output ROOT file and TTree
output_root_file = TFile(output_file, "RECREATE")
sorted_tree = TTree("wfm", "Skimmed tree containing waveform information ordered by event number")

# Define variables to hold data
event = array('i', [0])
w = std.vector('double')()
t = std.vector('double')()

# Create branches in the TTree
sorted_tree.Branch("event", event, "event/I")
sorted_tree.Branch("w", w)
sorted_tree.Branch("t", t)

# Fill the sorted TTree with entries
for entry in entries:
    event[0] = entry[0]
    w.clear()
    t.clear()
    for val in entry[1]:
        w.push_back(val)
    for val in entry[2]:
        t.push_back(val)
    sorted_tree.Fill()

# Write the sorted TTree to the output file
sorted_tree.Write()
output_root_file.Close()

print(f"Combined and sorted TTree written to {output_file}")
