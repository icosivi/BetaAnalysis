import ROOT
from ROOT import TFile, TTree, std
from array import array
import os
import re
import argparse
import sys

#path = "SPS_TestBeam_ROOT_files_Transcend_230V_comb_events/"
#output_file = "JustATest/SPS_TestBeam-2024_230V.root"


def main():
  parser = argparse.ArgumentParser(description='[L2 pre-proc]: merge ROOT files across channel into a single ROOT file.')
  parser.add_argument('-i', '--input', required=True, help='Path to output directory (do not need to state the files, just the path)')
  parser.add_argument('-o', '--output', required=True, help='Path of output directory (will create it if it doesn\'t exist)')

  args = parser.parse_args()

  file_pattern = re.compile(rf"SPS_TestBeam-2024_C[0-9]{{1}}\.root")
  all_files = os.listdir(args.input)
  input_files = [f for f in all_files if file_pattern.match(f)]

  entries = []

  input_trees = std.vector('TTree*')()
  input_files_handles = []

  # Open the input files and retrieve the TTrees
  for file_name in input_files:
    file = TFile.Open(args.input+"/"+file_name)
    input_files_handles.append(file)  # Keep file handles open
    tree = file.Get("wfm")
    if tree:
      input_trees.push_back(tree)
    else:
      print(f"[L2 pre-proc]: Error retrieving TTree 'wfm' from file: {args.input + file_name}")

  print("[L2 pre-proc]: Number of TTrees in input_trees:", input_trees.size())
  for tree in input_trees:
    tree.Print()

  # Determine the maximum number of events among all input files
  max_num_entries = 0
  for input_tree in input_trees:
    num_entries = input_tree.GetEntries()
    if num_entries > max_num_entries:
      max_num_entries = num_entries

  # Create the output file and TTree
  if not os.path.exists(args.output):
    print(f"[L2 pre-proc]: Output directory '{args.output}' does not exist. Creating it.")
    os.makedirs(args.output)
  output_root_file = TFile(args.output+"/SPS_TestBeam-2024.root", "RECREATE")
  combined_tree = TTree("wfm", "Combined tree with branches from all input files")

  # Define variables to hold data
  event = array('i', [0])
  w1 = std.vector('double')()
  w2 = std.vector('double')()
  w3 = std.vector('double')()
  w4 = std.vector('double')()
  t1 = std.vector('double')()
  t2 = std.vector('double')()
  t3 = std.vector('double')()
  t4 = std.vector('double')()

  # Create branches in the combined TTree
  combined_tree.Branch("event", event, "event/I")
  combined_tree.Branch("w1", w1)
  combined_tree.Branch("w2", w2)
  combined_tree.Branch("w3", w3)
  combined_tree.Branch("w4", w4)
  combined_tree.Branch("t1", t1)
  combined_tree.Branch("t2", t2)
  combined_tree.Branch("t3", t3)
  combined_tree.Branch("t4", t4)

  # Loop over the maximum number of entries and fill the combined TTree
  for i in range(max_num_entries):
    w1.clear()
    w2.clear()
    w3.clear()
    w4.clear()
    t1.clear()
    t2.clear()
    t3.clear()
    t4.clear()
    
    for idx, input_tree in enumerate(input_trees):
      if i < input_tree.GetEntries():
        input_tree.GetEntry(i)
        event[0] = input_tree.event
            
        if idx == 0:
          for val in input_tree.w:
            w1.push_back(val)
          for val in input_tree.t:
            t1.push_back(val)
        elif idx == 1:
          for val in input_tree.w:
            w2.push_back(val)
          for val in input_tree.t:
            t2.push_back(val)
        elif idx == 2:
          for val in input_tree.w:
            w3.push_back(val)
          for val in input_tree.t:
            t3.push_back(val)
        elif idx == 3:
          for val in input_tree.w:
            w4.push_back(val)
          for val in input_tree.t:
            t4.push_back(val)
      else:
        # If i is greater than the number of entries in this file, fill with dummy values
        event[0] = -999  # Choose a suitable dummy event value
        if idx == 0:
          w1.push_back(0.0)
          t1.push_back(0.0)
        elif idx == 1:
          w2.push_back(0.0)
          t2.push_back(0.0)
        elif idx == 2:
          w3.push_back(0.0)
          t3.push_back(0.0)
        elif idx == 3:
          w4.push_back(0.0)
          t4.push_back(0.0)
    
    combined_tree.Fill()

  # Write the combined TTree to the output file
  combined_tree.Write()
  output_root_file.Close()

  # Close the input files
  for file in input_files_handles:
    file.Close()

  print(f"[L2 pre-proc]: Combined TTree written to {args.output}")

if __name__ == "__main__":
  main()
