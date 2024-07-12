import ROOT
from ROOT import TFile, TTree, std
from array import array  # Import the array module
import os
import re
import argparse
import sys

def main():
  parser = argparse.ArgumentParser(description='[L1 pre-proc]: merge ROOT files across a channel into a single ROOT file for said channel.')
  parser.add_argument('-i', '--input', required=True, help='Path to output directory (do not need to state the files, just the path)')
  parser.add_argument('-o', '--output', required=True, help='Path of output directory (will create it if it doesn\'t exist)')
  parser.add_argument('-ch', '--channel', type=int, required=True, help='Channel number for the event files for that channel to be merged from L0')

  args = parser.parse_args()

  CVAL = "C"+str(args.channel)
  file_pattern = re.compile(rf"SPS_TestBeam-2024_{CVAL}_[0-9]{{5}}\.root")
  all_files = os.listdir(args.input)
  input_files = [f for f in all_files if file_pattern.match(f)]

  output_file = args.output + "/SPS_TestBeam-2024_"+CVAL+".root"
  if not os.path.exists(args.output):
    print(f"[L1 pre-proc]: Output directory '{args.output}' does not exist. Creating it.")
    os.makedirs(args.output)

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
    file = TFile.Open(args.input+"/"+file_name)
    tree = file.Get("wfm")
    if file_name == input_files[0]:
      total_number_entries = tree.GetEntries()
    collect_entries(tree)
    print(f"[L1 pre-proc]: Collected entries in {file_name}")
    file.Close()

  # Sort entries by 'event'
  entries.sort(key=lambda x: x[0])
  print(f"[L1 pre-proc]: Sorting entries in all files")

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
  entry_counter = 0
  for entry in entries:
    event[0] = entry[0]
    w.clear()
    t.clear()
    for val in entry[1]:
      w.push_back(val)
    for val in entry[2]:
      t.push_back(val)
    sorted_tree.Fill()
    entry_counter += 1
    if entry_counter % total_number_entries == 0:
      print(f"[L1 pre-proc]: {entry_counter} events processed of total {total_number_entries*len(input_files)}")

  # Write the sorted TTree to the output file
  sorted_tree.Write()
  output_root_file.Close()

  print(f"[L1 pre-proc]: Combined and sorted TTree written to {output_file}")

if __name__ == "__main__":
  main()
