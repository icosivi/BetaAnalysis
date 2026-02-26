#!/usr/bin/env python3
"""Utility to produce stats files more comfortably."""
import argparse
import os
import re
# import sys
import ROOT

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
DEFAULT_CFG = os.path.join(SCRIPT_DIR, "beta_config.ini")
SRC_FILES = [
    os.path.join(SCRIPT_DIR, "src/general.cpp"),
    os.path.join(SCRIPT_DIR, "src/Chameleon.cpp"),
    os.path.join(SCRIPT_DIR, "src/ConfigFile.cpp"),
    os.path.join(SCRIPT_DIR, "src/Analyzer.cpp"),
    os.path.join(SCRIPT_DIR, "analisi.C"),
]

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("input_files", nargs="+", metavar="input_file",
                    help="Input root files from BetaDAQ")
parser.add_argument("-c", "--config-file", default=DEFAULT_CFG,
                    help="Config file (in&out files will be overridden); default: %(default)s")
parser.add_argument("--out-prefix", default="stats_",
                    help="Prefix for the output file names; default: %(default)s")
args = parser.parse_args()

ROOT.gROOT.SetBatch(True)

for src in SRC_FILES:
    print(f"Compiling/loading {src}")
    ROOT.gSystem.CompileMacro(src, "kg")
print()

for ifp in args.input_files:
    # # Use fork to process each file in an isolated process
    # pid = os.fork()
    # if pid != 0:
    #     print("[PARENT] Waiting child", pid)
    #     os.waitpid(pid, 0)
    #     print("[PARENT] Child done")
    #     continue
    # print("[CHILD] Running analysis")

    # Pick the output file name and write a modified config file
    ifd, ifn = os.path.split(ifp)
    ofp = os.path.join(ifd, f"{args.out_prefix}{ifn}")

    # Run the analysis
    ROOT.analisi(args.config_file, ifp, ofp)

    # # Exit the child (forked) process and return control to main process
    # print("[CHILD] Done")
    # sys.exit()
