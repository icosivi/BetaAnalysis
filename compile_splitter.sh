#!/bin/bash
# Compile rawfile_splitter against the local EUDAQ2 and ROOT installations.

EUDAQ=/home/tb_pc/Desktop/TestBeam/Tracking/eudaq
ROOT_INST=/home/tb_pc/Desktop/DigiDAQ/ROOT/install

g++ -std=c++17 -O2 \
    -I${EUDAQ}/main/lib/core/include \
    -I${EUDAQ}/extern/include \
    -I${EUDAQ}/build/include \
    -I${ROOT_INST}/include \
    rawfile_splitter.cc \
    -L${EUDAQ}/lib \
    -leudaq_core \
    -leudaq_module_tlu \
    -leudaq_module_eudet \
    -leudaq_module_adenium \
    -Wl,-rpath,${EUDAQ}/lib \
    -L${ROOT_INST}/lib \
    -lCore -lTree -lRIO -lNet -lThread \
    -Wl,-rpath,${ROOT_INST}/lib \
    -o rawfile_splitter \
    && echo "OK: ./rawfile_splitter ready" \
    || echo "FAILED"
