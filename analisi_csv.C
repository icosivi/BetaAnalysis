// analisi_csv.C
// Extracts aligned waveform windows around the signal peak to a CSV file.
//
// Usage:
//   root -l -q 'analisi_csv.C(x_min, x_max, y_min, y_max)'
//   root -l -q 'analisi_csv.C(-3.0, 3.0, -3.0, 3.0, "output.csv")'
//
// Output CSV columns:
//   x_pos1, y_pos1, x_pos2, y_pos2,
//   w1[imax-20], ..., w1[imax], ..., w1[imax+20],
//   w2[imax-20], ..., w2[imax], ..., w2[imax+20],
//   ...
//   w17[imax-20], ..., w17[imax], ..., w17[imax+20]
//
// imax = index of the max sample in the channel with the highest pmax.

#include <TFile.h>
#include <TTree.h>
#include <TBranch.h>

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <algorithm>
#include <iomanip>

void analisi_csv(float x_min   = -5.0f,
                 float x_max   =  5.0f,
                 float y_min   = -5.0f,
                 float y_max   =  5.0f,
                 const char *outfile = "waveforms.csv")
{
    // ── Configuration ─────────────────────────────────────────────────────
    const char *RAW_FILE   = "/media/SSD_4TB/TB11/raw/Run5_DCsquare_500um_new_board1_tracker5_240V_20C_HG.root";
    const char *STATS_FILE = "/media/SSD_4TB/TB11/stats/stats_Run5_tracks.root";

    const int N_CH   = 17;          // channels: w1 … w17  /  pmax[1] … pmax[17]
    const int N_SIDE = 20;          // samples on each side of the peak
    const int N_WIN  = 2*N_SIDE+1;  // 41 samples per channel per event
    const int N_SAMP = 1024;        // samples per waveform in raw file

    // ── Open files ────────────────────────────────────────────────────────
    TFile *fraw   = TFile::Open(RAW_FILE,   "READ");
    TFile *fstats = TFile::Open(STATS_FILE, "READ");

    if (!fraw   || fraw->IsZombie())   { std::cerr << "Cannot open raw file\n";   return; }
    if (!fstats || fstats->IsZombie()) { std::cerr << "Cannot open stats file\n"; return; }

    TTree *traw   = (TTree*)fraw->Get("wfm");
    TTree *tstats = (TTree*)fstats->Get("Analysis");

    if (!traw)   { std::cerr << "Tree 'wfm' not found\n";      return; }
    if (!tstats) { std::cerr << "Tree 'Analysis' not found\n"; return; }

    Long64_t nEntries = tstats->GetEntries();
    std::cout << "Stats entries : " << nEntries << "\n"
              << "Raw entries   : " << traw->GetEntries() << "\n"
              << "Spatial cut   : x=[" << x_min << "," << x_max << "]"
              << "  y=[" << y_min << "," << y_max << "]\n"
              << "Output        : " << outfile << "\n\n";

    // ── Stats branches ────────────────────────────────────────────────────
    Long64_t          ev_raw = 0;
    float             x1=0, y1=0, x2=0, y2=0;
    std::vector<float> *pmax = nullptr;

    tstats->SetBranchStatus("*", 0);
    tstats->SetBranchStatus("event",  1);
    tstats->SetBranchStatus("x_pos1", 1);
    tstats->SetBranchStatus("y_pos1", 1);
    tstats->SetBranchStatus("x_pos2", 1);
    tstats->SetBranchStatus("y_pos2", 1);
    tstats->SetBranchStatus("pmax",   1);

    tstats->SetBranchAddress("event",  &ev_raw);
    tstats->SetBranchAddress("x_pos1", &x1);
    tstats->SetBranchAddress("y_pos1", &y1);
    tstats->SetBranchAddress("x_pos2", &x2);
    tstats->SetBranchAddress("y_pos2", &y2);
    tstats->SetBranchAddress("pmax",   &pmax);

    // ── Raw waveform branches (fixed C-arrays w1[1024] … w17[1024]) ──────
    // Disable all branches first, then enable only the ones we need.
    traw->SetBranchStatus("*", 0);
    for (int ch = 1; ch <= N_CH; ch++)
        traw->SetBranchStatus(Form("w%d", ch), 1);

    // One buffer per channel.
    float wbuf[N_CH+1][N_SAMP];   // wbuf[ch][sample], ch index 1..N_CH
    for (int ch = 1; ch <= N_CH; ch++) {
        traw->SetBranchAddress(Form("w%d", ch), wbuf[ch]);
    }

    // ── CSV output ────────────────────────────────────────────────────────
    std::ofstream csv(outfile);
    if (!csv.is_open()) { std::cerr << "Cannot open output: " << outfile << "\n"; return; }

    // Header
    csv << "x_pos1,y_pos1,x_pos2,y_pos2";
    for (int ch = 1; ch <= N_CH; ch++)
        for (int s = -N_SIDE; s <= N_SIDE; s++)
            csv << ",w" << ch << "_" << (s >= 0 ? "p" : "m") << std::abs(s);
    csv << "\n";

    // ── Main loop ─────────────────────────────────────────────────────────
    Long64_t n_selected = 0;

    for (Long64_t i = 0; i < nEntries; i++) {
        tstats->GetEntry(i);

        // Spatial cut
        if (x1 < x_min || x1 > x_max || y1 < y_min || y1 > y_max) continue;
        if (!pmax || (int)pmax->size() < N_CH+1) continue;

        // Channel with max pmax among [1..N_CH]
        int   best_ch  = 1;
        float best_val = (*pmax)[1];
        for (int ch = 2; ch <= N_CH; ch++) {
            if ((*pmax)[ch] > best_val) {
                best_val = (*pmax)[ch];
                best_ch  = ch;
            }
        }

        // Get corresponding raw event
        traw->GetEntry(ev_raw);

        // Find index of max sample in best channel
        const float *wbest = wbuf[best_ch];
        int imax = (int)(std::max_element(wbest, wbest + N_SAMP) - wbest);

        // Require enough margin on both sides
        if (imax < N_SIDE || imax + N_SIDE >= N_SAMP) continue;

        // Write CSV row
        csv << std::setprecision(6)
            << x1 << "," << y1 << "," << x2 << "," << y2;

        for (int ch = 1; ch <= N_CH; ch++) {
            for (int s = -N_SIDE; s <= N_SIDE; s++)
                csv << "," << wbuf[ch][imax + s];
        }
        csv << "\n";
        ++n_selected;

        if (n_selected % 5000 == 0)
            std::cout << "  Selected: " << n_selected
                      << "  (processed " << i+1 << "/" << nEntries << ")\r"
                      << std::flush;
    }

    csv.close();
    fraw->Close();
    fstats->Close();

    std::cout << "\n\nDone.\n"
              << "  Events selected : " << n_selected << "\n"
              << "  Output file     : " << outfile << "\n";
}
