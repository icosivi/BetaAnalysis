// rawfile_splitter.cc
// Splits an EUDAQ2 native .raw file (and optionally a ROOT digitizer file)
// at a user-specified list of cumulative event counts.
//
// Compile with compile_splitter.sh, then:
//
//   ./rawfile_splitter input.raw output_base cut1 [cut2 cut3 ...]
//   ./rawfile_splitter input.raw output_base cut1 [cut2 ...] --root input.root
//
// cut_N = cumulative data-event count after which a new chunk is opened.
// Example: cuts 100037 200500 → chunks [0,100037) [100037,200500) [200500,end)
//
// Output: output_base_000000.raw [+.root], _000001.raw [+.root], ...

#include <iostream>
#include <string>
#include <memory>
#include <deque>
#include <algorithm>
#include <cstdio>
#include <cstdlib>

// EUDAQ2
#include "eudaq/FileDeserializer.hh"
#include "eudaq/FileSerializer.hh"
#include "eudaq/Factory.hh"
#include "eudaq/Event.hh"

// ROOT
#include <TFile.h>
#include <TTree.h>

int main(int argc, char* argv[]) {
    if (argc < 4) {
        std::cerr << "Usage: " << argv[0]
                  << " input.raw output_base cut1 [cut2 ...] [--root input.root]\n\n"
                  << "  input.raw    – EUDAQ2 native raw file\n"
                  << "  output_base  – prefix for output files\n"
                  << "  cut1 cut2 …  – cumulative event counts at which to split\n"
                  << "  --root FILE  – (optional) ROOT digitizer file, split in sync\n\n"
                  << "Example:\n"
                  << "  " << argv[0]
                  << " run7.raw /out/run7 100037 200500 300001 --root run7.root\n"
                  << "  → run7_000000.raw+.root (100037 ev)\n"
                  << "    run7_000001.raw+.root (100463 ev)\n"
                  << "    run7_000002.raw+.root (99501 ev)\n"
                  << "    run7_000003.raw+.root (rest)\n";
        return 1;
    }

    const std::string input_raw  = argv[1];
    const std::string output_base = argv[2];

    // Parse remaining args: numbers → cut points, --root FILE → ROOT input
    std::deque<size_t> cut_points;
    bool        has_root   = false;
    std::string input_root;

    for (int i = 3; i < argc; i++) {
        std::string arg = argv[i];
        if (arg == "--root" && i + 1 < argc) {
            input_root = argv[++i];
            has_root   = true;
        } else {
            cut_points.push_back(std::stoul(arg));
        }
    }
    std::sort(cut_points.begin(), cut_points.end());

    if (cut_points.empty()) {
        std::cerr << "Error: at least one cut point is required.\n";
        return 1;
    }

    std::cout << "RAW input  : " << input_raw << "\n"
              << "Base       : " << output_base << "\n"
              << "Cut points :";
    for (size_t c : cut_points) std::cout << " " << c;
    std::cout << "\n";
    if (has_root) std::cout << "ROOT input : " << input_root << "\n";
    std::cout << "\n";

    // ── Open ROOT input file ───────────────────────────────────────────────
    TFile    *fin_root       = nullptr;
    TTree    *tin_root       = nullptr;
    Long64_t  n_root_entries = 0;

    if (has_root) {
        fin_root = TFile::Open(input_root.c_str(), "READ");
        if (!fin_root || fin_root->IsZombie()) {
            std::cerr << "Cannot open ROOT file: " << input_root << "\n";
            return 1;
        }
        tin_root = (TTree*)fin_root->Get("wfm");
        if (!tin_root) {
            std::cerr << "Tree 'wfm' not found in " << input_root << "\n";
            return 1;
        }
        n_root_entries = tin_root->GetEntries();
        std::cout << "ROOT entries : " << n_root_entries << "\n\n";
    }

    // ── EUDAQ2 deserializer ────────────────────────────────────────────────
    eudaq::FileDeserializer des(input_raw, false);
    eudaq::EventSP bore_ev;

    size_t   total_data_events = 0;
    size_t   file_idx          = 0;
    Long64_t root_entry        = 0;

    std::unique_ptr<eudaq::FileSerializer> ser;
    TFile *fout_root = nullptr;
    TTree *tout_root = nullptr;

    auto make_suffix = [](size_t idx) -> std::string {
        char buf[16];
        snprintf(buf, sizeof(buf), "_%06zu", idx);
        return buf;
    };

    auto open_next_chunk = [&]() {
        const std::string sfx      = make_suffix(file_idx++);
        const std::string raw_out  = output_base + sfx + ".raw";
        std::cout << "  → " << raw_out << "\n" << std::flush;
        ser = std::make_unique<eudaq::FileSerializer>(raw_out, /*overwrite=*/true);
        if (bore_ev) ser->write(*bore_ev);

        if (has_root) {
            const std::string root_out = output_base + sfx + ".root";
            std::cout << "  → " << root_out << "\n" << std::flush;
            fout_root = TFile::Open(root_out.c_str(), "RECREATE");
            fout_root->cd();
            tout_root = tin_root->CloneTree(0);
            tout_root->SetAutoSave(0);
        }
    };

    auto close_root_chunk = [&]() {
        if (fout_root) {
            fout_root->cd();
            tout_root->Write("", TObject::kOverwrite);
            fout_root->Close();
            delete fout_root;
            fout_root = nullptr;
            tout_root = nullptr;
        }
    };

    // ── Main loop ──────────────────────────────────────────────────────────
    while (des.HasData()) {
        uint32_t id;
        des.PreRead(id);

        eudaq::EventUP ev =
            eudaq::Factory<eudaq::Event>::Create<eudaq::Deserializer&>(id, des);
        if (!ev) {
            std::cerr << "Warning: null event (type=0x" << std::hex << id
                      << std::dec << "), stopping.\n";
            break;
        }

        if (ev->IsBORE()) {
            bore_ev = std::move(ev);
            open_next_chunk();
            continue;
        }

        if (ev->IsEORE()) {
            if (ser) { ser->write(*ev); ser.reset(); }
            close_root_chunk();
            break;
        }

        // Data event
        if (!ser) open_next_chunk();

        ser->write(*ev);

        if (has_root && root_entry < n_root_entries) {
            tin_root->GetEntry(root_entry++);
            tout_root->Fill();
        }

        ++total_data_events;

        if (total_data_events % 100000 == 0)
            std::cout << "  " << total_data_events << " events processed\n" << std::flush;

        // Split if we've reached the next cut point
        if (!cut_points.empty() && total_data_events >= cut_points.front()) {
            cut_points.pop_front();
            ser.reset();
            close_root_chunk();
            open_next_chunk();
        }
    }

    if (ser) ser.reset();
    close_root_chunk();
    if (fin_root) { fin_root->Close(); delete fin_root; }

    std::cout << "\nDone.\n"
              << "  Total data events  : " << total_data_events << "\n"
              << "  Output chunks      : " << file_idx << "\n";
    if (has_root)
        std::cout << "  ROOT entries copied: " << root_entry << "\n";
    return 0;
}
