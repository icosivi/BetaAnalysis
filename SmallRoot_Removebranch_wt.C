#include "TFile.h"
#include "TTree.h"

void create_subset_root_file(int max_events = 200000) {
    // Open the input ROOT file
    TFile *infile = TFile::Open("stats_DC_RSD_A1_FNAL_DC1_275V_.root", "READ");
    if (!infile || infile->IsZombie()) {
        printf("Error: Cannot open input file\n");
        return;
    }

    // Get the TTree from the file
    TTree *intree = (TTree*)infile->Get("Analysis"); // Replace "TreeName" with your tree's name
    if (!intree) {
        printf("Error: Cannot find TTree in file\n");
        return;
    }

    // Create the output ROOT file
    TFile *outfile = TFile::Open("stats_DC_RSD_A1_FNAL_DC1_275V_small.root", "RECREATE");
    if (!outfile || outfile->IsZombie()) {
        printf("Error: Cannot create output file\n");
        return;
    }

    // Clone the tree structure but not the data
    TTree *outtree = intree->CloneTree(0);

    // Copy the first max_events entries
    int nentries = std::min(max_events, (int)intree->GetEntries());
    for (int i = nentries; i < nentries+100; ++i) {
        intree->GetEntry(i);
        outtree->Fill();
    }

    // Write the new tree and close the files
    outfile->Write();
    outfile->Close();
    infile->Close();

    //  printf("Subset file with %d events saved to %s\n", nentries);
}



void remove_branch() {
    // Open the input ROOT file
  const char* remove_branch_name = "t";
  const char* remove_branch_name2 = "w";
  TFile *infile = TFile::Open("TB7_run_1960/stats_small_1960.root", "READ");
  TFile *outfile = TFile::Open("nicolo.root", "RECREATE");
  //    TFile *infile = TFile::Open(input_file, "READ");
  if (!infile || infile->IsZombie()) {
    printf("Error: Cannot open input file\n");
    return;
    }

    // Get the TTree from the file
    TTree *intree = (TTree*)infile->Get("Analysis"); // Replace "TreeName" with your tree's name
    if (!intree) {
        printf("Error: Cannot find TTree in file\n");
        infile->Close();
        return;
    }

    // Create the output ROOT file
    //   TFile *outfile = TFile::Open(output_file, "RECREATE");
    if (!outfile || outfile->IsZombie()) {
        printf("Error: Cannot create output file\n");
        infile->Close();
        return;
    }

    // Disable the branch to be removed
    intree->SetBranchStatus(remove_branch_name, 0);
    intree->SetBranchStatus(remove_branch_name2, 0);

    // Clone the tree structure excluding the branch
    TTree *outtree = intree->CloneTree(-1);

    // Write the new tree and close the files
    outtree->Write();
    outfile->Close();
    infile->Close();

    //   printf("Branch '%s' removed and file saved to %s\n", remove_branch_name, output_file);
}
