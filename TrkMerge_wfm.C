void TrkMerge_wfm(TString originalFileName = "/media/tb_pc/320935fb-07d3-4c77-8da1-ced50467a779/DC_RSD/TB8/raw/W3_C44_DC19_220V_Run304.root", TString treeName= "wfm",TString  newFileName= "/media/tb_pc/320935fb-07d3-4c77-8da1-ced50467a779/DC_RSD/TB8/raw/with_tracks/W3_C44_DC19_220V_Run304_Tracks_DUT.root", TString  TrkFile= "/home/tb_pc/Desktop/TestBeam/Tracking/data8/data/tracks_304_withDUT.csv") {

  #include <iostream>
  double fxtrk1, fytrk1, fxtrk2, fytrk2, fchi2trk1;
  int fntrk, pippo;
  double xtrkm1, ytrkm1, xtrkm2, ytrkm2;
  int NShift, Skip;
  //  Run31_160V_gigi2_trackerN31.root

  TTreeReader     fReader;  //!the tree reader

    // Open the original ROOT file
    TFile *file = new TFile(originalFileName, "READ");

    // Get the TTree from the original ROOT file
    TTree *tree = (TTree*) file->Get(treeName);


    // Create a new ROOT file to write the TTree to
    TFile *newFile = new TFile(newFileName, "RECREATE");

    // Create a new TTree in the new ROOT file with the same name and structure as the original TTree
    TTree *newTree = tree->CloneTree(0);

    Long64_t nEntries = tree->GetEntries();
    fReader.SetTree(tree);

    std::cout << "Total number of events = " << tree->GetEntries() << "\n";

    float xtrk1;
    float xtrk2;
    float ytrk1;
    float ytrk2;
    float chi2trk;
    TBranch *branch0 = newTree->Branch("xtrk1", &xtrk1, "xtrk1/F");
    TBranch *branch1 = newTree->Branch("xtrk2", &xtrk2, "xtrk2/F");
    TBranch *branch2 = newTree->Branch("ytrk1", &ytrk1, "ytrk1/F");
    TBranch *branch3 = newTree->Branch("ytrk2", &ytrk2, "ytrk2/F");
    TBranch *branch4 = newTree->Branch("chi2trk", &chi2trk, "chi2trk/F");


    std::ifstream Filein;

    Filein.open(TrkFile, std::ios::in);

    if(!Filein.is_open()) std::cout << "It failed" << std::endl;

    else std::cout << "Opened file " << TrkFile << std::endl;

    
    //    TTreeReaderValue<Double_t> x_pos = {fReader, "x_pos"};
    //  TTreeReaderValue<Double_t> y_pos = {fReader, "y_pos"};
    

    std::cout << "Total number of events = " << tree->GetEntries() << "\n";
    // Loop over the events in the original TTree and fill the new TTree in the new file with each event

    //    Filein >>  fntrk >> fxtrk1 >> fytrk1 >> fchi2trk1 >> fxtrk2 >> fytrk2 >> fchi2trk1  ;   

    // needs to start from 0, otherwise it is offsink with the gigi2 file. 
    for (Long64_t i=0; i<tree->GetEntries(); i++) 
    //for (Long64_t i=1; i<1000; i++)
      {	    
	tree->GetEntry(i);
	fReader.SetLocalEntry(i);
	// skip one events to be in sink with the tracker file
	  
	if (i>0) 
	  //Filein >>  fntrk >> fxtrk1 >> fytrk1 >> fchi2trk1 >> fxtrk2 >> fytrk2 >> fchi2trk1  ; 
    Filein >>  fntrk >> fxtrk1 >> fytrk1 >> fxtrk2 >> fytrk2 >> fchi2trk1 ;   

	if (fntrk % 10000 == 0)
	    {
	      std::cout << "Event = " << i-1 << " " << fntrk << "\n";	   
	    }
	  xtrk1 = fxtrk1;
	  ytrk1 = fytrk1;
	  xtrk2 = fxtrk2;
	  ytrk2 = fytrk2;
	  chi2trk = fchi2trk1;
	  newTree->Fill();

      }
    
	
	// Write the changes to the file and close it
    newFile->cd();
    newTree->Write();
    newFile->Close();
    
    return;
    
}
