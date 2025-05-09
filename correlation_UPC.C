#include "read_tree.h" // read the TChains
#include "tracking_correction.h" // tracking correction
#include "define_histograms_mult.h" // histogram definition
#define pi_mass 0.1396

/*
--> Arguments
input_file: text file with a list of root input files: Forest or Skims
ouputfile: just a counting number to run on Condor
isMC: 0 for false --> data and > 0 for true --> MC
doquicktest: 0 for false and > 0 for true --> tests with 1000 events
domixing: 0 with mixing and 1 without mixing
Nmixevents: number of events to mix
mincentormult: minimum centrality of multiplicity in the mixing
minvz: minimum vertez between events in the mixing
hbt3d: 0 with 3D and 1 without 3D
gamov: 0 means no GAMOV Coulomb correction and > 0 means GAMOV is added
cent_bool: 0 to use centrality and different than 0 for multiplicity. For multiplicity we may need to change the bins in define_histograms.h
syst:systematic uncertainties
	--> 0: nominal
	--> 1: |vz| < 3
	--> 2: 3 < |vz| < 15
	--> 3: tracking tight
	--> 4: tracking loose
	--> 5: centrality up
	--> 6: centrality down
	--> 7: removal of duplication cut
	--> 8: removal of Npixel cut
	--> 9: gamov + 15% weight
	--> 10: gamov - 15% weight
	--> 11: Number of mix events from 10 to 20
	--> 12: Number of mix events from 10 to 5
	--> 13: Number of vz separation btw events from 2cm to 3cm
	--> 14: Number of vz separation btw events from 2cm to 1cm
	--> 15: Number of multiplicity/centrality separation btw events from 5 to 10 units
	--> 16: Number of multiplicity/centrality separation btw events from 5 to 3 units
*/
void correlation_UPC(TString input_file, TString ouputfile, int doquicktest, int domixing, int Nmixevents, int mincentormult, float minvz, int hbt3d, int gamov, int syst){

	clock_t sec_start, sec_end;
	sec_start = clock(); // start timing measurement
	TDatime* date = new TDatime(); // to add date in the output file

	bool do_quicktest; if(doquicktest == 0){do_quicktest = false;}else{do_quicktest = true;}
	bool do_mixing; if(domixing == 0){do_mixing = true;}else{do_mixing = false;}
	bool do_hbt3d; if(hbt3d == 0){do_hbt3d = true;}else{do_hbt3d = false;}
	bool do_gamov; if(gamov == 0){do_gamov = false;}else{do_gamov = true;}
	
	if(syst == 4 || syst == 5) do_gamov = true;
	
	bool dosplit = false;
	if(syst != 3) dosplit = true; 
	
	TString systematics = "nonapplied_nominal";
	if(syst == 0) systematics =  "nominal";
	if(syst == 1) systematics =  "vznarrow";
	if(syst == 2) systematics =  "vzwide";
	if(syst == 3) systematics =  "removeduplicatedcut";
	if(syst == 4) systematics =  "gamovplus15";
	if(syst == 5) systematics =  "gamovminus15";
	if(syst == 6){ systematics =  "Nmix20"; Nmixevents = Nmixevents + 10; }
	if(syst == 7){ systematics =  "Nmix05"; Nmixevents = Nmixevents - 5; }
	if(syst == 8){ systematics =  "minvz3"; minvz = minvz + 1.0; }
	if(syst == 9){ systematics =  "minvz1"; minvz = minvz - 1.0; }
	if(syst == 10){ systematics =  "mincentormult10"; mincentormult = mincentormult + 5; }
	if(syst == 11){ systematics =  "mincentormult3"; mincentormult = mincentormult - 2; }
	if(syst == 12){ systematics =  "hfthres4"; }
	if(syst == 13){ systematics =  "hfthres6"; }
	if(syst == 14){ systematics =  "hfthres8"; }
	if(syst == 15){ systematics =  "lowpfthres"; }
	if(syst == 16){ systematics =  "highpfthres"; }
	
	// Track or particle efficiency file
	/*
	TFile *fileeff;
	fileeff = TFile::Open("efftables/XeXe_eff_table_94x_cent.root");
	if(syst != 1 && syst != 2 && syst != 3 && syst != 4) fileeff = TFile::Open("efftables/XeXe_eff_table_94x_cent.root");
	if(syst == 1) fileeff = TFile::Open("efftables/XeXe_eff_narrow_table_94x_cent.root");
	if(syst == 2) fileeff = TFile::Open("efftables/XeXe_eff_wide_table_94x_cent.root");
	if(syst == 3) fileeff = TFile::Open("efftables/XeXe_eff_tight_table_94x_cent.root");
	if(syst == 4) fileeff = TFile::Open("efftables/XeXe_eff_loose_table_94x_cent.root");
	
	TH2 *trkeff_file010 = nullptr; 
	fileeff->GetObject("rTotalEff3D_0_10", trkeff_file010);
	TH2 *trkeff_file1030 = nullptr; 
	fileeff->GetObject("rTotalEff3D_10_30", trkeff_file1030);
	TH2 *trkeff_file3050 = nullptr; 
	fileeff->GetObject("rTotalEff3D_30_50", trkeff_file3050);
	TH2 *trkeff_file5070 = nullptr; 
	fileeff->GetObject("rTotalEff3D_50_70", trkeff_file5070);
	TH2 *trkeff_file70100 = nullptr; 
	fileeff->GetObject("rTotalEff3D_70_100", trkeff_file70100);	
	*/
	// Read the list of input file(s)
	fstream inputfile;
	inputfile.open(Form("%s",input_file.Data()), ios::in);
	if(!inputfile.is_open()){cout << "List of input files not founded!" << endl; return;}{cout << "List of input files founded! --> " << input_file.Data() << endl;}

	// Make a chain and a vector of file names
	std::vector<TString> file_name_vector;
	string file_chain;
	while(getline(inputfile, file_chain)){file_name_vector.push_back(file_chain.c_str());}
	inputfile.close();

	// Read the trees to be added in the Chain
	TChain *hea_tree = new TChain("analyzer/eventTree"); // event quantities
	// loop to add all the trees to the chain
	for (std::vector<TString>::iterator listIterator = file_name_vector.begin(); listIterator != file_name_vector.end(); listIterator++){
        TFile *testfile = TFile::Open(*listIterator,"READ");
		if(testfile && !testfile->IsZombie() && !testfile->TestBit(TFile::kRecovered)){ // safety against corrupted files
			cout << "Adding file " << *listIterator << " to the chains" << endl; // adding files to the chains for each step
			hea_tree->Add(*listIterator);
		}else{cout << "File: " << *listIterator << " failed!" << endl;}		
	}
	file_name_vector.clear();	
	
    // Read the desired branchs in the trees
	read_tree(hea_tree); // access the tree informations
	
    // Use sumw2() to make sure about histogram uncertainties in ROOT
	sw2(); 

	int nevents = hea_tree->GetEntries(); // number of events
	cout << endl;
	cout << "Total number of events in those files: "<< nevents << endl;
	cout << endl;
	cout << "-------------------------------------------------" << endl;

	//define vectors for mixing
	std::vector<int> multiplicity_vector;
	std::vector<double> vz_vector;
	std::vector<int> phpos_vector;
	std::vector<std::vector<ROOT::Math::PtEtaPhiMVector>> track_4vector;
	std::vector<std::vector<double>> track_weights_vector;
	std::vector<std::vector<int>> track_charge_vector;


	double nev = (double)nevents;
	
	for(int i = 0; i < nevents; i++){

		hea_tree->GetEntry(i); // get events from ttree

		if(i != 0 && (i % 100000) == 0){double alpha = (double)i; cout << " Running -> percentage: " << std::setprecision(3) << ((alpha / nev) * 100) << "%" << endl;} // % processed
		if(do_quicktest){if(i != 0 && i % 1000 == 0 ) break;} // just for quick tests

		int Ntroff = mult;
		double vertexz = (double) vtx[2];

		// fill quantities before any cut
		vzhist_beforecuts->Fill(vertexz);
		multiplicity_beforecuts->Fill(Ntroff);
				
		Nevents->Fill(0); // Filled after each event cut -> 0 means no cuts	
		// vertex selection 
		
		if(syst == 1){ if(fabs(vertexz) > 3.0) continue; 
		} else if(syst == 2){ if( fabs(vertexz) > 15.0 || fabs(vertexz) < 3.0 ) continue; 
		} else{ if(fabs(vertexz) > 15.0) continue; }
		Nevents->Fill(1); // Vertex Z cut	

      		if( zdcSumPlus <= -50000 || zdcSumPlus > 300000 || zdcSumMinus <= -50000 || zdcSumMinus > 300000 ) continue; // clean up XZDC
      		Nevents->Fill(2); // ZDC cut
        	if( triggers[4] || triggers[5] ) { hZDCPlusZDCOr_beforecuts->Fill( zdcSumPlus ); hZDCMinusZDCOr_beforecuts->Fill( zdcSumMinus );}
        	if( triggers[6] || triggers[7] ) { hZDCPlusZeroBias_beforecuts->Fill( zdcSumPlus ); hZDCMinusZeroBias_beforecuts->Fill( zdcSumMinus ); }

		float ZDC_MINUS_CUT = 1000;
		float ZDC_PLUS_CUT = 1000;
		//If needed for systematics we can add it here using the systematic flag , e.g.
		//if( syst == XX ) {ZDC_MINUS_CUT = XXXX; ZDC_PLUS_CUT = XXXX;} 

        	bool posPhoton;
      		if( zdcSumMinus > ZDC_MINUS_CUT && zdcSumPlus < ZDC_PLUS_CUT ){ //photon traveling in + direction
	      		posPhoton = true;
      		} else if( zdcSumMinus < ZDC_MINUS_CUT && zdcSumPlus > ZDC_PLUS_CUT ){ //photon traveling in - direction
        		posPhoton = false;
        	} else { continue; }
      		Nevents->Fill(3); // Photon travel side cuts

		// HF threshold
		double HFThreshold = 0.0; // have to decide the nominal
		if(syst == 12 ) HFThreshold = 4.0;
		if(syst == 13 ) HFThreshold = 6.0;
		if(syst == 14 ) HFThreshold = 8.0;

		//Find Eta Gaps
        	double etaGapPos = 0;
        	double etaGapNeg = 0;
        	double ETSum = 0;
        	int numNegHFClusters=0, numPosHFClusters=0;
		etagaps(syst, HFThreshold, ETSum, numNegHFClusters, numPosHFClusters, pfPt, pfEta, pfE, pfID, isPrimary, etaGapPos, etaGapNeg, h_h, h_e, h_mu, h_gamma, h_h0, h_HFhad, h_HFem, posPhoton);

		if( HFThreshold > 0.0 && (numNegHFClusters < 1 || numPosHFClusters < 1) ) continue; // HF Coincidence
		Nevents->Fill(4);
		double sumGapCut = 3.5;
        	if( Ntroff > 30 ) sumGapCut = 3.0;
        	// can be also modified as systematics (if needed)
		//cout << "posPhoton: "<< (int) posPhoton << endl;
     		if( posPhoton && etaGapPos < sumGapCut ) continue;
      		if( !posPhoton && etaGapNeg < sumGapCut ) continue;
      		Nevents->Fill(5); // Multiplicity cut	

		if( Ntroff < 0 ) continue; // remove events with multiplicity < 0
		if( Ntroff > 250 ) continue; // remove events with multiplicity > 250
      		Nevents->Fill(6); // Multiplicity cut	
      	
        	if( Ntroff <= 20 && !triggers[4] && !triggers[5] ) continue;
      		Nevents->Fill(7); // Multiplicity cut	
        	if( Ntroff > 20 && !triggers[0] && !triggers[1] && !triggers[2] && !triggers[3] ) continue;
       		Nevents->Fill(8); // Multiplicity cut	

		// Fill event histograms after all cuts
		vzhist->Fill(vertexz);
		multiplicity->Fill( Ntroff );
        	if(posPhoton) multposplus->Fill( Ntroff ); 
        	if(!posPhoton) multposminus->Fill( Ntroff ); 
        	if( triggers[4] || triggers[5] ) { hZDCPlusZDCOr->Fill( zdcSumPlus ); hZDCMinusZDCOr->Fill( zdcSumMinus );}
        	if( triggers[6] || triggers[7] ) { hZDCPlusZeroBias->Fill( zdcSumPlus ); hZDCMinusZeroBias->Fill( zdcSumMinus ); }	

		// Vectors used for objects
		std::vector<ROOT::Math::PtEtaPhiMVector> tracks_reco;
		std::vector<double> track_weight_reco;
		std::vector<int> track_charge_reco;
		
		// ------------------- Reconstruction level (Data and MC) ----------------------------
		// Start loop over reco tracks (trksize is number of reco tracks)
	    	int numPF = pfID->size();
    		if( numPF == 0 ) continue;
		for (int j = 0; j < numPF; j++){  // Loop over PF Candidates

 	        	double trkEta = pfEta->at(j);
     	        	//if( posPhoton ) trkEta = -trkEta;
 	        	double trkPt = pfPt->at(j);
 	        	double trkPhi = pfPhi->at(j);
 	        	int trkCharge = 1; // for future usage 

		   	if( trkPt <= 0.0 ) continue; // remove negative and 0 pT tracks
   		   	if( pfID->at(j) != 1 ) continue; // keep only charged hadrons
		   	if( !isPrimary->at(j) ) continue; // remove non-primary tracks
		   	if( fabs( trkEta ) > 2.4 ) continue; // tracker acceptance

		    	// Track QA array
			double x_reco_trk[6]={trkPt,trkEta,trkPhi,(double) trkCharge,(double) Ntroff, (double) posPhoton}; 
			hist_reco_trk->Fill(x_reco_trk);

			// Track efficiency correction -> future usage
			double trk_weight = 1.0;
			// trk_weight = trk_weight*getTrkCorrWeight(trkeff_file, trkpt->at(j), trketa->at(j));
			// hist_reco_trk_corr->Fill(x_reco_trk,trk_weight);
			
			ROOT::Math::PtEtaPhiMVector TrackFourVector;
      			TrackFourVector.SetM(pi_mass);
      			TrackFourVector.SetEta(trkEta);
     			TrackFourVector.SetPhi(trkPhi);
     			TrackFourVector.SetPt(trkPt);  
     		
     			tracks_reco.push_back(TrackFourVector);
			track_charge_reco.push_back(trkCharge); 
			track_weight_reco.push_back(trk_weight); 

		} // End loop over tracks
		
		if(tracks_reco.size() > 1){
			Nevents->Fill(9); // filled after each event cut
			twoparticlecorrelation(tracks_reco, track_charge_reco, track_weight_reco, hist_pairSS_Mass, hist_dpt_cos_SS, hist_qinv_SS, hist_qinv_SS_INV, hist_qinv_SS_ROT, hist_q3D_SS, hist_q3D_SS_INV, hist_q3D_SS_ROT, hist_pairOS_Mass, hist_dpt_cos_OS, hist_qinv_OS, hist_qinv_OS_INV, hist_qinv_OS_ROT, hist_q3D_OS, hist_q3D_OS_INV, hist_q3D_OS_ROT, Ntroff, (int) posPhoton, dosplit, do_hbt3d, do_gamov, syst); // HBT correlations done at this step
			track_4vector.push_back(tracks_reco); // save 4 vector for mixing
			track_charge_vector.push_back(track_charge_reco); // save charge vector for mixing
			track_weights_vector.push_back(track_weight_reco); // save eff weight vector for mixing
			multiplicity_vector.push_back(Ntroff); // save multiplicity vector for mixing
			vz_vector.push_back(vertexz); // save vz vector for mixing
			phpos_vector.push_back((int) posPhoton);// save photon side vector for mixing
		}
		
		tracks_reco.clear();
		track_weight_reco.clear();
		track_charge_reco.clear();
		
	} // End of event loop
	
	// do the mixing after the event selections
	if(do_mixing) cout << "Time for mixing" << endl;
	if(do_mixing) MixEvents(mincentormult, Nmixevents, multiplicity_vector, vz_vector, phpos_vector, minvz, track_4vector, track_charge_vector, track_weights_vector, hist_qinv_SS_MIX, hist_q3D_SS_MIX, hist_qinv_OS_MIX, hist_q3D_OS_MIX, dosplit, do_hbt3d, do_gamov, syst, NeventsAss);
	
	// Output file name
	cout << endl;
	cout << "Writing histograms on ";
	cout << endl;

	// Make an output file
	string file_output = Form("%s_syst_%s_Nmix_%i_Mixint_%i_Vzint_%.f_on_%i", ouputfile.Data(), systematics.Data(), Nmixevents, mincentormult, minvz,date->GetDate()); // output file
	std::replace(file_output.begin(), file_output.end(), '.', 'p'); // replace . to p	

	// Open, write and close the output file
	TFile *MyFile = new TFile(Form("%s.root", file_output.c_str()), "RECREATE");
	if(MyFile->IsOpen()) cout << "output file: " << file_output.c_str() << ".root" << endl;
	MyFile->cd(); 

	// Write in different folders (see histogram_definition.h)
	// Control plots
	MyFile->mkdir("EventQA_histograms"); 
	MyFile->cd("EventQA_histograms"); 
	write_eventQA();

	MyFile->mkdir("HBT_1D"); 
	MyFile->cd("HBT_1D"); 
	write_HBT1D();
	if(do_hbt3d){
		MyFile->mkdir("HBT_3D"); 
		MyFile->cd("HBT_3D"); 
		write_HBT3D();
	}
	

	MyFile->Close();

	cout << endl;
	cout << "------------------------------------- DONE --------------------------------------" << endl;
	cout << endl;


	sec_end = clock(); // stop time counting
	cout << "========================================" << endl;
	cout << "Total running time: " << (double)(sec_end - sec_start) / CLOCKS_PER_SEC << " [s]" << endl;
	cout << "========================================" << endl;

}
