#pragma once
#include "call_libraries.h"  // call libraries from ROOT and C++
#include "function_definition.h" // function definition
#include <algorithm>
#include <iostream>
#include <random>
#define coscutmix 0.99996
#define dptcutmix 0.04

void MixEvents(int centrality_or_ntrkoff_int, int nEvt_to_mix, std::vector<int> ev_multiplicity, std::vector<double> vtx_z_vec, std::vector<int> ev_phpos, double vzcut, std::vector<std::vector<ROOT::Math::PtEtaPhiMVector>> Track_Vector, std::vector<std::vector<int>> Track_Chg_Vector, std::vector<std::vector<double>> Track_Eff_Vector, THnSparseD *histo_SS, THnSparseD *histo_SSLCMS, THnSparseD *histo_SS3D, THnSparseD *histo_OS, THnSparseD *histo_OSLCMS, THnSparseD *histo_OS3D, bool docostdptcut, bool do_hbt3d, bool dogamovcorrection, int systematic, TH1I* NeventsAss){

    const int aux_n_evts = static_cast<int>(vtx_z_vec.size());
	
	cout << "Running ... " << endl;
	cout << "Total # of events: " << aux_n_evts << endl; 
   
	std::mt19937 rng(12345); // fixed seed for reproducibility; replace with std::random_device{}() for true randomness
   
	for(int nevt_trg = 0; nevt_trg < aux_n_evts; nevt_trg++){ // first loop over all events

		if (nevt_trg != 0 && (nevt_trg % (aux_n_evts / 10)) == 0) std::cout << nevt_trg << " out of " << aux_n_evts << std::endl;
		auto& Trk_nevt_trg_vec = Track_Vector[nevt_trg]; // track 4-vector for each trigger event
		auto& Trk_chg_nevt_trg_vec = Track_Chg_Vector[nevt_trg]; // track charge vector for each trigger event
		auto& Trk_eff_nevt_trg_vec = Track_Eff_Vector[nevt_trg]; // track efficiency vector for each trigger event
		const int nMix_nevt_trg = static_cast<int>(Trk_nevt_trg_vec.size());
		if (nMix_nevt_trg == 0) continue;
		int mult_trg = ev_multiplicity[nevt_trg];
		int ev_phside = ev_phpos[nevt_trg];
		double vz_trg = vtx_z_vec[nevt_trg];
        // do not modify original vzcut; use a local copy
        double vzcut_local = vzcut;
        if (std::fabs(vz_trg) > 7.0 && std::fabs(vz_trg) < 10.0) vzcut_local = 2.0 * vzcut_local;
        if (std::fabs(vz_trg) > 10.0) vzcut_local = 5.0 * vzcut_local;
		// Build list of candidate events with similar multiplicity
		std::vector<int> assocCandidates;
		int totalevents = 0;
		for (int iev = 0; iev < aux_n_evts; iev++) {
			if (iev == nevt_trg) continue;
			if ( fabs(ev_multiplicity[iev] - mult_trg) > centrality_or_ntrkoff_int ) continue; 
		    if ( fabs(vtx_z_vec[iev] - vz_trg) > vzcut_local ) continue;
			if ( ev_phpos[iev] != ev_phside) ) continue;
			totalevents = totalevents + 1;
			if( totalevents > int(1000 * nEvt_to_mix) ) break;
			assocCandidates.push_back(iev);
		}
		// If no candidates, skip
		if (assocCandidates.empty()) continue;
		// Shuffle candidates
		std::shuffle(assocCandidates.begin(), assocCandidates.end(), rng);
		// Take first NEventtoMix candidates
		int nAssocToUse = std::min(nEvt_to_mix, (int)assocCandidates.size());
		NeventsAss->Fill(nAssocToUse);   
		for (int iass = 0; iass < nAssocToUse; iass++) {
			int nevt_assoc = assocCandidates[iass];
			auto& Track_nevt_ass_vec = Track_Vector[nevt_assoc];
			auto& Trk_chg_nevt_ass_vec = Track_Chg_Vector[nevt_assoc]; // track charge vector for each associate event
			auto& Trk_eff_nevt_ass_vec = Track_Eff_Vector[nevt_assoc]; // track efficiency vector for each associate event

			int nMix_nevt_ass = (int)Track_nevt_ass_vec.size();
			if (nMix_nevt_ass == 0) continue;
			// loop and fill correlation histograms
			for (int imix = 0; imix < nMix_nevt_trg; imix++) {
				for (int iimix = 0; iimix < nMix_nevt_ass; iimix++) {    	    
					double eff_trk_imix = Trk_eff_nevt_trg_vec[imix];
					double eff_trk_iimix = Trk_eff_nevt_ass_vec[iimix];
					double tot_eff = eff_trk_imix * eff_trk_iimix;
					if(docostdptcut){ if(splitcomb(Trk_nevt_trg_vec[imix],Track_nevt_ass_vec[iimix],coscutmix,dptcutmix)) continue; }
					ROOT::Math::PtEtaPhiMVector psum2 = Trk_nevt_trg_vec[imix] + Track_nevt_ass_vec[iimix];
					double kt = (psum2.Pt())/2.;
					double qinv = GetQ(Trk_nevt_trg_vec[imix],Track_nevt_ass_vec[iimix]);
					double qlcms = GetQLCMS(Trk_nevt_trg_vec[imix],Track_nevt_ass_vec[iimix]);
					double qlong = GetQlongLCMS(Trk_nevt_trg_vec[imix],Track_nevt_ass_vec[iimix]);
					double qout = GetQout(Trk_nevt_trg_vec[imix],Track_nevt_ass_vec[iimix]);
					double qside = GetQside(Trk_nevt_trg_vec[imix],Track_nevt_ass_vec[iimix]);
					double x_2pc_hbt[4]={qinv, kt, (double)ev_multiplicity[nevt_trg], (double) ev_phpos[nevt_trg]}; 
					double x_2pc_hbt_lcms[4]={qlcms, kt, (double)ev_multiplicity[nevt_trg], (double) ev_phpos[nevt_trg]};
					double x_2pc_hbt_3D[6]={qlong, qout, qside, kt, (double)ev_multiplicity[nevt_trg], (double) ev_phpos[nevt_trg]}; 
					double coulomb_ss = 1.0;
					double coulomb_os = 1.0;
					if(dogamovcorrection){ coulomb_ss = CoulombSS(qinv,systematic); coulomb_os = CoulombOS(qinv,systematic); }
					if(Trk_chg_nevt_trg_vec[imix]*Trk_chg_nevt_ass_vec[iimix] > 0){
						histo_SS->Fill(x_2pc_hbt,coulomb_ss*tot_eff);
						histo_SSLCMS->Fill(x_2pc_hbt_lcms,coulomb_ss*tot_eff);
						if(do_hbt3d) histo_SS3D->Fill(x_2pc_hbt_3D,coulomb_ss*tot_eff);
					}else{
						histo_OS->Fill(x_2pc_hbt,coulomb_os*tot_eff);			
						histo_OSLCMS->Fill(x_2pc_hbt_lcms,coulomb_os*tot_eff);
						if(do_hbt3d) histo_OS3D->Fill(x_2pc_hbt_3D,coulomb_os*tot_eff);			
					}
				} // end of iimix loop						
			} // end of imix loop
		} // end of associate loop
	} // end of trigger events loop
} // end of function
