#include "call_libraries.h"  // call libraries from ROOT and C++
#include "random_mixing.h" // random mixing

// define the bins
// qinv
const int nQBins = 400;   // number of qinv bins
const double minQ = 0.0;  // minimum qinv
const double maxQ = 4.0;  // maximumm qinv

// q3D
const int nQBins3D = 100;   // number of q3D bins
const double minQ3D = 0.0;  // minimum q3D
const double maxQ3D = 2.0;  // maximumm q3D

//kT
const int nKtBins = 6; // number of average transverse momentum bins
double KtBins[nKtBins+1] = {0.1,0.3,0.4,0.5,0.6,0.7,1.0}; 

// multiplicity
//const int nCentBins = 5; // number of multiplicity bins
//double CentBins[nCentBins+1] = {0.0, 20.0, 60.0, 100.0, 140.0, 200.0}; // 0-10%, 10-30%, 30-50%, 50-70%, 70-90% (0.5 in 0.5%)

// multiplicity
const int nCentBins = 11; // number of multiplicity bins
double CentBins[nCentBins+1] = {0.0, 5.0, 10.0, 20.0, 30.0, 40.0, 50.0, 70.0, 100.0, 150.0, 200.0, 250.0};  // multiplicity bins to match PbPb

// Event histograms
TH1I *Nevents = new TH1I("Nevents", "Nevents", 10, 0, 10);
TH1D *vzhist_beforecuts = new TH1D("vzhist_beforecuts", "vzhist_beforecuts", 80, -20., 20.);
TH1D *vzhist = new TH1D("vzhist", "vzhist", 80, -20., 20.);
TH1D *multiplicity_beforecuts = new TH1D("multiplicity_beforecuts", "multiplicity_beforecuts", 400, 0.0, 400.0);
TH1D *multiplicity = new TH1D("multiplicity", "multiplicity", 400, 0.0, 400.0);
TH1D *multposplus = new TH1D("multposplus", "multposplus", 400, 0.0, 400.0);
TH1D *multposminus = new TH1D("multposminus", "multposminus", 400, 0.0, 400.0);
TH1I *NeventsAss = new TH1I("NeventsAss", "NeventsAss", 40, 0, 40);

// ZDC
TH1D *hZDCPlusZeroBias_beforecuts  = new TH1D( "hZDCPlusZeroBias_beforecuts",  "", 150, 0, 15000 );
TH1D *hZDCMinusZeroBias_beforecuts = new TH1D( "hZDCMinusZeroBias_beforecuts", "", 150, 0, 15000 );
TH1D *hZDCPlusZDCOr_beforecuts  = new TH1D( "hZDCPlusZDCOr_beforecuts", "", 150, 0, 15000 );
TH1D *hZDCMinusZDCOr_beforecuts = new TH1D( "hZDCMinusZDCOr_beforecuts", "", 150, 0, 15000 );
TH1D *hZDCPlusZeroBias  = new TH1D( "hZDCPlusZeroBias",  "Zero Bias Trigger ZDC Plus Distribution;ZDC Offline Energy Sum (GeV);", 150, 0, 15000 );
TH1D *hZDCMinusZeroBias = new TH1D( "hZDCMinusZeroBias", "Zero Bias Trigger ZDC Minus Distribution;ZDC Offline Energy Sum (GeV);", 150, 0, 15000 );
TH1D *hZDCPlusZDCOr  = new TH1D( "hZDCPlusZDCOr", "ZDCOr Trigger ZDC Plus Distribution;ZDC Offline Energy Sum (GeV);", 150, 0, 15000 );
TH1D *hZDCMinusZDCOr = new TH1D( "hZDCMinusZDCOr", "ZDCOr Trigger ZDC Minus Distribution;ZDC Offline Energy Sum (GeV);", 150, 0, 15000 );

// PF Candidates
TH2D *h_h      = new TH2D( "h_h",     "Charged Hadrons; #eta; p_{T}", 1100, -5.5, 5.5, 500, 0, 50 );
TH2D *h_e      = new TH2D( "h_e",     "Electrons; #eta; p_{T}", 1100, -5.5, 5.5, 500, 0, 50 );
TH2D *h_mu     = new TH2D( "h_mu",    "Muons; #eta; p_{T}", 1100, -5.5, 5.5, 500, 0, 50 );
TH2D *h_gamma  = new TH2D( "h_gamma", "Photons; #eta; E_{T}", 1100, -5.5, 5.5, 500, 0, 50 );
TH2D *h_h0     = new TH2D( "h_h0",    "Neutral Hadrons; #eta; E_{T}", 1100, -5.5, 5.5, 500, 0, 50 );
TH2D *h_HFhad  = new TH2D( "h_HFhad", "HF Hadronic; #eta; E", 1100, -5.5, 5.5, 500, 0, 50 );
TH2D *h_HFem   = new TH2D( "h_HFem",  "HF Electromagnetic; #eta; E", 1100, -5.5, 5.5, 500, 0, 50 );

// Track/Particle histograms
// Axis : 0 -> track pT, 1 -> trk eta, 2 -> trk phi, 3 -> trk charge, 4 -> multiplicity bin, 5 -> photon side
int	bins_trk[6]      =   { 200   ,  24  ,   30				      , 3   , nCentBins  			, 3 };
double xmin_trk[6]   =   { 0.0   , -2.4 ,   -TMath::Pi()  		  , -1.5, CentBins[0]			, -1.5};
double xmax_trk[6]   =   { 50.0  ,  2.4 ,   TMath::Pi()  		  ,  1.5, CentBins[nCentBins]	, 1.5};
// --> Reco
THnSparseD *hist_reco_trk = new THnSparseD("hist_reco_trk", "hist_reco_trk", 6, bins_trk, xmin_trk, xmax_trk);


// HBT histograms for tests
// Axis : 0 -> qinv, 1 -> kT, 2 -> multiplicity bin
int	bins_qinv[4]      =   { nQBins , nKtBins  		   ,   nCentBins, 3};
double xmin_qinv[4]   =   { minQ   , KtBins[0] 		   ,   CentBins[0], -1.5};
double xmax_qinv[4]   =   { maxQ   , KtBins[nKtBins] ,   CentBins[nCentBins], 1.5};
THnSparseD *hist_qinv_SS = new THnSparseD("hist_qinv_SS", "hist_qinv_SS", 4, bins_qinv, xmin_qinv, xmax_qinv);
THnSparseD *hist_qinv_SS_INV = new THnSparseD("hist_qinv_SS_INV", "hist_qinv_SS_INV", 4, bins_qinv, xmin_qinv, xmax_qinv);
THnSparseD *hist_qinv_SS_ROT = new THnSparseD("hist_qinv_SS_ROT", "hist_qinv_SS_ROT", 4, bins_qinv, xmin_qinv, xmax_qinv);
THnSparseD *hist_qinv_SS_MIX = new THnSparseD("hist_qinv_SS_MIX", "hist_qinv_SS_MIX", 4, bins_qinv, xmin_qinv, xmax_qinv);
THnSparseD *hist_qinv_OS = new THnSparseD("hist_qinv_OS", "hist_qinv_OS", 4, bins_qinv, xmin_qinv, xmax_qinv);
THnSparseD *hist_qinv_OS_INV = new THnSparseD("hist_qinv_OS_INV", "hist_qinv_OS_INV", 4, bins_qinv, xmin_qinv, xmax_qinv);
THnSparseD *hist_qinv_OS_ROT = new THnSparseD("hist_qinv_OS_ROT", "hist_qinv_OS_ROT", 4, bins_qinv, xmin_qinv, xmax_qinv);
THnSparseD *hist_qinv_OS_MIX = new THnSparseD("hist_qinv_OS_MIX", "hist_qinv_OS_MIX", 4, bins_qinv, xmin_qinv, xmax_qinv);
THnSparseD *hist_qinv_SS_gen = new THnSparseD("hist_qinv_SS_gen", "hist_qinv_SS_gen", 4, bins_qinv, xmin_qinv, xmax_qinv);
THnSparseD *hist_qinv_SS_gen_INV = new THnSparseD("hist_qinv_SS_gen_INV", "hist_qinv_SS_gen_INV", 4, bins_qinv, xmin_qinv, xmax_qinv);
THnSparseD *hist_qinv_SS_gen_ROT = new THnSparseD("hist_qinv_SS_gen_ROT", "hist_qinv_SS_gen_ROT", 4, bins_qinv, xmin_qinv, xmax_qinv);
THnSparseD *hist_qinv_SS_gen_MIX = new THnSparseD("hist_qinv_SS_gen_MIX", "hist_qinv_SS_gen_MIX", 4, bins_qinv, xmin_qinv, xmax_qinv);
THnSparseD *hist_qinv_OS_gen = new THnSparseD("hist_qinv_OS_gen", "hist_qinv_OS_gen", 4, bins_qinv, xmin_qinv, xmax_qinv);
THnSparseD *hist_qinv_OS_gen_INV = new THnSparseD("hist_qinv_OS_gen_INV", "hist_qinv_OS_gen_INV", 4, bins_qinv, xmin_qinv, xmax_qinv);
THnSparseD *hist_qinv_OS_gen_ROT = new THnSparseD("hist_qinv_OS_gen_ROT", "hist_qinv_OS_gen_ROT", 4, bins_qinv, xmin_qinv, xmax_qinv);
THnSparseD *hist_qinv_OS_gen_MIX = new THnSparseD("hist_qinv_OS_gen_MIX", "hist_qinv_OS_gen_MIX", 4, bins_qinv, xmin_qinv, xmax_qinv);

// HBT histograms for tests
// Axis : 0 -> qlong, 1 -> qout, 2 -> qside, 3 -> kT, 4 -> multiplicity bin
int	bins_q3D[6]      =   { nQBins3D, nQBins3D, nQBins3D,  nKtBins 		      , nCentBins            , 3};
double xmin_q3D[6]   =   { minQ3D  , minQ3D  , minQ3D  ,  KtBins[0] 		  , CentBins[0]          , -1.5};
double xmax_q3D[6]   =   { maxQ3D  , maxQ3D  , maxQ3D  ,  KtBins[nKtBins]     , CentBins[nCentBins]  , 1.5};
THnSparseD *hist_q3D_SS = new THnSparseD("hist_q3D_SS", "hist_q3D_SS", 6, bins_q3D, xmin_q3D, xmax_q3D);
THnSparseD *hist_q3D_SS_INV = new THnSparseD("hist_q3D_SS_INV", "hist_q3D_SS_INV", 6, bins_q3D, xmin_q3D, xmax_q3D);
THnSparseD *hist_q3D_SS_ROT = new THnSparseD("hist_q3D_SS_ROT", "hist_q3D_SS_ROT", 6, bins_q3D, xmin_q3D, xmax_q3D);
THnSparseD *hist_q3D_SS_MIX = new THnSparseD("hist_q3D_SS_MIX", "hist_q3D_SS_MIX", 6, bins_q3D, xmin_q3D, xmax_q3D);
THnSparseD *hist_q3D_OS = new THnSparseD("hist_q3D_OS", "hist_q3D_OS", 6, bins_q3D, xmin_q3D, xmax_q3D);
THnSparseD *hist_q3D_OS_INV = new THnSparseD("hist_q3D_OS_INV", "hist_q3D_OS_INV", 6, bins_q3D, xmin_q3D, xmax_q3D);
THnSparseD *hist_q3D_OS_ROT = new THnSparseD("hist_q3D_OS_ROT", "hist_q3D_OS_ROT", 6, bins_q3D, xmin_q3D, xmax_q3D);
THnSparseD *hist_q3D_OS_MIX = new THnSparseD("hist_q3D_OS_MIX", "hist_q3D_OS_MIX", 6, bins_q3D, xmin_q3D, xmax_q3D);

TH2D *hist_dpt_cos_SS = new TH2D("hist_dpt_cos_SS", "hist_dpt_cos_SS",1000, 0.99910, 1.0001, 100, 0, 0.25);
TH2D *hist_dpt_cos_OS = new TH2D("hist_dpt_cos_OS", "hist_dpt_cos_OS",1000, 0.99910, 1.0001, 100, 0, 0.25);
TH1D *hist_pairSS_Mass = new TH1D("hist_pairSS_Mass", "Invariant mass same-sign tracks", 500, 0, 1.0);
TH1D *hist_pairOS_Mass = new TH1D("hist_pairOS_Mass", "Invariant mass opposite-sign tracks", 500, 0, 1.0);

void sw2(){

	// Event histograms
	Nevents->Sumw2();
	vzhist_beforecuts->Sumw2();
	vzhist->Sumw2();
	multiplicity_beforecuts->Sumw2();
	multiplicity->Sumw2();
	multposplus->Sumw2();
	multposminus->Sumw2();
	NeventsAss->Sumw2();
	// ZDC
	hZDCPlusZeroBias_beforecuts->Sumw2();
	hZDCMinusZeroBias_beforecuts->Sumw2();
	hZDCPlusZDCOr_beforecuts->Sumw2();
	hZDCMinusZDCOr_beforecuts->Sumw2();
	hZDCPlusZeroBias->Sumw2();
	hZDCMinusZeroBias->Sumw2();
	hZDCPlusZDCOr->Sumw2();
	hZDCMinusZDCOr->Sumw2();
	// PF Candidates
	h_h->Sumw2();
	h_e->Sumw2();
	h_mu->Sumw2();
	h_gamma->Sumw2();
	h_h0 ->Sumw2();
	h_HFhad->Sumw2();
	h_HFem ->Sumw2();

	hist_reco_trk->Sumw2();
	hist_qinv_SS->Sumw2();
	hist_qinv_SS_INV->Sumw2();
	hist_qinv_SS_ROT->Sumw2();
	hist_qinv_SS_MIX->Sumw2();
	hist_qinv_OS->Sumw2();
	hist_qinv_OS_INV->Sumw2();
	hist_qinv_OS_ROT->Sumw2();
	hist_qinv_OS_MIX->Sumw2();
	hist_q3D_SS->Sumw2();
	hist_q3D_SS_INV->Sumw2();
	hist_q3D_SS_ROT->Sumw2();
	hist_q3D_SS_MIX->Sumw2();
	hist_q3D_OS->Sumw2();
	hist_q3D_OS_INV->Sumw2();
	hist_q3D_OS_ROT->Sumw2();
	hist_q3D_OS_MIX->Sumw2();
	hist_pairSS_Mass->Sumw2();
	hist_pairOS_Mass->Sumw2();
	
	hist_reco_trk->GetAxis(4)->Set(bins_trk[4],CentBins);	
	hist_qinv_SS->GetAxis(1)->Set(bins_qinv[1],KtBins);
	hist_qinv_SS_INV->GetAxis(1)->Set(bins_qinv[1],KtBins);
	hist_qinv_SS_ROT->GetAxis(1)->Set(bins_qinv[1],KtBins);
	hist_qinv_SS_MIX->GetAxis(1)->Set(bins_qinv[1],KtBins);
	hist_qinv_OS->GetAxis(1)->Set(bins_qinv[1],KtBins);
	hist_qinv_OS_INV->GetAxis(1)->Set(bins_qinv[1],KtBins);
	hist_qinv_OS_ROT->GetAxis(1)->Set(bins_qinv[1],KtBins);
	hist_qinv_OS_MIX->GetAxis(1)->Set(bins_qinv[1],KtBins);
	hist_qinv_SS->GetAxis(2)->Set(bins_qinv[2],CentBins);
	hist_qinv_SS_INV->GetAxis(2)->Set(bins_qinv[2],CentBins);
	hist_qinv_SS_ROT->GetAxis(2)->Set(bins_qinv[2],CentBins);
	hist_qinv_SS_MIX->GetAxis(2)->Set(bins_qinv[2],CentBins);
	hist_qinv_OS->GetAxis(2)->Set(bins_qinv[2],CentBins);
	hist_qinv_OS_INV->GetAxis(2)->Set(bins_qinv[2],CentBins);
	hist_qinv_OS_ROT->GetAxis(2)->Set(bins_qinv[2],CentBins);
	hist_qinv_OS_MIX->GetAxis(2)->Set(bins_qinv[2],CentBins);
	hist_q3D_SS->GetAxis(3)->Set(bins_q3D[3],KtBins);
	hist_q3D_SS_INV->GetAxis(3)->Set(bins_q3D[3],KtBins);
	hist_q3D_SS_ROT->GetAxis(3)->Set(bins_q3D[3],KtBins);
	hist_q3D_SS_MIX->GetAxis(3)->Set(bins_q3D[3],KtBins);
	hist_q3D_OS->GetAxis(3)->Set(bins_q3D[3],KtBins);
	hist_q3D_OS_INV->GetAxis(3)->Set(bins_q3D[3],KtBins);
	hist_q3D_OS_ROT->GetAxis(3)->Set(bins_q3D[3],KtBins);
	hist_q3D_OS_MIX->GetAxis(3)->Set(bins_q3D[3],KtBins);
	hist_q3D_SS->GetAxis(4)->Set(bins_q3D[4],CentBins);
	hist_q3D_SS_INV->GetAxis(4)->Set(bins_q3D[4],CentBins);
	hist_q3D_SS_ROT->GetAxis(4)->Set(bins_q3D[4],CentBins);
	hist_q3D_SS_MIX->GetAxis(4)->Set(bins_q3D[4],CentBins);
	hist_q3D_OS->GetAxis(4)->Set(bins_q3D[4],CentBins);
	hist_q3D_OS_INV->GetAxis(4)->Set(bins_q3D[4],CentBins);
	hist_q3D_OS_ROT->GetAxis(4)->Set(bins_q3D[4],CentBins);
	hist_q3D_OS_MIX->GetAxis(4)->Set(bins_q3D[4],CentBins);

}

// histogram writting

void write_eventQA(){
	// Event histograms
	Nevents->Write();
	vzhist_beforecuts->Write();
	vzhist->Write();
	multiplicity_beforecuts->Write();
	multiplicity->Write();
	multposplus->Write();
	multposminus->Write();
	NeventsAss->Write();
	// ZDC
	hZDCPlusZeroBias_beforecuts->Write();
	hZDCMinusZeroBias_beforecuts->Write();
	hZDCPlusZDCOr_beforecuts->Write();
	hZDCMinusZDCOr_beforecuts->Write();
	hZDCPlusZeroBias->Write();
	hZDCMinusZeroBias->Write();
	hZDCPlusZDCOr->Write();
	hZDCMinusZDCOr->Write();
	// PF Candidates
	h_h->Write();
	h_e->Write();
	h_mu->Write();
	h_gamma->Write();
	h_h0 ->Write();
	h_HFhad->Write();
	h_HFem ->Write();
	// Tracks
	hist_reco_trk->Write();
}

void write_HBT1D(){
	hist_qinv_SS->Write();
	hist_qinv_SS_INV->Write();
	hist_qinv_SS_ROT->Write();
	hist_qinv_SS_MIX->Write();
	hist_qinv_OS->Write();
	hist_qinv_OS_INV->Write();
	hist_qinv_OS_ROT->Write();
	hist_qinv_OS_MIX->Write();
	hist_dpt_cos_SS->Write();
	hist_dpt_cos_OS->Write();
	hist_pairSS_Mass->Write();
	hist_pairOS_Mass->Write();
}

void write_HBT3D(){
	hist_q3D_SS->Write();
	hist_q3D_SS_INV->Write();
	hist_q3D_SS_ROT->Write();
	hist_q3D_SS_MIX->Write();
	hist_q3D_OS->Write();
	hist_q3D_OS_INV->Write();
	hist_q3D_OS_ROT->Write();
	hist_q3D_OS_MIX->Write();
}

