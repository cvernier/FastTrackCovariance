#include <TMath.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TH1.h>
#include <TAxis.h>
#include <TLegend.h>
#include <TCanvas.h>
#include <TGraph.h>
#include <TString.h>
#include <THStack.h>
#include <TFile.h>
#include <iostream>
#include "geometry_scripts/SolGeom.h"
#include "trkcovariance_scripts/SolTrack.h"
#include "trkcovariance_scripts/SolGridCov.h"
//
void plot_materialSiD()
{
	SolGeom *G;				// Initialize geometry
	G = new SolGeom("geometry_files/GeoSiD.txt",5.0);	// Geometry with selected detectors
	
	// *****************************
	// Make plot of material       *
	//******************************
	// 
	THStack *hMat = new THStack("hMat", "SiD: Material vs. cos(#theta)");

	TH1D *hPipe = new TH1D("hPisa", "SiD: Pipe material vs. Theta", 100, 0., 1.);
	hPipe->SetFillColor(kRed);
	hPipe->SetLineColor(kBlack);
	hPipe->GetYaxis()->SetRange(0.0, 100.);
	TH1D *hVtx = new TH1D("hVtx", "SiD: Silicon VTX material vs. Theta", 100, 0., 1.);
	hVtx->SetFillColor(kGreen);
	hVtx->SetLineColor(kBlack);
	hVtx->GetYaxis()->SetRange(0.0, 100.);
	TH1D *hTrk = new TH1D("hTrk", "SiD: Silicon TRK material vs. Theta", 100, 0., 1.);
	hTrk->SetFillColor(kCyan);
	hTrk->SetLineColor(kBlack);
	hTrk->GetYaxis()->SetRange(0.0, 100.);

	//hPipe->GetXaxis()->SetTitle("cos(#theta)");
	//hPipe->GetYaxis()->SetTitle("% radiation length");


	Int_t nStep = hPipe->GetNbinsX();
	for (Int_t i = 1; i <= nStep; i++)
	{
		Double_t CosTh = hPipe->GetBinCenter(i);
		Double_t th = TMath::ACos(CosTh);
		Double_t *mat = new Double_t[7]; //7 subsetectors for SiD, inc. magnet, which is not included here
		//0: Pipe, 1: VTX, 2: TRK, 3: VTXDSK, 4: VTXDSK_FWD, 5: TRKDSK, 6: MAG,

		mat = G->FracX0(th);
		Double_t mPipe = 100.*mat[0]; hPipe->SetBinContent(i, mPipe);
		Double_t mVtx = 100 * (mat[1] + mat[3] + mat[4]); hVtx->SetBinContent(i, mVtx);
		Double_t mTrk = 100.*(mat[2] + mat[5]); hTrk->SetBinContent(i, mTrk);
	}
	hMat->Add(hPipe);
	hMat->Add(hVtx);
	hMat->Add(hTrk);


	TCanvas *cmat = new TCanvas("cmat", "SiD material in tracking", 100, 100, 500, 500);
	cmat->cd(1);
	//TLegend *lg = new TLegend(0.5, 0.6, 0.9, 0.9);
	TLegend *lg = new TLegend(0.1, 0.9, 0.5, 0.7);
	lg->AddEntry(hPipe, "Beam pipe", "f");
	lg->AddEntry(hVtx, "Vertex", "f");
	lg->AddEntry(hTrk, "Tracker", "f");
	hMat->SetMaximum(15.);
	hMat->Draw(); lg->Draw();


	// Set axis titles on the stack (not individual histograms)
	hMat->GetXaxis()->SetTitle("cos(#theta)");
	hMat->GetYaxis()->SetTitle("% radiation length");


	cmat->Update();


	cmat->SaveAs("SiD_material_tracking.pdf");
	//
	
}

