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
#include <TSystem.h>
#include <iostream>
#include "geometry_scripts/SolGeom.h"
#include "trkcovariance_scripts/SolTrack.h"
#include "trkcovariance_scripts/SolGridCov.h"
//
void CompGeom(Double_t Ang, 
              const char* geom_file1, const char* geom_file2, const char* geom_file3,
              const char* det_name1, const char* det_name2, const char* det_name3,
              Double_t B_field1, Double_t B_field2, Double_t B_field3)
{

	Int_t iang = TMath::Nint(Ang); // Angle to integer

	//
	// Track to draw
	Double_t x[3] = { 0.0, 0.0, 0.0 };		// Track starting point
	Double_t ppt = 5.0;						// Track pt
	Double_t th = Ang * TMath::Pi() / 180.;
	Double_t ppz = ppt / TMath::Tan(th);	// Track pz
	Double_t p[3] = { ppt, 0.0, ppz };		// Track momentum


	//
	//**************************
	//	Initialize geometry    *
	//**************************
    //
    // First detector geometry (e.g. IDEA)
    TString filename1;
    filename1.Form("geometry_files/%s.txt", geom_file1);
    SolGeom *Gdet1;            
    char* fname1 = const_cast<char*>(filename1.Data());
    Gdet1 = new SolGeom(fname1, B_field1);     
    Gdet1->Draw();            
    TCanvas *cc_det1 = Gdet1->cnv();                    
    cc_det1->cd(1);
    SolTrack *trk_det1 = new SolTrack(x, p, Gdet1);    
    TGraph *gr_det1 = trk_det1->TrkPlot();            
    gr_det1->Draw("PLSAME");                        
    TCanvas *cnv_det1 = new TCanvas("cnv_det1", Form("%s Geometry", det_name1));
    cnv_det1->Divide(1, 1); cnv_det1->cd(1); cc_det1->DrawClonePad();
    TString filename_det1_track;
    filename_det1_track.Form("plot_dump/%s_geometry_track_%d.pdf", det_name1, iang);
    cnv_det1->SaveAs(filename_det1_track);

    //
    // Second detector geometry (e.g. CLD)
    TString filename2;
    filename2.Form("geometry_files/%s.txt", geom_file2);
    SolGeom *Gdet2;            
    char* fname2 = const_cast<char*>(filename2.Data());
    Gdet2 = new SolGeom(fname2, B_field2); 
    Gdet2->Draw();            
    TCanvas *cc_det2 = Gdet2->cnv();
    cc_det2->cd(1);
    SolTrack *trk_det2 = new SolTrack(x, p, Gdet2);    
    TGraph *gr_det2 = trk_det2->TrkPlot();            
    gr_det2->Draw("PLSAME");                        
    TCanvas *cnv_det2 = new TCanvas("cnv_det2", Form("%s Geometry", det_name2));
    cnv_det2->Divide(1, 1); cnv_det2->cd(1); cc_det2->DrawClonePad();
    TString filename_det2_track;
    filename_det2_track.Form("plot_dump/%s_geometry_track_%d.pdf", det_name2, iang);
    cnv_det2->SaveAs(filename_det2_track);
    cc_det2->Close(); gSystem->ProcessEvents(); delete cc_det2; cc_det2 = 0;

    //
    // Third detector geometry (e.g. SiD)
    TString filename3;
    filename3.Form("geometry_files/%s.txt", geom_file3);
    SolGeom *Gdet3;            
    char* fname3 = const_cast<char*>(filename3.Data());
    Gdet3 = new SolGeom(fname3, B_field3);   
    Gdet3->Draw();            
    TCanvas *cc_det3 = Gdet3->cnv();
    cc_det3->cd(1);
    SolTrack *trk_det3 = new SolTrack(x, p, Gdet3);    
    TGraph *gr_det3 = trk_det3->TrkPlot();            
    gr_det3->Draw("PLSAME");                        
    TCanvas *cnv_det3 = new TCanvas("cnv_det3", Form("%s Geometry", det_name3));
    cnv_det3->Divide(1, 1); cnv_det3->cd(1); cc_det3->DrawClonePad();
    TString filename_det3_track;
    filename_det3_track.Form("plot_dump/%s_geometry_track_%d.pdf", det_name3, iang);
    cnv_det3->SaveAs(filename_det3_track);
    cnv_det3->Close(); gSystem->ProcessEvents(); delete cc_det3; cc_det3 = 0;


	//********************************************************
	//	
	//******************************************************
	// Compare track parameter resolutions vs pt and theta *
	//******************************************************
	//
	TCanvas *resol = new TCanvas("resol", "Comparison of resolutions", 100, 100, 500, 500);
	TString CnvTitle; 
	CnvTitle.Form("Comparison of resolutions - Track angle %d deg.",iang);
	resol->SetTitle(CnvTitle);
	
	resol->Divide(2, 2);
	// Define graphs for first detector
    TGraph *grpt_det1;                // pt resolution graphs
    TGraph *grptms_det1;             // pt resolution graphs MS only
    TGraph *grd0_det1;               // D resolution graphs
    TGraph *grz0_det1;               // z0 resolution graphs
    TGraph *grth_det1;               // theta resolution
    
	// Define graphs for second detector
    TGraph *grpt_det2;                // pt resolution graphs
    TGraph *grptms_det2;             // pt resolution graphs MS only
    TGraph *grd0_det2;               // D resolution graphs
    TGraph *grz0_det2;               // z0 resolution graphs
    TGraph *grth_det2;               // theta resolution

	// Define graphs for first detector
    TGraph *grpt_det3;                // pt resolution graphs
    TGraph *grptms_det3;             // pt resolution graphs MS only
    TGraph *grd0_det3;               // D resolution graphs
    TGraph *grz0_det3;               // z0 resolution graphs
    TGraph *grth_det3;               // theta resolution

	// Setup graph arrays
	Int_t Npt = 200;			// Nr. of points per graph
	Double_t * pt = new Double_t[Npt];
	Double_t * pp = new Double_t[Npt];

	Double_t *spt_det1 = new Double_t[Npt];
	Double_t *sptms_det1 = new Double_t[Npt];
	Double_t *sd0_det1 = new Double_t[Npt];
	Double_t *sz0_det1 = new Double_t[Npt];
	Double_t *sth_det1 = new Double_t[Npt];
	
	Double_t *spt_det2 = new Double_t[Npt];
	Double_t *sptms_det2 = new Double_t[Npt];
	Double_t *sd0_det2 = new Double_t[Npt];
	Double_t *sz0_det2 = new Double_t[Npt];
	Double_t *sth_det2 = new Double_t[Npt];
	
	Double_t *spt_det3 = new Double_t[Npt];
	Double_t *sptms_det3 = new Double_t[Npt];
	Double_t *sd0_det3 = new Double_t[Npt];
	Double_t *sz0_det3 = new Double_t[Npt];
	Double_t *sth_det3 = new Double_t[Npt];

	// Fill graph arrays
	Double_t ptmin = 2.0;
	Double_t ptmax = 100;
	Double_t pts = (ptmax - ptmin) / (Double_t)(Npt-1);
	for (Int_t k = 0; k < Npt; k++)	// Loop on pt
	{
		Double_t x[3]; Double_t p[3];
		x[0] = 0; x[1] = 0; x[2] = 0;			// Set origin
		pt[k] = ptmin+k*pts;					// Set transverse momentum
		p[0] = pt[k]; p[1] = 0;	p[2] = pt[k] / TMath::Tan(th);
		pp[k] = pt[k]/TMath::Sin(th);			// Set momentum
		
		
		// Fill det1 arrays
		SolTrack *tr_det1 = new SolTrack(x, p, Gdet1);	// Initialize track
		Bool_t Res = kTRUE;		// Enable detector resolution effects
		Bool_t MS  = kTRUE;		// Enable multiple scattering
		tr_det1->CovCalc(Res,MS);					// Calculate covariance
		spt_det1[k] = tr_det1->s_pt();							// Dpt/pt
		//spt_det1[k] = tr_det1->s_pt()/pt[k];							// Dpt/pt^2 -Jim: change this to sigma pt/pt^2 to match what is shown in ILC TDR - 2*sigma(pt)/pt^2
		sd0_det1[k] = tr_det1->s_D()*1e6;							// D  res. - change to microns
		sz0_det1[k] = tr_det1->s_z0()*1e6;						// z0 res. - change to microns
		sth_det1[k] = tr_det1->s_ct() / (1 + pow(tr_det1->ct(), 2));	// theta resolution
		//
		Res = kFALSE;
		MS = kTRUE;
		tr_det1->CovCalc(Res,MS);					// Calculate covariance with only MS
		sptms_det1[k] = tr_det1->s_pt();							// Dpt/pt
		//sptms_det1[k] = tr_det1->s_pt()/pt[k];							// Dpt/pt^2 -Jim: change this to sigma pt/pt^2 to match what is shown in ILC TDR - 2*sigma(pt)/pt^2

		// Fill det2 arrays
		SolTrack *tr_det2 = new SolTrack(x, p, Gdet2);	// Initialize track
		Res = kTRUE;		// Enable detector resolution effects
		MS  = kTRUE;		// Enable multiple scattering
		tr_det2->CovCalc(Res,MS);					// Calculate covariance
		spt_det2[k] = tr_det2->s_pt();							// Dpt/pt
		//spt_det1[k] = tr_det1->s_pt()/pt[k];							// Dpt/pt^2 -Jim: change this to sigma pt/pt^2 to match what is shown in ILC TDR - 2*sigma(pt)/pt^2
		sd0_det2[k] = tr_det2->s_D()*1e6;							// D  res. - change to microns
		sz0_det2[k] = tr_det2->s_z0()*1e6;						// z0 res. - change to microns
		sth_det2[k] = tr_det2->s_ct() / (1 + pow(tr_det2->ct(), 2));	// theta resolution
		//
		Res = kFALSE;
		MS = kTRUE;
		tr_det2->CovCalc(Res,MS);					// Calculate covariance with only MS
		sptms_det2[k] = tr_det2->s_pt();							// Dpt/pt
		//sptms_det1[k] = tr_det1->s_pt()/pt[k];							// Dpt/pt^2 -Jim: change this to sigma pt/pt^2 to match what is shown in ILC TDR - 2*sigma(pt)/pt^2

		// Fill det3 arrays
		SolTrack *tr_det3 = new SolTrack(x, p, Gdet3);	// Initialize track
		Res = kTRUE;		// Enable detector resolution effects
		MS  = kTRUE;		// Enable multiple scattering
		tr_det3->CovCalc(Res,MS);					// Calculate covariance
		spt_det3[k] = tr_det3->s_pt();							// Dpt/pt
		//spt_det1[k] = tr_det1->s_pt()/pt[k];							// Dpt/pt^2 -Jim: change this to sigma pt/pt^2 to match what is shown in ILC TDR - 2*sigma(pt)/pt^2
		sd0_det3[k] = tr_det3->s_D()*1e6;							// D  res. - change to microns
		sz0_det3[k] = tr_det3->s_z0()*1e6;						// z0 res. - change to microns
		sth_det3[k] = tr_det3->s_ct() / (1 + pow(tr_det3->ct(), 2));	// theta resolution
		//
		Res = kFALSE;
		MS = kTRUE;
		tr_det3->CovCalc(Res,MS);					// Calculate covariance with only MS
		sptms_det3[k] = tr_det3->s_pt();							// Dpt/pt
		//sptms_det1[k] = tr_det1->s_pt()/pt[k];							// Dpt/pt^2 -Jim: change this to sigma pt/pt^2 to match what is shown in ILC TDR - 2*sigma(pt)/pt^2


	}
	//
	// Compare pt resolution
	resol->cd(1);
	// det1
	grpt_det1 = new TGraph(Npt, pt, spt_det1);			// pt resolution
	grpt_det1->SetLineColor(kRed);
	grpt_det1->SetMarkerColor(kRed);
	grpt_det1->SetTitle("#sigma_{pt}/pt");
	//grpt_det1->SetTitle("#sigma_{pt}/pt^2 (GeV^{-1})");
	grpt_det1->SetMinimum(0.0);
	grpt_det1->GetXaxis()->SetTitle("pt (GeV)");
	grpt_det1->Draw("APL");
	grptms_det1 = new TGraph(Npt, pt, sptms_det1);			// pt resolution MS only
	grptms_det1->SetLineColor(kRed);
	grptms_det1->SetMarkerColor(kRed);
	grptms_det1->SetLineStyle(7);
	grptms_det1->SetTitle("#sigma_{pt}/pt");
	//grptms_det1->SetTitle("#sigma_{pt}/pt^2 (GeV^{-1})");
	grptms_det1->SetMinimum(0.0);
	grptms_det1->GetXaxis()->SetTitle("pt (GeV)");
	grptms_det1->Draw("SAME");
	// det2
	grpt_det2 = new TGraph(Npt, pt, spt_det2);			// pt resolution
	grpt_det2->SetLineColor(kBlue);
	grpt_det2->SetMarkerColor(kBlue);
	grpt_det2->SetTitle("#sigma_{pt}/pt");
	//grpt_det2->SetTitle("#sigma_{pt}/pt^2 (GeV^{-1})");
	grpt_det2->SetMinimum(0.0);
	grpt_det2->GetXaxis()->SetTitle("pt (GeV)");
	grpt_det2->Draw("SAME");
	grptms_det2 = new TGraph(Npt, pt, sptms_det2);			// pt resolution MS only
	grptms_det2->SetLineColor(kBlue);
	grptms_det2->SetMarkerColor(kBlue);
	grptms_det2->SetLineStyle(7);
	grptms_det2->SetTitle("#sigma_{pt}/pt");
	//grptms_det2->SetTitle("#sigma_{pt}/pt^2 (GeV^{-1})");
	grptms_det2->SetMinimum(0.0);
	grptms_det2->GetXaxis()->SetTitle("pt (GeV)");
	grptms_det2->Draw("SAME");
	//SiD
	grpt_det3 = new TGraph(Npt, pt, spt_det3);			// pt resolution
	grpt_det3->SetLineColor(kGreen);
	grpt_det3->SetMarkerColor(kGreen);
	grpt_det3->SetTitle("#sigma_{pt}/pt");
	//grpt_det3->SetTitle("#sigma_{pt}/pt^2 (GeV^{-1})");
	grpt_det3->SetMinimum(0.0);
	grpt_det3->GetXaxis()->SetTitle("pt (GeV)");
	grpt_det3->Draw("SAME");
	grptms_det3 = new TGraph(Npt, pt, sptms_det3);			// pt resolution MS only
	grptms_det3->SetLineColor(kGreen);
	grptms_det3->SetMarkerColor(kGreen);
	grptms_det3->SetLineStyle(7);
	grptms_det3->SetTitle("#sigma_{pt}/pt");
	//grptms_det3->SetTitle("#sigma_{pt}/pt^2 (GeV^{-1})");
	grptms_det3->SetMinimum(0.0);
	grptms_det3->GetXaxis()->SetTitle("pt (GeV)");
	grptms_det3->Draw("SAME");

	
	TLegend *lgpt = new TLegend(0.2, 0.9, 0.6, 0.70);
	TString LgTitle; 
    LgTitle.Form("Track angle %d deg.", iang);
    lgpt->SetHeader(LgTitle);
    lgpt->AddEntry(grpt_det1, det_name1, "L");
    lgpt->AddEntry(grptms_det1, Form("%s MS only", det_name1), "L");
    lgpt->AddEntry(grpt_det2, det_name2, "L");
    lgpt->AddEntry(grptms_det2, Form("%s MS only", det_name2), "L");
    lgpt->AddEntry(grpt_det3, det_name3, "L");
    lgpt->AddEntry(grptms_det3, Form("%s MS only", det_name3), "L");
    lgpt->Draw();

	// Compare d0 resolution
	resol->cd(2);
	grd0_det1 = new TGraph(Npt, pp, sd0_det1);			// D resolution
	grd0_det1->SetLineColor(kBlue);
	grd0_det1->SetMarkerColor(kBlue);
	grd0_det1->SetTitle("#sigma_{D_{0}} (#mum)");
	grd0_det1->SetMinimum(0.0);
	grd0_det1->GetXaxis()->SetTitle("p (GeV)");
	grd0_det1->Draw("APL");
	grd0_det2 = new TGraph(Npt, pp, sd0_det2);			// D resolution
	grd0_det2->SetLineColor(kRed);
	grd0_det2->SetMarkerColor(kRed);
	grd0_det2->SetTitle("#sigma_{D_{0}} (#mum)");
	grd0_det2->SetMinimum(0.0);
	grd0_det2->GetXaxis()->SetTitle("p (GeV)");
	grd0_det2->Draw("SAME"); 
	grd0_det3 = new TGraph(Npt, pp, sd0_det3);			// D resolution
	grd0_det3->SetLineColor(kGreen);
	grd0_det3->SetMarkerColor(kGreen);
	grd0_det3->SetTitle("#sigma_{D_{0}} (#mum)");
	grd0_det3->SetMinimum(0.0);
	grd0_det3->GetXaxis()->SetTitle("p (GeV)");
	grd0_det3->Draw("SAME");


	TLegend *lgd0 = new TLegend(0.2, 0.9, 0.6, 0.70);
	lgd0->SetHeader(LgTitle);
	lgd0->AddEntry(grpt_det1, det_name1, "L");
	lgd0->AddEntry(grpt_det2, det_name2, "L");
	lgd0->AddEntry(grpt_det3, det_name3, "L");
	lgd0->Draw();






	// Compare z0 resolution
	resol->cd(3);
	//IDEA
	grz0_det1 = new TGraph(Npt, pp, sz0_det1);			// z0 resolution
	grz0_det1->SetLineColor(kBlue);
	grz0_det1->SetMarkerColor(kBlue);
	grz0_det1->SetTitle("#sigma_{Z_{0}} (#mum)");
	grz0_det1->GetXaxis()->SetTitle("p (GeV)");
	grz0_det1->Draw("APL");
	//CLD
	grz0_det2 = new TGraph(Npt, pp, sz0_det2);			// z0 resolution
	grz0_det2->SetLineColor(kRed);
	grz0_det2->SetMarkerColor(kRed);
	grz0_det2->SetTitle("#sigma_{Z_{0}} (#mum)");
	grz0_det2->GetXaxis()->SetTitle("p (GeV)");
	grz0_det2->Draw("SAME");			
	//SiD
	grz0_det3 = new TGraph(Npt, pp, sz0_det3);			// z0 resolution
	grz0_det3->SetLineColor(kGreen);
	grz0_det3->SetMarkerColor(kGreen);
	grz0_det3->SetTitle("#sigma_{Z_{0}} (#mum)");
	grz0_det3->GetXaxis()->SetTitle("p (GeV)");
	grz0_det3->Draw("SAME");

	TLegend *lgz0 = new TLegend(0.2, 0.9, 0.6, 0.70);
	lgz0->SetHeader(LgTitle);
	lgz0->AddEntry(grpt_det1, det_name1, "L");
	lgz0->AddEntry(grpt_det2, det_name2, "L");
	lgz0->AddEntry(grpt_det3, det_name3, "L");
	lgz0->Draw();



	resol->cd(4);

	// Compare theta resolution
	//IDEA
	grth_det1 = new TGraph(Npt, pp, sth_det1);			// theta resolution
	grth_det1->SetLineColor(kBlue);
	grth_det1->SetMarkerColor(kBlue);
	grth_det1->SetTitle("#sigma_{#theta} (rad)");
	grth_det1->SetMinimum(0.0);
	grth_det1->GetXaxis()->SetTitle("p (GeV)");
	grth_det1->Draw("APL");
	//CLD
	grth_det2 = new TGraph(Npt, pp, sth_det2);			// theta resolution
	grth_det2->SetLineColor(kRed);
	grth_det2->SetMarkerColor(kRed);
	grth_det2->SetTitle("#sigma_{#theta} (rad)");
	grth_det2->SetMinimum(0.0);
	grth_det2->GetXaxis()->SetTitle("p (GeV)");
	grth_det2->Draw("SAME");
	//SiD
	grth_det3 = new TGraph(Npt, pp, sth_det3);			// theta resolution
	grth_det3->SetLineColor(kGreen);
	grth_det3->SetMarkerColor(kGreen);
	grth_det3->SetTitle("#sigma_{#theta} (rad)");
	grth_det3->SetMinimum(0.0);
	grth_det3->GetXaxis()->SetTitle("p (GeV)");
	grth_det3->Draw("SAME");

	TLegend *lgth = new TLegend(0.2, 0.9, 0.6, 0.70);
	lgth->SetHeader(LgTitle);
	lgth->AddEntry(grpt_det1, det_name1, "L");
	lgth->AddEntry(grpt_det2, det_name2, "L");
	lgth->AddEntry(grpt_det3, det_name3, "L");
	lgth->Draw();


    TString filename;
    filename.Form("plot_dump/%s_%s_%s_resolution_tracking_%d.pdf", 
                 det_name1, det_name2, det_name3, iang);
    resol->SaveAs(filename);
	//


	TCanvas *resolp = new TCanvas("resolp", "Comparison of pt resolutions", 50, 50, 500, 500);
	resolp->Divide(1, 1);
	// Compare pt resolution
	// CLD
	resolp->cd(1); 
	grpt_det2->SetMaximum(0.005);
	grpt_det2->Draw("APL");
	grptms_det2->Draw("SAME");
	// SiD
	grpt_det3->Draw("SAME");
	grptms_det3->Draw("SAME");
	// IDEA
	grpt_det1->SetLineColor(kBlue);
	grpt_det1->SetLineStyle(2);
	grpt_det1->SetMarkerColor(kBlue);
	grpt_det1->SetTitle("#sigma_{pt}/pt");
	grpt_det1->SetMinimum(0.0);
	grpt_det1->GetXaxis()->SetTitle("pt (GeV)");
	grpt_det1->Draw("SAME");
	grptms_det1->Draw("SAME");




	TLegend *lgpt1 = new TLegend(0.2, 0.9, 0.6, 0.70);
	lgpt1->SetHeader(LgTitle);
	lgpt1->AddEntry(grpt_det1, det_name1, "L");
	lgpt1->AddEntry(grptms_det1, det_name1, "L");
	lgpt1->AddEntry(grpt_det2, det_name2, "L");
	lgpt1->AddEntry(grptms_det2, det_name2, "L");
	lgpt1->AddEntry(grpt_det3, det_name3, "L");
	lgpt1->AddEntry(grptms_det3, det_name3, "L");
	lgpt1->Draw();

    TCanvas *resolp = new TCanvas("resolp", "Comparison of pt resolutions", 50, 50, 500, 500);
    // [Rest of the pt resolution plotting code remains similar, just update names]
    
    TString filenamep;
    filenamep.Form("plot_dump/%s_%s_%s_resolution_tracking_pt_%d.pdf", 
                   det_name1, det_name2, det_name3, iang);
    resolp->SaveAs(filenamep);


}

