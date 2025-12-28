#include <iostream>
#include "geometry_scripts/SolGeom.h"
#include "trkcovariance_scripts/SolGridCov.h"
#include "trkcovariance_scripts/ObsTrk.h"

#include "TVector3.h"
#include "TLorentzVector.h"
#include "TH1D.h"
#include "TF1.h"
#include "TVectorD.h"
#include "TMatrixDSym.h"
#include "TMatrixD.h"
#include "TRandom.h"

using namespace std;

//
// Regularized symmetric positive definite  matrix inversion
//
TMatrixDSym SymRegInv(TMatrixDSym &Smat)
{
	//
	Int_t N = Smat.GetNrows();
	TMatrixDSym D(N); D.Zero();
	for (Int_t i = 0; i < N; i++) D(i, i) = 1.0/TMath::Sqrt(Smat(i, i));
	TMatrixDSym RegMat = Smat.Similarity(D); 
	TDecompChol rChl(RegMat);
	Bool_t OK = rChl.Decompose();
	if (!OK)
	{
		cout << "RegMat: input matrix not positive definite"; RegMat.Print();
	}
	RegMat = rChl.Invert(OK);
	TMatrixDSym RegOut = RegMat.Similarity(D);
	//
	return RegOut;
}
//
TVectorD Fillf(TVectorD par, TVectorD xin)
{
	//
	// Decode input arrays
	//
	Double_t D = par(0);
	Double_t p0 = par(1);
	Double_t C = par(2);
	Double_t z0 = par(3);
	Double_t ct = par(4);
	//
	Double_t x = xin(0);
	Double_t y = xin(1);
	Double_t z = xin(2);
	Double_t R2 = x * x + y * y;
	Double_t D2 = D*D;
	Double_t UpCD = 1.0 + C*D;
	Double_t Up2CD = 1.0 + 2 * C*D;
	Double_t Az = C*(z - z0) / ct;
	//
	// Calculate constraints
	//
	TVectorD f(2);
	f(0) = (R2 * C + UpCD*D)/Up2CD - y * TMath::Cos(p0) + x * TMath::Sin(p0);
	f(1) = TMath::Sin(Az)*TMath::Sin(Az) - C*C*(R2 - D2) / Up2CD;
	//
	return f;
}
//

TMatrixD FillD(TVectorD par, TVectorD xin)
{
	//
	// Decode input arrays
	//
	Double_t D = par(0);
	Double_t p0 = par(1);
	Double_t C = par(2);
	Double_t z0 = par(3);
	Double_t ct = par(4);
	//
	Double_t x = xin(0);
	Double_t y = xin(1);
	Double_t z = xin(2);
	Double_t R2 = x * x + y * y;
	Double_t D2 = D*D;
	Double_t UpCD = 1.0 + C*D;
	Double_t Up2CD = 1.0 + 2 * C*D;
	Double_t Az = C*(z - z0) / ct;
	//
	// Calculate matrix elements
	//
	TMatrixD Do(2, 5); Do.Zero();
	Do(0, 0) = 1.0 - 2 * C*(D + C*(R2 - D2)) / (Up2CD*Up2CD);	// df(0)/dD
	Do(0, 1) = x * TMath::Cos(p0) + y * TMath::Sin(p0);			// df(0)/dphi0
	Do(0, 2) = (R2 - D2) / (Up2CD*Up2CD);						// df(0)/dC
	Do(0, 3) = 0.0;												// df(0)/dz0
	Do(0, 4) = 0.0;												// df(0)/dct
	Do(1, 0) = 2 * C*C*(C*R2 + D*UpCD) / (Up2CD*Up2CD);			// df(1)/dD
	Do(1, 1) = 0.0;												// df(1)/dphi0
	Do(1, 2) = TMath::Sin(2 * Az)*Az / C - 2 * C*(R2 - D2)*UpCD / (Up2CD*Up2CD); // df(1)/dC
	Do(1, 3) = -TMath::Sin(2 * Az)*C / ct;						// df(1)/dz0
	Do(1, 4) = -TMath::Sin(2 * Az)*Az / ct;						// df(1)/dct
	//
	return Do;
}
//
TMatrixD FillB(TVectorD par, TVectorD xin)
{
	//
	// Decode input arrays
	//
	Double_t D = par(0);
	Double_t p0 = par(1);
	Double_t C = par(2);
	Double_t z0 = par(3);
	Double_t ct = par(4);
	//
	Double_t x = xin(0);
	Double_t y = xin(1);
	Double_t z = xin(2);
	Double_t R2 = x * x + y * y;
	Double_t D2 = D*D;
	Double_t UpCD = 1.0 + C*D;
	Double_t Up2CD = 1.0 + 2 * C*D;
	Double_t Az = C*(z - z0) / ct;
	//
	// Calculate constraints
	//
	TMatrixD B(2, 3); B.Zero();
	B(0, 0) = 2 * C*x/Up2CD + TMath::Sin(p0);
	B(0, 1) = 2 * C*y/Up2CD - TMath::Cos(p0);
	B(0, 2) = 0.0;
	B(1, 0) = -2 * x*C*C / Up2CD;
	B(1, 1) = -2 * y*C*C / Up2CD;
	B(1, 2) = TMath::Sin(2 * Az)*C / ct;
	//
	return B;
}
//
TVectorD Vertex0(Int_t Ntr, ObsTrk **tracks)
{
	//
	// Preliminary estimate of the vertex position
	// based on transformation of track into points
	// and vertices into lines
	// No steering of track parameters
	// No error calculation
	//
	TVectorD xv(3);		// returned vertex position
	//
	TMatrixDSym H(2); 
	TVectorD xvt(2);
	TVectorD cxy(2); 
	//
	// Loop on tracks for transverse fit
	//
	TVectorD x0(2); x0.Zero();
	Double_t Rv = 0.0;	    // Radius of first iteration
	Int_t Ntry = 0;
	Int_t TryMax = 10;
	Double_t epsi = 1000.;		// Starting stability
	Double_t eps = 0.0001;		// vertex stability required
	while (epsi > eps && Ntry < TryMax)
	{
		H.Zero(); cxy.Zero();
		for (Int_t i = 0; i < Ntr; i++)
		{
			// Get track helix parameters and their covariance matrix 
			ObsTrk *t = tracks[i];
			TVectorD par = t->GetObsPar();
			TMatrixDSym C = t->GetCov();
			//
			// Transverse fit
			Double_t D0i = par(0);
			Double_t phi = par(1);
			Double_t Ci = par(2);
			Double_t Di = (D0i*(1. + Ci*D0i) + Rv*Rv*Ci) / (1. + 2 * Ci*D0i);
			Double_t sDi2 = C(0, 0);
			//
			TVectorD ni(2);
			ni(0) = -TMath::Sin(phi);
			ni(1) =  TMath::Cos(phi);
			TMatrixDSym Hadd(2);
			Hadd.Rank1Update(ni, 1);		// Tensor product of vector ni with itself
			H += (1.0 / sDi2)*Hadd;
			cxy += (Di / sDi2)*ni;
			
		}
		//
		TMatrixDSym Cov = SymRegInv(H);
		xvt = Cov*cxy;
		xv.SetSub(0, xvt);	// Store x,y of vertex
		Rv = TMath::Sqrt(xv(0)*xv(0) + xv(1)*xv(1));
		TVectorD dx = xvt - x0;
		epsi = H.Similarity(dx);
		x0 = xvt;
		Ntry++;
		//cout << "Vtx0: Iteration #" << Ntry << ", eps = " << epsi << ", x = " << xv(0) << ", y = " << xv(1);
	}
	//
	// Longitudinal fit
	Double_t Rv2 = Rv*Rv;
	//
	// Loop on tracks for longitudinal fit
	Double_t hz = 0.0;
	Double_t cz = 0.0;
	for (Int_t i = 0; i < Ntr; i++)
	{
		// Get track helix parameters and their covariance matrix 
		ObsTrk *t = tracks[i];
		TVectorD par = t->GetObsPar();
		TMatrixDSym C = t->GetCov();
		//
		// Longitudinal fit
		Double_t zi = par(3);
		Double_t cti = par(4);
		Double_t Di = par(0);
		Double_t Ci = par(2);
		Double_t sZi2 = C(3, 3);
		//
		hz += 1 / sZi2;
		Double_t arg = TMath::Sqrt(TMath::Max(0.0, Rv*Rv - Di*Di)/(1.+2*Ci*Di));
		cz += (cti*arg + zi)/sZi2;
	}
	xv(2) = cz / hz;
	//cout << ", z = " << xv(2) << endl;
	//
	return xv;
}
//
Double_t Vertex(Int_t Ntr, ObsTrk **tracks, TVectorD &x, TMatrixDSym &covX)
{
	//
	// Get approximate vertex evaluation
	//
	TVectorD x0 = Vertex0(Ntr, tracks);
	//cout << "Preliminary vertex" << endl; x0.Print();
	TVectorD dx(3);		// Solution x variation
	//
	TVectorD f(2);		// Constraints
	TMatrixD D(2, 5);	// df/d alf (constraint over parameters)
	TMatrixD B(2, 3);	// df/dx    (constraint over x) 
	// Stored quantities
	TVectorD **fi = new TVectorD*[Ntr];
	TVectorD **pi = new TVectorD*[Ntr];
	TMatrixD **Di = new TMatrixD*[Ntr];
	TMatrixD **Bi = new TMatrixD*[Ntr];
	TMatrixDSym **Wi = new TMatrixDSym*[Ntr];
	TMatrixDSym **Ci = new TMatrixDSym*[Ntr];
	//
	// Loop on tracks to calculate everything
	//
	Int_t Ntry = 0;
	Int_t TryMax = 100;
	Double_t eps = 0.001; // vertex stability
	Double_t epsi = 1000.;
	x = x0;
	// Protect for vertices close to 0
	/*
	Double_t pvx = 1.0e-8;
	if (x(0)*x(0) + x(1)*x(1) < pvx)
	{
		Double_t rn = gRandom->Rndm();
		Double_t phrn = TMath::TwoPi()*rn;
		x(0) += 2.*pvx*TMath::Cos(phrn);
		x(1) += 2.*pvx*TMath::Sin(phrn);
	}
	*/
	while (epsi > eps && Ntry < TryMax)		// Iterate until found vertex is stable
	{
		TVectorD BtWf(3); BtWf.Zero();
		covX.Zero();		// Reset vertex covariance
		// 
		for (Int_t i = 0; i < Ntr; i++)
		{
			// Get track helix parameters and their covariance matrix 
			ObsTrk *t = tracks[i];
			TVectorD par0 = t->GetObsPar();
			TMatrixDSym C = t->GetCov();
			Ci[i] = new TMatrixDSym(C);
			TVectorD par(5);
			if (Ntry > 0) par = *pi[i];
			else
			{
				par = par0;
				pi[i] = new TVectorD(par0);
			}
			// Fill D
			D = FillD(par, x);
			Di[i] = new TMatrixD(D);
			// Fill B
			B = FillB(par, x);
			Bi[i] = new TMatrixD(B);
			//cout << "Bi" << endl; Bi[i]->Print();
			// Fill constraints
			f = Fillf(par, x);
			fi[i] = new TVectorD(f);
			//cout << "fi" << endl; fi[i]->Print();
			//
			TMatrixDSym W = C.Similarity(D);
			//cout << "W: "; W.Print();
			//W.Invert();
			W = SymRegInv(W);
			Wi[i] = new TMatrixDSym(W);
			//cout << "Wi" << endl; Wi[i]->Print();
			TMatrixD Bt(TMatrixD::kTransposed, B);
			TMatrixDSym W1(W);
			TMatrixDSym BtWB = W1.Similarity(Bt);
			covX += BtWB;
			BtWf += Bt * (W*f);
		}
		// Update vertex covariance
		TMatrixDSym Hess = covX;
		//cout << "Hesse: "; Hess.Print();
		//covX.Invert();
		covX = SymRegInv(Hess);
		//covX.Print();
		// update vertex position
		dx = (-1.0*covX) * BtWf;
		//dx.Print();
		x += dx;
		// Update track parameters
		for (Int_t i = 0; i < Ntr; i++)
		{
			TVectorD lambda = *fi[i] + (*Bi[i]) * dx;
			TMatrixD Dt(TMatrixD::kTransposed, *Di[i]);
			*pi[i] = *pi[i] - ((*Ci[i])*Dt) * lambda;
		}
		// update vertex stability
		epsi = Hess.Similarity(dx);
		Ntry++;
		//if (epsi >10)
		//cout << "Vtx:  Iteration #"<<Ntry<<", eps = "<<epsi<<", x = " << x(0) << ", " << x(1) << ", " << x(2) << endl;
	}
	//
	// Calculate Chi2
	//
	Double_t Chi2 = 0.0;
	for (Int_t i = 0; i < Ntr; i++)
	{
		TVectorD lambda = *fi[i] + (*Bi[i]) * dx;
		TMatrixDSym Wp = *Wi[i];
		Chi2 += Wp.Similarity(lambda);
	}
	//
	return Chi2;
}
Double_t fchi2(Double_t *x, Double_t *p)
{
	Double_t Ndof = p[0];
	Double_t Norm = p[1];
	Double_t z = x[0];
	//
	Double_t ex = Ndof / 2.0 - 1.0;
	Double_t chi2 = pow(z, ex)*TMath::Exp(-z / 2.0);
	Double_t bin = 5.*Ndof / 100.;
	chi2 = Norm * chi2 * bin / (pow(2, Ndof / 2.0)*TMath::Gamma(Ndof / 2.0));
	return chi2;
}






// needed input:


// - input Covariance Matrix in root file form
// 
// - provide vector position at origin, and vector momentum

int main() {

  Int_t Nvtx = 100;
  Int_t Ntr = 2;
  
  
  // initialize Geometry
  SolGeom *G = new SolGeom();
  Double_t Bfield = G->B();

  // initialize tracking resolution
  SolGridCov *GC = new SolGridCov();

  GC->Write("test.root",G);
  GC->Read("test.root");
  
  // Ranges
  //
  Double_t ThDegMin = 40.0;
  Double_t ThDegMax = 140.0;
  Double_t Lmin = 0.001;
  Double_t Lmax = 0.01;
  Double_t dTheta = 0.10;
  Double_t dPhi = 0.20;
  Double_t Pmin = 1.0;
  Double_t Pmax = 10.0;
  //
  // Histograms
  TH1D *hXpull0 = new TH1D("hXpull0", "Pull X vertex0 component", 100, -5., 5.);
  TH1D *hYpull0 = new TH1D("hYpull0", "Pull Y vertex0 component", 100, -5., 5.);
  TH1D *hZpull0 = new TH1D("hZpull0", "Pull Z vertex0 component", 100, -5., 5.);
  TH1D *hXpull = new TH1D("hXpull", "Pull X vertex component", 100, -5., 5.);
  TH1D *hYpull = new TH1D("hYpull", "Pull Y vertex component", 100, -5., 5.);
  TH1D *hZpull = new TH1D("hZpull", "Pull Z vertex component", 100, -5., 5.);
  Double_t Ndof = 2.0 * Ntr - 3.0;
  TH1D *hChi2 = new TH1D("hChi2", "Vertex #chi^{2}", 100, 0., 5 * Ndof);
  TF1 *fch = new TF1("fch", fchi2, 0., 5 * Ndof, 2);
  fch->SetParameter(0, Ndof);
  fch->SetParameter(1, Nvtx);
  //
  // Loop on # vertices
  //
  for (Int_t n = 0; n < Nvtx; n++)
  {
	  //
	  // Extract a direction
	  Double_t ThMin = ThDegMin * TMath::Pi() / 180.;
	  Double_t ThMax = ThDegMax * TMath::Pi() / 180.;
	  Double_t rnTh = gRandom->Rndm();
	  Double_t Th = ThMin + rnTh * (ThMax - ThMin);
	  Double_t rnPh = gRandom->Rndm();
	  Double_t Ph = 2 * TMath::Pi()*rnPh;
	  //
	  // Extract a flight distance (m)
	  Double_t rnL = gRandom->Rndm();
	  Double_t Lvtx = Lmin + rnL * (Lmax - Lmin);
	  //
	  TVector3 x;
	  Double_t zz0 = 0.0;
	  x(0) = Lvtx * TMath::Sin(Th)*TMath::Cos(Ph);
	  x(1) = Lvtx * TMath::Sin(Th)*TMath::Sin(Ph);
	  x(2) = Lvtx * TMath::Cos(Th) + zz0;;
	  //
	  cout << "True vertex: x = " << x(0) << ", y = " << x(1) << ", z = " << x(2) << endl;
	  //
	  // Loop on tracks
	  ObsTrk **tracks = new ObsTrk*[Ntr];	// Smear tracks according to covariance matrix
	  for (Int_t i = 0; i < Ntr; i++)
	  {
		  Double_t rnP = gRandom->Rndm();
		  Double_t Ptot = Pmin + rnP * (Pmax - Pmin);
		  Double_t ThP = gRandom->Gaus(Th, dTheta);
		  Double_t PhP = gRandom->Gaus(Ph, dPhi);
		  //
		  TVector3 P;
		  P(0) = Ptot * TMath::Sin(ThP)*TMath::Cos(PhP);
		  P(1) = Ptot * TMath::Sin(ThP)*TMath::Sin(PhP);
		  P(2) = Ptot * TMath::Cos(ThP);
		  Double_t Q = 1.0;
		  if (gRandom->Rndm() > 0.5)Q = -1.0;
		  tracks[i] = new ObsTrk(x, P, Q, Bfield, GC);
		  TVectorD par = tracks[i]->GetObsPar();
		  cout << "i = " << i << ", loading par = " << endl; par.Print();
	  }
	  //Double_t xa[3];
	  //x.GetXYZ(xa);
	  TVectorD xvtx(3); xvtx.Zero();
	  cout << "xvtx = " << endl; xvtx.Print();
	  TMatrixDSym covX(3);
	  Double_t Chi2 = Vertex(Ntr, tracks, xvtx, covX);
	  //
	  cout << "Fit vertex: x = " << xvtx(0) << ", y = " << xvtx(1) << ", z = " << xvtx(2) << endl;
	  Double_t PullX = (xvtx(0) - x(0)) / TMath::Sqrt(covX(0, 0));
	  hXpull->Fill(PullX);
	  Double_t PullY = (xvtx(1) - x(1)) / TMath::Sqrt(covX(1, 1));
	  hYpull->Fill(PullY);
	  Double_t PullZ = (xvtx(2) - x(2)) / TMath::Sqrt(covX(2, 2));
	  hZpull->Fill(PullZ);
	  hChi2->Fill(Chi2);
  }
  //
  // Plots
  TCanvas *cnv = new TCanvas("cnv", "Normalized residuals", 20, 20, 500, 500);
  cnv->Divide(2, 2);
  cnv->cd(1);
  hXpull->Fit("gaus");
  hXpull->Draw();
  cnv->cd(2);
  hYpull->Fit("gaus");
  hYpull->Draw();
  cnv->cd(3);
  hZpull->Fit("gaus");
  hZpull->Draw();
  cnv->cd(4);
  hChi2->Fit("fch");
  hChi2->Draw();
  cnv->SaveAs("vtx_residuals.pdf"); //Jim
}
