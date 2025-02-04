#include <TString.h>
#include <TROOT.h>

void LoadAll(TString dname)
{
gROOT->Reset();
//TString Action = ".L geometry_scripts/SolGeom" + dname + ".cxx+";
TString Action = ".L geometry_scripts/SolGeom.cxx+";
gROOT->ProcessLine(Action);
gROOT->ProcessLine(".L trkcovariance_scripts/SolTrack.cxx+");
gROOT->ProcessLine(".L trkcovariance_scripts/SolGridCov.cxx+");
gROOT->ProcessLine(".L trkcovariance_scripts/ObsTrk.cxx+");
}
