#ifndef initClasses_h
#define initClasses_h

#include "TSystem.h"
#include "TROOT.h"
#include "TDirectory.h"
#include "TSystem.h"

#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TF1.h"
#include "TProfile.h"
#include "TVectorD.h"
#include "TCanvas.h"
#include "TPaveText.h"
#include "TPave.h"
#include "TPad.h"
#include "TFrame.h"
#include "TAxis.h"
#include "TLegend.h"

#include "RooWorkspace.h"
#include "RooChi2Var.h"
#include "RooDataSet.h"
#include "RooDataHist.h"
#include "RooHist.h"
#include "RooPlot.h"
#include "RooPlotable.h"
#include "RooFitResult.h"
#include "RooCategory.h"
#include "RooSimultaneous.h"
#include "RooArgSet.h"
#include "RooArgList.h"
#include "RooConstVar.h"
#include "RooRealVar.h"
#include "RooFormulaVar.h"
#include "RooGenericPdf.h"
#include "RooAddPdf.h"
#include "RooProdPdf.h"
#include "RooAbsPdf.h"
#include "RooCBShape.h"
#include "RooGaussian.h"
#include "RooChebychev.h"
#include "RooPolynomial.h"
#include "RooExponential.h"

#include "../CMS/tdrstyle.C"
#include "../CMS/CMS_lumi.C"

#include <ctime>
#include <iostream>
#include <sstream>

using namespace std;
using namespace RooFit;

typedef struct StartEnd {
  int Start, End;
} StartEnd;

typedef struct MinMax {
  double Min, Max;
} MinMax;

typedef struct EvtPar {
  StartEnd RunNb;
  int      TriggerBit;
} EvtPar;

typedef struct DiMuonPar {
  MinMax ctau, ctauErr, M, Pt, AbsRap;
} DiMuonPar;

typedef struct SiMuonPar {
  MinMax Pt, Eta;
} SiMuonPar;

typedef struct InputOpt {
  int        oniaMode;
  bool       isData, do2DFit, inExcStat, doSimulFit;
  EvtPar     PbPb, pp;
} InputOpt;

typedef struct KinCuts {
  StartEnd   Centrality;
  SiMuonPar  sMuon;
  DiMuonPar  dMuon;
} KinCuts;


struct ParticleMass { double JPsi, Psi2S, Y1S, Y2S, Y3S, Z; };
ParticleMass Mass = {3.096, 3.686, 9.460, 10.023, 10.355, 91.188};


enum class MassModel 
  { 
    SingleGaussian, DoubleGaussian, SingleCrystalBall, DoubleCrystalBall, GaussianAndCrystalBall, 
    FirstOrderPolynomial, SecondOrderPolynomial, ThirdOrderPolynomial, FourthOrderPolynomial, 
    FirstOrderChebychev, SecondOrderChebychev, ThirdOrderChebychev, FourthOrderChebychev,
    Exponential
  };

enum class CtauModel 
  {     
    DoubleGaussianResolution, SingleGaussianResolution,
    TripleDecay, SingleSidedDecay, Delta
  };

typedef struct CtauPNP {
  CtauModel    Prompt, NonPrompt;
} CtauPNP;

typedef struct CtauMassModel {
  MassModel  Mass;
  CtauPNP    Ctau;
} FitModel;

typedef struct CharmModel {
  CtauMassModel  Jpsi, Psi2S, Bkg;
  CtauModel CtauRes;
} CharmModel;

typedef struct OniaModel {
  CharmModel  PbPb, PP;
} OniaModel;

typedef struct SignDataSets {
  RooDataSet* SS, OS;
} SignDataSets; 

const char* massModelName(MassModel massmodel) {
   if (massmodel == MassModel::SingleGaussian) return "SingleGaussian";
   if (massmodel == MassModel::DoubleGaussian) return "DoubleGaussian";
   if (massmodel == MassModel::SingleCrystalBall) return "SingleCrystalBall";
   if (massmodel == MassModel::DoubleCrystalBall) return "DoubleCrystalBall";
   if (massmodel == MassModel::GaussianAndCrystalBall) return "GaussianAndCrystalBall";
   if (massmodel == MassModel::FirstOrderPolynomial) return "FirstOrderPolynomial";
   if (massmodel == MassModel::SecondOrderPolynomial) return "SecondOrderPolynomial";
   if (massmodel == MassModel::ThirdOrderPolynomial) return "ThirdOrderPolynomial";
   if (massmodel == MassModel::FourthOrderPolynomial) return "FourthOrderPolynomial";
   if (massmodel == MassModel::FirstOrderChebychev) return "FirstOrderChebychev";
   if (massmodel == MassModel::SecondOrderChebychev) return "SecondOrderChebychev";
   if (massmodel == MassModel::ThirdOrderChebychev) return "ThirdOrderChebychev";
   if (massmodel == MassModel::FourthOrderChebychev) return "FourthOrderChebychev";
   if (massmodel == MassModel::Exponential) return "Exponential";
   return "undefined";
};

int massModelNpar(MassModel massmodel) {
   if (massmodel == MassModel::SingleGaussian) return 2;
   if (massmodel == MassModel::DoubleGaussian) return 4;
   if (massmodel == MassModel::SingleCrystalBall) return 4;
   if (massmodel == MassModel::DoubleCrystalBall) return 8;
   if (massmodel == MassModel::GaussianAndCrystalBall) return 6;
   if (massmodel == MassModel::FirstOrderPolynomial) return 1;
   if (massmodel == MassModel::SecondOrderPolynomial) return 2;
   if (massmodel == MassModel::ThirdOrderPolynomial) return 3;
   if (massmodel == MassModel::FourthOrderPolynomial) return 4;
   if (massmodel == MassModel::FirstOrderChebychev) return 1;
   if (massmodel == MassModel::SecondOrderChebychev) return 2;
   if (massmodel == MassModel::ThirdOrderChebychev) return 3;
   if (massmodel == MassModel::FourthOrderChebychev) return 4;
   if (massmodel == MassModel::Exponential) return 1;
   return 0;
};


#endif // #ifndef initClasses_h
