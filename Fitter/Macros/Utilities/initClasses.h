#ifndef initClasses_h
#define initClasses_h

#include "TSystem.h"
#include "TROOT.h"
#include "TDirectory.h"
#include "TSystem.h"

#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
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
  MinMax M, Pt, AbsRap;
} DiMuonPar;

typedef struct SiMuonPar {
  MinMax Pt, Eta;
} SiMuonPar;

typedef struct InputOpt {
  int        oniaMode;
  bool       isData, doFit, inExcStat, doSimulFit;
  EvtPar     PbPb, pp;
} InputOpt;

typedef struct KinCuts {
  StartEnd   Centrality;
  SiMuonPar  sMuon;
  DiMuonPar  dMuon;
} KinCuts;


struct ParticleMass { double JPsi, Psi2S, Y1S, Y2S, Y3S, Z; };
ParticleMass Mass = {3.096, 3.686, 9.460, 10.023, 10.355, 91.188};

#endif // #ifndef initClasses_h
