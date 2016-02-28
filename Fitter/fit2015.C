#ifndef __CINT__
#include "RooGlobalFunc.h"
#endif

#include "Macros/Utilities/initClasses.h"
#include "Macros/makeWorkspace2015.C"
#include "Macros/buildModelJpsi2015.C"
#include "Macros/fitMeanPTJpsi2015.C"
#include "Macros/drawPlot.C"

void SetOptions(struct InputOpt* opt, bool isData = true, int oniamode = 1, bool doFit = false, bool inExcStat = false, bool doSimulFit = false);

float rangeY_PbPb = 300;  float rangeY_PP   = 300; int nbins = 74;

void SetKinCuts(struct KinCuts* cut, int oniamode = 1, bool inExcStat= false) 
{
  cut->sMuon.Pt.Min  =  0.0;  cut->sMuon.Pt.Max  = 100000.0;
  cut->sMuon.Eta.Min = -2.4;  cut->sMuon.Eta.Max = 2.4;

  if (oniamode==1){        // Charmonia
    cut->dMuon.M.Min = inExcStat ? 2.2 : 2.6; 
    cut->dMuon.M.Max = inExcStat ? 4.2 : 3.5;  
 
    cut->dMuon.AbsRap.Min = 1.6;
    cut->dMuon.AbsRap.Max = 2.4;
    cut->dMuon.Pt.Min  =  6.5;
    cut->dMuon.Pt.Max  =  30.0;
    cut->Centrality.Start = 80;
    cut->Centrality.End = 200;

    rangeY_PbPb = 40000; //pp: 300000 PbPb: 10000
    rangeY_PP   = 40000; //pp: 300000 PbPb: 10000
    
  }
  
  return;
};


void fit2015(
             int  oniamode   = 1,        // oniamode-> 3: Z,  2: Upsilon and 1: J/Psi
             bool isData     = true,     // isData = false for MC, true for Data
             bool isPbPb     = true,     // isPbPb = false for pp, true for PbPb
	     bool doFit      = true,     // 
	     bool zoomPsi    = false,    // Zoom Psi(2S) peak on extra pad
	     bool incSS      = false,    // Include Same Sign data
	     bool getMeanPT  = true,     // Compute the mean PT
             bool inExcStat  = true,     // if inExcStat is true, then the excited states are fitted
	     bool doSimulFit = true      // Do simultaneous fit
             ) 
{

  if (inExcStat==false) { doSimulFit = false; }
  TreeName = TString("myTree");
  vector<TString> FileNames_PbPb; vector<TString> FileNames_PP;
  FileNames_PbPb.push_back(TString("/home/andre/MCSAMPLES/OniaTree_HIOniaL1DoubleMu0ABCD_HIRun2015-PromptReco-v1_Run_262548_263757_withOLDNEWCUT_withTRIG_MASSCUT.root"));
  FileNames_PP.push_back(TString("/home/andre/MCSAMPLES/OniaTree_DoubleMu_Run2015E-PromptReco-v1_Run_262157_262328_with2011CUT_withTRIG_MASSCUT.root"));
  
  struct InputOpt opt; SetOptions(&opt, isData, oniamode, doFit, inExcStat, doSimulFit);
  struct KinCuts  cut; SetKinCuts(&cut, oniamode, inExcStat);
 
  RooWorkspace myws;
  if (doSimulFit || isPbPb) {
    makeWorkspace2015(myws, FileNames_PbPb, &opt, &cut, true, nbins);
  }
  if (doSimulFit || !isPbPb) {
    makeWorkspace2015(myws, FileNames_PP, &opt, &cut, false, nbins);
  }
  
  if (doFit) {
    int sigModel=0, bkgModel=0;
    if (isData) {
      if (oniamode==1){
        sigModel = inExcStat ? 11 : 5;
        bkgModel = 4;
      }      
    } else {
      if (oniamode==1){
        sigModel = inExcStat ? 2 : 3; // gaussian   
        bkgModel = 2;
      }
    }

    if (doSimulFit) {

      if (opt.oniaMode==1) { 
	buildModelJpsi2015(myws, opt, sigModel, bkgModel, false);
      	buildModelJpsi2015(myws, opt, sigModel, bkgModel, true);
      }
  
      RooCategory* sample = new RooCategory("sample","sample");
      sample->defineType("PbPb"); sample->defineType("PP");
      RooDataSet* combData = new RooDataSet("combData","combined data", *myws.var("invariantMass"), Index(*sample), Import("PbPb", *((RooDataSet*)myws.data("dataOS_PbPb"))), Import("PP", *((RooDataSet*)myws.data("dataOS_PP"))));
      RooSimultaneous* simPdf = new RooSimultaneous("simPdf", "simultaneous pdf", *sample);
      simPdf->addPdf(*myws.pdf("pdf_PbPb"), "PbPb"); simPdf->addPdf(*myws.pdf("pdf_PP"), "PP"); 
      RooFitResult* fitMass = simPdf->fitTo(*combData, SumW2Error(kTRUE), Extended(kTRUE), Save(), NumCPU(8), Range("MassWindow")); 

      if (getMeanPT && opt.oniaMode==1) { fitMeanPTJpsi2015(myws, opt, cut, true, nbins); }
      drawPlot(myws, opt, cut, true, zoomPsi, incSS, getMeanPT, rangeY_PbPb, nbins);
      if (getMeanPT && opt.oniaMode==1) { fitMeanPTJpsi2015(myws, opt, cut, false, nbins); }
      drawPlot(myws, opt, cut, false, zoomPsi, incSS, getMeanPT, rangeY_PP, nbins);
 
    } else {
      if (opt.oniaMode==1) { buildModelJpsi2015(myws, opt, sigModel, bkgModel, isPbPb); }

      if (isPbPb) {
       	RooFitResult* fitMass = myws.pdf("pdf_PbPb")->fitTo(*myws.data("dataOS_PbPb"), SumW2Error(kTRUE), Extended(kTRUE), Save(), NumCPU(8));
	if (getMeanPT && opt.oniaMode==1) { fitMeanPTJpsi2015(myws, opt, cut, isPbPb, nbins); }
	drawPlot(myws, opt, cut, true, zoomPsi, incSS, getMeanPT, rangeY_PbPb, nbins);
      } else {
	RooFitResult* fitObject = myws.pdf("pdf_PP")->fitTo(*myws.data("dataOS_PP"), SumW2Error(kTRUE), Extended(kTRUE), Save(), NumCPU(8), Range("MassWindow"));
	if (getMeanPT && opt.oniaMode==1) { fitMeanPTJpsi2015(myws, opt, cut, isPbPb, nbins); }
	drawPlot(myws, opt, cut, false, zoomPsi, incSS, getMeanPT, rangeY_PP, nbins);
      }
    }	
  }

};

void SetOptions(struct InputOpt* opt, bool isData, int oniamode, bool doFit, bool inExcStat, bool doSimulFit) 
{
  opt->isData     = isData; 
  opt->oniaMode   = oniamode;
  opt->doFit      = doFit; 
  opt->inExcStat  = inExcStat; 
  opt->doSimulFit = doSimulFit;
  opt->pp.RunNb.Start = 262157; opt->PbPb.RunNb.Start = 262620;
  opt->pp.RunNb.End   = 262328; opt->PbPb.RunNb.End   = 263757;
  opt->pp.TriggerBit    = (int) PP::HLT_HIL1DoubleMu0_v1; 
  opt->PbPb.TriggerBit  = (int) HI::HLT_HIL1DoubleMu0_v1; 
  return;
};
