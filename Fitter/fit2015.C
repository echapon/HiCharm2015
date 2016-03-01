#ifndef __CINT__
#include "RooGlobalFunc.h"
#endif

#include "Macros/Utilities/initClasses.h"
#include "Macros/makeWorkspace2015.C"
#include "Macros/buildCharmoniaCtauMassModel.C"
#include "Macros/fitMeanPTJpsi2015.C"
#include "Macros/drawMassPlot.C"

void SetOptions(struct InputOpt* opt, bool isData = true, int oniamode = 1, bool doFit = false, bool inExcStat = false, bool doSimulFit = false);

float rangeY_PbPb = 300;  float rangeY_PP   = 300; int nbins = 74;

void SetKinCuts(struct KinCuts* cut, int oniamode = 1, bool inExcStat= false) 
{
  cut->sMuon.Pt.Min  =  0.0;  cut->sMuon.Pt.Max  = 100000.0;
  cut->sMuon.Eta.Min = -2.4;  cut->sMuon.Eta.Max = 2.4;
  cut->dMuon.ctau.Min = -3.0; cut->dMuon.ctau.Max = 5.0;  
  cut->dMuon.ctauErr.Min = 0.0; cut->dMuon.ctauErr.Max = 0.5;  

  if (oniamode==1){        // Charmonia
    cut->dMuon.M.Min = inExcStat ? 2.2 : 2.6; 
    cut->dMuon.M.Max = inExcStat ? 4.5 : 3.5;  
 
    cut->dMuon.AbsRap.Min = 1.6;
    cut->dMuon.AbsRap.Max = 2.4;
    cut->dMuon.Pt.Min  =  6.5;
    cut->dMuon.Pt.Max  =  30.0;
    cut->Centrality.Start = 0;
    cut->Centrality.End = 200;

    rangeY_PbPb = 40000; //pp: 300000 PbPb: 10000
    rangeY_PP   = 40000; //pp: 300000 PbPb: 10000
    
  }
  
  return;
};


void fit2015(
             int  oniamode    = 1,        // oniamode-> 3: Z,  2: Upsilon and 1: J/Psi
             bool isData      = true,     // isData = false for MC, true for Data
             bool isPbPb      = false,     // isPbPb = false for pp, true for PbPb
	     bool zoomPsi     = false,    // Zoom Psi(2S) peak on extra pad
	     bool setLogScale = true,    // Draw plot with log scale
	     bool incSS       = false,    // Include Same Sign data
	     bool do2DFit     = false,     // 
	     bool getMeanPT   = false,     // Compute the mean PT
             bool inExcStat   = false,    // if inExcStat is true, then the excited states are fitted
	     bool doSimulFit  = true     // Do simultaneous fit
             ) 
{

  if (inExcStat==false) { doSimulFit = false; }
  struct map< string, vector<TString> > FileNames;
  FileNames["DATA_PbPb"].push_back(TString("root://eoscms//eos/cms/store/group/phys_heavyions/dileptons/Data2015/PbPb502TeV/TTrees/PromptAOD/OniaTree_HIOniaL1DoubleMu0_HIRun2015-PromptReco-v1_Run_262620_263757.root"));
  FileNames["DATA_PbPb"].push_back(TString("root://eoscms//eos/cms/store/group/phys_heavyions/dileptons/Data2015/PbPb502TeV/TTrees/PromptAOD/OniaTree_HIOniaL1DoubleMu0B_HIRun2015-PromptReco-v1_Run_263322_263757.root"));
  FileNames["DATA_PbPb"].push_back(TString("root://eoscms//eos/cms/store/group/phys_heavyions/dileptons/Data2015/PbPb502TeV/TTrees/PromptAOD/OniaTree_HIOniaL1DoubleMu0C_HIRun2015-PromptReco-v1_Run_263322_263757.root"));
  FileNames["DATA_PbPb"].push_back(TString("root://eoscms//eos/cms/store/group/phys_heavyions/dileptons/Data2015/PbPb502TeV/TTrees/PromptAOD/OniaTree_HIOniaL1DoubleMu0D_HIRun2015-PromptReco-v1_Run_263322_263757.root"));
  FileNames["DATA_PP"].push_back(TString("root://eoscms//eos/cms/store/group/phys_heavyions/dileptons/Data2015/pp502TeV/TTrees/PromptAOD/OniaTree_DoubleMu_Run2015E-PromptReco-v1_Run_262157_262328.root"));
  
  struct InputOpt opt; SetOptions(&opt, isData, oniamode, do2DFit, inExcStat, doSimulFit);
  struct KinCuts  cut; SetKinCuts(&cut, oniamode, inExcStat);
 
  RooWorkspace myws;
  if (doSimulFit || isPbPb) {
    makeWorkspace2015(myws, FileNames["DATA_PbPb"], &opt, &cut, nbins, true, true);
  }
  if (doSimulFit || !isPbPb) {
    makeWorkspace2015(myws, FileNames["DATA_PP"], &opt, &cut, nbins, false, true);
  }
  myws.var("invMass")->setRange("MassWindow", cut.dMuon.M.Min, cut.dMuon.M.Max);
  
  struct OniaModel model;
  if (oniamode==1){
    model.PbPb.Jpsi.Mass   = MassModel::SingleCrystalBall;
    model.PbPb.Psi2S.Mass  = MassModel::SingleCrystalBall;
    model.PbPb.Bkg.Mass    = MassModel::SecondOrderChebychev;
    model.PP.Jpsi.Mass     = MassModel::DoubleCrystalBall;
    model.PP.Psi2S.Mass    = MassModel::DoubleCrystalBall;
    model.PP.Bkg.Mass      = MassModel::SecondOrderChebychev;
    if (do2DFit) {
      model.PbPb.CtauRes              = CtauModel::DoubleGaussianResolution;
      model.PbPb.Jpsi.Ctau.Prompt     = CtauModel::Delta;  
      model.PbPb.Psi2S.Ctau.Prompt    = CtauModel::Delta;
      model.PbPb.Bkg.Ctau.Prompt      = CtauModel::Delta;
      model.PbPb.Jpsi.Ctau.NonPrompt  = CtauModel::SingleSidedDecay;
      model.PbPb.Psi2S.Ctau.NonPrompt = CtauModel::SingleSidedDecay;
      model.PbPb.Bkg.Ctau.NonPrompt   = CtauModel::TripleDecay;
      model.PP.CtauRes                = CtauModel::DoubleGaussianResolution;
      model.PP.Jpsi.Ctau.Prompt       = CtauModel::Delta;  
      model.PP.Psi2S.Ctau.Prompt      = CtauModel::Delta;
      model.PP.Bkg.Ctau.Prompt        = CtauModel::Delta;
      model.PP.Jpsi.Ctau.NonPrompt    = CtauModel::SingleSidedDecay;
      model.PP.Psi2S.Ctau.NonPrompt   = CtauModel::SingleSidedDecay;
      model.PP.Bkg.Ctau.NonPrompt     = CtauModel::TripleDecay;
    }
  }      

  if (doSimulFit) { 
    // Build the Fit Model
    if (opt.oniaMode==1) { 
      if (!buildCharmoniaCtauMassModel(myws, opt, model.PP, false, do2DFit)) { return; }
      if (!buildCharmoniaCtauMassModel(myws, opt, model.PbPb, true, do2DFit))  { return; }
    }
    // Create the combided datasets and models
    RooCategory* sample = new RooCategory("sample","sample"); sample->defineType("PbPb"); sample->defineType("PP");
    RooDataSet*  combData = new RooDataSet("combData","combined data", *myws.var("invMass"), Index(*sample), 
					   Import("PbPb", *((RooDataSet*)myws.data("dOS_DATA_PbPb"))), 
					   Import("PP",   *((RooDataSet*)myws.data("dOS_DATA_PP")))
					   );
    RooSimultaneous* simPdf = new RooSimultaneous("simPdf", "simultaneous pdf", *sample);
    simPdf->addPdf(*myws.pdf("pdfMASS_Tot_PbPb"), "PbPb"); simPdf->addPdf(*myws.pdf("pdfMASS_Tot_PP"), "PP"); 
    // Do the simultaneous fit
    RooFitResult* fitMass = simPdf->fitTo(*combData, SumW2Error(kTRUE), Extended(kTRUE), Save(), NumCPU(8), Range("MassWindow"),NormRange("MassWindow")); 
    // Compute Mean PT 
    if (getMeanPT && opt.oniaMode==1) { fitMeanPTJpsi2015(myws, opt, cut, true, nbins); }
    if (getMeanPT && opt.oniaMode==1) { fitMeanPTJpsi2015(myws, opt, cut, false, nbins); }
    // Draw the mass plots
    drawMassPlot(myws, opt, cut, true, zoomPsi, setLogScale, incSS, getMeanPT, rangeY_PbPb, nbins);
    drawMassPlot(myws, opt, cut, false, zoomPsi, setLogScale, incSS, getMeanPT, rangeY_PP, nbins);
 
  } else {
    if (isPbPb) {
      // Build the Fit Model
      if (opt.oniaMode==1) { if (!buildCharmoniaCtauMassModel(myws, opt, model.PbPb, true, do2DFit)) { return; } }
      // Fit the Datasets
      myws.pdf("pdfMASS_Tot_PbPb")->fitTo(*myws.data("dOS_DATA_PbPb"), SumW2Error(kTRUE), Extended(kTRUE), Range("MassWindow"),NormRange("MassWindow"), NumCPU(8));
      if (getMeanPT && opt.oniaMode==1) { fitMeanPTJpsi2015(myws, opt, cut, isPbPb, nbins); }
      drawMassPlot(myws, opt, cut, true, zoomPsi, setLogScale, incSS, getMeanPT, rangeY_PbPb, nbins);
    } else {
      // Build the Fit Model
      if (opt.oniaMode==1) { if (!buildCharmoniaCtauMassModel(myws, opt, model.PP, false, do2DFit)) { return; } }
      // Fit the Datasets
      myws.pdf("pdfMASS_Tot_PP")->fitTo(*myws.data("dOS_DATA_PP"), SumW2Error(kTRUE), Extended(kTRUE), Save(), NumCPU(8), Range("MassWindow"),NormRange("MassWindow"));
      if (getMeanPT && opt.oniaMode==1) { fitMeanPTJpsi2015(myws, opt, cut, isPbPb, nbins); }
      drawMassPlot(myws, opt, cut, false, zoomPsi, setLogScale, incSS, getMeanPT, rangeY_PP, nbins);
    }
  }	
};

void SetOptions(struct InputOpt* opt, bool isData, int oniamode, bool do2DFit, bool inExcStat, bool doSimulFit) 
{
  opt->isData     = isData; 
  opt->oniaMode   = oniamode;
  opt->inExcStat  = inExcStat; 
  opt->do2DFit    = do2DFit;
  opt->doSimulFit = doSimulFit;
  opt->pp.RunNb.Start = 262157; opt->PbPb.RunNb.Start = 262620;
  opt->pp.RunNb.End   = 262328; opt->PbPb.RunNb.End   = 263757;
  opt->pp.TriggerBit    = (int) PP::HLT_HIL1DoubleMu0_v1; 
  opt->PbPb.TriggerBit  = (int) HI::HLT_HIL1DoubleMu0_v1; 
  return;
};
