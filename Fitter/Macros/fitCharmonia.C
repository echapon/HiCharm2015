
#include "Utilities/initClasses.h"
#include "buildCharmoniaMassModel.C"
#include "drawMassPlot.C"

#include <algorithm>

bool importDataset(RooWorkspace& myws, RooWorkspace& inputWS, struct KinCuts cut, string label);
bool setModel( struct OniaModel& model, map<string, string>  parIni, bool isPbPb, bool inExcStat);

void setOptions(struct InputOpt* opt, bool inExcStat = false, bool doSimulFit = false);

bool isFitAlreadyFound(RooArgSet *newpars, string outputDir, string plotLabel, string TAG, struct KinCuts cut, bool isPbPb);
bool compareSnapshots(RooArgSet *pars1, const RooArgSet *pars2);

bool fitCharmonia( RooWorkspace&  inputWorkspace, // Workspace with all the input RooDatasets
		   struct KinCuts cut,            // Variable containing all kinematic cuts
		   map<string, string>  parIni,   // Variable containing all initial parameters
		   string outputDir,              // Path to output directory
		   string TAG,                    // Specifies the type of datasets: i.e, DATA, MCJPSINP, ...
		   bool isPbPb      = false,      // isPbPb = false for pp, true for PbPb
		   bool zoomPsi     = false,      // Zoom Psi(2S) peak on extra pad
		   bool setLogScale = true,       // Draw plot with log scale
		   bool incSS       = false,      // Include Same Sign data
		   bool getMeanPT   = false,      // Compute the mean PT
		   bool inExcStat   = false,      // if inExcStat is true, then the excited states are fitted
		   bool doSimulFit  = true,       // Do simultaneous fit
		   int nbins        = 74,         // Number of bins for fitting
                   int numCores     = 2           // Number of cores used for fitting
		   )  
{
  if (inExcStat==false) { 
    doSimulFit = false; 
    cut.dMuon.M.Min = 2.6;
    cut.dMuon.M.Max = 3.5;
  }

  struct InputOpt opt; setOptions(&opt, inExcStat, doSimulFit);
  
  int    numEntriesPbPb, numEntriesPP;
  string plotLabelPbPb,  plotLabelPP;

  struct OniaModel model;
  RooWorkspace     myws("workspace", "local workspace");

  if (doSimulFit || isPbPb) {
    // Set models based on initial parameters
    if (!setModel(model, parIni, true, inExcStat)) { return false; }

    // Import the local datasets
    string label = Form("%s_%s", TAG.c_str(), "PbPb");
    if (!importDataset(myws, inputWorkspace, cut, label)) { return false; }

    numEntriesPbPb = myws.data(Form("dOS_%s_PbPb", TAG.c_str()))->sumEntries();
    plotLabelPbPb = Form("Jpsi_%s_Bkg_%s", parIni["Model_Jpsi_PbPb"].c_str(), parIni["Model_Bkg_PbPb"].c_str());
  }
  if (doSimulFit || !isPbPb) {
    // Set models based on initial parameters
    if (!setModel(model, parIni, false, inExcStat)) { return false; }

    // Import the local datasets
    string label = Form("%s_%s", TAG.c_str(), "PP");
    if (!importDataset(myws, inputWorkspace, cut, label)) { return false; }

    numEntriesPP = myws.data(Form("dOS_%s_PP", TAG.c_str()))->sumEntries();
    plotLabelPP = Form("Jpsi_%s_Bkg_%s", parIni["Model_Jpsi_PP"].c_str(), parIni["Model_Bkg_PP"].c_str());
  }

 
  if (doSimulFit) { 
    // Build the Fit Model
    if (!buildCharmoniaMassModel(myws, opt, model.PP, parIni, false, numEntriesPP))  { return false; }
    if (!buildCharmoniaMassModel(myws, opt, model.PbPb, parIni, true, numEntriesPbPb)) { return false; }

    // check if we have already done this fit. If yes, do nothing and return true.
    bool found = true;
    RooArgSet *newpars = myws.pdf("pdfMASS_Tot_PbPb")->getParameters(*(myws.var("invMass")));
    found = found && isFitAlreadyFound(newpars, outputDir, plotLabelPbPb, TAG, cut, true);
    newpars = myws.pdf("pdfMASS_Tot_PP")->getParameters(*(myws.var("invMass")));
    found = found && isFitAlreadyFound(newpars, outputDir, plotLabelPP, TAG, cut, false);
    if (found) {
       cout << "[INFO] This fit was already done, so I'll just go to the next one." << endl;
       return true;
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
    RooFitResult* fitMass = simPdf->fitTo(*combData, SumW2Error(kTRUE), Extended(kTRUE), Save(), NumCPU(numCores), Range("MassWindow")); 

    // Draw the mass plots
    drawMassPlot(myws, outputDir, plotLabelPbPb, TAG, opt, cut, true, zoomPsi, setLogScale, incSS, getMeanPT, nbins);
    drawMassPlot(myws, outputDir, plotLabelPP, TAG, opt, cut, false, zoomPsi, setLogScale, incSS, getMeanPT, nbins);
    
  } else {
     if (isPbPb) {
       // Build the Fit Model
       if (!buildCharmoniaMassModel(myws, opt, model.PbPb, parIni, true, numEntriesPbPb)) { return false; }

       // check if we have already done this fit. If yes, do nothing and return true.
       RooArgSet *newpars = myws.pdf("pdfMASS_Tot_PbPb")->getParameters(*(myws.var("invMass")));
       bool found =  isFitAlreadyFound(newpars, outputDir, plotLabelPbPb, TAG, cut, true);
       if (found) {
         cout << "[INFO] This fit was already done, so I'll just go to the next one." << endl;
         return true;
       }

       // Fit the Datasets
       myws.pdf("pdfMASS_Tot_PbPb")->fitTo(*myws.data("dOS_DATA_PbPb"), SumW2Error(kTRUE), Extended(kTRUE), Range("MassWindow"), NumCPU(numCores));

       // Draw the mass plot
       drawMassPlot(myws, outputDir, plotLabelPbPb, TAG,  opt, cut, true, zoomPsi, setLogScale, incSS, getMeanPT, nbins);

     } else {
       // Build the Fit Model
       if (!buildCharmoniaMassModel(myws, opt, model.PP, parIni, false, numEntriesPP)) { return false; }

       // check if we have already done this fit. If yes, do nothing and return true.
       RooArgSet *newpars = myws.pdf("pdfMASS_Tot_PP")->getParameters(*(myws.var("invMass")));
       bool found =  isFitAlreadyFound(newpars, outputDir, plotLabelPP, TAG, cut, false);
       if (found) {
         cout << "[INFO] This fit was already done, so I'll just go to the next one." << endl;
         return true;
       }

       // Fit the Datasets
       myws.pdf("pdfMASS_Tot_PP")->fitTo(*myws.data("dOS_DATA_PP"), SumW2Error(kTRUE), Extended(kTRUE), Save(), NumCPU(numCores), Range("MassWindow"));

       // Draw the mass plot
       drawMassPlot(myws, outputDir, plotLabelPP, TAG, opt, cut, false, zoomPsi, setLogScale, incSS, getMeanPT, nbins);
     }
  }

  return true;
};

bool setModel( struct OniaModel& model, map<string, string>  parIni, bool isPbPb, bool inExcStat)
{
  if (isPbPb) {
    if (parIni.count("Model_Bkg_PbPb")>0) {
      model.PbPb.Bkg.Mass = MassModelDictionary[parIni["Model_Bkg_PbPb"]];
    } else { 
      cout << "[ERROR] Background mass model for PbPb was not found in the initial parameters!" << endl; return false;
    }
    if (parIni.count("Model_Jpsi_PbPb")>0) {
      model.PbPb.Jpsi.Mass = MassModelDictionary[parIni["Model_Jpsi_PbPb"]];
    } else { 
      cout << "[ERROR] Jpsi mass model for PbPb was not found in the initial parameters!" << endl; return false;
    }
    if (inExcStat) {
      if (parIni.count("Model_Psi2S_PbPb")>0) {
	model.PbPb.Psi2S.Mass = MassModelDictionary[parIni["Model_Psi2S_PbPb"]];
      } else { 
	cout << "[ERROR] psi(2S) mass model for PbPb was not found in the initial parameters!" << endl; return false;
      }
    }
  } else {
    if (parIni.count("Model_Bkg_PP")>0) {
      model.PP.Bkg.Mass = MassModelDictionary[parIni["Model_Bkg_PP"]];
    } else { 
      cout << "[ERROR] Background mass model for PP was not found in the initial parameters!" << endl; return false;
    }
    if (parIni.count("Model_Jpsi_PP")>0) {
      model.PP.Jpsi.Mass = MassModelDictionary[parIni["Model_Jpsi_PP"]];
    } else { 
      cout << "[ERROR] Jpsi mass model for PP was not found in the initial parameters!" << endl; return false;
    }
    if (inExcStat) {
      if (parIni.count("Model_Psi2S_PP")>0) {
	model.PP.Psi2S.Mass = MassModelDictionary[parIni["Model_Psi2S_PP"]];
      } else { 
	cout << "[ERROR] psi(2S) mass model for PP was not found in the initial parameters!" << endl; return false;
      }
    }
  }

  return true;
}   

bool importDataset(RooWorkspace& myws, RooWorkspace& inputWS, struct KinCuts cut, string label)
{
  string indMuonMass    = Form("(%.6f < invMass && invMass < %.6f)",       cut.dMuon.M.Min,       cut.dMuon.M.Max);
  string indMuonRap     = Form("(%.6f < abs(rap) && abs(rap) < %.6f)",     cut.dMuon.AbsRap.Min,  cut.dMuon.AbsRap.Max);
  string indMuonPt      = Form("(%.6f < pt && pt < %.6f)",                 cut.dMuon.Pt.Min,      cut.dMuon.Pt.Max);
  string indMuonCtau    = Form("(%.6f < ctau && ctau < %.6f)",             cut.dMuon.ctau.Min,    cut.dMuon.ctau.Max);
  string indMuonCtauErr = Form("(%.6f < ctauErr && ctauErr < %.6f)",       cut.dMuon.ctauErr.Min, cut.dMuon.ctauErr.Max);
  string inCentrality   = Form("(%d <= cent && cent <= %d)",               cut.Centrality.Start,  cut.Centrality.End);

  string strCut         = indMuonMass +"&&"+ indMuonRap +"&&"+ indMuonPt +"&&"+ indMuonCtau +"&&"+ indMuonCtauErr;
  if (label.find("PbPb")!=std::string::npos){ strCut = strCut +"&&"+ inCentrality; } 

  // Reduce and import the datasets
  if (!(inputWS.data(Form("dOS_%s", label.c_str())))){ 
    cout << "[ERROR] The dataset " <<  Form("dOS_%s", label.c_str()) << " was not found!" << endl;
    return false;
  }
  RooDataSet* dataOS = (RooDataSet*)inputWS.data(Form("dOS_%s", label.c_str()))->reduce(strCut.c_str());
  if (dataOS->sumEntries()==0){ 
    cout << "[ERROR] No events from dataset " <<  Form("dOS_%s", label.c_str()) << " passed the kinematic cuts!" << endl;
    return false;
  }  
  myws.import(*dataOS);

  if (!(inputWS.data(Form("dSS_%s", label.c_str())))){ 
    cout << "[ERROR] The dataset " <<  Form("dSS_%s", label.c_str()) << " was not found!" << endl;
    return false;
  } 
  RooDataSet* dataSS = (RooDataSet*)inputWS.data(Form("dSS_%s", label.c_str()))->reduce(strCut.c_str());
  if (dataSS->sumEntries()==0){ 
    cout << "[WARNING] No events from dataset " <<  Form("dSS_%s", label.c_str()) << " passed the kinematic cuts!" << endl;
  }
  myws.import(*dataSS);

  // Set the range of each global parameter in the local workspace
  myws.var("invMass")->setMin(cut.dMuon.M.Min);        
  myws.var("invMass")->setMax(cut.dMuon.M.Max);
  myws.var("pt")->setMin(cut.dMuon.Pt.Min);            
  myws.var("pt")->setMax(cut.dMuon.Pt.Max);
  myws.var("rap")->setMin(cut.dMuon.AbsRap.Min);       
  myws.var("rap")->setMax(cut.dMuon.AbsRap.Max);
  myws.var("ctau")->setMin(cut.dMuon.ctau.Min);        
  myws.var("ctau")->setMax(cut.dMuon.ctau.Max);
  myws.var("ctauErr")->setMin(cut.dMuon.ctauErr.Min);  
  myws.var("ctauErr")->setMax(cut.dMuon.ctauErr.Max);
  if (label.find("PbPb")!=std::string::npos){
    myws.var("cent")->setMin(cut.Centrality.Start);      
    myws.var("cent")->setMax(cut.Centrality.End);
  }

  myws.var("invMass")->setRange("MassWindow", cut.dMuon.M.Min, cut.dMuon.M.Max);

  return true;
};

void setOptions(struct InputOpt* opt, bool inExcStat, bool doSimulFit) 
{
  opt->inExcStat        = inExcStat; 
  opt->doSimulFit       = doSimulFit;
  opt->pp.RunNb.Start   = 262157; opt->PbPb.RunNb.Start = 262620;
  opt->pp.RunNb.End     = 262328; opt->PbPb.RunNb.End   = 263757;
  opt->pp.TriggerBit    = (int) PP::HLT_HIL1DoubleMu0_v1; 
  opt->PbPb.TriggerBit  = (int) HI::HLT_HIL1DoubleMu0_v1; 
  return;
};

bool isFitAlreadyFound(RooArgSet *newpars, string outputDir, string plotLabel, string TAG, struct KinCuts cut, bool isPbPb) {
  TFile *file = new TFile(Form("%sresult/%s/FIT_%s_%s_%s_%s_pt%.0f%.0f_rap%.0f%.0f_cent%d%d.root", outputDir.c_str(), TAG.c_str(), TAG.c_str(), "Psi2SJpsi", (isPbPb?"PbPb":"PP"), plotLabel.c_str(), (cut.dMuon.Pt.Min*10.0), (cut.dMuon.Pt.Max*10.0), (cut.dMuon.AbsRap.Min*10.0), (cut.dMuon.AbsRap.Max*10.0), cut.Centrality.Start, cut.Centrality.End));
  if (!file) return false;

  RooWorkspace *ws = (RooWorkspace*) file->Get("workspace");
  if (!ws) return false;

  const RooArgSet *params = ws->getSnapshot(Form("pdfMASS_Tot_%s_parIni", (isPbPb?"PbPb":"PP")));
  if (!params) return false;
  
  return compareSnapshots(newpars, params);
}

bool compareSnapshots(RooArgSet *pars1, const RooArgSet *pars2) {
  TIterator* parIt = pars1->createIterator(); 

  for (RooRealVar* it = (RooRealVar*)parIt->Next(); it!=NULL; it = (RooRealVar*)parIt->Next() ) {
    double val = pars2->getRealValue(it->GetName(),-1e99);
    if (val==-1e99) return false;          // the parameter was not found!
    if (val != it->getVal()) return false; // the parameter was found, but with a different value!
    if ( ((RooRealVar&)(*pars2)[it->GetName()]).getMin() != it->getMin() ) return false; // the parameter has different lower limit
    if ( ((RooRealVar&)(*pars2)[it->GetName()]).getMax() != it->getMax() ) return false; // the parameter has different upper limit
  }

  return true;
}
