
#include "Utilities/initClasses.h"
#include "buildCharmoniaMassModel.C"
#include "drawMassPlot.C"

#include <algorithm>

void setCtauCuts(struct KinCuts& cut, bool isPbPb);
int importDataset(RooWorkspace& myws, RooWorkspace& inputWS, struct KinCuts cut, string label);
bool setModel( struct OniaModel& model, map<string, string>  parIni, bool isPbPb, bool incJpsi, bool incPsi2S, bool incBkg);

void setOptions(struct InputOpt* opt);

bool isFitAlreadyFound(RooArgSet *newpars, string outputDir, string plotLabel, string TAG, struct KinCuts cut, bool isPbPb, bool doSimulFit);
bool compareSnapshots(RooArgSet *pars1, const RooArgSet *pars2);

bool fitCharmonia( RooWorkspace&  inputWorkspace, // Workspace with all the input RooDatasets
		   struct KinCuts cut,            // Variable containing all kinematic cuts
		   map<string, string>  parIni,   // Variable containing all initial parameters
		   string outputDir,              // Path to output directory
                   // Select the type of datasets to fit
		   string DSTAG,                  // Specifies the type of datasets: i.e, DATA, MCJPSINP, ...
		   bool isPbPb      = false,      // isPbPb = false for pp, true for PbPb
                   // Select the type of object to fit
                   bool incJpsi     = true,       // Includes Jpsi model
                   bool incPsi2S    = true,       // Includes Psi(2S) model
                   bool incBkg      = true,       // Includes Background model
                   // Select the fitting options
                   bool cutCtau     = false,      // Apply prompt ctau cuts
                   bool doSimulFit  = false,      // Do simultaneous fit
                   bool wantPureSMC = false,      // Flag to indicate if we want to fit pure signal MC
                   int  numCores    = 2,          // Number of cores used for fitting
                   // Select the drawing options
                   bool setLogScale = true,       // Draw plot with log scale
                   bool incSS       = false,      // Include Same Sign data
                   bool zoomPsi     = false,      // Zoom Psi(2S) peak on extra pad
                   int  nBins       = 74,         // Number of bins used for plotting
                   bool getMeanPT   = false       // Compute the mean PT (NEED TO FIX)
		   )  
{

  // Define the mass range
  if (cut.dMuon.M.Max==5 && cut.dMuon.M.Min==2) { 
    // Default mass values, means that the user did not specify a mass range
    if ( incJpsi && !incPsi2S) {
      cut.dMuon.M.Min = 2.6;
      cut.dMuon.M.Max = 3.5;
    }
    else if ( !incJpsi && incPsi2S) {
      cut.dMuon.M.Min = 3.0;
      cut.dMuon.M.Max = 4.1;
    }
    else {
      cut.dMuon.M.Min = 2.2;
      cut.dMuon.M.Max = 4.5;
    }
  }
  parIni["invMassNorm"] = Form("RooFormulaVar::%s('( -1.0 + 2.0*( @0 - @1 )/( @2 - @1) )', {%s, mMin[%.6f], mMax[%.6f]})", "invMassNorm", "invMass", cut.dMuon.M.Min, cut.dMuon.M.Max );
  // Apply the ctau cuts to reject non-prompt charmonia
  if (cutCtau) { setCtauCuts(cut, isPbPb); }  

  


  // Check if input dataset is MC
  bool isMC = false;
  if (DSTAG.find("MC")!=std::string::npos) {
    if (incJpsi && incPsi2S) { 
      cout << "[ERROR] We can only fit one type of signal using MC" << endl; return false; 
    }
    isMC = true;
  }
  if (isMC && wantPureSMC) wantPureSMC=true;
  else wantPureSMC=false;

  struct InputOpt opt; setOptions(&opt);
  
  string plotLabelPbPb,  plotLabelPP;

  struct OniaModel model;
  RooWorkspace     myws("workspace", "local workspace");

  bool doFit = true;
  if (doSimulFit || !isPbPb) {
    
    // Set models based on initial parameters
    if (!setModel(model, parIni, false, incJpsi, incPsi2S, incBkg)) { return false; }
    
    // Import the local datasets
    string label = Form("%s_%s", DSTAG.c_str(), "PP");
    if (wantPureSMC) label = Form("%s_%s_NoBkg", DSTAG.c_str(), "PP");
    string dsName = Form("dOS_%s", label.c_str());
    int importID = importDataset(myws, inputWorkspace, cut, label);
    if (importID<0) { return false; }
    else if (importID==0) { doFit = false; }
    
    // Build the Fit Model    
    double numEntries = myws.data(dsName.c_str())->sumEntries();
    if (!buildCharmoniaMassModel(myws, model.PP, parIni, false, doSimulFit, incBkg, incJpsi, incPsi2S, numEntries))  { return false; }

    if (incJpsi)  { plotLabelPP = plotLabelPP + Form("_Jpsi_%s", parIni["Model_Jpsi_PP"].c_str());   } 
    if (incPsi2S) { plotLabelPP = plotLabelPP + Form("_Psi2S_%s", parIni["Model_Psi2S_PP"].c_str()); }
    if (incBkg)   { plotLabelPP = plotLabelPP + Form("_Bkg_%s", parIni["Model_Bkg_PP"].c_str());     }
  }
  if (doSimulFit || isPbPb) {
    
    // Set models based on initial parameters
    if (!setModel(model, parIni, true, incJpsi, incPsi2S, incBkg)) { return false; }
    
    // Import the local datasets
    string label = Form("%s_%s", DSTAG.c_str(), "PbPb");
    if (wantPureSMC) label = Form("%s_%s_NoBkg", DSTAG.c_str(), "PbPb");
    string dsName = Form("dOS_%s", label.c_str());
    int importID = importDataset(myws, inputWorkspace, cut, label);
    if (importID<0) { return false; }
    else if (importID==0) { doFit = false; }
    
    // Build the Fit Model
    double    numEntries = myws.data(dsName.c_str())->sumEntries();
    if (!buildCharmoniaMassModel(myws, model.PbPb, parIni, true, doSimulFit, incBkg, incJpsi, incPsi2S, numEntries)) { return false; }

    if (incJpsi)  { plotLabelPbPb = plotLabelPbPb + Form("_Jpsi_%s", parIni["Model_Jpsi_PbPb"].c_str());   } 
    if (incPsi2S) { plotLabelPbPb = plotLabelPbPb + Form("_Psi2S_%s", parIni["Model_Psi2S_PbPb"].c_str()); }
    if (incBkg)   { plotLabelPbPb = plotLabelPbPb + Form("_Bkg_%s", parIni["Model_Bkg_PbPb"].c_str());     }
  }
  
  if (doFit)
  {
    if (doSimulFit) {
      // check if we have already done this fit. If yes, do nothing and return true.
      bool found = true;
      RooArgSet *newpars = myws.pdf("pdfMASS_Tot_PbPb")->getParameters(*(myws.var("invMass")));
      found = found && isFitAlreadyFound(newpars, outputDir, plotLabelPbPb, DSTAG, cut, true, true);
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
      RooFitResult* fitResult = simPdf->fitTo(*combData, Offset(kTRUE), Extended(kTRUE), NumCPU(numCores), Range("MassWindow"), Save(), Minimizer("Minuit2","Migrad"));
      fitResult->Print();

      // Create the output files
      drawMassPlot(myws, outputDir, opt, cut, plotLabelPbPb, DSTAG, true, incJpsi, incPsi2S, incBkg, cutCtau, doSimulFit, false, setLogScale, incSS, zoomPsi, nBins, getMeanPT);
      drawMassPlot(myws, outputDir, opt, cut, plotLabelPP, DSTAG, false, incJpsi, incPsi2S, incBkg, cutCtau, doSimulFit, false, setLogScale, incSS, zoomPsi, nBins, getMeanPT);
      
      // Delete the objects used during the simultaneous fit
      delete sample; delete combData; delete simPdf;
      
    }
    else {
      if (isPbPb) {
        
        string pdfName = "pdfMASS_Tot_PbPb";
        string dsName = Form("dOS_%s_PbPb", DSTAG.c_str());
        if (wantPureSMC) dsName = Form("dOS_%s_PbPb_NoBkg", DSTAG.c_str());

        // check if we have already done this fit. If yes, do nothing and return true.
        RooArgSet *newpars = myws.pdf(pdfName.c_str())->getParameters(*(myws.var("invMass")));
        bool found =  isFitAlreadyFound(newpars, outputDir, wantPureSMC ? (plotLabelPbPb+"_NoBkg") : plotLabelPbPb, DSTAG, cut, true, false);
        if (found) {
          cout << "[INFO] This fit was already done, so I'll just go to the next one." << endl;
          return true;
        }

        bool isWeighted = myws.data(dsName.c_str())->isWeighted();
        
        // Fit the Datasets
        if (incJpsi || incPsi2S) {
          if (isWeighted) {
            RooFitResult* fitResult = myws.pdf(pdfName.c_str())->fitTo(*myws.data(dsName.c_str()), Extended(kTRUE), SumW2Error(kTRUE), Range("MassWindow"), NumCPU(numCores), Save());
            fitResult->Print();
          } else {
            RooFitResult* fitResult = myws.pdf(pdfName.c_str())->fitTo(*myws.data(dsName.c_str()), Extended(kTRUE), Range("MassWindow"), NumCPU(numCores), Save());
            fitResult->Print();
          }  
        } else {
          RooFitResult* fitResult = myws.pdf(pdfName.c_str())->fitTo(*myws.data(dsName.c_str()), Extended(kTRUE), Range("SideBand1,SideBand2"), NumCPU(numCores), Save());
          fitResult->Print();
        }
        
        // Create the output files
        drawMassPlot(myws, outputDir, opt, cut, wantPureSMC ? (plotLabelPbPb+"_NoBkg") : plotLabelPbPb, DSTAG, true, incJpsi, incPsi2S, incBkg, cutCtau, doSimulFit, wantPureSMC, setLogScale, incSS, zoomPsi, nBins, getMeanPT);
      }
      else {
        cut.Centrality.Start = 0;
        cut.Centrality.End = 200;
        
        string pdfName = "pdfMASS_Tot_PP";
        string dsName = Form("dOS_%s_PP", DSTAG.c_str());
        if (wantPureSMC) dsName = Form("dOS_%s_PP_NoBkg", DSTAG.c_str());
        
        // check if we have already done this fit. If yes, do nothing and return true.
        RooArgSet *newpars = myws.pdf(pdfName.c_str())->getParameters(*(myws.var("invMass")));
        bool found =  isFitAlreadyFound(newpars, outputDir, wantPureSMC ? (plotLabelPP+"_NoBkg") : plotLabelPP, DSTAG, cut, false, false);
        if (found) {
          cout << "[INFO] This fit was already done, so I'll just go to the next one." << endl;
          return true;
        }

        bool isWeighted = myws.data(dsName.c_str())->isWeighted();
        
        // Fit the Datasets
        if (incJpsi || incPsi2S) {
          RooFitResult* fitResult = myws.pdf(pdfName.c_str())->fitTo(*myws.data(dsName.c_str()), Extended(kTRUE), Range("MassWindow"), NumCPU(numCores), Save());
          fitResult->Print(); 
        } else {
          RooFitResult* fitResult = myws.pdf(pdfName.c_str())->fitTo(*myws.data(dsName.c_str()), Extended(kTRUE), Range("SideBand1,SideBand2"), NumCPU(numCores), Save());
          fitResult->Print();
        }
        
        // Draw the mass plot
        drawMassPlot(myws, outputDir, opt, cut, wantPureSMC ? (plotLabelPP+"_NoBkg") : plotLabelPP, DSTAG, false, incJpsi, incPsi2S, incBkg, cutCtau, doSimulFit, wantPureSMC, setLogScale, incSS, zoomPsi, nBins, getMeanPT);
      }
    }
  }

  return true;
};


void setCtauCuts(struct KinCuts& cut, bool isPbPb) 
{
  if (cut.dMuon.AbsRap.Max<=1.6 && isPbPb) {
    cut.dMuon.ctau.Max = 0.03;
  }
  if (cut.dMuon.AbsRap.Min>=1.6 && isPbPb) {
    cut.dMuon.ctau.Max = 0.05;
  }
  if (cut.dMuon.AbsRap.Max<=1.6 && !isPbPb) {
    cut.dMuon.ctau.Max = 0.03;
  }
  if (cut.dMuon.AbsRap.Min>=1.6 && !isPbPb) {
    cut.dMuon.ctau.Max = 0.05;
  }
};


bool setModel( struct OniaModel& model, map<string, string>  parIni, bool isPbPb, bool incJpsi, bool incPsi2S, bool incBkg)
{
  if (isPbPb && incBkg) {
    if (parIni.count("Model_Bkg_PbPb")>0) {
      model.PbPb.Bkg.Mass = MassModelDictionary[parIni["Model_Bkg_PbPb"]];
      if (model.PbPb.Bkg.Mass==MassModel(0)) {
        cout << "[ERROR] The background model: " << parIni["Model_Bkg_PbPb"] << " is invalid" << endl; return false;
      }
    } else { 
      cout << "[ERROR] Background mass model for PbPb was not found in the initial parameters!" << endl; return false;
    }
  }
  if (isPbPb && incJpsi) {
    if (parIni.count("Model_Jpsi_PbPb")>0) {
      model.PbPb.Jpsi.Mass = MassModelDictionary[parIni["Model_Jpsi_PbPb"]];
      if (model.PbPb.Jpsi.Mass==MassModel(0)) {
        cout << "[ERROR] The Jpsi model: " << parIni["Model_Jpsi_PbPb"] << " is invalid" << endl; return false;
      }
    } else { 
      cout << "[ERROR] Jpsi mass model for PbPb was not found in the initial parameters!" << endl; return false;
    }
  }
  if (isPbPb && incPsi2S) {
    if (parIni.count("Model_Psi2S_PbPb")>0) {
      model.PbPb.Psi2S.Mass = MassModelDictionary[parIni["Model_Psi2S_PbPb"]];
      if (model.PbPb.Psi2S.Mass==MassModel(0)) {
        cout << "[ERROR] The psi2S model: " << parIni["Model_Psi2S_PbPb"] << " is invalid" << endl; return false;
      }
    } else { 
      cout << "[ERROR] psi(2S) mass model for PbPb was not found in the initial parameters!" << endl; return false;
    }
  }
  if (!isPbPb && incBkg) {
    if (parIni.count("Model_Bkg_PP")>0) {
      model.PP.Bkg.Mass = MassModelDictionary[parIni["Model_Bkg_PP"]];
      if (model.PP.Bkg.Mass==MassModel(0)) {
        cout << "[ERROR] The background model: " << parIni["Model_Bkg_PP"] << " is invalid" << endl; return false;
      }
    } else { 
      cout << "[ERROR] Background mass model for PP was not found in the initial parameters!" << endl; return false;
    }
  }
  if (!isPbPb && incJpsi) {
    if (parIni.count("Model_Jpsi_PP")>0) {
      model.PP.Jpsi.Mass = MassModelDictionary[parIni["Model_Jpsi_PP"]];
      if (model.PP.Jpsi.Mass==MassModel(0)) {
        cout << "[ERROR] The Jpsi model: " << parIni["Model_Jpsi_PbPb"] << " is invalid" << endl; return false;
      }
    } else { 
      cout << "[ERROR] Jpsi mass model for PP was not found in the initial parameters!" << endl; return false;
    }
  }
  if (!isPbPb && incPsi2S) {
    if (parIni.count("Model_Psi2S_PP")>0) {
      model.PP.Psi2S.Mass = MassModelDictionary[parIni["Model_Psi2S_PP"]];
      if (model.PP.Psi2S.Mass==MassModel(0)) {
        cout << "[ERROR] The psi2S model: " << parIni["Model_Psi2S_PbPb"] << " is invalid" << endl; return false;
      }
    } else { 
      cout << "[ERROR] psi(2S) mass model for PP was not found in the initial parameters!" << endl; return false;
    }
  }

  return true;
};
   

int importDataset(RooWorkspace& myws, RooWorkspace& inputWS, struct KinCuts cut, string label)
{
  string indMuonMass    = Form("(%.6f < invMass && invMass < %.6f)",       cut.dMuon.M.Min,       cut.dMuon.M.Max);
  string indMuonRap     = Form("(%.6f <= abs(rap) && abs(rap) < %.6f)",    cut.dMuon.AbsRap.Min,  cut.dMuon.AbsRap.Max);
  string indMuonPt      = Form("(%.6f <= pt && pt < %.6f)",                cut.dMuon.Pt.Min,      cut.dMuon.Pt.Max);
  string indMuonCtau    = Form("(%.6f < ctau && ctau < %.6f)",             cut.dMuon.ctau.Min,    cut.dMuon.ctau.Max);
  string indMuonCtauErr = Form("(%.6f < ctauErr && ctauErr < %.6f)",       cut.dMuon.ctauErr.Min, cut.dMuon.ctauErr.Max);
  string inCentrality   = Form("(%d <= cent && cent < %d)",                cut.Centrality.Start,  cut.Centrality.End);

  string strCut         = indMuonMass +"&&"+ indMuonRap +"&&"+ indMuonPt +"&&"+ indMuonCtau +"&&"+ indMuonCtauErr;
  if (label.find("PbPb")!=std::string::npos){ strCut = strCut +"&&"+ inCentrality; } 

  // Reduce and import the datasets
  if (!(inputWS.data(Form("dOS_%s", label.c_str())))){ 
    cout << "[ERROR] The dataset " <<  Form("dOS_%s", label.c_str()) << " was not found!" << endl;
    return -1;
  }
  RooDataSet* dataOS = (RooDataSet*)inputWS.data(Form("dOS_%s", label.c_str()))->reduce(strCut.c_str());
  if (dataOS->sumEntries()==0){ 
    cout << "[WARNING] No events from dataset " <<  Form("dOS_%s", label.c_str()) << " passed the kinematic cuts!" << endl;
    return 0;
  }  
  myws.import(*dataOS);
  delete dataOS;

  if (label.find("NoBkg")==std::string::npos) // Don't try to find SS dataset if label contais NoBkg
  {
    if (!(inputWS.data(Form("dSS_%s", label.c_str())))){
      cout << "[ERROR] The dataset " <<  Form("dSS_%s", label.c_str()) << " was not found!" << endl;
      return -1;
    }
    RooDataSet* dataSS = (RooDataSet*)inputWS.data(Form("dSS_%s", label.c_str()))->reduce(strCut.c_str());
    if (dataSS->sumEntries()==0){
      cout << "[WARNING] No events from dataset " <<  Form("dSS_%s", label.c_str()) << " passed the kinematic cuts!" << endl;
    }
    myws.import(*dataSS);
    delete dataSS;
  }
  
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

  if (label.find("MC")!=std::string::npos)
  {
    if (label.find("PSI2S")!=std::string::npos)
    {
      if (cut.dMuon.AbsRap.Min >= 1.6) myws.var("invMass")->setRange("MassWindow", cut.dMuon.M.Min, 3.95);
      else myws.var("invMass")->setRange("MassWindow", cut.dMuon.M.Min, 3.85);
    }
    else
    {
      if (cut.dMuon.AbsRap.Min >= 1.6) myws.var("invMass")->setRange("MassWindow", cut.dMuon.M.Min, 3.32);
      else myws.var("invMass")->setRange("MassWindow", cut.dMuon.M.Min, 3.26);
    }
   
  }
  else myws.var("invMass")->setRange("MassWindow", cut.dMuon.M.Min, cut.dMuon.M.Max);
  
  if (cut.dMuon.M.Min<2.8) { myws.var("invMass")->setRange("SideBand1",  cut.dMuon.M.Min, 2.8); }
  if (cut.dMuon.M.Max>4.0) { myws.var("invMass")->setRange("SideBand2",  4.0, cut.dMuon.M.Max); }

  return 1;
};

void setOptions(struct InputOpt* opt) 
{
  opt->pp.RunNb.Start   = 262157; opt->PbPb.RunNb.Start = 262620;
  opt->pp.RunNb.End     = 262328; opt->PbPb.RunNb.End   = 263757;
  opt->pp.TriggerBit    = (int) PP::HLT_HIL1DoubleMu0_v1; 
  opt->PbPb.TriggerBit  = (int) HI::HLT_HIL1DoubleMu0_v1; 
  return;
};

bool isFitAlreadyFound(RooArgSet *newpars, string outputDir, string plotLabel, string TAG, struct KinCuts cut, bool isPbPb, bool doSimulFit) {
  string FileName = "";
  if (doSimulFit) {
   FileName = Form("%sresult/%s/FIT_%s_%s_%s%s_pt%.0f%.0f_rap%.0f%.0f_cent%d%d.root", outputDir.c_str(), TAG.c_str(), TAG.c_str(), "Psi2SJpsi", "COMB", plotLabel.c_str(), (cut.dMuon.Pt.Min*10.0), (cut.dMuon.Pt.Max*10.0), (cut.dMuon.AbsRap.Min*10.0), (cut.dMuon.AbsRap.Max*10.0), cut.Centrality.Start, cut.Centrality.End);
  } else {
    FileName = Form("%sresult/%s/FIT_%s_%s_%s%s_pt%.0f%.0f_rap%.0f%.0f_cent%d%d.root", outputDir.c_str(), TAG.c_str(), TAG.c_str(), "Psi2SJpsi", (isPbPb?"PbPb":"PP"), plotLabel.c_str(), (cut.dMuon.Pt.Min*10.0), (cut.dMuon.Pt.Max*10.0), (cut.dMuon.AbsRap.Min*10.0), (cut.dMuon.AbsRap.Max*10.0), cut.Centrality.Start, cut.Centrality.End);
  }

  if (gSystem->AccessPathName(FileName.c_str())) return false; // File was not found

  TFile *file = new TFile(FileName.c_str());

  if (!file) return false;

  RooWorkspace *ws = (RooWorkspace*) file->Get("workspace");
  if (!ws) {
    file->Close(); delete file;
    return false;
  }

  const RooArgSet *params = ws->getSnapshot(Form("pdfMASS_Tot_%s_parIni", (isPbPb?"PbPb":"PP")));
  if (!params) {
    delete ws;
    file->Close(); delete file;
    return false;
  }

  bool result = compareSnapshots(newpars, params);

  delete ws;
  file->Close(); delete file; 

  return result;
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
