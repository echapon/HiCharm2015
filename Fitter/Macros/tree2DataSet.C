// -*- C++ -*-
//
// Package:    Fitter
// 
/*
 Description: TTree to RooDataSet converter.
 Implementation:
     This program creates two RooDataSets (opposite-sign and same-sign dimuons) from an Onia Tree.
*/
// Original Author:  Andre Stahl,
//         Created:  Feb 26 19:08 CET 2016
//
//

#include "Utilities/initOniaTree.C"
#include "Utilities/EVENTUTILS.h"
#include "Utilities/initClasses.h"


string findMyTree(string FileName);
bool    getTChain(TChain* fChain, vector<string> FileNames);
void    iniBranch(TChain* fChain);  

bool tree2DataSet(RooWorkspace& Workspace, vector<string> InputFileNames, string DSName, string OutputFileName, bool UpdateDS=false)
{
  RooDataSet* dataOS = NULL; RooDataSet* dataSS = NULL;

  if (gSystem->AccessPathName(OutputFileName.c_str()) || UpdateDS) {
    cout << "[INFO] Creating RooDataSet for " << DSName << endl;
    TreeName = findMyTree(InputFileNames[0]); if(TreeName==""){return false;}

    TChain* theTree = new TChain(TreeName.c_str(),"");
    if(!getTChain(theTree, InputFileNames)){ return false; }     // Import files to TChain
    initOniaTree(theTree);                                       // Initialize the Onia Tree
    iniBranch(theTree);                                          // Initialize the Branches

    RooRealVar* mass    = new RooRealVar("invMass","#mu#mu mass", 2.0, 5.0, "GeV/c^{2}");
    RooRealVar* ctau    = new RooRealVar("ctau","c_{#tau}", -10.0, 10.0, "cm");
    RooRealVar* ctauErr = new RooRealVar("ctauErr","#sigma_{c#tau}", -1.0, 1.0, "cm");	
    RooRealVar* ptQQ    = new RooRealVar("pt","#mu#mu p_{T}", 0.0, 50.0, "GeV/c");
    RooRealVar* rapQQ   = new RooRealVar("rap","#mu#mu y", -2.4, 2.4, "");
    RooRealVar* cent    = new RooRealVar("cent","centrality", 0.0, 200.0, "");  
    RooArgSet cols(*mass, *ctau, *ctauErr, *ptQQ, *rapQQ, *cent);

    dataOS = new RooDataSet(Form("dOS_%s", DSName.c_str()), "dOS", cols);
    dataSS = new RooDataSet(Form("dSS_%s", DSName.c_str()), "dSS", cols);
    
    Long64_t nentries = theTree->GetEntries();
    cout << "[INFO] Starting to process " << nentries << " nentries" << endl;
    for (Long64_t jentry=0; jentry<nentries;jentry++) {

      if (jentry%1000000==0) cout << jentry << "/" << nentries << endl;
       
      if (theTree->LoadTree(jentry)<0) break;
      if (theTree->GetTreeNumber()!=fCurrent) {
	fCurrent = theTree->GetTreeNumber();
	cout << "[INFO] Processing Root File: " << InputFileNames[fCurrent] << endl;
      }
      Reco_QQ_4mom->Clear();
      Reco_QQ_mumi_4mom->Clear();
      Reco_QQ_mupl_4mom->Clear();
      theTree->GetEntry(jentry);  
  
      for (int iQQ=0; iQQ<Reco_QQ_size; iQQ++) {
	TLorentzVector *RecoQQ4mom = (TLorentzVector*) Reco_QQ_4mom->At(iQQ);
	mass->setVal(RecoQQ4mom->M());
        ctau->setVal(Reco_QQ_ctau3D[iQQ]);
        ctauErr->setVal(Reco_QQ_ctauErr3D[iQQ]);
	ptQQ->setVal(RecoQQ4mom->Pt());
	rapQQ->setVal(RecoQQ4mom->Rapidity());
	cent->setVal(Centrality);
	
	if ( 
	    ( RecoQQ::areMuonsInAcceptance2015(iQQ) ) &&  // 2015 Global Muon Acceptance Cuts
	    ( RecoQQ::passQualityCuts2015(iQQ)      ) &&  // 2015 Soft Global Muon Quality Cuts
	    ( RecoQQ::isTriggerMatch(iQQ, 0)        )     // HLT_HIL1DoubleMu0_v1
	    )
	  {
	    if (Reco_QQ_sign[iQQ]==0) { // Opposite-Sign dimuons
	      dataOS->add(cols);
	    } else {                    // Like-Sign dimuons 
	      dataSS->add(cols);
	    }
	  }
      }
    }
    TFile *DBFile = TFile::Open(OutputFileName.c_str(),"RECREATE");
    DBFile->cd();
    dataOS->Write(Form("dOS_%s", DSName.c_str())); 
    dataSS->Write(Form("dSS_%s", DSName.c_str())); 
    DBFile->Write(); DBFile->Close(); delete DBFile;

  } else {

    cout << "[INFO] Loading RooDataSet from " << OutputFileName << endl;
    
    TFile *DBFile = TFile::Open(OutputFileName.c_str(),"READ");
    dataOS = (RooDataSet*)DBFile->Get(Form("dOS_%s", DSName.c_str()));  
    dataSS = (RooDataSet*)DBFile->Get(Form("dSS_%s", DSName.c_str()));
  }

  if(!dataOS || !dataSS){ cout << "[ERROR] " << DSName << " was not found" << endl; return false; } 
  Workspace.import(*dataOS); Workspace.import(*dataSS);
						   
  return true;
};

string findMyTree(string FileName)
{
  TFile *f = TFile::Open(FileName.c_str(), "READ");
  if(f->GetListOfKeys()->Contains("hionia")){ return "hionia/myTree"; }
  else if(f->GetListOfKeys()->Contains("myTree")){ return "myTree"; }
  else { cout << "[ERROR] myTree was not found in: " << FileName << endl; }
  return ""; 
};
  
bool getTChain(TChain *fChain, vector<string> FileNames) 
{
  cout << "[INFO] Extrating TTree " << TreeName.c_str() << endl;
  for (vector<string>::iterator FileName = FileNames.begin() ; FileName != FileNames.end(); ++FileName){
    cout << "[INFO] Adding TFile " << FileName->c_str() << endl;
    fChain->Add(Form("%s/%s", FileName->c_str(),  TreeName.c_str()));
  } 
  if (!fChain) { cout << "[ERROR] fChain was not created, some input files are missing" << endl; return false; } 
  return true;
};

void iniBranch(TChain* fChain)
{
  cout << "[INFO] Initializing Branches of " << TreeName.c_str() << endl;
  fChain->GetBranch("Reco_QQ_4mom")->SetAutoDelete(kFALSE);   
  fChain->GetBranch("Reco_QQ_mupl_4mom")->SetAutoDelete(kFALSE); 
  fChain->GetBranch("Reco_QQ_mumi_4mom")->SetAutoDelete(kFALSE); 
  fChain->SetBranchStatus("*",0);
  RecoQQ::iniBranches(fChain); 
  fChain->SetBranchStatus("Centrality",1); 
  fChain->SetBranchStatus("Reco_QQ_size",1); 
  fChain->SetBranchStatus("Reco_QQ_sign",1); 
  fChain->SetBranchStatus("Reco_QQ_4mom",1); 
  fChain->SetBranchStatus("Reco_QQ_mupl_4mom",1); 
  fChain->SetBranchStatus("Reco_QQ_mumi_4mom",1);
  fChain->SetBranchStatus("Reco_QQ_ctau",1); 
  fChain->SetBranchStatus("Reco_QQ_ctau3D",1); 
  fChain->SetBranchStatus("Reco_QQ_ctauErr",1); 
  fChain->SetBranchStatus("Reco_QQ_ctauErr3D",1);
};
