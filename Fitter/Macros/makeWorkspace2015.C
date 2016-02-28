#include "Utilities/initOniaTree.C"
#include "Utilities/EVENTUTILS.h"
#include "Utilities/initClasses.h"

bool useDBFile = true;

bool getTChain(TChain* fChain, vector<TString> FileNames);
void iniBranch(TChain* fChain);
Bool_t passKinematicCuts(Int_t iRecoQQ, struct KinCuts* cut, Bool_t isPbPb);
  

void makeWorkspace2015(RooWorkspace& ws, vector<TString> FileNames, struct InputOpt* opt, struct KinCuts* cut, bool isPbPb, int nbins = 54)
{

  TString FILEDBNAME;
  gSystem->mkdir("./DBFiles/", kTRUE);
  if (isPbPb) {
    FILEDBNAME = Form("./DBFiles/DBFile_%s_%s_Run_%d_%d_%s_pt_%.0f_%.0f_rap_%.0f_%.0f_cent_%d_%d.root", (opt->isData ? "DATA" : "MC"), "PbPb", opt->PbPb.RunNb.Start, opt->PbPb.RunNb.End, (opt->oniaMode==1 ? "Jpsi" : "Upsilon"), (cut->dMuon.Pt.Min*10.0), (cut->dMuon.Pt.Max*10.0), (cut->dMuon.AbsRap.Min*10.0), (cut->dMuon.AbsRap.Max*10.0), cut->Centrality.Start, cut->Centrality.End);      
  } else {
    FILEDBNAME = Form("./DBFiles/DBFile_%s_%s_Run_%d_%d_%s_pt_%.0f_%.0f_rap_%.0f_%.0f.root", (opt->isData ? "DATA" : "MC"), "pp", opt->pp.RunNb.Start, opt->pp.RunNb.End, (opt->oniaMode==1 ? "Jpsi" : "Upsilon"), (cut->dMuon.Pt.Min*10.0), (cut->dMuon.Pt.Max*10.0), (cut->dMuon.AbsRap.Min*10.0), (cut->dMuon.AbsRap.Max*10.0));     
  }
  TFile *DBFile = new TFile(FILEDBNAME,"READ");

  TProfile*   dProf  = NULL; RooDataSet* dataOS = NULL; RooDataSet* dataSS = NULL; RooDataHist* DHProf = NULL;
  if (!(DBFile->IsOpen()&&useDBFile)) {
    if (DBFile->IsOpen()){DBFile->Close();}  
    DBFile = new TFile(FILEDBNAME,"RECREATE");

    TChain* theTree = new TChain(TreeName.Data(),"");
    if(!getTChain(theTree, FileNames)) return;
    initOniaTree(theTree);
    iniBranch(theTree);

    RooRealVar* mass = new RooRealVar("invariantMass","#mu#mu mass", cut->dMuon.M.Min, cut->dMuon.M.Max, "GeV/c^{2}");	  
    RooArgSet cols(*mass);

    dProf  = new TProfile(Form("dProf_%s", (isPbPb?"PbPb":"PP")), "Profile of pt versus mass", nbins, cut->dMuon.M.Min, cut->dMuon.M.Max); dProf->Sumw2();  
    dataOS = new RooDataSet(Form("dataOS_%s", (isPbPb?"PbPb":"PP")), Form("dataOS_%s", (isPbPb?"PbPb":"PP")), cols);
    dataSS = new RooDataSet(Form("dataSS_%s", (isPbPb?"PbPb":"PP")), Form("dataSS_%s", (isPbPb?"PbPb":"PP")), cols);
    
    Long64_t nentries = theTree->GetEntries();
    cout << "[INFO] Starting to process " << nentries << " nentries" << endl;
    for (Long64_t jentry=0; jentry<nentries;jentry++) {

      if (jentry%1000000==0) cout << jentry << "/" << nentries << endl;
       
      if (theTree->LoadTree(jentry)<0) break;
      if (theTree->GetTreeNumber()!=fCurrent) {
	fCurrent = theTree->GetTreeNumber();
	cout << "[INFO] Processing Tree number: " << fCurrent << endl;
      }
      Reco_QQ_4mom->Clear();
      Reco_QQ_mumi_4mom->Clear();
      Reco_QQ_mupl_4mom->Clear();
      theTree->GetEntry(jentry);  
  
      for (int iQQ=0; iQQ<Reco_QQ_size; iQQ++) {
	TLorentzVector *RecoQQ4mom = (TLorentzVector*) Reco_QQ_4mom->At(iQQ);
	mass->setVal(RecoQQ4mom->M());
	
	if ( 
	    ( passKinematicCuts(iQQ, cut, isPbPb)   ) &&
	    ( RecoQQ::areMuonsInAcceptance2011(iQQ) ) &&
	    ( RecoQQ::passQualityCuts2015(iQQ)      ) &&
	    ( RecoQQ::isTriggerMatch(iQQ, (isPbPb ? opt->PbPb.TriggerBit : opt->pp.TriggerBit)) )
	    )
	  {
	    if (Reco_QQ_sign[iQQ]==0) { 
	      dataOS->add(cols);
	      dProf->Fill(RecoQQ4mom->M(), RecoQQ4mom->Pt()); 
	    } else { 
	      dataSS->add(cols);
	    }
	  }
      }
    }
    DHProf = new RooDataHist(Form("DHProf_%s", (isPbPb?"PbPb":"PP")), "DHProf", *mass, Import(*dProf) );
    DBFile->cd();
    dataOS->Write(Form("dataOS_%s", (isPbPb?"PbPb":"PP"))); dataSS->Write(Form("dataSS_%s", (isPbPb?"PbPb":"PP"))); 
    DHProf->Write(Form("DHProf_%s", (isPbPb?"PbPb":"PP"))); dProf->Write(Form("dProf_%s", (isPbPb?"PbPb":"PP")));
    DBFile->Write(); DBFile->Close(); delete DBFile;
  } else {
    dataOS = (RooDataSet*)DBFile->Get(Form("dataOS_%s", (isPbPb?"PbPb":"PP")));  dataSS = (RooDataSet*)DBFile->Get(Form("dataSS_%s", (isPbPb?"PbPb":"PP")));
    DHProf = (RooDataHist*)DBFile->Get(Form("DHProf_%s", (isPbPb?"PbPb":"PP")));
  }

  ws.import(*dataSS);
  ws.import(*dataOS);
  ws.import(*DHProf);

};

  
bool getTChain(TChain* fChain, vector<TString> FileNames) 
{
  cout << "[INFO] Extrating TTree " << TreeName.Data() << endl;
  for (vector<TString>::iterator FileName = FileNames.begin() ; FileName != FileNames.end(); ++FileName){
    cout << "[INFO] Adding TFile " << FileName->Data() << endl;
    fChain->Add(Form("%s/%s", FileName->Data(),  TreeName.Data()));
  } 
  if (!fChain) { cout << "[ERROR] myTree was not found" << endl; return false; } 
  else {return true;}
};

void iniBranch(TChain* fChain)
{
  cout << "[INFO] Initializing Branches of " << TreeName.Data() << endl;
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
};

Bool_t passKinematicCuts(Int_t iRecoQQ, struct KinCuts* cut, Bool_t isPbPb)
{
  TLorentzVector *RecoQQ4mom = (TLorentzVector*) Reco_QQ_4mom->At(iRecoQQ);
  TLorentzVector *RecoQQmupl = (TLorentzVector*) Reco_QQ_mupl_4mom->At(iRecoQQ);
  TLorentzVector *RecoQQmumi = (TLorentzVector*) Reco_QQ_mumi_4mom->At(iRecoQQ);

  Bool_t indMuonMass = (cut->dMuon.M.Min   < RecoQQ4mom->M()   && RecoQQ4mom->M()   < cut->dMuon.M.Max  );
  Bool_t indMuonEta = (cut->sMuon.Eta.Min < RecoQQmupl->Eta() && RecoQQmupl->Eta() < cut->sMuon.Eta.Max);
  Bool_t insMuonEta = (cut->sMuon.Eta.Min < RecoQQmumi->Eta() && RecoQQmumi->Eta() < cut->sMuon.Eta.Max);
  Bool_t indMuonRap = (cut->dMuon.AbsRap.Min < abs(RecoQQ4mom->Rapidity()) && abs(RecoQQ4mom->Rapidity()) < cut->dMuon.AbsRap.Max);
  Bool_t indMuonPt = (cut->dMuon.Pt.Min  < RecoQQ4mom->Pt()  && RecoQQ4mom->Pt()  < cut->dMuon.Pt.Max );
  Bool_t insMuonPt = (RecoQQmupl->Pt()   > cut->sMuon.Pt.Min && RecoQQmumi->Pt()  > cut->sMuon.Pt.Min );
  Bool_t inCentrality = (isPbPb ? (cut->Centrality.Start <= Centrality && Centrality <= cut->Centrality.End) : true);

  return ( indMuonMass && indMuonEta && insMuonEta && indMuonRap && indMuonPt && insMuonPt && inCentrality);
};
