#include "Utilities/initOniaTree.C"
#include "Utilities/EVENTUTILS.h"
#include "Utilities/initClasses.h"

bool useDBFile = true;

TString findMyTree(TString FileName);
bool    getTChain(TChain* fChain, vector<TString> FileNames);
void    iniBranch(TChain* fChain);
Bool_t  passKinematicCuts(Int_t iRecoQQ, struct KinCuts* cut, Bool_t isPbPb);
  

void makeWorkspace2015(RooWorkspace& ws, vector<TString> FileNames, struct InputOpt* opt, struct KinCuts* cut, int nbins = 54, bool isPbPb=true, bool isData=true, string MCTYPE="")
{

  TString FILEDBNAME;
  if (isData) {
    gSystem->mkdir("./DBFiles/DATA/", kTRUE);
    if (isPbPb) {
      FILEDBNAME = Form("./DBFiles/DATA/DBFile_DATA_PbPb_Run_%d_%d_%s_pt_%.0f_%.0f_rap_%.0f_%.0f_cent_%d_%d.root", opt->PbPb.RunNb.Start, opt->PbPb.RunNb.End, (opt->oniaMode==1 ? "Charm" : "Bottom"), (cut->dMuon.Pt.Min*10.0), (cut->dMuon.Pt.Max*10.0), (cut->dMuon.AbsRap.Min*10.0), (cut->dMuon.AbsRap.Max*10.0), cut->Centrality.Start, cut->Centrality.End);      
    } else {
      FILEDBNAME = Form("./DBFiles/DATA/DBFile_DATA_PP_Run_%d_%d_%s_pt_%.0f_%.0f_rap_%.0f_%.0f.root", opt->pp.RunNb.Start, opt->pp.RunNb.End, (opt->oniaMode==1 ? "Charm" : "Bottom"), (cut->dMuon.Pt.Min*10.0), (cut->dMuon.Pt.Max*10.0), (cut->dMuon.AbsRap.Min*10.0), (cut->dMuon.AbsRap.Max*10.0));     
    }
  } else {
    gSystem->mkdir("./DBFiles/MC/", kTRUE);
    if (isPbPb) {
      FILEDBNAME = Form("./DBFiles/MC/DBFile_MC_PbPb %s_Run_%d_%d_pt_%.0f_%.0f_rap_%.0f_%.0f_cent_%d_%d.root", MCTYPE.c_str(), opt->PbPb.RunNb.Start, opt->PbPb.RunNb.End, (cut->dMuon.Pt.Min*10.0), (cut->dMuon.Pt.Max*10.0), (cut->dMuon.AbsRap.Min*10.0), (cut->dMuon.AbsRap.Max*10.0), cut->Centrality.Start, cut->Centrality.End);      
    } else {
      FILEDBNAME = Form("./DBFiles/MC/DBFile_MC_PP_%s_Run_%d_%d_pt_%.0f_%.0f_rap_%.0f_%.0f.root", MCTYPE.c_str(), opt->pp.RunNb.Start, opt->pp.RunNb.End, (cut->dMuon.Pt.Min*10.0), (cut->dMuon.Pt.Max*10.0), (cut->dMuon.AbsRap.Min*10.0), (cut->dMuon.AbsRap.Max*10.0));     
    }
  }
  TFile *DBFile = new TFile(FILEDBNAME,"READ");

  TProfile*   dProf  = NULL; RooDataSet* dataOS = NULL; RooDataSet* dataSS = NULL; RooDataHist* DHProf = NULL;
  if (!(DBFile->IsOpen()&&useDBFile)) {
    if (DBFile->IsOpen()){DBFile->Close();}  
    DBFile = new TFile(FILEDBNAME,"RECREATE");
    
    TreeName = findMyTree(FileNames[0]);
    TChain* theTree = new TChain(TreeName.c_str(),"");
    if(!getTChain(theTree, FileNames)) return;
    initOniaTree(theTree);
    iniBranch(theTree);

    RooRealVar* mass = new RooRealVar("invMass","#mu#mu mass", cut->dMuon.M.Min, cut->dMuon.M.Max, "GeV/c^{2}");
    RooRealVar* ctau = new RooRealVar("ctau","c_{#tau}", cut->dMuon.ctau.Min, cut->dMuon.ctau.Max, "cm");
    RooRealVar* ctauErr = new RooRealVar("ctauErr","#sigma_{c#tau}", cut->dMuon.ctauErr.Min, cut->dMuon.ctauErr.Max, "cm");	  
    RooArgSet cols(*mass, *ctau, *ctauErr);

    dProf  = new TProfile(Form("dProf_%s_%s", (isData?"DATA":Form("MC%s", MCTYPE.c_str())), (isPbPb?"PbPb":"PP")), "Mean PT", nbins, cut->dMuon.M.Min, cut->dMuon.M.Max); dProf->Sumw2();  
    dataOS = new RooDataSet(Form("dOS_%s_%s", (isData?"DATA":Form("MC%s", MCTYPE.c_str())), (isPbPb?"PbPb":"PP")), "dOS", cols);
    dataSS = new RooDataSet(Form("dSS_%s_%s", (isData?"DATA":Form("MC%s", MCTYPE.c_str())), (isPbPb?"PbPb":"PP")), "dSS", cols);
    
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
        ctau->setVal(Reco_QQ_ctau3D[iQQ]);
        ctauErr->setVal(Reco_QQ_ctauErr3D[iQQ]);
	
	if ( 
	    ( passKinematicCuts(iQQ, cut, isPbPb)   ) &&
	    ( RecoQQ::areMuonsInAcceptance2015(iQQ) ) &&
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
    DHProf = new RooDataHist(Form("dHProf_%s_%s", (isData?"DATA":Form("MC%s", MCTYPE.c_str())), (isPbPb?"PbPb":"PP")), "DHProf", *mass, Import(*dProf) );
    DBFile->cd();
    dataOS->Write(Form("dOS_%s_%s", (isData?"DATA":Form("MC%s", MCTYPE.c_str())), (isPbPb?"PbPb":"PP"))); 
    dataSS->Write(Form("dSS_%s_%s", (isData?"DATA":Form("MC%s", MCTYPE.c_str())), (isPbPb?"PbPb":"PP"))); 
    DHProf->Write(Form("dHProf_%s_%s", (isData?"DATA":Form("MC%s", MCTYPE.c_str())), (isPbPb?"PbPb":"PP"))); 
    dProf->Write(Form("dProf_%s_%s", (isData?"DATA":Form("MC%s", MCTYPE.c_str())), (isPbPb?"PbPb":"PP")));
    DBFile->Write(); DBFile->Close(); delete DBFile;
  } else {
    dataOS = (RooDataSet*)DBFile->Get(Form("dOS_%s_%s", (isData?"DATA":Form("MC%s", MCTYPE.c_str())), (isPbPb?"PbPb":"PP")));  
    dataOS = (RooDataSet*)dataOS->reduce(Form("(invMass>%.2f) && (invMass<%.2f)",cut->dMuon.M.Min, cut->dMuon.M.Max));
    dataSS = (RooDataSet*)DBFile->Get(Form("dSS_%s_%s", (isData?"DATA":Form("MC%s", MCTYPE.c_str())), (isPbPb?"PbPb":"PP")));
    dataSS = (RooDataSet*)dataSS->reduce(Form("(invMass>%.2f) && (invMass<%.2f)",cut->dMuon.M.Min, cut->dMuon.M.Max));
    DHProf = (RooDataHist*)DBFile->Get(Form("dHProf_%s_%s", (isData?"DATA":Form("MC%s", MCTYPE.c_str())), (isPbPb?"PbPb":"PP")));
  }

  ws.import(*dataSS);
  ws.import(*dataOS);
  ws.import(*DHProf);

  ws.var("invMass")->setMin(cut->dMuon.M.Min);  ws.var("invMass")->setMax(cut->dMuon.M.Max);

};

TString findMyTree(TString FileName)
{
  TFile *f = TFile::Open(FileName.Data(), "READ");
  if(f->GetListOfKeys()->Contains("hionia")){ return TString("hionia/myTree"); }
  else if(f->GetListOfKeys()->Contains("myTree")){ return TString("myTree"); }
  return TString(); 
}
  
bool getTChain(TChain* fChain, vector<TString> FileNames) 
{
  ///// ADD OPTION TO SEARCH ON EITHER MYTREE OR HIONIA/MYTREE
  cout << "[INFO] Extrating TTree " << TreeName.c_str() << endl;
  for (vector<TString>::iterator FileName = FileNames.begin() ; FileName != FileNames.end(); ++FileName){
    cout << "[INFO] Adding TFile " << FileName->Data() << endl;
    fChain->Add(Form("%s/%s", FileName->Data(),  TreeName.c_str()));
  } 
  if (!fChain) { cout << "[ERROR] myTree was not found" << endl; return false; } 
  else {return true;}
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
