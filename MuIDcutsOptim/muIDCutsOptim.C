#define muIDCutsOptim_cxx
#include "muIDCutsOptim.h"
#include "Riostream.h"
#include "TH2.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TLorentzVector.h"
#include "TString.h"
#include "TLatex.h"
#include "TObjArray.h"
#include "TObjString.h"

#include <iostream>
#include <map>

using namespace std;


//_____________________________
void muIDCutsOptim::Loop()
{
  //  In a ROOT session, you can do:
  //      root> .L muIDCutsOptim.C
  //      root> muIDCutsOptim t
  //      root> t.GetEntry(12); // Fill t data members with entry number 12
  //      root> t.Show();       // Show values of entry 12
  //      root> t.Show(16);     // Read and show values of entry 16
  //      root> t.Loop();       // Loop on all entries
  //
  
  //  This is the loop skeleton where:
  //    jentry is the global entry number in the chain
  //    ientry is the entry number in the current Tree
  //  Note that the argument to GetEntry must be:
  //    jentry for TChain::GetEntry
  //    ientry for TTree::GetEntry and TBranch::GetEntry
  //

  if (fChain == 0) return;
  
  
  //===== Definition of file type, particle and if we want the high purity cut included in the soft muon ID cuts
  const char* fileType = "MC";
  const char* particle = "JPsi";
  fIspp = kFALSE;
  fIncludeHighPurity = kFALSE;
  //=====
  
  
  Long64_t nentries(0);
  if ( !strcmp(fileType,"MC") ) nentries = fChain->GetEntries();
  else nentries = 1500000;
  
  std::cout << "# Events = " << nentries << std::endl;
  
  
  //===== Definition of some cuts
  Double_t leMinvSig = 3.0; // Signal range
  Double_t ueMinvSig = 3.2;
  
  Double_t leMinvBkg1 = 2.7; // Bkg 1 range (sideband)
  Double_t ueMinvBkg1 = 2.9;
  
  Double_t leMinvBkg2 = 3.3; // Bkg 1 range (sideband)
  Double_t ueMinvBkg2 = 3.5;

  Double_t Ptmin = 3.0; // Pt cut (6.5) (3.0)
  Double_t Ptmax = 6.5; // (12.0) (6.5)
  
  Double_t Ymin = 2.0; // Y cut (0.0) (2.0)
  Double_t Ymax = 2.4; // (2.4) (2.4)
  //=====
  
  
  //===== Definition of MuID variables and binning
  const int nvar = 19;
  const char* varname[19] = {"isGoodMuon", "highPurity", "TrkMuArb", "TMOneStaTight", "nPixValHits",
    "nMuValHits", "nTrkHits", "normChi2_inner", "normChi2_global", "nPixWMea",
    "nTrkWMea", "StationsMatched", "dxy", "dxyErr", "dz",
    "dzErr", "ptErr_inner", "ptErr_global","VtxProb"};
  const int nbins[19] = {2, 2, 2, 2, 11,
    56, 36, 100, 100, 6,
    19, 7, 100, 100, 100,
    100, 100, 100, 200};
  const double lowedge[19] = {0, 0, 0, 0, 0,
    0, 0, 0, 0, 0,
    0, 0, 0, 0, 0,
    0, 0, 0, 0};
  const double hiedge[19] = {2, 2, 2, 2, 11,
    56, 36, 10, 50, 6,
    19, 7, 2, 0.25, 20,
    5, 0.15, 0.15, 1.};
  //=====
  
  
  SelectBranches(nvar,varname); // Activates only desired branches (General and MuID ones in this case)
  
  
  //====== Create histograms and arrays to store them ========
  TObjArray* aDistr = new TObjArray();
  TObjArray* aSig = new TObjArray();
  TObjArray* aBkg = new TObjArray();
  TObjArray* aSigC = new TObjArray();
  TObjArray* aBkgC = new TObjArray();
  TObjArray* aGlob = new TObjArray();
  TObjArray* aMinvCuts = new TObjArray();
  TObjArray* aSigSoftMuCutsButOne = new TObjArray();
  TObjArray* aBkgSoftMuCutsButOne = new TObjArray();
  TObjArray* aSigSoftMuCuts = new TObjArray();
  TObjArray* aBkgSoftMuCuts = new TObjArray();
  
  map<TString, TH1F*> hists_sig, hists_bkg, hists_sig_distr, hists_bkg_distr, hists_distr, hists_sig_distr_CutTest, hists_bkg_distr_CutTest, hists_sig_SofMu, hists_bkg_SofMu;
  
  TString sHPName("");
  if ( fIncludeHighPurity ) sHPName += "highPurity"; // This is for the histograms creation
  
  TObjArray* dummyA(0x0);
  TString partType("");
  for (int i=0; i<nvar; i++)
  {
    // Single cuts histos
    hists_sig[Form("%s",varname[i])] = new TH1F(Form("h%s_sig",varname[i]),Form("Signal di-#mu %s cut distribution",varname[i]),nbins[i],lowedge[i],hiedge[i]);
    hists_bkg[Form("%s",varname[i])] = new TH1F(Form("h%s_bkg",varname[i]),Form("Background di-#mu %s cut distribution",varname[i]),nbins[i],lowedge[i],hiedge[i]);
    
    // Variable distributions histos
    hists_sig_distr[Form("%s",varname[i])] = new TH1F(Form("h%s_sig_distr",varname[i]),Form("Signal single-#mu %s distribution",varname[i]),nbins[i],lowedge[i],hiedge[i]);
    hists_bkg_distr[Form("%s",varname[i])] = new TH1F(Form("h%s_bkg_distr",varname[i]),Form("Background single-#mu %s distribution",varname[i]),nbins[i],lowedge[i],hiedge[i]);
    
    partType.Clear();
    if ( !strcmp(varname[i],"VtxProb") ) partType += "di-#mu";
    else partType += "di-#mu";
    hists_distr[Form("%s",varname[i])] = new TH1F(Form("h%s_distr",varname[i]),Form("Sig+Bkg %s %s distribution",partType.Data(),varname[i]),nbins[i],lowedge[i],hiedge[i]);
    
    dummyA = new TObjArray();
    dummyA->SetName(Form("mInv_%s",varname[i]));
    aMinvCuts->Add(dummyA);
    
    if ( !strcmp(varname[i],"isGoodMuon") || !strcmp(varname[i],sHPName.Data()) || !strcmp(varname[i],"nPixWMea") || !strcmp(varname[i],"nTrkWMea") || !strcmp(varname[i],"dxy") || !strcmp(varname[i],"dz") )
    {
      // Histos for all the soft muon cuts applied but the one on the histo
      hists_sig_distr_CutTest[Form("%s",varname[i])] = new TH1F(Form("hSoftMuNoCut_%s_sig",varname[i]),Form("Signal single-#mu %s distribution (w/ rest of soft MuID cuts)",varname[i]),nbins[i],lowedge[i],hiedge[i]);
      hists_bkg_distr_CutTest[Form("%s",varname[i])] = new TH1F(Form("hSoftMuNoCut_%s_bkg",varname[i]),Form("Background single-#mu %s distribution (w/ rest of soft MuID cuts)",varname[i]),nbins[i],lowedge[i],hiedge[i]);
      
      hists_sig_SofMu[Form("%s",varname[i])] = new TH1F(Form("hSoftMuCut_%s_sig",varname[i]),Form("Signal di-#mu %s cut distribution (w/ all Soft MuID cuts except %s)",varname[i],varname[i]),nbins[i],lowedge[i],hiedge[i]);
      hists_bkg_SofMu[Form("%s",varname[i])] = new TH1F(Form("hSoftMuCut_%s_bkg",varname[i]),Form("Background di-#mu %s cut distribution (w/ all Soft MuID cuts except %s)",varname[i],varname[i]),nbins[i],lowedge[i],hiedge[i]);
    }
    else
    {
      // Histos for all the soft muon cuts applied
      hists_sig_distr_CutTest[Form("%s",varname[i])] = new TH1F(Form("hSoftMuFull_%s_sig",varname[i]),Form("Signal %s %s distribution (w/ all Soft MuID cuts)",partType.Data(),varname[i]),nbins[i],lowedge[i],hiedge[i]);
      hists_bkg_distr_CutTest[Form("%s",varname[i])] = new TH1F(Form("hSoftMuFull_%s_bkg",varname[i]),Form("Signal %s %s distribution (w/ all Soft MuID cuts)",partType.Data(),varname[i]),nbins[i],lowedge[i],hiedge[i]);
      
      hists_sig_SofMu[Form("%s",varname[i])] = new TH1F(Form("hSoftMuCut_%s_sig",varname[i]),Form("Signal di-#mu %s cut distribution (w/ all Soft MuID cuts)",varname[i]),nbins[i],lowedge[i],hiedge[i]);
      hists_bkg_SofMu[Form("%s",varname[i])] = new TH1F(Form("hSoftMuCut_%s_bkg",varname[i]),Form("Background di-#mu %s cut distribution (w/ all Soft MuID cuts)",varname[i]),nbins[i],lowedge[i],hiedge[i]);
    }
  }


  // Histos for global distributions
  TH1* hpt = new TH1D("hPtSig",Form("Dimuon p_{T} distribution (%2.1f < m_{#mu^{+}#mu^{-}} < %2.1f)",leMinvSig,ueMinvSig),300,0.0,15.0);
  TH1* hptw = new TH1D("hPtSigW",Form("Dimuon p_{T} distribution weighted (%2.1f < m_{#mu^{+}#mu^{-}} < %2.1f)",leMinvSig,ueMinvSig),300,0.0,15.0);
  
  TH1* hMinvw = new TH1D("hmInv","Dimuon m_{inv} distribution",600,0.,15.0);
  TH1* hRapw = new TH1D("hRap","Dimuon rapidity distribution",180,-4.5,4.5);
  TH1* hPtw = new TH1D("hPt","Dimuon p_{T} distribution",300,0.0,15.0);
  TH1* hCent = new TH1D("hCent","Centrality distribution",200,0.,100);
  TH1* hCentw = new TH1D("hCentw","Weighted centrality distribution",200,0.,100);
  //======================
  
  
  TLorentzVector *tlvmupl(0x0);
  TLorentzVector *tlvmumi(0x0);
  TLorentzVector *tlvqq(0x0);

  Long64_t nbytes = 0, nb = 0;
  //*************** Event loop
  for (Long64_t jentry=0; jentry<nentries;jentry++)
  {
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);   nbytes += nb;
    
    if ( jentry % 1000 == 0 ) std::cout << jentry << "/" << nentries << std::endl;
    
    Double_t weight = fChain->GetWeight()*GetNColl(Centrality); // This is the weight for the Pt and centrality. The weight stored in the chain is filterEficiency/nentries

    //Fill centrality histo
    hCentw->Fill(Centrality/2.,weight);
    hCent->Fill(Centrality/2.);
    
    
    //*************** Loop on reco dimuons
//    for (int i=0; i<Reco_QQ_size; i++)
//    {
    Int_t i(0);
    Bool_t kQQfound(kFALSE);
    while ( i < Reco_QQ_size && !kQQfound ) // This is just to take only the first dimuon passing the basic conditions, to not double count single muons
    {
      if ( !isTriggerSelected(i) || Reco_QQ_sign[i]!=0 )
      {
        i++;
        continue; // Checks if the event is selected by the trigger
      }
      
//      if (Reco_QQ_sign[i]!=0) continue;
      tlvmupl = (TLorentzVector*) Reco_QQ_mupl_4mom->At(i);
      tlvmumi = (TLorentzVector*) Reco_QQ_mumi_4mom->At(i);
      
      if (!IsAccept(tlvmupl->Pt(), tlvmupl->Eta()) || !IsAccept(tlvmumi->Pt(), tlvmumi->Eta()))
      {
        i++;
        continue; // Basic single muon kinematic cuts
      }
      
      tlvqq = (TLorentzVector*) Reco_QQ_4mom->At(i);
      Double_t mass = tlvqq->M();
      Double_t QQPt = tlvqq->Pt();
      Double_t QQY = TMath::Abs(tlvqq->Rapidity());
      
      if ( QQPt < Ptmin || QQPt > Ptmax || QQY < Ymin || QQY > Ymax )
      {
        i++;
        continue; // The simulation is restricted to pt in [Ptmin,Ptmax], so this rejects background dmuons outside the range
      }
      
      bool issig = (mass>leMinvSig && mass<ueMinvSig);
      bool isbkg = ((mass>leMinvBkg1 && mass<ueMinvBkg1) || (mass>leMinvBkg2 && mass<ueMinvBkg2));
      
      if ( !issig && ! isbkg )
      {
        i++;
        continue; // We keep only dimuons in the sidebands and signal ranges
      }
      
      kQQfound = kTRUE; // If at this point the dimuon has passed the basic conditions we keep only this dimuon from the event
      
      const double varValP[19] = {Reco_QQ_mupl_isGoodMuon[i],Reco_QQ_mupl_highPurity[i],Reco_QQ_mupl_TrkMuArb[i],Reco_QQ_mupl_TMOneStaTight[i],Reco_QQ_mupl_nPixValHits[i],Reco_QQ_mupl_nMuValHits[i],Reco_QQ_mupl_nTrkHits[i],Reco_QQ_mupl_normChi2_inner[i],Reco_QQ_mupl_normChi2_global[i],Reco_QQ_mupl_nPixWMea[i],Reco_QQ_mupl_nTrkWMea[i],Reco_QQ_mupl_StationsMatched[i],Reco_QQ_mupl_dxy[i],Reco_QQ_mupl_dxyErr[i],Reco_QQ_mupl_dz[i],Reco_QQ_mupl_dzErr[i],Reco_QQ_mupl_ptErr_inner[i],Reco_QQ_mupl_ptErr_global[i],Reco_QQ_VtxProb[i]};
      const double varValM[19] = {Reco_QQ_mumi_isGoodMuon[i],Reco_QQ_mumi_highPurity[i],Reco_QQ_mumi_TrkMuArb[i],Reco_QQ_mumi_TMOneStaTight[i],Reco_QQ_mumi_nPixValHits[i],Reco_QQ_mumi_nMuValHits[i],Reco_QQ_mumi_nTrkHits[i],Reco_QQ_mumi_normChi2_inner[i],Reco_QQ_mumi_normChi2_global[i],Reco_QQ_mumi_nPixWMea[i],Reco_QQ_mumi_nTrkWMea[i],Reco_QQ_mumi_StationsMatched[i],Reco_QQ_mumi_dxy[i],Reco_QQ_mumi_dxyErr[i],Reco_QQ_mumi_dz[i],Reco_QQ_mumi_dzErr[i],Reco_QQ_mumi_ptErr_inner[i],Reco_QQ_mumi_ptErr_global[i],0.}; //Reco_QQ_VtxProb set to 0. since it is a dimuon variable
      

      //**************************************************
      //*************** Fill histos **********************
      //**************************************************
      
      //===== General distributions
      hMinvw->Fill(mass,weight);
      hRapw->Fill(tlvqq->Rapidity(),weight);
      hPtw->Fill(QQPt,weight);
      
      if (issig)
      {
        hpt->Fill(QQPt);
        hptw->Fill(QQPt,weight);
      }
      //=====
      
      
      //===== MuID distributions
      for (int ihist=0; ihist<nvar; ihist++)
      {
        Bool_t isQQVar = kFALSE;
        if ( !strcmp(varname[ihist],"VtxProb") ) isQQVar = kTRUE; // Just to avoid filling the histos twice if testing a variable corresponding to a dimuon
          
        hists_distr[varname[ihist]]->Fill(varValP[ihist],weight); // Signal+ Background variable distribution
        if ( !isQQVar ) hists_distr[varname[ihist]]->Fill(varValM[ihist],weight);
        
        if (issig)
        {
          hists_sig_distr[varname[ihist]]->Fill(varValP[ihist],weight);  // Signal variable distribution
          if ( !isQQVar ) hists_sig_distr[varname[ihist]]->Fill(varValM[ihist],weight);
        }
        if (isbkg)
        {
          hists_bkg_distr[varname[ihist]]->Fill(varValP[ihist],weight);  // Background variable distribution
          if ( !isQQVar ) hists_bkg_distr[varname[ihist]]->Fill(varValM[ihist],weight);
        }

        
        ////_______ Soft muons cut testing MuID distributions
        Int_t nTest = whichMuIDVar(varname[ihist]); // If the variable is not on the soft muon cuts the full cut is applied
        
        if (isSoftMuon(varValP,nTest))
        {
          if (issig) hists_sig_distr_CutTest[varname[ihist]]->Fill(varValP[ihist],weight); // Variable distribution with the rest of soft MuID cuts applied (all the soft muon cuts if the variable does't belong to the soft cut)
          if (isbkg) hists_bkg_distr_CutTest[varname[ihist]]->Fill(varValP[ihist],weight);
        }
        if ( !isQQVar && isSoftMuon(varValM,nTest) )
        {
          if (issig) hists_sig_distr_CutTest[varname[ihist]]->Fill(varValM[ihist],weight);
          if (isbkg) hists_bkg_distr_CutTest[varname[ihist]]->Fill(varValM[ihist],weight);
        }
        ////_______
        
        
        dummyA = static_cast<TObjArray*>(aMinvCuts->FindObject(Form("mInv_%s",varname[ihist])));
        
        TH1* dummyH(0x0);
        Double_t nBins = hists_sig[varname[ihist]]->GetNbinsX();
        for ( Int_t j=1 ; j<=nBins ; j++)
        {
          //// Cuts on MuID variables
          Double_t cutVal = hists_sig[varname[ihist]]->GetBinLowEdge(j);
          if (Cut(varname[ihist],i,cutVal))
          {
            if (issig) hists_sig[varname[ihist]]->Fill(hists_sig[varname[ihist]]->GetBinCenter(j),weight); // Variable cut distribution
            if (isbkg) hists_bkg[varname[ihist]]->Fill(hists_bkg[varname[ihist]]->GetBinCenter(j),weight);
            
            if ( !(dummyH = static_cast<TH1*>(dummyA->FindObject(Form("mInv_%s_Cut%2.4f",varname[ihist],cutVal)))) ) // Creates the minv histogram and store in array (only once)
            {
              dummyH = new TH1D(Form("mInv_%s_Cut%2.4f",varname[ihist],cutVal),Form("Dimuon m_{inv} distribution (cut on %s: %2.4f)",varname[ihist],cutVal),600,0.,15.0);
              dummyA->Add(dummyH);
            }
            
            dummyH->Fill(mass,weight); // Invariant mass histo for each variable cut
            
            if ( isSoftMuon(varValP,nTest) && isSoftMuon(varValM,nTest) ) // If the variable is not in the soft MuID cut, the full soft cut is applied and the variable cut distribution is filled
            {
              if (issig) hists_sig_SofMu[varname[ihist]]->Fill(hists_sig[varname[ihist]]->GetBinCenter(j),weight);
              if (isbkg) hists_bkg_SofMu[varname[ihist]]->Fill(hists_bkg[varname[ihist]]->GetBinCenter(j),weight);
            }
          }
        }
      }
      //=======
      i++;
    }
  }
  
  //====== Save histos ======
  TFile *f = new TFile(Form("histos_%s_%s_%s_Pt%2.1f_%2.1f_Y%2.1f_%2.1f.root",fileType,particle,fIncludeHighPurity ? "HPincl" : "NoHPincl",Ptmin,Ptmax,Ymin,Ymax),"RECREATE");
  
  for (int i=0; i<nvar; i++)
  {
    aSigC->Add(hists_sig[varname[i]]);
    aBkgC->Add(hists_bkg[varname[i]]);
    
    aSig->Add(hists_sig_distr[varname[i]]);
    aBkg->Add(hists_bkg_distr[varname[i]]);
  
    aDistr->Add(hists_distr[varname[i]]);
    
    aSigSoftMuCutsButOne->Add(hists_sig_distr_CutTest[varname[i]]);
    aBkgSoftMuCutsButOne->Add(hists_bkg_distr_CutTest[varname[i]]);
    
    aSigSoftMuCuts->Add(hists_sig_SofMu[varname[i]]);
    aBkgSoftMuCuts->Add(hists_bkg_SofMu[varname[i]]);
  }
  
  aSigC->Write("SigSingleCut", TObject::kOverwrite | TObject::kSingleKey);
  aBkgC->Write("BkgSingleCut", TObject::kOverwrite | TObject::kSingleKey);
  
  aSig->Write("SigDistr", TObject::kOverwrite | TObject::kSingleKey);
  aBkg->Write("BkgDistr", TObject::kOverwrite | TObject::kSingleKey);
  
  aDistr->Write("Distr", TObject::kOverwrite | TObject::kSingleKey);
  
  aSigSoftMuCutsButOne->Write("SigDistrSoftMuCutButOne", TObject::kOverwrite | TObject::kSingleKey);
  aBkgSoftMuCutsButOne->Write("BkgDistrSoftMuCutButOne", TObject::kOverwrite | TObject::kSingleKey);
  
  aSigSoftMuCuts->Write("SigSoftMuCut", TObject::kOverwrite | TObject::kSingleKey);
  aBkgSoftMuCuts->Write("BkgSoftMuCut", TObject::kOverwrite | TObject::kSingleKey);
  
  
  aGlob->Add(hMinvw);
  aGlob->Add(hRapw);
  aGlob->Add(hPtw);
  aGlob->Add(hpt);
  aGlob->Add(hptw);
  aGlob->Add(hCent);
  aGlob->Add(hCentw);

  aGlob->Write("Global", TObject::kOverwrite | TObject::kSingleKey);
  
  
  aMinvCuts->Write("MinvWCuts", TObject::kOverwrite | TObject::kSingleKey);
  
  f->Close();
  //=========
  
}


//_____________________________
void muIDCutsOptim::SelectBranches(Int_t nvar,const char* varList [ ])
{
  //       To read only selected branches, Insert statements like:
  // METHOD1:
  //    fChain->SetBranchStatus("*",0);  // disable all branches
  //    fChain->SetBranchStatus("branchname",1);  // activate branchname
  // METHOD2: replace line
  //    fChain->GetEntry(jentry);       //read all branches
  //by  b_branchname->GetEntry(ientry); //read only this branch
  
  fChain->SetBranchStatus("*",0);  // disable all branches
  
  //====== Enable basic branches ======
  fChain->SetBranchStatus("Centrality",1);
  fChain->SetBranchStatus("Reco_QQ_size",1);
  fChain->SetBranchStatus("Reco_QQ_sign",1);
  fChain->SetBranchStatus("Reco_QQ_mupl_4mom",1);
  fChain->SetBranchStatus("Reco_QQ_mumi_4mom",1);
  fChain->SetBranchStatus("Reco_QQ_4mom",1);
  
  fChain->SetBranchStatus("Reco_QQ_mupl_pt_inner",1);
  fChain->SetBranchStatus("Reco_QQ_mumi_pt_inner",1);
  
  fChain->SetBranchStatus("Reco_QQ_mupl_pt_global",1);
  fChain->SetBranchStatus("Reco_QQ_mumi_pt_global",1);
  
  fChain->SetBranchStatus("HLTriggers",1);
  fChain->SetBranchStatus("Reco_QQ_trig",1);
  
  fChain->SetBranchStatus("Reco_QQ_VtxProb",1);
  //=========
  
  //====== Enable MuID related branches ======
  for (int i=0; i<nvar; i++)
  {
    if ( strcmp(varList[i],"VtxProb") ) // VtxProb only exists for the dimuon, so mupl and mumi branches don't exist
    {
      fChain->SetBranchStatus(Form("Reco_QQ_mupl_%s",varList[i]),1);
      fChain->SetBranchStatus(Form("Reco_QQ_mumi_%s",varList[i]),1);
    }
  }
  //=========
}


//_____________________________
Double_t muIDCutsOptim::GetNColl(Int_t centr)
{
  // Returns the corresponding Ncoll value to the "centr" centrality bin
  
  if ( fIspp ) return 1.;
    
  Int_t normCent = TMath::Nint(centr/2.);
  
  Int_t lcent = 0;
  Int_t ucent = 0;
  for ( Int_t i = 0 ; i < fCentBins ; i++ )
  {
    ucent = fCentBinning[i];
    if ( (normCent >= lcent) && (normCent < ucent) ) return fCentMap[ucent];
    else lcent = ucent;
  }
  return 1.;
}


//_____________________________
void muIDCutsOptim::SetCentralityMap(const char* file)
{
  // Creates a mapping between centrality and Ncoll, based on a text file (taken from: https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideHeavyIonCentrality)
  
  if ( strlen(file) > 0 )
  {
    char line[1024];
    ifstream in(file);
    float lcent;
    float ucent;
    float Ncoll;
    
    fCentBins = 0;
    while ( in.getline(line,1024,'\n'))
    {
      sscanf(line,"%f %f %f",&lcent,&ucent,&Ncoll);
      
      fCentMap[ucent] = Ncoll;
      fCentBinning[fCentBins++] = ucent;
    }
    if ( fCentBins == 0 ) std::cout << "No centrality map could be defined: The file provided is empty" << std::endl;
  }
  else
  {
    fCentBins = 0;
    std::cout << "No centrality map could be defined: No file provided" << std::endl;
  }
}


//_____________________________
Bool_t muIDCutsOptim::isSoftMuon(const double varVal [ ], Int_t cutType)
{
  // Soft MuID cuts (https://github.com/CMS-HIN-dilepton/cmssw/blob/Onia_HI_75X/HiAnalysis/HiOnia/plugins/HiOniaAnalyzer.cc#L1398):
  //     - varVal[0] = isGoodMuon > 0
  //     - varVal[1] = highPurity > 0
  //     - varVal[9] = nPixWMea > 0
  //     - valVal[10] = nTrkWMea > 5
  //     - varVal[12] = dxy < 0.3
  //     - varVal[14] = dz < 20
  //
  // Depending on cutType variable, the corresponding cut is removed but the others are applied
  
  Double_t cutHP(0.);
  if ( !fIncludeHighPurity ) cutHP = -1.; // If the HighPurity cut is not included we set the cut to -1. (equivalent to no cut)
  
  switch (cutType)
  {
    case 1: // No isGoodMuon requirement
      if ( (varVal[1] > cutHP) && (varVal[9] > 0.) && (varVal[10] > 5.) && (varVal[12] < 0.3)  && (varVal[14] < 20.) ) return kTRUE;
      else return kFALSE;
      break;
      
    case 2: // No highPurity requirement
      if ( (varVal[0] > 0.) && (varVal[9] > 0.) && (varVal[10] > 5.) && (varVal[12] < 0.3)  && (varVal[14] < 20.) ) return kTRUE;
      else return kFALSE;
      break;
      
    case 3: // No nPixWMea requirement
      if ( (varVal[0] > 0.) && (varVal[1] > cutHP) && (varVal[10] > 5.) && (varVal[12] < 0.3)  && (varVal[14] < 20.) ) return kTRUE;
      else return kFALSE;
      break;
      
    case 4: // No nTrkWMea requirement
      if ( (varVal[0] > 0.) && (varVal[1] > cutHP) && (varVal[9] > 0.) && (varVal[12] < 0.3)  && (varVal[14] < 20.) ) return kTRUE;
      else return kFALSE;
      break;
      
    case 5: // No dxy requirement
      if ( (varVal[0] > 0.) && (varVal[1] > cutHP) && (varVal[9] > 0.) && (varVal[10] > 5.) && (varVal[14] < 20.) ) return kTRUE;
      else return kFALSE;
      break;
      
    case 6: // No dz requirement
      if ( (varVal[0] > 0.) && (varVal[1] > cutHP) && (varVal[9] > 0.) && (varVal[10] > 5.) && (varVal[12] < 0.3) ) return kTRUE;
      else return kFALSE;
      break;
      
    default: // Full soft muon cut
      if ( (varVal[0] > 0.) && (varVal[1] > cutHP) && (varVal[9] > 0.) && (varVal[10] > 5.) && (varVal[12] < 0.3)  && (varVal[14] < 20.) ) return kTRUE;
      else return kFALSE;
      break;
  }

}


//_____________________________
Bool_t muIDCutsOptim::isTriggerSelected(Int_t diMuN)
{
  // Checks if the event and dimuon are selected by the trigger
  
  return ( ((HLTriggers&(ULong64_t)pow(2,0))==pow(2,0)) && ((Reco_QQ_trig[diMuN]&(ULong64_t)pow(2,0))==pow(2,0)) );
}


//_____________________________
Int_t muIDCutsOptim::whichMuIDVar(const char* varName)
{
  // Checks if the event and dimuon are selected by the trigger
  
  if ( !strcmp(varName,"isGoodMuon") ) return 1; // Otherwise all the cuts except one are aplied
  if ( fIncludeHighPurity && !strcmp(varName,"highPurity") ) return 2; // If fIncludeHighPurity =kFALSE, it returns 0 for varName = "highPurity"
  if ( !strcmp(varName,"nPixWMea") ) return 3;
  if ( !strcmp(varName,"nTrkWMea") ) return 4;
  if ( !strcmp(varName,"dxy") ) return 5;
  if ( !strcmp(varName,"dz") ) return 6;
  else return 0;
}

