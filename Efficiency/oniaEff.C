#define oniaEff_cxx
#include "oniaEff.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include "TFile.h"

#include <iostream>

// define bins
const int nbins_centmid = 6;
const float bins_centmid[nbins_centmid+1] = {0, 20, 40, 60, 80, 100, 200};
const int nbins_centfwd = 3;
const float bins_centfwd[nbins_centfwd+1] = {0, 40, 80, 200};
const int nbins_ptmid = 5;
const float bins_ptmid[nbins_ptmid+1] = {6.5, 9, 12, 15, 20, 30};
const int nbins_ptfwd = 3;
const float bins_ptfwd[nbins_ptfwd+1] = {3, 6.5, 12, 30};

// NP rejection cuts
const double ctaucutmid_pp = 0.03;
const double ctaucutfwd_pp = 0.05;
const double ctaucutmid_pbpb = 0.03;
const double ctaucutfwd_pbpb = 0.05;

// other settings
const double maxdr = 0.03;
const double massjpsi = 3.096;
const double masspsip = 3.686;
const double massdown = 0.25;
const double massup = 0.15;

using namespace HI;
using namespace std;

void oniaEff::Loop(const char* fname, bool ispbpb, bool isPsip)
{
//   In a ROOT session, you can do:
//      root> .L oniaEff.C
//      root> oniaEff t
//      root> t.GetEntry(12); // Fill t data members with entry number 12
//      root> t.Show();       // Show values of entry 12
//      root> t.Show(16);     // Read and show values of entry 16
//      root> t.Loop();       // Loop on all entries
//

//     This is the loop skeleton where:
//    jentry is the global entry number in the chain
//    ientry is the entry number in the current Tree
//  Note that the argument to GetEntry must be:
//    jentry for TChain::GetEntry
//    ientry for TTree::GetEntry and TBranch::GetEntry
//
//       To read only selected branches, Insert statements like:
// METHOD1:
   fChain->SetBranchStatus("*",0);  // disable all branches
   fChain->SetBranchStatus("Centrality",1);  
   fChain->SetBranchStatus("HLTriggers",1);  
   fChain->SetBranchStatus("Reco_QQ_*",1);  
   fChain->SetBranchStatus("Gen_QQ_*",1);  
// METHOD2: replace line
//    fChain->GetEntry(jentry);       //read all branches
//by  b_branchname->GetEntry(ientry); //read only this branch
   if (fChain == 0) return;

   Long64_t nentries = fChain->GetEntriesFast();

   TFile *f = new TFile(fname, "RECREATE");
   f->cd();

   // define the histos
   // we need 2 numerators (with and without ctau3d cut), 1 denominator. Let's make them in arrays
   TH1F *hnum_centmid = new TH1F("hnum_centmid","hnum_centmid",nbins_centmid,bins_centmid);
   TH1F *hnum_centfwd = new TH1F("hnum_centfwd","hnum_centfwd",nbins_centfwd,bins_centfwd);
   TH1F *hnum_ptmid = new TH1F("hnum_ptmid","hnum_ptmid",nbins_ptmid,bins_ptmid);
   TH1F *hnum_ptfwd = new TH1F("hnum_ptfwd","hnum_ptfwd",nbins_ptfwd,bins_ptfwd);
   TH1F *hnumcut_centmid = new TH1F("hnumcut_centmid","hnumcut_centmid",nbins_centmid,bins_centmid);
   TH1F *hnumcut_centfwd = new TH1F("hnumcut_centfwd","hnumcut_centfwd",nbins_centfwd,bins_centfwd);
   TH1F *hnumcut_ptmid = new TH1F("hnumcut_ptmid","hnumcut_ptmid",nbins_ptmid,bins_ptmid);
   TH1F *hnumcut_ptfwd = new TH1F("hnumcut_ptfwd","hnumcut_ptfwd",nbins_ptfwd,bins_ptfwd);
   TH1F *hden_centmid = new TH1F("hden_centmid","hden_centmid",nbins_centmid,bins_centmid);
   TH1F *hden_centfwd = new TH1F("hden_centfwd","hden_centfwd",nbins_centfwd,bins_centfwd);
   TH1F *hden_ptmid = new TH1F("hden_ptmid","hden_ptmid",nbins_ptmid,bins_ptmid);
   TH1F *hden_ptfwd = new TH1F("hden_ptfwd","hden_ptfwd",nbins_ptfwd,bins_ptfwd);

   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      // if (ientry>1000) break;

      // skip the event if Gen_QQ_size != 1 for now
      b_Gen_QQ_size->GetEntry(ientry); //read only this branch
      if (Gen_QQ_size !=1) continue;

      // apply gen acceptance cuts
      b_Gen_QQ_mupl_4mom->GetEntry(ientry); //read only this branch
      b_Gen_QQ_mumi_4mom->GetEntry(ientry); //read only this branch
      TLorentzVector *tlvgenpl = (TLorentzVector*) Gen_QQ_mupl_4mom->At(0);
      TLorentzVector *tlvgenmi = (TLorentzVector*) Gen_QQ_mumi_4mom->At(0);
      if (!isGlobalMuonInAccept2015(tlvgenpl) || !isGlobalMuonInAccept2015(tlvgenmi)) continue;

      // fill denominators
      b_Gen_QQ_4mom->GetEntry(ientry); //read only this branch
      TLorentzVector *tlvgenqq = (TLorentzVector*) Gen_QQ_4mom->At(0);
      double genpt = tlvgenqq->Pt();
      b_Centrality->GetEntry(ientry); //read only this branch

      double weight = ispbpb ? fChain->GetWeight()*findNcoll(Centrality) : 1.;

      if (fabs(tlvgenqq->Rapidity()) < 1.6) {
         hden_centmid->Fill(Centrality,weight);
         hden_ptmid->Fill(genpt,weight);
      } else if (fabs(tlvgenqq->Rapidity()) < 2.4) {
         hden_centfwd->Fill(Centrality,weight);
         hden_ptfwd->Fill(genpt,weight);
      } else continue;

      // is there at least one reco dimuon?
      b_Reco_QQ_size->GetEntry(ientry);
      if (Reco_QQ_size==0) continue;

      // ok, now read all other branches
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;

      // loop on the reconstructed dimuons to find the one matched to the gen one
      double mindr=999.;
      int ibestqq=-1;
      for (int i=0; i<Reco_QQ_size; i++) {
         // acceptance
         if (!areMuonsInAcceptance2015(i)) continue;
         // sign
         if (Reco_QQ_sign[i]!=0) continue;
         // quality cuts
         if (!passQualityCuts2015(i)) continue;
         // trigger matching
         if (!isTriggerMatch(i,0)) continue; // HLT_HIL1DoubleMu0_v1
         // mass cut
         double mass = ((TLorentzVector*) Reco_QQ_4mom->At(i))->M();
         double mass0 = isPsip ? masspsip : massjpsi;
         if (mass<(mass0-massdown) || mass>(mass0+massup)) continue;
         // gen-reco matching
         TLorentzVector *tlvrecpl = (TLorentzVector*) Reco_QQ_mupl_4mom->At(i);
         TLorentzVector *tlvrecmi = (TLorentzVector*) Reco_QQ_mumi_4mom->At(i);
         double dr = max(tlvrecpl->DeltaR(*tlvgenpl),tlvrecmi->DeltaR(*tlvgenmi));
         if (dr<mindr) {
            mindr = dr;
            ibestqq = i;
         }
      } // Reco_QQ loop

      if (ibestqq<0) continue;

      // fill the numerators
      if (fabs(tlvgenqq->Rapidity()) < 1.6) {
         hnum_centmid->Fill(Centrality,weight);
         hnum_ptmid->Fill(genpt,weight);
      } else {
         hnum_centfwd->Fill(Centrality,weight);
         hnum_ptfwd->Fill(genpt,weight);
      }

      // apply ctau3D cut
      bool ctaucutok = false;
      TLorentzVector *tlvrecqq = (TLorentzVector*) Reco_QQ_4mom->At(ibestqq);
      if (ispbpb) {
         if (fabs(tlvrecqq->Rapidity()) < 1.6) {
            ctaucutok = (Reco_QQ_ctau3D[ibestqq] < ctaucutmid_pbpb);
         } else {
            ctaucutok = (Reco_QQ_ctau3D[ibestqq] < ctaucutfwd_pbpb);
         }
      } else {
         if (fabs(tlvrecqq->Rapidity()) < 1.6) {
            ctaucutok = (Reco_QQ_ctau3D[ibestqq] < ctaucutmid_pp);
         } else {
            ctaucutok = (Reco_QQ_ctau3D[ibestqq] < ctaucutfwd_pp);
         }
      }

      if (!ctaucutok) continue;
      if (fabs(tlvgenqq->Rapidity()) < 1.6) {
         hnumcut_centmid->Fill(Centrality,weight);
         hnumcut_ptmid->Fill(genpt,weight);
      } else {
         hnumcut_centfwd->Fill(Centrality,weight);
         hnumcut_ptfwd->Fill(genpt,weight);
      }

   } // event loop

   f->Write();
   f->Close();
}
