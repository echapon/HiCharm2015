#include "TFile.h"
#include "TH1F.h"
#include "TGraphAsymmErrors.h"
#include "TString.h"

#include "../Fitter/Macros/CMS/CMS_lumi.C"
#include "../Fitter/Macros/CMS/tdrstyle.C"

#include <iostream>

using namespace std;

void plotEffs() {
   TFile *fjpsi_pp = new TFile("histos_jpsi_pp.root");
   TFile *fpsi2s_pp = new TFile("histos_psi2s_pp.root");
   TFile *fnpjpsi_pp = new TFile("histos_npjpsi_pp.root");
   TFile *fjpsi_pbpb = new TFile("histos_jpsi_pbpb.root");
   TFile *fpsi2s_pbpb = new TFile("histos_psi2s_pbpb.root");
   TFile *fnpjpsi_pbpb = new TFile("histos_npjpsi_pbpb.root");

   // first, let's draw simple efficiencies
   // we'll draw on the same plot the efficiencies for prompt and non-prompt J/psi, and psi(2S)
   // both for pp and pbpb, and both for the pp and centrality dependences, and both for midrapidity and forward rapidities, and both with and without ctau cut (2*2*2*2 = 16 plots)

   TString colltag, deptag, raptag, cuttag;
   for (int icoll=0; icoll<2; icoll++) {
      colltag = (icoll==0) ? "pp" : "pbpb";
      TString name = "histos_jpsi_" + colltag + ".root";
      TFile *fjpsi = new TFile(name);
      name = "histos_psi2s_" + colltag + ".root";
      TFile *fpsi2s = new TFile(name);
      name = "histos_npjpsi_" + colltag + ".root";
      TFile *fnpjpsi = new TFile(name);

      for (int idep=0; idep<=icoll; idep++) { // do not do centrality for pp
         deptag = (idep==0) ? "pt" : "cent";

         for (int irap=0; irap<2; irap++) {
            raptag = (irap==0) ? "mid" : "fwd";

            for (int icut=0; icut<2; icut++) {
               cuttag = (icut==0) ? "_" : "cut_";

               setTDRStyle();
               TCanvas *c1 = new TCanvas();

               TString hname;
               hname = "hnum" + cuttag + deptag + raptag;
               TH1F *hjpsinum = (TH1F*) fjpsi->Get(hname);
               TH1F *hpsi2snum = (TH1F*) fpsi2s->Get(hname);
               TH1F *hnpjpsinum = (TH1F*) fnpjpsi->Get(hname);
               hname = "hden_" + deptag + raptag;
               TH1F *hjpsiden = (TH1F*) fjpsi->Get(hname);
               TH1F *hpsi2sden = (TH1F*) fpsi2s->Get(hname);
               TH1F *hnpjpsiden = (TH1F*) fnpjpsi->Get(hname);

               TGraphAsymmErrors *tg_jpsi = new TGraphAsymmErrors(hjpsinum,hjpsiden,(icoll==0) ? "" : "norm");
               tg_jpsi->SetMarkerColor(kBlack);
               tg_jpsi->SetLineColor(kBlack);
               TGraphAsymmErrors *tg_psi2s = new TGraphAsymmErrors(hpsi2snum,hpsi2sden,(icoll==0) ? "" : "norm");
               tg_psi2s->SetMarkerColor(kRed);
               tg_psi2s->SetLineColor(kRed);
               TGraphAsymmErrors *tg_npjpsi = new TGraphAsymmErrors(hnpjpsinum,hnpjpsiden,(icoll==0) ? "" : "norm");
               tg_npjpsi->SetMarkerColor(kBlue);
               tg_npjpsi->SetLineColor(kBlue);

               TH1F *haxes = new TH1F("haxes","haxes",1,0,(idep==1) ? 200 : 30);
               haxes->GetYaxis()->SetTitle("Efficiency");
               haxes->GetXaxis()->SetTitle((idep==1) ? "Centrality bin" : "p_{T}");

               haxes->Draw();
               tg_jpsi->Draw("P");
               tg_psi2s->Draw("P");
               tg_npjpsi->Draw("P");

               double yshift =0; if (icoll==1 && irap==1) yshift=0.3;
               TLegend *tleg = new TLegend(0.5,0.26+yshift,0.88,0.46+yshift);
               tleg->SetBorderSize(0);
               tleg->AddEntry(tg_jpsi,"J/#psi (prompt)","lp");
               tleg->AddEntry(tg_psi2s,"#psi(2S)","lp");
               tleg->AddEntry(tg_npjpsi,"J/#psi (non-prompt)","lp");
               tleg->Draw();

               TLatex tl; 
               tl.DrawLatex((idep==0) ? 1.5 : 10, 0.9, colltag + TString(", ") 
                     + ((irap==0) ? "|y|<1.6" : "|y|>1.6") + TString(", ") 
                     + ((icut==0) ? "no cut" : "ctau3D cut"));

               TString cname = "singleff_" + colltag + "_" + deptag + "_" + raptag + "_" + cuttag;

               c1->SaveAs(cname + ".root");
               c1->SaveAs(cname + ".png");
               c1->SaveAs(cname + ".pdf");


               // now, let's draw simple ratios of efficiencies: psi(2S)/J/psi
               // but do it only once
               if (icoll>0) continue;
               hname = "hnum" + cuttag + deptag + raptag;
               TH1F *hjpsipp = (TH1F*) fjpsi_pp->Get(hname); hjpsipp->SetName(TString(hjpsipp->GetName()) + "_jpsipp");
               TH1F *hpsi2spp = (TH1F*) fpsi2s_pp->Get(hname); hpsi2spp->SetName(TString(hpsi2spp->GetName()) + "_psi2spp");
               TH1F *hjpsipbpb = (TH1F*) fjpsi_pbpb->Get(hname); hjpsipbpb->SetName(TString(hjpsipbpb->GetName()) + "_jpsipbpb");
               TH1F *hpsi2spbpb = (TH1F*) fpsi2s_pbpb->Get(hname); hpsi2spbpb->SetName(TString(hpsi2spbpb->GetName()) + "_psi2spbpb");
               hname = "hden_" + deptag + raptag;
               TH1F *hjpsidenpp = (TH1F*) fjpsi_pp->Get(hname); hjpsidenpp->SetName(TString(hjpsidenpp->GetName()) + "_jpsidenpp");
               TH1F *hpsi2sdenpp = (TH1F*) fpsi2s_pp->Get(hname); hpsi2sdenpp->SetName(TString(hpsi2sdenpp->GetName()) + "_psi2sdenpp");
               TH1F *hjpsidenpbpb = (TH1F*) fjpsi_pbpb->Get(hname); hjpsipbpb->SetName(TString(hjpsidenpbpb->GetName()) + "_jpsidenpbpb");
               TH1F *hpsi2sdenpbpb = (TH1F*) fpsi2s_pbpb->Get(hname); hpsi2sdenpbpb->SetName(TString(hpsi2sdenpbpb->GetName()) + "_psi2sdenpbpb");
               hjpsidenpp->Sumw2(true); hjpsipp->Divide(hjpsidenpp);
               hpsi2sdenpp->Sumw2(true); hpsi2spp->Divide(hpsi2sdenpp);
               hjpsipbpb->Divide(hjpsidenpbpb);
               hpsi2spbpb->Divide(hpsi2sdenpbpb);
               cout << hpsi2spp->GetBinContent(2) << " / " << hjpsipp->GetBinContent(2) << " = ";
               hpsi2spp->Divide(hjpsipp);
               cout << hpsi2spp->GetBinContent(2) << endl;
               hpsi2spp->SetMarkerColor(kBlack);
               hpsi2spp->SetLineColor(kBlack);
               cout << hpsi2spbpb->GetBinError(2) << " / " << hjpsipbpb->GetBinError(2) << " = ";
               hpsi2spbpb->Divide(hjpsipbpb);
               cout << hpsi2spbpb->GetBinError(2) << endl;
               hpsi2spbpb->SetMarkerColor(kRed);
               hpsi2spbpb->SetLineColor(kRed);

               haxes->GetYaxis()->SetTitle("Eff(#psi(2S)) / Eff(J/#psi)");
               haxes->GetYaxis()->SetRangeUser(0.5,1.5);
               haxes->SetBinContent(1,1);
               haxes->Draw();
               hpsi2spp->Draw("E same");
               hpsi2spbpb->Draw("E same");

               TLegend *tleg2 = new TLegend(0.7,0.17,0.89,0.31);
               tleg2->SetBorderSize(0);
               tleg2->AddEntry(hpsi2spp,"pp","lp");
               tleg2->AddEntry(hpsi2spbpb,"pbpb","lp");
               tleg2->Draw();

               tl.DrawLatex((idep==0) ? 1.5 : 10, 1.4, ((irap==0) ? "|y|<1.6" : "|y|>1.6") + TString(", ") 
                     + ((icut==0) ? "no cut" : "ctau3D cut"));

               cname = "simpleratio_" + deptag + "_" + raptag + "_" + cuttag;

               c1->SaveAs(cname + ".root");
               c1->SaveAs(cname + ".png");
               c1->SaveAs(cname + ".pdf");


               // at last, the double raio
               hpsi2spbpb->Divide(hpsi2spp);
               haxes->GetYaxis()->SetTitle("[Eff(#psi(2S)) / Eff(J/#psi)]_{PbPb} / [Eff(#psi(2S)) / Eff(J/#psi)]_{pp}");
               haxes->GetYaxis()->SetTitleSize(0.04);
               haxes->GetYaxis()->SetTitleOffset(2);
               haxes->GetYaxis()->SetRangeUser(0.5,1.5);
               haxes->Draw();
               hpsi2spbpb->Draw("same");
               tl.DrawLatex((idep==0) ? 1.5 : 10, 1.4, ((irap==0) ? "|y|<1.6" : "|y|>1.6") + TString(", ") 
                     + ((icut==0) ? "no cut" : "ctau3D cut"));
               cname = "doubleratio_" + deptag + "_" + raptag + "_" + cuttag;
               c1->SaveAs(cname + ".root");
               c1->SaveAs(cname + ".png");
               c1->SaveAs(cname + ".pdf");

               // clean behind ourselves
               delete c1;
               delete tg_jpsi, tg_psi2s, tg_npjpsi;
               delete haxes;
               delete tleg, tleg2;
            } // icut loop (with / without ctau cut)
         } // irap loop (mid / fwd)
      } // idep loop (pt / centrality)
   } // icoll loop (pp / pbpb)
}
