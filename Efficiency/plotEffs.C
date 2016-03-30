#include "TFile.h"
#include "TH1F.h"
#include "TGraphAsymmErrors.h"
#include "TString.h"
#include "TCanvas.h"
#include "TLegend.h"

#include "../Fitter/Macros/CMS/CMS_lumi.C"
#include "../Fitter/Macros/CMS/tdrstyle.C"

#include <iostream>
#include <fstream>

using namespace std;

void setErr(TH1F *hist);
void fixCentPp(TH1F *hist);
TH1F* integrateHist(TH1F *hist);
void printHist(vector<TH1F*> hist, const char* filename);
void printGraph(vector<TGraphAsymmErrors*> tg, const char* filename);
void inittex(const char* filename, string xname, vector<string> yname);
void closetex(const char* filename);

void plotEffs() {
   TFile *fjpsi_pp = new TFile("files/histos_jpsi_pp.root");
   TFile *fpsi2s_pp = new TFile("files/histos_psi2s_pp.root");
   TFile *fnpjpsi_pp = new TFile("files/histos_npjpsi_pp.root");
   TFile *fjpsi_pbpb = new TFile("files/histos_jpsi_pbpb.root");
   TFile *fpsi2s_pbpb = new TFile("files/histos_psi2s_pbpb.root");
   TFile *fnpjpsi_pbpb = new TFile("files/histos_npjpsi_pbpb.root");

   ofstream file_nocut("files/syst_PbPb_eff_MCstat_nocut.csv");
   ofstream file_ctaucut("files/syst_PbPb_eff_MCstat_ctaucut.csv");
   file_nocut << "MC statistics in efficiency (no ctau cut)" << endl;
   file_ctaucut << "MC statistics in efficiency (with ctau cut)" << endl;

   // first, let's draw simple efficiencies
   // we'll draw on the same plot the efficiencies for prompt and non-prompt J/psi, and psi(2S)
   // both for pp and pbpb, and both for the pp and centrality dependences, and both for midrapidity and forward rapidities, and both with and without ctau cut (2*2*2*2 = 16 plots)

   TString colltag, deptag, raptag, cuttag;
   for (int icoll=0; icoll<2; icoll++) {
      colltag = (icoll==0) ? "pp" : "pbpb";
      TString name = "files/histos_jpsi_" + colltag + ".root";
      TFile *fjpsi = new TFile(name);
      name = "files/histos_psi2s_" + colltag + ".root";
      TFile *fpsi2s = new TFile(name);
      name = "files/histos_npjpsi_" + colltag + ".root";
      TFile *fnpjpsi = new TFile(name);

      for (int idep=0; idep<3; idep++) { // idep=2 -> integrated
         deptag = (idep==0) ? "pt" : "cent";

         for (int irap=0; irap<2; irap++) {
            raptag = (irap==0) ? "mid" : "fwd";

            for (int icut=0; icut<2; icut++) {
               cuttag = (icut==0) ? "_" : "cut_";

               setTDRStyle();
               gStyle->SetEndErrorSize(3);
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

               if (idep==2) {
                  integrateHist(hjpsinum);
                  integrateHist(hpsi2snum);
                  integrateHist(hnpjpsinum);
                  integrateHist(hjpsiden);
                  integrateHist(hpsi2sden);
                  integrateHist(hnpjpsiden);
               }

               if (icoll==0 && idep==1) { // centrality for pp... fill all the bins
                  fixCentPp(hjpsinum);
                  fixCentPp(hpsi2snum);
                  fixCentPp(hnpjpsinum);
                  fixCentPp(hjpsiden);
                  fixCentPp(hpsi2sden);
                  fixCentPp(hnpjpsiden);
               }

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
               TLatex tl; TString cname;

               if (idep<2) {
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

                  tl.DrawLatex((idep==0) ? 1.5 : 10, 0.9, colltag + TString(", ") 
                        + ((irap==0) ? "|y|<1.6" : "|y|>1.6") + TString(", ") 
                        + ((icut==0) ? "no cut" : "ctau3D cut"));

                  cname = "files/singleff_" + colltag + "_" + deptag + "_" + raptag + "_" + cuttag;

                  c1->SaveAs(cname + ".root");
                  c1->SaveAs(cname + ".png");
                  c1->SaveAs(cname + ".pdf");
               }


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

               if (icoll==0 && idep==1) { // centrality for pp... fill all the bins
                  fixCentPp(hjpsipp);
                  fixCentPp(hpsi2spp);
                  fixCentPp(hjpsidenpp);
                  fixCentPp(hpsi2sdenpp);
               }

               // hjpsipp->Sumw2(true); hjpsidenpp->Sumw2(true); hjpsipp->Divide(hjpsidenpp);
               // hpsi2spp->Sumw2(true); hpsi2sdenpp->Sumw2(true); hpsi2spp->Divide(hpsi2sdenpp);
               setErr(hjpsipp); setErr(hjpsidenpp); hjpsipp->Divide(hjpsidenpp);
               setErr(hpsi2spp); setErr(hpsi2sdenpp); hpsi2spp->Divide(hpsi2sdenpp);
               hjpsipbpb->Divide(hjpsidenpbpb);
               hpsi2spbpb->Divide(hpsi2sdenpbpb);
               hpsi2spp->Divide(hjpsipp);
               hpsi2spp->SetMarkerColor(kBlack);
               hpsi2spp->SetLineColor(kBlack);
               hpsi2spbpb->Divide(hjpsipbpb);
               hpsi2spbpb->SetMarkerColor(kRed);
               hpsi2spbpb->SetLineColor(kRed);

               if (idep<2) {
                  haxes->GetYaxis()->SetTitle("Eff(#psi(2S)) / Eff(J/#psi)");
                  haxes->GetYaxis()->SetRangeUser(0.5,1.5);
                  haxes->SetBinContent(1,1);
                  haxes->Draw();
                  hpsi2spp->Draw("E1 same");
                  hpsi2spbpb->Draw("E1 same");

                  TLegend *tleg2 = new TLegend(0.7,0.17,0.89,0.31);
                  tleg2->SetBorderSize(0);
                  tleg2->AddEntry(hpsi2spp,"pp","lp");
                  tleg2->AddEntry(hpsi2spbpb,"pbpb","lp");
                  tleg2->Draw();

                  tl.DrawLatex((idep==0) ? 1.5 : 10, 1.4, ((irap==0) ? "|y|<1.6" : "|y|>1.6") + TString(", ") 
                        + ((icut==0) ? "no cut" : "ctau3D cut"));

                  cname = "files/simpleratio_" + deptag + "_" + raptag + "_" + cuttag;

                  c1->SaveAs(cname + ".root");
                  c1->SaveAs(cname + ".png");
                  c1->SaveAs(cname + ".pdf");
               }


               // at last, the double ratio
               hpsi2spbpb->Divide(hpsi2spp);

               if (idep<2) {
                  haxes->GetYaxis()->SetTitle("[Eff(#psi(2S)) / Eff(J/#psi)]_{PbPb} / [Eff(#psi(2S)) / Eff(J/#psi)]_{pp}");
                  haxes->GetYaxis()->SetTitleSize(0.04);
                  haxes->GetYaxis()->SetTitleOffset(2);
                  haxes->GetYaxis()->SetRangeUser(0.5,1.5);
                  haxes->Draw();
                  hpsi2spbpb->Draw("E1 same");
                  tl.DrawLatex((idep==0) ? 1.5 : 10, 1.4, ((irap==0) ? "|y|<1.6" : "|y|>1.6") + TString(", ") 
                        + ((icut==0) ? "no cut" : "ctau3D cut"));
                  cname = "files/doubleratio_" + deptag + "_" + raptag + "_" + cuttag;
                  c1->SaveAs(cname + ".root");
                  c1->SaveAs(cname + ".png");
                  c1->SaveAs(cname + ".pdf");
               }

               // print the uncertainty values to the csv
               ofstream *file = (icut==0) ? &file_nocut : &file_ctaucut;
               double rapmin, rapmax, ptmin, ptmax, centmin, centmax, value;
               rapmin = (irap==0) ? 0 : 1.6;
               rapmax = (irap==0) ? 1.6 : 2.4;
               if (idep==0) {
                  centmin = 0;
                  centmax = 200;
                  for (int ibin=1; ibin<hpsi2spbpb->GetNbinsX()+1; ibin++) {
                     ptmin = hpsi2spbpb->GetXaxis()->GetBinLowEdge(ibin);
                     ptmax = hpsi2spbpb->GetXaxis()->GetBinUpEdge(ibin);
                     value = hpsi2spbpb->GetBinError(ibin);
                     *file << rapmin << ", " << rapmax << ", " << ptmin << ", " << ptmax << ", " << centmin << ", " << centmax << ", " << value << endl;
                  }
               } else if (idep==1) {
                  ptmin = (irap==0) ? 6.5 : 3;
                  ptmax = 30;
                  for (int ibin=1; ibin<hpsi2spbpb->GetNbinsX()+1; ibin++) {
                     centmin = hpsi2spbpb->GetXaxis()->GetBinLowEdge(ibin);
                     centmax = hpsi2spbpb->GetXaxis()->GetBinUpEdge(ibin);
                     value = hpsi2spbpb->GetBinError(ibin);
                     *file << rapmin << ", " << rapmax << ", " << ptmin << ", " << ptmax << ", " << centmin << ", " << centmax << ", " << value << endl;
                  }
               }
               else {
                  ptmin = (irap==0) ? 6.5 : 3;
                  ptmax = 30;
                  centmin = 0;
                  centmax = 200;
                  value = hpsi2spbpb->GetBinError(1);
                  *file << rapmin << ", " << rapmax << ", " << ptmin << ", " << ptmax << ", " << centmin << ", " << centmax << ", " << value << endl;
               }

               // clean behind ourselves
               delete c1;
               delete tg_jpsi; delete tg_psi2s; delete tg_npjpsi;
               delete haxes;
            } // icut loop (without / with ctau cut)
         } // irap loop (mid / fwd)
      } // idep loop (pt / centrality)
   } // icoll loop (pp / pbpb)

   file_nocut.close();
   file_ctaucut.close();
}

void setErr(TH1F *hist) {
   int nbins = hist->GetNbinsX();
   for (int i=1; i<nbins+1; i++) {
      hist->SetBinError(i,sqrt(hist->GetBinContent(i)));
   }
}

void fixCentPp(TH1F *hist) {
   int nbins = hist->GetNbinsX();
   float y = hist->GetBinContent(1);
   float dy = hist->GetBinError(1);
   for (int i=2; i<nbins+1; i++) {
      hist->SetBinContent(i,y);
      hist->SetBinError(i,dy);
   }
}

TH1F* integrateHist(TH1F *hist) {
   TString name = hist->GetName(); name = name + "_int";
   TString title = hist->GetTitle(); title = title + " integrated";
   double integral, integralerror;
   int nbins = hist->GetNbinsX();
   integral = hist->IntegralAndError(1,nbins,integralerror);
   TH1F *ans = new TH1F(name, title, 1, hist->GetXaxis()->GetBinLowEdge(1), hist->GetXaxis()->GetBinUpEdge(nbins));
   ans->SetBinContent(1,integral);
   ans->SetBinError(1,integralerror);
   return ans;
}

void printHist(vector<TH1F*> hist, const char* filename) {
   ofstream file(filename, ofstream::app);
   file.setprecision(3);
   if (hist.size()==0) return;
   int nbins = hist[0]->GetNbinsX();
   for (int ibin=1; ibin<nbins+1; ibin++) {
      file << "[" << hist[0]->GetXaxis()->GetBinLowEdge(ibin) << "-" << hist[0]->GetXaxis()->GetBinUpEdge(ibin) << "]";
      for (vector<TH1F*>::const_iterator ith=hist.begin(); ith!=hist.end(); ith++) {
         file << " & $" << ith->GetBinContent(ibin) << " \\pm " << ith->GetBinError(ibin) << "$";
      }
      file << "\\\\" << endl;
   }
   delete hist_int;
   file.close();
}

void printGraph(vector<TGraphAsymmErrors*> tg, const char* filename) {
   ofstream file(filename, ofstream::app);
   file.setprecision(3);
   if (hist.size()==0) return;
   int nbins = tg[0]->GetN();
   for (int ibin=0; ibin<nbins; ibin++) {
      file << "[" << tg[0]->GetX()[ibin]-tg[0]->GetErrorXlow(ibin) << "-" << tg[0]->GetX()[ibin]-tg[0]->GetErrorXhigh(ibin);
      for (vector<TGraphAsymmErrors*>::const_iterator itg=tg.begin(); itg!=tg.end(); itg++) {
         file << " & $" << itg->GetY()[ibin] << "_{-" << itg->GetErrorYlow(ibin) << "}^{+" << itg->GetErrorYhigh(ibin) << "} $";
      }
      file << "\\\\" << endl;
   }
   delete hist_int;
   file.close();
}

void inittex(const char* filename, string xname, vector<string> yname) {
   ofstream file(filename);
   file << "\\begin{tabular}{c"; 
   for (int i=0; i<yname.size; i++) file << "c";
   file << "}" << endl;
   file << "\\hline" << endl;
   file << xname;
   for (int i=0; i<yname.size(); i++) file << " & " << yname[i];
   file<< "\\\\" << endl;
   file << "\\hline" << endl;
   file.close();
}

void closetex(const char* filename) {
   ofstream file(filename, ofstream::app);
   file << "\\hline" << endl;
   file << "\\end{tabular}" << endl;
   file.close();
}
