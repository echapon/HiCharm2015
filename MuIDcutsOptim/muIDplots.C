#include <Riostream.h>
#include "TFile.h"
#include "TH1.h"
#include "TH1F.h"
#include "TCanvas.h"
#include "TGaxis.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TMath.h"
#include "TLatex.h"
#include "TLine.h"
#include "TString.h"
#include "TObjString.h"

Double_t bkgNorm = 0.5;

Bool_t isHPincluded(kFALSE);
TString fKinRange("");

//___________________________________________________________
Double_t GetSoftMuCutValue(const char* varName)
{
  // Returns the cut value for the desired variable in soft muon cut set
  
  if ( !strcmp(varName,"isGoodMuon") ) return 1.;
  else if  ( isHPincluded && !strcmp(varName,"highPurity") ) return 1.;
  else if  ( !strcmp(varName,"nPixWMea") ) return 1.;
  else if  ( !strcmp(varName,"nTrkWMea") ) return 6.;
  else if  ( !strcmp(varName,"dxy") ) return -0.3;
  else if  ( !strcmp(varName,"dz") ) return -20;
  else
  {
    std::cout << "The variable " << varName << " does not belong to the soft muon cut set" << std::endl;
    return 999.;
  }
}


//___________________________________________________________
void plotSingleCutEfficiency(TFile* fsig, TFile* fbkg, const char* varname [ ])
{
  // Computes the single cut efficiencies of signal and background, and the significance.
  
  TObjArray* aSig = static_cast<TObjArray*>(fsig->FindObjectAny("SigSingleCut"));
  TObjArray* aBSig = static_cast<TObjArray*>(fsig->FindObjectAny("BkgSingleCut"));
  
  TObjArray* aSBkg = static_cast<TObjArray*>(fbkg->FindObjectAny("SigSingleCut"));
  TObjArray* aBkg = static_cast<TObjArray*>(fbkg->FindObjectAny("BkgSingleCut"));

  TString sign("");
  for (int i=0; i<19; i++) {
    TCanvas *c1 = new TCanvas();
    TH1F *hsig = (TH1F*) aSig->FindObject(Form("h%s_sig",varname[i])); // Signal from MC
    TH1F *hSbkg = (TH1F*) aSBkg->FindObject(Form("h%s_sig",varname[i])); // Signal from Data
    TH1F *hbkg = (TH1F*) aBkg->FindObject(Form("h%s_bkg",varname[i])); // Background from Data
    TH1F *hbkg2 = (TH1F*) aBSig->FindObject(Form("h%s_bkg",varname[i]));  // Background from MC
    
    
    if (!hsig || !hSbkg || !hbkg || !hbkg2 )
    {
      std::cout << "No histos found" << std::endl;
      return;
    }
    
    double nbkg = hbkg->GetMaximum(); // Maximum data bkg (total)
    double nbkg2 = hbkg2->GetMaximum(); // Maximum MC bkg (total)
    double nsig = hsig->GetMaximum()- bkgNorm*nbkg2; // MC signal (total) substracted from MC bkg
    double nSbkg = hSbkg->GetMaximum()- bkgNorm*nbkg; // Data signal (total) substracted from MC bkg

    Double_t scale_sig = nSbkg/nsig;
    nsig *= scale_sig;
    
    int nbins = hsig->GetNbinsX();
    double binmin = hsig->GetXaxis()->GetXmin();
    double binmax = hsig->GetXaxis()->GetXmax();
    
    TH1F *hsigeff = new TH1F(Form("hsigeff_%s",varname[i]),Form("Signal efficiency, bkg rejection and significance;%s;Signal efficiency",varname[i]),nbins,binmin,binmax);
    hsigeff->SetLineColor(kRed);
    TH1F *hbkgeff = new TH1F(Form("hbkgeff_%s",varname[i]),Form("hbkgeff;%s;1 - Background efficiency",varname[i]),nbins,binmin,binmax);
    hbkgeff->SetLineColor(kBlue);
    TH1F *hsignif = new TH1F(Form("hsignif_%s",varname[i]),Form("hsignif;%s;S / #sqrt{S+B}",varname[i]),nbins,binmin,binmax);
    hsignif->SetLineColor(kGreen);

    for (int j=1; j<nbins+1; j++) {
      double nbkgj = hbkg->GetBinContent(j);
      double nbkgj2 = hbkg2->GetBinContent(j);
      double nsigj = scale_sig*(hsig->GetBinContent(j) - bkgNorm*nbkgj2);

      double effSig = nsigj/nsig;
      if ( effSig > 1.0 ) // In pp MC, since there is no injected bkg, the bkg substraction can lead to eff slightly bigger than 1
      {
        if ( effSig < 1.01 ) effSig = 1.0;
        else cout << "Error: the efficiency is bigger than 100%" << endl;
      }
      hsigeff->SetBinContent(j,effSig);;
      hbkgeff->SetBinContent(j,1.-nbkgj/nbkg);
      hsignif->SetBinContent(j,nsigj+nbkgj>0 ? nsigj/sqrt(nsigj+nbkgj) : 0);
    }
    
    
    // Draw the values of the efficiencies at the soft muon cut value
    Double_t cutVal = GetSoftMuCutValue(varname[i]);
    sign.Clear();
    if ( cutVal < 0. )
    {
      cutVal = TMath::Abs(cutVal) - 1E-6;
      sign += "<";
    }
    else sign += ">";
    
    Bool_t printEff(kFALSE);
    if ( cutVal < 100. ) printEff = kTRUE;
    
    TLatex *  text(0x0);
    if ( printEff )
    {
      Double_t sigEff = hsigeff->GetBinContent(hsigeff->FindBin(cutVal));
      Double_t bkgEff = hbkgeff->GetBinContent(hbkgeff->FindBin(cutVal));
      
      text = new TLatex(0.1451149,0.5835095,Form("#varepsilon_{sig} [%%] =  %2.2f ; Rej_{bkg} [%%] =  %2.2f (%s %s %2.1f)",sigEff*100.,bkgEff*100.,varname[i],sign.Data(),cutVal));
      text->SetNDC();
      text->SetTextFont(42);
      text->SetTextSize(0.0422833);
      text->SetLineWidth(2);
    }
    
    TLatex *textKin = new TLatex(0.1465517,0.653277,fKinRange.Data());
    textKin->SetNDC();
    textKin->SetTextFont(42);
    textKin->SetTextSize(0.0422833);
    textKin->SetLineWidth(2);
  
    
    hsigeff->GetYaxis()->SetRangeUser(0,1.05);
    hsigeff->GetXaxis()->SetLabelSize(0.05);
    hsigeff->GetXaxis()->SetTitleSize(0.05);
    hsigeff->GetXaxis()->SetTitleOffset(0.93);
    hsigeff->GetYaxis()->SetLabelSize(0.05);
    hsigeff->GetYaxis()->SetTitleSize(0.05);
    hsigeff->GetYaxis()->SetTitleOffset(0.9);
    
    hsigeff->Draw("");
    hbkgeff->Draw("same");
    c1->Update();
    
    // scale hsignif to the pad coordinates
    Float_t rightmax = 1.1*hsignif->GetMaximum();
    Float_t rightmin = 0.98*hsignif->GetMaximum();
    Float_t scale = gPad->GetUymax()/rightmax;
    hsignif->Scale(scale);
    hsignif->Draw("same");
    // draw an axis on the right side
    TGaxis *axis = new TGaxis(gPad->GetUxmax(),gPad->GetUymin(),
                              gPad->GetUxmax(), gPad->GetUymax(),rightmin,rightmax,510,"+L");
    axis->SetLineColor(kGreen);
    axis->SetLabelColor(kGreen);
    axis->Draw();
    
    TLegend *leg = new TLegend(0.454023,0.3002114,0.8721264,0.4904863);
    leg->SetBorderSize(0);
    leg->SetLineColor(1);
    leg->SetLineStyle(1);
    leg->SetLineWidth(1);
    leg->SetFillColor(0);
    leg->SetFillStyle(0);
    leg->AddEntry(hsigeff,"Signal efficiency","lp");
    leg->AddEntry(hbkgeff,"1 - Background efficiency","lp");
    leg->AddEntry(hsignif,"S / #sqrt{S+B}","lp");
    leg->Draw();
    
    if ( text ) text->Draw("same");
    textKin->Draw("same");
    
    c1->Update();
    c1->SaveAs(Form("SingleCutsEfficiency_%s.pdf",varname[i]));
  }
}


//___________________________________________________________
void plotSingleCutEfficiencyWSoftMuCut(TFile* fsig, TFile* fbkg, const char* varname [ ])
{
  // Computes the single cut efficiencies (after aplying the soft muon cuts) of signal and background, and the significance.
  // If the variable belongs to the soft muon cuts, all the soft muon cuts except the variable are applied
  
  TObjArray* aSig = static_cast<TObjArray*>(fsig->FindObjectAny("SigSoftMuCut"));
  TObjArray* aBSig = static_cast<TObjArray*>(fsig->FindObjectAny("BkgSoftMuCut"));
  
  TObjArray* aSBkg = static_cast<TObjArray*>(fbkg->FindObjectAny("SigSoftMuCut"));
  TObjArray* aBkg = static_cast<TObjArray*>(fbkg->FindObjectAny("BkgSoftMuCut"));
  
  TString sign("");
  for (int i=0; i<19; i++) {
    TCanvas *c1 = new TCanvas();
    TH1F *hsig = (TH1F*) aSig->FindObject(Form("hSoftMuCut_%s_sig",varname[i]));
    TH1F *hSbkg = (TH1F*) aSBkg->FindObject(Form("hSoftMuCut_%s_sig",varname[i]));
    TH1F *hbkg = (TH1F*) aBkg->FindObject(Form("hSoftMuCut_%s_bkg",varname[i]));
    TH1F *hbkg2 = (TH1F*) aBSig->FindObject(Form("hSoftMuCut_%s_bkg",varname[i]));
    
    
    if (!hsig || !hSbkg || !hbkg || !hbkg2 )
    {
      std::cout << "No histos found" << std::endl;
      return;
    }
    
    double nbkg = hbkg->GetMaximum();
    double nbkg2 = hbkg2->GetMaximum();
    double nsig = hsig->GetMaximum()- bkgNorm*nbkg2;
    double nSbkg = hSbkg->GetMaximum()- bkgNorm*nbkg;
    
    Double_t scale_sig = nSbkg/nsig;
    nsig *= scale_sig;
    
    int nbins = hsig->GetNbinsX();
    double binmin = hsig->GetXaxis()->GetXmin();
    double binmax = hsig->GetXaxis()->GetXmax();
    
    TH1F *hsigeff = new TH1F(Form("hSoftMuCut_sigeff_%s",varname[i]),Form("Signal efficiency, bkg rejection and significance (soft #mu-ID cuts applied);%s;Signal efficiency",varname[i]),nbins,binmin,binmax);
    hsigeff->SetLineColor(kRed);
    TH1F *hbkgeff = new TH1F(Form("hSoftMuCut_bkgeff_%s",varname[i]),Form("hSoftMuCut_bkgeff;%s;1 - Background efficiency",varname[i]),nbins,binmin,binmax);
    hbkgeff->SetLineColor(kBlue);
    TH1F *hsignif = new TH1F(Form("hSoftMuCut_signif_%s",varname[i]),Form("hSoftMuCut_signif;%s;S / #sqrt{S+B}",varname[i]),nbins,binmin,binmax);
    hsignif->SetLineColor(kGreen);
    
    for (int j=1; j<nbins+1; j++) {
      double nbkgj = hbkg->GetBinContent(j);
      double nbkgj2 = hbkg2->GetBinContent(j);
      double nsigj = scale_sig*(hsig->GetBinContent(j) - bkgNorm*nbkgj2);
      
      double effSig = nsigj/nsig;
      if ( effSig > 1.0 ) // In pp MC, since there is no injected bkg, the bkg substraction can lead to eff slightly bigger than 1
      {
        if ( effSig < 1.01 ) effSig = 1.0;
        else cout << "Error: the efficiency is bigger than 100%" << endl;
      }
      hsigeff->SetBinContent(j,effSig);;
      hbkgeff->SetBinContent(j,1.-nbkgj/nbkg);
      hsignif->SetBinContent(j,nsigj+nbkgj>0 ? nsigj/sqrt(nsigj+nbkgj) : 0);
    }
    
    
    // Draw the values of the efficiencies at the soft muon cut value
    Double_t cutVal = GetSoftMuCutValue(varname[i]);
    sign.Clear();
    if ( cutVal < 0. )
    {
      cutVal = TMath::Abs(cutVal) - 1E-6;
      sign += "<";
    }
    else sign += ">";
    
    Bool_t printEff(kFALSE);
    if ( cutVal < 100. ) printEff = kTRUE;
    
    TLatex *  text(0x0);
    if ( printEff )
    {
      Double_t sigEff = hsigeff->GetBinContent(hsigeff->FindBin(cutVal));
      Double_t bkgEff = hbkgeff->GetBinContent(hbkgeff->FindBin(cutVal));
      
      text = new TLatex(0.1451149,0.5835095,Form("#varepsilon_{sig} [%%] =  %2.2f ; Rej_{bkg} [%%] =  %2.2f (%s %s %2.1f)",sigEff*100.,bkgEff*100.,varname[i],sign.Data(),cutVal));
      text->SetNDC();
      text->SetTextFont(42);
      text->SetTextSize(0.0422833);
      text->SetLineWidth(2);
    }
    
    TLatex *textKin = new TLatex(0.1465517,0.653277,fKinRange.Data());
    textKin->SetNDC();
    textKin->SetTextFont(42);
    textKin->SetTextSize(0.0422833);
    textKin->SetLineWidth(2);
    
    hsigeff->GetYaxis()->SetRangeUser(0,1.05);
    hsigeff->GetXaxis()->SetLabelSize(0.05);
    hsigeff->GetXaxis()->SetTitleSize(0.05);
    hsigeff->GetXaxis()->SetTitleOffset(0.93);
    hsigeff->GetYaxis()->SetLabelSize(0.05);
    hsigeff->GetYaxis()->SetTitleSize(0.05);
    hsigeff->GetYaxis()->SetTitleOffset(0.9);
    
    hsigeff->Draw("");
    hbkgeff->Draw("same");
    c1->Update();
    
    // scale hsignif to the pad coordinates
    Float_t rightmax = 1.1*hsignif->GetMaximum();
    Float_t rightmin = 0.98*hsignif->GetMaximum();
    Float_t scale = gPad->GetUymax()/rightmax;
    hsignif->Scale(scale);
    hsignif->Draw("same");
    // draw an axis on the right side
    TGaxis *axis = new TGaxis(gPad->GetUxmax(),gPad->GetUymin(),
                              gPad->GetUxmax(), gPad->GetUymax(),rightmin,rightmax,510,"+L");
    axis->SetLineColor(kGreen);
    axis->SetLabelColor(kGreen);
    axis->Draw();
    
    TLegend *leg = new TLegend(0.454023,0.3002114,0.8721264,0.4904863);
    leg->SetBorderSize(0);
    leg->SetLineColor(1);
    leg->SetLineStyle(1);
    leg->SetLineWidth(1);
    leg->SetFillColor(0);
    leg->SetFillStyle(0);
    leg->AddEntry(hsigeff,"Signal efficiency","lp");
    leg->AddEntry(hbkgeff,"1 - Background efficiency","lp");
    leg->AddEntry(hsignif,"S / #sqrt{S+B}","lp");
    leg->Draw();
    
    if ( text ) text->Draw("same");
    textKin->Draw("same");
    
    c1->Update();
    c1->SaveAs(Form("SingleCutsEfficiencyWSoftMuCut_%s.pdf",varname[i]));
  }
}


//___________________________________________________________
void plotSigMuIDDistrWSoftMuCutButOne(TFile* fsig, TFile* fbkg, const char* varname [ ])
{
  // Plot the muon ID variable distributions after aplying all the muon cuts but the variable itself.
  
  TObjArray* aSig = static_cast<TObjArray*>(fsig->FindObjectAny("SigDistrSoftMuCutButOne"));
  TObjArray* aBSig = static_cast<TObjArray*>(fsig->FindObjectAny("BkgDistrSoftMuCutButOne"));
  
  TObjArray* aSBkg = static_cast<TObjArray*>(fbkg->FindObjectAny("SigDistrSoftMuCutButOne"));
  TObjArray* aBkg = static_cast<TObjArray*>(fbkg->FindObjectAny("BkgDistrSoftMuCutButOne"));
  
  TH1F *hsig(0x0);
  TH1F *hBsig(0x0);
  TH1F *hSbkg(0x0);
  TH1F *hbkg(0x0);
  
  Int_t nVar = 6;
  if ( !isHPincluded ) nVar = 5;
  
  for (int i=0; i<nVar; i++)
  {
    TCanvas *c1 = new TCanvas();
    hsig = (TH1F*) aSig->FindObject(Form("hSoftMuNoCut_%s_sig",varname[i]));
    hBsig = (TH1F*) aBSig->FindObject(Form("hSoftMuNoCut_%s_bkg",varname[i]));
    
    hSbkg = (TH1F*) aSBkg->FindObject(Form("hSoftMuNoCut_%s_sig",varname[i]));
    hbkg = (TH1F*) aBkg->FindObject(Form("hSoftMuNoCut_%s_bkg",varname[i]));
    
    
    if (!hsig || !hBsig || !hSbkg || !hbkg )
    {
      std::cout << "No MuID test histos found" << std::endl;
      return;
    }
    
    //======= Creation of histograms
    int nbins = hsig->GetNbinsX();
    double binmin = hsig->GetXaxis()->GetXmin();
    double binmax = hsig->GetXaxis()->GetXmax();
    
    TH1F *hsigSubs = new TH1F(Form("hsigSubsDistr_%s",varname[i]),Form("Comparison data-MC signal %s distributions (rest of soft #mu-ID cuts applied);%s;",varname[i],varname[i]),nbins,binmin,binmax);
    hsigSubs->SetLineColor(kRed);
    
    TH1F *hSbkgSubs = new TH1F(Form("hSbkgSubsDistr_%s",varname[i]),Form("%s distribution (rest of soft #mu-ID cuts applied);%s;",varname[i],varname[i]),nbins,binmin,binmax);
    hSbkgSubs->SetLineColor(kBlue);
    hSbkgSubs->SetMarkerStyle(20);
    hSbkgSubs->SetMarkerColor(kBlue);
    //=======
    
    
    //======= Background substraction
    for (Int_t j = 1 ; j < hsig->GetNbinsX() + 1 ; j++ )
    {
      Double_t sigSubs = hsig->GetBinContent(j) - bkgNorm*hBsig->GetBinContent(j);
      if ( sigSubs < 0. ) sigSubs = 0.;
      Double_t SbkgSubs =  hSbkg->GetBinContent(j) - bkgNorm*hbkg->GetBinContent(j);
      if ( SbkgSubs < 0. ) SbkgSubs = 0.;
        
      hsigSubs->SetBinContent(j,sigSubs);
      hsigSubs->SetBinError(j,sqrt(sigSubs));
      
      hSbkgSubs->SetBinContent(j,SbkgSubs);
      hSbkgSubs->SetBinError(j,sqrt(SbkgSubs));
    }
    //=======
    
    
    hsigSubs->GetYaxis()->SetRangeUser(0.,1.);
    hSbkgSubs->GetYaxis()->SetRangeUser(0.,1.);
    
    
    //======= Histogram normalization
    Double_t intSig = hsigSubs->Integral();
    Double_t intSbkg = hSbkgSubs->Integral();
    for (Int_t j = 1 ; j < hsig->GetNbinsX() + 1 ; j++ )
    {
      Double_t sig = hsigSubs->GetBinContent(j);
      hsigSubs->SetBinContent(j,sig/intSig);
      
      Double_t sigError(0.);
      if ( sig > 0. ) sigError = (sig/intSig)*TMath::Sqrt( TMath::Power(1./TMath::Sqrt(sig),2.) + TMath::Power(1./TMath::Sqrt(intSig),2.) );
      hsigSubs->SetBinError(j,sigError);
      
      Double_t bkg = hSbkgSubs->GetBinContent(j);
      hSbkgSubs->SetBinContent(j,bkg/intSbkg);
      
      Double_t bkgError(0.);
      if ( bkg > 0. ) bkgError = (bkg/intSbkg)*TMath::Sqrt( TMath::Power(1./TMath::Sqrt(bkg),2.) + TMath::Power(1./TMath::Sqrt(intSbkg),2.) );
      hSbkgSubs->SetBinError(j,bkgError);
    }
    //=======
    
    
    //======= Integral of the distributions from the cut value
    Double_t cutVal = GetSoftMuCutValue(varname[i]);
    Int_t lBin(1);
    Int_t uBin(1);
    
    if ( cutVal > 0. )
    {
      lBin = 1;
      uBin = hsigSubs->FindBin(cutVal-1E-6);
    }
    else
    {
      lBin = hsigSubs->FindBin(abs(cutVal-1E-6));
      uBin = hsigSubs->GetNbinsX();
    }
    
    Double_t sigInt = hsigSubs->Integral(lBin,uBin);
    Double_t bkgInt = hSbkgSubs->Integral(lBin,uBin);
    
    Double_t ratio = sigInt/bkgInt;
    if ( sigInt == bkgInt ) ratio = 1.;
    TLatex *   text = new TLatex(0.1707901,0.8203883,Form("%% Evts removed (MC) / %% Evts removed (Data) = %2.2f/%2.2f = %2.2f",sigInt*100,bkgInt*100,ratio));
    text->SetNDC();
    text->SetTextFont(42);
    text->SetLineWidth(2);
    text->SetTextSize(0.03794038);
    
    TLatex* textKin = new TLatex(0.1465517,0.653277,fKinRange.Data());
    textKin->SetNDC();
    textKin->SetTextFont(42);
    textKin->SetTextSize(0.0422833);
    textKin->SetLineWidth(2);
    //=======
    
    hsigSubs->GetXaxis()->SetLabelSize(0.05);
    hsigSubs->GetXaxis()->SetTitleSize(0.05);
    hsigSubs->GetXaxis()->SetTitleOffset(0.93);
    hsigSubs->GetYaxis()->SetLabelSize(0.05);
    hsigSubs->GetYaxis()->SetTitleSize(0.05);
    hsigSubs->GetYaxis()->SetTitleOffset(0.9);
    
    
    // Draw line at the cut value
    TLine* l = new TLine(abs(cutVal),0.,abs(cutVal),1.);
    l->SetLineStyle(2);
    l->SetLineWidth(4);
    
    // Draw the histos
    hsigSubs->Draw("hist");
    hSbkgSubs->Draw("same");
    c1->Update();
    
    text->Draw("same");
    textKin->Draw("same");
    l->Draw("same");
    
    TLegend *leg = new TLegend(0.454023,0.3002114,0.8721264,0.4904863);
    leg->SetBorderSize(0);
    leg->SetLineColor(1);
    leg->SetLineStyle(1);
    leg->SetLineWidth(1);
    leg->SetFillColor(0);
    leg->SetFillStyle(0);
    leg->AddEntry(hsigSubs,"MC signal (bkg substracted)","lp");
    leg->AddEntry(hSbkgSubs,"Data signal (bkg substracted)","lp");
    leg->Draw("same");
    
    c1->Update();
    c1->SaveAs(Form("MuIDDistrWSoftMuCutButOne_%s.pdf",varname[i]));
  }
}


//___________________________________________________________
void plotSigMuIDDistrWSoftMuCut(TFile* fsig, TFile* fbkg, const char* varname [ ])
{
  // Plot the muon ID variables distributions after aplying the full soft muon ID cut
  
  TObjArray* aSig = static_cast<TObjArray*>(fsig->FindObjectAny("SigDistrSoftMuCutButOne"));
  TObjArray* aBSig = static_cast<TObjArray*>(fsig->FindObjectAny("BkgDistrSoftMuCutButOne"));
  
  TObjArray* aSBkg = static_cast<TObjArray*>(fbkg->FindObjectAny("SigDistrSoftMuCutButOne"));
  TObjArray* aBkg = static_cast<TObjArray*>(fbkg->FindObjectAny("BkgDistrSoftMuCutButOne"));
  
  TH1F *hsig(0x0);
  TH1F *hBsig(0x0);
  TH1F *hSbkg(0x0);
  TH1F *hbkg(0x0);
  
  Int_t nVar = 13;
  if ( !isHPincluded ) nVar = 14;
  
  for (int i=0; i<nVar; i++) {
    TCanvas *c1 = new TCanvas();
    hsig = (TH1F*) aSig->FindObject(Form("hSoftMuFull_%s_sig",varname[i]));
    hBsig = (TH1F*) aBSig->FindObject(Form("hSoftMuFull_%s_bkg",varname[i]));
    
    hSbkg = (TH1F*) aSBkg->FindObject(Form("hSoftMuFull_%s_sig",varname[i]));
    hbkg = (TH1F*) aBkg->FindObject(Form("hSoftMuFull_%s_bkg",varname[i]));
    
    
    if (!hsig || !hBsig || !hSbkg || !hbkg )
    {
      std::cout << "No MuID test histos found" << std::endl;
      return;
    }
    
    //======= Creation of histograms
    int nbins = hsig->GetNbinsX();
    double binmin = hsig->GetXaxis()->GetXmin();
    double binmax = hsig->GetXaxis()->GetXmax();
    
    TH1F *hsigSubs = new TH1F(Form("hsigSubs_SoftMuCut_%s",varname[i]),Form("Comparison data-MC signal %s distributions (soft #mu-ID cuts applied);%s;",varname[i],varname[i]),nbins,binmin,binmax);
    hsigSubs->SetLineColor(kRed);
    
    TH1F *hSbkgSubs = new TH1F(Form("hSbkgSubs_SoftMuCut_%s",varname[i]),Form("hSbkgSubs_DistrSoftMuCut;%s;",varname[i]),nbins,binmin,binmax);
    hSbkgSubs->SetLineColor(kBlue);
    hSbkgSubs->SetMarkerStyle(20);
    hSbkgSubs->SetMarkerColor(kBlue);
    //=======
    
    
    //======= Background substraction
    for (Int_t j = 1 ; j < hsig->GetNbinsX() + 1 ; j++ )
    {
      Double_t sigSubs = hsig->GetBinContent(j) - bkgNorm*hBsig->GetBinContent(j);
      Double_t SbkgSubs =  hSbkg->GetBinContent(j) - bkgNorm*hbkg->GetBinContent(j);
      
      hsigSubs->SetBinContent(j,sigSubs);
      hsigSubs->SetBinError(j,sqrt(sigSubs));
      
      hSbkgSubs->SetBinContent(j,SbkgSubs);
      hSbkgSubs->SetBinError(j,sqrt(SbkgSubs));
    }
    //=======
    
    
    //======= Histogram normalization
    Double_t intSig = hsigSubs->Integral();
    Double_t intSbkg = hSbkgSubs->Integral();
    for (Int_t j = 1 ; j < hsig->GetNbinsX() + 1 ; j++ )
    {
      Double_t sig = hsigSubs->GetBinContent(j);
      hsigSubs->SetBinContent(j,sig/intSig);
      
      Double_t sigError(0.);
      if ( sig > 0. ) sigError = (sig/intSig)*TMath::Sqrt( TMath::Power(1./TMath::Sqrt(sig),2.) + TMath::Power(1./TMath::Sqrt(intSig),2.) );
      hsigSubs->SetBinError(j,sigError);
      
      Double_t bkg = hSbkgSubs->GetBinContent(j);
      hSbkgSubs->SetBinContent(j,bkg/intSbkg);
      
      Double_t bkgError(0.);
      if ( bkg > 0. ) bkgError = (bkg/intSbkg)*TMath::Sqrt( TMath::Power(1./TMath::Sqrt(bkg),2.) + TMath::Power(1./TMath::Sqrt(intSbkg),2.) );
      hSbkgSubs->SetBinError(j,bkgError);
    }
    
    Double_t max = TMath::Max(hsigSubs->GetMaximum(),hSbkgSubs->GetMaximum());
    max = max + max*0.1;
    hsigSubs->GetYaxis()->SetRangeUser(0.,max);
    hSbkgSubs->GetYaxis()->SetRangeUser(0.,max);
    //=======
    
    //=======
    TLatex* textKin = new TLatex(0.1465517,0.653277,fKinRange.Data());
    textKin->SetNDC();
    textKin->SetTextFont(42);
    textKin->SetTextSize(0.0422833);
    textKin->SetLineWidth(2);
    //=======
    
    hsigSubs->GetXaxis()->SetLabelSize(0.05);
    hsigSubs->GetXaxis()->SetTitleSize(0.05);
    hsigSubs->GetXaxis()->SetTitleOffset(0.93);
    hsigSubs->GetYaxis()->SetLabelSize(0.05);
    hsigSubs->GetYaxis()->SetTitleSize(0.05);
    hsigSubs->GetYaxis()->SetTitleOffset(0.9);
    
    
    hsigSubs->Draw("hist");
    hSbkgSubs->Draw("same");
    textKin->Draw("same");
    c1->Update();
    
    
    TLegend *leg = new TLegend(0.454023,0.3002114,0.8721264,0.4904863);
    leg->SetBorderSize(0);
    leg->SetLineColor(1);
    leg->SetLineStyle(1);
    leg->SetLineWidth(1);
    leg->SetFillColor(0);
    leg->SetFillStyle(0);
    leg->AddEntry(hsigSubs,"MC signal (bkg substracted)","lp");
    leg->AddEntry(hSbkgSubs,"Data signal (bkg substracted)","lp");
    leg->Draw("same");
    
    c1->Update();
    c1->SaveAs(Form("MuIDDistrWSoftMuCut_%s.pdf",varname[i]));
  }
}



//___________________________________________________________
void plotSigBkgCompMuIDDistr(TFile* fsig, TFile* fbkg, const char* varname [ ])
{
  // Plot the data-MC comparison of muon ID variable distributions
  
  TObjArray* aSig = static_cast<TObjArray*>(fsig->FindObjectAny("SigDistr"));
  TObjArray* aBSig = static_cast<TObjArray*>(fsig->FindObjectAny("BkgDistr"));
  
  TObjArray* aBkg = static_cast<TObjArray*>(fbkg->FindObjectAny("BkgDistr"));
  
  TH1F *hsig(0x0);
  TH1F *hBsig(0x0);
  TH1F *hbkg(0x0);
  
  for (int i=0; i<19; i++)
  {
    TCanvas *c1 = new TCanvas();
    hsig = (TH1F*) aSig->FindObject(Form("h%s_sig_distr",varname[i]));
    hBsig = (TH1F*) aBSig->FindObject(Form("h%s_bkg_distr",varname[i]));
    
    hbkg = (TH1F*) aBkg->FindObject(Form("h%s_bkg_distr",varname[i]));
    
    
    if (!hsig || !hBsig || !hbkg )
    {
      std::cout << "No MuID test histos found" << std::endl;
      return;
    }
    
    //======= Creation of histograms
    int nbins = hsig->GetNbinsX();
    double binmin = hsig->GetXaxis()->GetXmin();
    double binmax = hsig->GetXaxis()->GetXmax();
    
    TH1F *hsigSubs = new TH1F(Form("hsigSubs_%s",varname[i]),Form("Comparison sig-bkg %s distribution;%s;",varname[i],varname[i]),nbins,binmin,binmax);
    hsigSubs->SetLineColor(kRed);
    hsigSubs->SetFillColor(kRed);
    hsigSubs->SetFillStyle(3005);
    
//    TH1F *hbkg = new TH1F(Form("hbkg_%s",varname[i]),Form("hbkg_DistrSoftMuCutBut;%s;",varname[i]),nbins,binmin,binmax);
    hbkg->SetLineColor(kBlue);
    hbkg->SetFillColor(kBlue);
    hbkg->SetFillStyle(3004);
    //=======
    
    
    //======= Background substraction
    for (Int_t j = 1 ; j < hsig->GetNbinsX() + 1 ; j++ )
    {
      Double_t sigSubs = hsig->GetBinContent(j) - bkgNorm*hBsig->GetBinContent(j);
      
      hsigSubs->SetBinContent(j,sigSubs);
      hsigSubs->SetBinError(j,sqrt(sigSubs));
    }
    //=======
    
    
    //======= Histogram normalization
    Double_t intSig = hsigSubs->Integral();
    Double_t intbkg = hbkg->Integral();
    for (Int_t j = 1 ; j < hsig->GetNbinsX() + 1 ; j++ )
    {
      Double_t sig = hsigSubs->GetBinContent(j);
      hsigSubs->SetBinContent(j,sig/intSig);
      
      Double_t sigError(0.);
      if ( sig > 0. ) sigError = (sig/intSig)*TMath::Sqrt( TMath::Power(1./TMath::Sqrt(sig),2.) + TMath::Power(1./TMath::Sqrt(intSig),2.) );
      hsigSubs->SetBinError(j,sigError);
      
      Double_t bkg = hbkg->GetBinContent(j);
      hbkg->SetBinContent(j,bkg/intbkg);
      
      Double_t bkgError(0.);
      if ( bkg > 0. ) bkgError = (bkg/intbkg)*TMath::Sqrt( TMath::Power(1./TMath::Sqrt(bkg),2.) + TMath::Power(1./TMath::Sqrt(intbkg),2.) );
      hbkg->SetBinError(j,bkgError);
    }
    
    Double_t max = TMath::Max(hsigSubs->GetMaximum(),hbkg->GetMaximum());
    max = max + max*0.1;
    hsigSubs->GetYaxis()->SetRangeUser(0.,max);
    hbkg->GetYaxis()->SetRangeUser(0.,max);
    //=======
  
    //=======
    TLatex* textKin = new TLatex(0.1465517,0.653277,fKinRange.Data());
    textKin->SetNDC();
    textKin->SetTextFont(42);
    textKin->SetTextSize(0.0422833);
    textKin->SetLineWidth(2);
    //=======
    
    hsigSubs->GetXaxis()->SetLabelSize(0.05);
    hsigSubs->GetXaxis()->SetTitleSize(0.05);
    hsigSubs->GetXaxis()->SetTitleOffset(0.93);
    hsigSubs->GetYaxis()->SetLabelSize(0.05);
    hsigSubs->GetYaxis()->SetTitleSize(0.05);
    hsigSubs->GetYaxis()->SetTitleOffset(0.9);
    
    // Draw the histos
    hsigSubs->Draw("hist");
    hbkg->Draw("samehist");
    textKin->Draw("same");
    c1->Update();
    
    TLegend *leg = new TLegend(0.454023,0.3002114,0.8721264,0.4904863);
    leg->SetBorderSize(0);
    leg->SetLineColor(1);
    leg->SetLineStyle(1);
    leg->SetLineWidth(1);
    leg->SetFillColor(0);
    leg->SetFillStyle(0);
    leg->AddEntry(hsigSubs,"MC signal (bkg substracted)","f");
    leg->AddEntry(hbkg,"Data background","f");
    leg->Draw("same");
    
    c1->Update();
    c1->SaveAs(Form("CompSigBkg_distr_%s.pdf",varname[i]));
  }
}


//___________________________________________________________
void plotMinvForMuIDCuts(TFile* fsig, TFile* fbkg, const char* varname [ ])
{
  TFile* file(0x0);
  TString stype("");
  if ( fsig && fbkg )
  {
    std::cout << "plotMinvForMuIDCuts only valid for plotting histos of one file at a time" << std::endl;
    return;
  }
  else if ( fsig )
  {
    file = fsig;
    stype += "MC";
  }
  else if ( fbkg )
  {
    file = fbkg;
    stype += "DATA";
  }
  else
  {
    std::cout << "Need a file to plot" << std::endl;
    return;
  }
  
  TObjArray* aSig = static_cast<TObjArray*>(file->FindObjectAny("MinvWCuts"));
  if ( !aSig )
  {
    std::cout << "No MinvWCuts folder found in file" << std::endl;
    return;
  }
  
  TObjArray* varMinv(0x0);
  TString hName("");
  
  Int_t nVar = 6;
  if ( !isHPincluded ) nVar = 5;
  
  for (int i=0; i<nVar; i++)
  {
    varMinv = static_cast<TObjArray*>(aSig->FindObject(Form("mInv_%s",varname[i])));
    if ( !varMinv )
    {
      std::cout << "No minv plots found" << std::endl;
      return;
    }
    
    TCanvas *c1 = new TCanvas();
    
    TLegend *leg = new TLegend(0.4770115,0.615222,0.8951149,0.8054968);
    leg->SetBorderSize(0);
    leg->SetLineColor(1);
    leg->SetLineStyle(1);
    leg->SetLineWidth(1);
    leg->SetFillColor(0);
    leg->SetFillStyle(0);
    
    //=======
    TLatex* textKin = new TLatex(0.1465517,0.653277,fKinRange.Data());
    textKin->SetNDC();
    textKin->SetTextFont(42);
    textKin->SetTextSize(0.0422833);
    textKin->SetLineWidth(2);
    //=======
  
    
    TIter histos(varMinv);
    TH1* h(0x0);
    Int_t color(1);
    while ( ( h = static_cast<TH1*>(histos.Next()) ) )
    {
      hName.Remove(0,hName.Sizeof());
      hName += h->GetName();
      hName.Remove(0,hName.Last('_')+1);
      hName.Remove(0,hName.First("t")+1);
      Double_t cutVal = hName.Atof();
      
      h->SetLineColor(color);
      h->SetMarkerColor(color++);
      
      h->GetXaxis()->SetLabelSize(0.05);
      h->GetXaxis()->SetTitleSize(0.05);
      h->GetXaxis()->SetTitleOffset(0.93);
      h->GetYaxis()->SetLabelSize(0.05);
      h->GetYaxis()->SetTitleSize(0.05);
      h->GetYaxis()->SetTitleOffset(0.9);
      
      h->SetTitle(Form("%s dimuon m_{inv} %s cuts comparison",stype.Data(),varname[i]));
      
      leg->AddEntry(h,Form("%s cut %2.4f",varname[i],cutVal),"pl");
      
      h->Draw("same");
      textKin->Draw("same");
      c1->Update();
    }
    
    leg->Draw("same");
    c1->Update();

    c1->SaveAs(Form("mInv_%s_CutsVariation_%s.pdf",stype.Data(),varname[i]));
  }
}


//___________________________________________________________
Bool_t getKinRange(TString sSig, TString sBkg)
{
  // Extract pt and y ranges from file names
  
  sSig.Remove(0,sSig.First('_')+1);
  sSig.Remove(0,sSig.First('_')+1);
  sSig.Remove(0,sSig.First('_')+1);
  
  sBkg.Remove(0,sBkg.First('_')+1);
  sBkg.Remove(0,sBkg.First('_')+1);
  sBkg.Remove(0,sBkg.First('_')+1);
  
  if ( sSig.CompareTo(sBkg.Data()) )
  {
    std::cout << "Signal and background are not computed in the same kinematical range or the soft muon cuts are different" << std::endl;
    return kFALSE;
  }
  
  sSig.Remove(0,sSig.First('_')+1);
  
  TString sSigY(sSig.Data());
  sSig.Remove(sSig.First('Y')-1,sSig.Sizeof());
  sSig.Remove(0,2); // Contains only the pt range
  TObjArray* aPt = static_cast<TObjArray*>(sSig.Tokenize("_"));
  
  TString sptMin = static_cast<TObjString*>(aPt->At(0))->GetString();
  Double_t ptMin = sptMin.Atof();
  
  TString sptMax = static_cast<TObjString*>(aPt->At(1))->GetString();
  Double_t ptMax = sptMax.Atof();
  
  sSigY.Remove(0,sSigY.First('Y'));
  sSigY.Remove(0,1); // Contains only the y range
  TObjArray* aY = static_cast<TObjArray*>(sSigY.Tokenize("_"));
  
  TString sYMin = static_cast<TObjString*>(aY->At(0))->GetString();
  Double_t YMin = sYMin.Atof();
  
  TString sYMax = static_cast<TObjString*>(aY->At(1))->GetString();
  Double_t YMax = sYMax.Atof();
  
  fKinRange += Form("%1.1f < p_{T} < %1.1f ; %1.1f < |y| < %1.1f",ptMin,ptMax,YMin,YMax);
  return kTRUE;
}


//___________________________________________________________
void muIDplots(const char* filename_sig, const char* filename_bkg)
{
  gStyle->SetOptStat(0);
  
  TFile *fsig = new TFile(filename_sig,"READ");
  TFile *fbkg = new TFile(filename_bkg,"READ");
  
  TString sSig(filename_sig);
  TString sBkg(filename_bkg);
  
  if ( sSig.Contains("NoHPincl") && sSig.Contains("NoHPincl") ) isHPincluded = kFALSE;

  if ( !getKinRange(sSig,sBkg) ) return;
  
  //============= Performing plots
  const char* varname[19] = {"isGoodMuon", "highPurity", "TrkMuArb", "TMOneStaTight", "nPixValHits",
    "nMuValHits", "nTrkHits", "normChi2_inner", "normChi2_global", "nPixWMea",
    "nTrkWMea", "StationsMatched", "dxy", "dxyErr", "dz",
    "dzErr", "ptErr_inner", "ptErr_global", "VtxProb"};
  
  plotSingleCutEfficiency(fsig,fbkg,varname);
  plotSingleCutEfficiencyWSoftMuCut(fsig,fbkg,varname);
  plotSigBkgCompMuIDDistr(fsig,fbkg,varname);
  
  
  if ( isHPincluded )
  {
    const char* varnameNoSoftMu[13] = {"TrkMuArb", "TMOneStaTight", "nPixValHits",
      "nMuValHits", "nTrkHits", "normChi2_inner", "normChi2_global", "StationsMatched", "dxyErr",
      "dzErr", "ptErr_inner", "ptErr_global", "VtxProb"};
    plotSigMuIDDistrWSoftMuCut(fsig,fbkg,varnameNoSoftMu);
    
    
    const char* varnameMuID[6] = {"isGoodMuon", "highPurity", "nPixWMea",
      "nTrkWMea", "dxy", "dz"};
    
    plotSigMuIDDistrWSoftMuCutButOne(fsig,fbkg,varnameMuID);
    plotMinvForMuIDCuts(fsig,0x0,varnameMuID);
    plotMinvForMuIDCuts(0x0,fbkg,varnameMuID);
  }
  else
  {
    const char* varnameNoSoftMu[14] = {"highPurity","TrkMuArb", "TMOneStaTight", "nPixValHits",
      "nMuValHits", "nTrkHits", "normChi2_inner", "normChi2_global", "StationsMatched", "dxyErr",
      "dzErr", "ptErr_inner", "ptErr_global", "VtxProb"};
    plotSigMuIDDistrWSoftMuCut(fsig,fbkg,varnameNoSoftMu);
    
    
    const char* varnameMuID[5] = {"isGoodMuon", "nPixWMea",
      "nTrkWMea", "dxy", "dz"};
    
    plotSigMuIDDistrWSoftMuCutButOne(fsig,fbkg,varnameMuID);
    plotMinvForMuIDCuts(fsig,0x0,varnameMuID);
    plotMinvForMuIDCuts(0x0,fbkg,varnameMuID);
  }

}



