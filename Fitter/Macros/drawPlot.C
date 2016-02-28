#ifndef drawPlot_C
#define drawPlot_C

#include "Utilities/initClasses.h"

void drawPlot(RooWorkspace& myws, struct InputOpt opt, struct KinCuts cut, bool isPbPb, bool zoomPsi, bool incSS, bool getMeanPT, float rangeY = 100, int nbins = 54) {

  RooPlot*   frame     = myws.var("invariantMass")->frame(Bins(nbins), Range(cut.dMuon.M.Min, cut.dMuon.M.Max)); RooPlot* frame2 = NULL;
  RooPlot*   framezoom = myws.var("invariantMass")->frame(Bins(19), Range(3.42,3.95));

  myws.data(Form("dataOS_%s", (isPbPb?"PbPb":"PP")))->plotOn(frame, Name("dataOS_FIT"), DataError(RooAbsData::SumW2), XErrorSize(0), 
							       MarkerColor(kBlack), LineColor(kBlack), MarkerSize(1.2));
  myws.pdf(Form("pdf_%s", (isPbPb?"PbPb":"PP")))->plotOn(frame,Name("theBKG"),Components(*myws.pdf(Form("bkgPDF_%s", (isPbPb?"PbPb":"PP")))),
							 Normalization(myws.data(Form("dataOS_%s", (isPbPb?"PbPb":"PP")))->sumEntries(),RooAbsReal::NumEvent), FillStyle(1001), FillColor(kAzure-9), 
							 VLines(), DrawOption("LCF"), LineColor(kBlue), LineStyle(kDashed),Precision(1e-4));
  if (incSS) { myws.data(Form("dataSS_%s", (isPbPb?"PbPb":"PP")))->plotOn(frame, Name("dataSS_FIT"), MarkerColor(kRed), LineColor(kRed), MarkerSize(1.2)); }
  myws.data(Form("dataOS_%s", (isPbPb?"PbPb":"PP")))->plotOn(frame, Name("dataOS_FIT"),DataError(RooAbsData::SumW2),XErrorSize(0), 
							     MarkerColor(kBlack), LineColor(kBlack), MarkerSize(1.2));
  myws.pdf(Form("pdf_%s", (isPbPb?"PbPb":"PP")))->plotOn(frame,Name("thePdf"),Normalization(myws.data(Form("dataOS_%s", (isPbPb?"PbPb":"PP")))->sumEntries(),RooAbsReal::NumEvent), 
							 LineColor(kBlack), LineStyle(1),Precision(1e-4));

  RooHist *hpull = frame->pullHist(0,0,true);
  hpull -> SetName("hpull");
  frame2 = myws.var("invariantMass")->frame(Title("Pull Distribution"),Bins(nbins),Range(cut.dMuon.M.Min,cut.dMuon.M.Max));
  frame2 -> addPlotable(hpull,"PX"); 

  if (zoomPsi) { 
    myws.data(Form("dataOS_%s", (isPbPb?"PbPb":"PP")))->plotOn(framezoom, Name("dataOS_FIT"), DataError(RooAbsData::SumW2), XErrorSize(0), 
									   MarkerColor(kBlack), LineColor(kBlack), MarkerSize(1.2)); 
    myws.pdf(Form("pdf_%s", (isPbPb?"PbPb":"PP")))->plotOn(framezoom, Name("theBKG"),Components(*myws.pdf(Form("bkgPDF_%s", (isPbPb?"PbPb":"PP")))),
							 Normalization(myws.data(Form("dataOS_%s", (isPbPb?"PbPb":"PP")))->sumEntries(),RooAbsReal::NumEvent), FillStyle(1001), FillColor(kAzure-9), 
							 VLines(), DrawOption("LCF"), LineColor(kBlue), LineStyle(kDashed),Precision(1e-4));
    myws.data(Form("dataOS_%s", (isPbPb?"PbPb":"PP")))->plotOn(framezoom, Name("dataOS_FIT"), DataError(RooAbsData::SumW2), XErrorSize(0), 
							     MarkerColor(kBlack), LineColor(kBlack), MarkerSize(1.2));
    myws.pdf(Form("pdf_%s", (isPbPb?"PbPb":"PP")))->plotOn(framezoom,Name("thePdf"),Normalization(myws.data(Form("dataOS_%s", (isPbPb?"PbPb":"PP")))->sumEntries(),RooAbsReal::NumEvent), 
							 LineColor(kBlack), LineStyle(1),Precision(1e-4));
  }			

  setTDRStyle();
  
  TCanvas* cFig = new TCanvas(Form("cFig_%s", (isPbPb?"PbPb":"PP")), "cFig",800,800);
  TPad *pad1 = new TPad(Form("pad1_%s", (isPbPb?"PbPb":"PP")),"",0,0.23,1,1);
  TPad *pad2 = new TPad(Form("pad2_%s", (isPbPb?"PbPb":"PP")),"",0,0,1,.228);
  TLine * pline = new TLine(cut.dMuon.M.Min,0.0,cut.dMuon.M.Max,0.0);
  
  TPad *pad4 = new TPad("pad4","This is pad4",0.55,0.46,0.97,0.87);
  pad4->SetFillStyle(0);
  pad4->SetLeftMargin(0.28);
  pad4->SetRightMargin(0.10);
  pad4->SetBottomMargin(0.21);
  pad4->SetTopMargin(0.072);

  float txtSize = opt.oniaMode==1 ? 0.032 : 0.028;
  if (opt.doFit) {
    float dx = opt.oniaMode==1 ? 0.63 : 0.61;
    if (!zoomPsi) {
      myws.pdf(Form("pdf_%s", (isPbPb?"PbPb":"PP")))->paramOn(frame,Layout(0.71,0.9285,0.7533),ShowConstants(kTRUE));
      frame->getAttText()->SetTextSize(0.022);
    }
    frame->SetTitle("");
    frame->GetXaxis()->SetTitle("");
    frame->GetXaxis()->CenterTitle(kTRUE);
    frame->GetXaxis()->SetTitleSize(0.045);
    frame->GetXaxis()->SetTitleFont(42);
    frame->GetXaxis()->SetTitleOffset(3);
    frame->GetXaxis()->SetLabelOffset(3);
    frame->GetYaxis()->SetLabelSize(0.04);
    frame->GetYaxis()->SetTitleSize(0.04);
    frame->GetYaxis()->SetTitleOffset(1.7);
    frame->GetYaxis()->SetTitleFont(42);
    frame->GetYaxis()->SetRangeUser(1, rangeY);
     
    cFig->cd();
    pad2->SetTopMargin(0.02);
    pad2->SetBottomMargin(0.4);
    pad2->SetFillStyle(4000); 
    pad2->SetFrameFillStyle(4000); 
    pad1->SetBottomMargin(0.015); 
    //plot fit
    pad1->Draw();
    pad1->cd(); 
    frame->Draw();
  }
  else{

    frame->SetTitle("");
    frame->GetXaxis()->SetTitle("#mu#mu mass (GeV/c^{2})");
    frame->GetXaxis()->CenterTitle(kTRUE);
    frame->GetXaxis()->SetTitleSize(0.045);
    frame->GetYaxis()->SetLabelSize(0.04);
    frame->GetXaxis()->SetTitleFont(42);

    frame->GetYaxis()->SetTitleSize(0.04);
    frame->GetYaxis()->SetTitleOffset(1.7);
    frame->GetYaxis()->SetTitleFont(42);

    cFig->cd();
 
    gPad->Update();
    frame->Draw();
  }

  TLatex *t = new TLatex(); t->SetNDC(); t->SetTextSize(txtSize);
  float dy = 0; float deltaY = 0.001;
  if (opt.oniaMode==1) { deltaY = 0.08; } 
  
  t->SetTextSize(0.03);
  t->DrawLatex(0.21, 0.86-dy, "2015 Soft Muon ID"); dy+=0.045;
  if (opt.oniaMode==1){
    if (isPbPb) {
      t->DrawLatex(0.21, 0.86-dy, "HLT_HIL1DoubleMu0_v1"); dy+=0.045;
    } else {
      t->DrawLatex(0.21, 0.86-dy, "HLT_HIL1DoubleMu0_v1"); dy+=0.045;
    } 
    if (isPbPb) {t->DrawLatex(0.21, 0.86-dy, Form("Cent. %d-%d%%", (int)(cut.Centrality.Start/2), (int)(cut.Centrality.End/2))); dy+=0.045;}
    t->DrawLatex(0.21, 0.86-dy, Form("%.1f < p_{T}^{#mu#mu} < %.1f GeV/c",cut.dMuon.Pt.Min,cut.dMuon.Pt.Max)); dy+=0.045;
    t->DrawLatex(0.21, 0.86-dy, Form("%.1f < |y^{#mu#mu}| < %.1f",cut.dMuon.AbsRap.Min,cut.dMuon.AbsRap.Max)); dy+=1.5*0.045;
    if (getMeanPT){
      t->DrawLatex(0.19, 0.86-dy, Form("<pt_{J/#psi}> = %.2f#pm%.2f GeV/c", myws.var(Form("<pt_{J/#psi}>^{%s}", (isPbPb?"PbPb":"PP")))->getValV(), myws.var(Form("<pt_{J/#psi}>^{%s}", (isPbPb?"PbPb":"PP")))->getError())); dy+=0.045;
      if (opt.inExcStat) {
	t->DrawLatex(0.19, 0.86-dy, Form("<pt_{#psi(2S)}> = %.2f#pm%.2f GeV/c", myws.var(Form("<pt_{#psi(2S)}>^{%s}", (isPbPb?"PbPb":"PP")))->getValV(), myws.var(Form("<pt_{#psi(2S)}>^{%s}", (isPbPb?"PbPb":"PP")))->getError())); dy+=0.045;
      }
      t->DrawLatex(0.19, 0.86-dy, Form("<pt_{bkg}> = %.2f#pm%.2f GeV/c", myws.var(Form("<pt_{bkg}>^{%s}", (isPbPb?"PbPb":"PP")))->getValV(), myws.var(Form("<pt_{bkg}>^{%s}", (isPbPb?"PbPb":"PP")))->getError())); dy+=0.045;
    }
  } 

  TLegend* leg = new TLegend(0.5175, 0.7802, 0.7180, 0.8809); leg->SetTextSize(0.03);
  leg->AddEntry(frame->findObject("dataOS_FIT"), (incSS?"Opposite Charge":"Data"),"pe");
  if (incSS) { leg->AddEntry(frame->findObject("dataSS_FIT"),"Same Charge","pe"); }
  leg->AddEntry(frame->findObject("thePdf"),"Total fit","l");
  leg->AddEntry(frame->findObject("theBKG"),"Background","fl");
  leg->Draw("same");

  //Drawing the title
  TString label;
  if (isPbPb) {
    if (opt.PbPb.RunNb.Start==opt.PbPb.RunNb.End){
      label = Form("PbPb Run %d", opt.PbPb.RunNb.Start);
    } else {
      label = Form("%s [%s %d-%d]", "PbPb", "HIOniaL1DoubleMu0", opt.PbPb.RunNb.Start, opt.PbPb.RunNb.End);
    }
  } else {
    if (opt.pp.RunNb.Start==opt.pp.RunNb.End){
      label = Form("PP Run %d", opt.pp.RunNb.Start);
    } else {
      label = Form("%s [%s %d-%d]", "PP", "DoubleMu0", opt.pp.RunNb.Start, opt.pp.RunNb.End);
    }
  }
  if(opt.doFit){ 
    CMS_lumi(pad1, isPbPb ? 105 : 104, 33, label);
        
    pad1->Update();
    cFig->cd(); 
     
    if (zoomPsi) {
      framezoom->SetName("zoom_frame_PbPb");
      framezoom->GetYaxis()->SetTitle(frame->GetYaxis()->GetTitle());
      framezoom->GetXaxis()->SetTitle(" ");
      framezoom->GetYaxis()->SetTitleOffset(1.3);
      framezoom->GetXaxis()->SetLabelOffset(0.012);
      framezoom->GetYaxis()->SetLabelSize(0.06);
      framezoom->GetXaxis()->SetLabelSize(0.06);
      framezoom->GetYaxis()->SetTitleSize(0.072);
      framezoom->GetXaxis()->SetTitleSize(0.072);
       
      pad4->Draw();
      pad4->cd();

      framezoom->Draw();
      pad4->Update();

      cFig->cd(); 
    }

    //---plot pull
    pad2->Draw();
    pad2->cd();
    
    frame2->GetYaxis()->CenterTitle(kTRUE);
    frame2->GetYaxis()->SetTitleOffset(0.4);
    frame2->GetYaxis()->SetTitleSize(0.1);
    frame2->GetYaxis()->SetLabelSize(0.1);
    frame2->GetYaxis()->SetTitle("Pull");
    frame2->GetXaxis()->CenterTitle(kTRUE);
    frame2->GetXaxis()->SetTitleOffset(1);
    frame2->GetXaxis()->SetTitleSize(0.12);
    frame2->GetXaxis()->SetLabelSize(0.1);
    frame2->GetXaxis()->SetTitle("m_{#mu^{+}#mu^{-}} (GeV/c^{2})");

    frame2->Draw(); 
    pline->Draw("same");
    pad2->Update();
  }
  else{
    CMS_lumi(cFig, isPbPb ? 105 : 104, 33, label);
    cFig->Update();
  }
 
  gSystem->mkdir("./Plots/root/", kTRUE); 
  cFig->SaveAs(Form("./Plots/root/%s_%sPrompt_pt%.0f%.0f_rap%.0f%.0f_cent%d%d_%d_%d.root", (opt.oniaMode==1?"Psi2SJpsi":"Upsilon"), (isPbPb?"PbPb":"PP"), (cut.dMuon.Pt.Min*10.0), (cut.dMuon.Pt.Max*10.0), (cut.dMuon.AbsRap.Min*10.0), (cut.dMuon.AbsRap.Max*10.0), cut.Centrality.Start, cut.Centrality.End, opt.PbPb.RunNb.Start, opt.PbPb.RunNb.End));
  gSystem->mkdir("./Plots/png/", kTRUE);
  cFig->SaveAs(Form("./Plots/png/%s_%sPrompt_pt%.0f%.0f_rap%.0f%.0f_cent%d%d_%d_%d.png", (opt.oniaMode==1?"Psi2SJpsi":"Upsilon"), (isPbPb?"PbPb":"PP"), (cut.dMuon.Pt.Min*10.0), (cut.dMuon.Pt.Max*10.0), (cut.dMuon.AbsRap.Min*10.0), (cut.dMuon.AbsRap.Max*10.0), cut.Centrality.Start, cut.Centrality.End, opt.PbPb.RunNb.Start, opt.PbPb.RunNb.End));
  gSystem->mkdir("./Plots/pdf/", kTRUE);
  //cFig->SaveAs(Form("./Plots/pdf/%s_%sPrompt_pt%.0f%.0f_rap%.0f%.0f_cent%d%d_%d_%d.pdf", (opt.oniaMode==1?"Psi2SJpsi":"Upsilon"), (isPbPb?"PbPb":"PP"), (cut.dMuon.Pt.Min*10.0), (cut.dMuon.Pt.Max*10.0), (cut.dMuon.AbsRap.Min*10.0), (cut.dMuon.AbsRap.Max*10.0), cut.Centrality.Start, cut.Centrality.End, opt.PbPb.RunNb.Start, opt.PbPb.RunNb.End));

}

#endif // #ifndef drawPlot_C
