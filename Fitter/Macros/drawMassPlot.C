#ifndef drawMassPlot_C
#define drawMassPlot_C

#include "Utilities/initClasses.h"

void printParameters(RooWorkspace myws, TPad* Pad, bool isPbPb);
void printChi2(RooWorkspace& myws, TPad* Pad, RooHist* hpull, string varLabel, string dataLabel, string pdfLabel); 

void drawMassPlot(RooWorkspace& myws, struct InputOpt opt, struct KinCuts cut, bool isPbPb, bool zoomPsi, bool setLogScale, 
		  bool incSS, bool getMeanPT, float rangeY = 100, int nbins = 54, bool isData=true, string MCTYPE="", const char* addLabel="") 
{

  if(zoomPsi) { setLogScale=false; }
  RooPlot*   frame     = myws.var("invMass")->frame(Bins(nbins), Range(cut.dMuon.M.Min, cut.dMuon.M.Max)); RooPlot* frame2 = NULL;
  RooPlot*   framezoom = myws.var("invMass")->frame(Bins(19), Range(Mass.Psi2S-0.265, Mass.Psi2S+0.265));

  myws.data(Form("dOS_%s_%s", (isData?"DATA":Form("MC%s", MCTYPE.c_str())), (isPbPb?"PbPb":"PP")))->plotOn(frame, Name("dOS"), DataError(RooAbsData::SumW2), XErrorSize(0), 
													   MarkerColor(kBlack), LineColor(kBlack), MarkerSize(1.2));
  myws.pdf(Form("pdfMASS_Tot_%s", (isPbPb?"PbPb":"PP")))->plotOn(frame,Name("BKG"),Components(*myws.pdf(Form("pdfMASS_Bkg_%s", (isPbPb?"PbPb":"PP")))),
								 Normalization(myws.data(Form("dOS_%s_%s", (isData?"DATA":Form("MC%s", MCTYPE.c_str())), (isPbPb?"PbPb":"PP")))->sumEntries(),
									       RooAbsReal::NumEvent), Range(cut.dMuon.M.Min, cut.dMuon.M.Max), 
								 FillStyle(1001), FillColor(kAzure-9), VLines(), DrawOption("LCF"), LineColor(kBlue), LineStyle(kDashed)
								 );
  if (incSS) { myws.data(Form("dSS_%s_%s", (isData?"DATA":Form("MC%s", MCTYPE.c_str())), (isPbPb?"PbPb":"PP")))->plotOn(frame, Name("dSS"), MarkerColor(kRed), LineColor(kRed), MarkerSize(1.2)); }
  myws.data(Form("dOS_%s_%s", (isData?"DATA":Form("MC%s", MCTYPE.c_str())), (isPbPb?"PbPb":"PP")))->plotOn(frame, Name("dOS"), DataError(RooAbsData::SumW2), XErrorSize(0), 
													   MarkerColor(kBlack), LineColor(kBlack), MarkerSize(1.2));
  myws.pdf(Form("pdfMASS_Tot_%s", (isPbPb?"PbPb":"PP")))->plotOn(frame,Name("PDF"), 
								 Normalization(myws.data(Form("dOS_%s_%s", (isData?"DATA":Form("MC%s", MCTYPE.c_str())), (isPbPb?"PbPb":"PP")))->sumEntries(),
									       RooAbsReal::NumEvent), 
								 LineColor(kBlack), LineStyle(1), Precision(1e-4), Range(cut.dMuon.M.Min, cut.dMuon.M.Max));


  RooHist *hpull = frame->pullHist(0, 0, true);
  hpull->SetName("hpull");
  frame2 = myws.var("invMass")->frame(Title("Pull Distribution"), Bins(nbins), Range(cut.dMuon.M.Min, cut.dMuon.M.Max));
  frame2->addPlotable(hpull, "PX"); 

  if (zoomPsi) { 
    myws.data(Form("dOS_%s_%s", (isData?"DATA":Form("MC%s", MCTYPE.c_str())), (isPbPb?"PbPb":"PP")))->plotOn(framezoom, Name("dOS"), DataError(RooAbsData::SumW2), XErrorSize(0), 
													     MarkerColor(kBlack), LineColor(kBlack), MarkerSize(1.2)); 
    myws.pdf(Form("pdfMASS_Tot_%s", (isPbPb?"PbPb":"PP")))->plotOn(framezoom, Name("BKG"),Components(*myws.pdf(Form("pdfMASS_Bkg_%s", (isPbPb?"PbPb":"PP")))),
								   Normalization(myws.data(Form("dOS_%s_%s", (isData?"DATA":Form("MC%s", MCTYPE.c_str())), (isPbPb?"PbPb":"PP")))->sumEntries(),
										 RooAbsReal::NumEvent), 
								   FillStyle(1001), FillColor(kAzure-9), VLines(), DrawOption("LCF"), LineColor(kBlue), LineStyle(kDashed),Precision(1e-4));
    myws.data(Form("dOS_%s_%s", (isData?"DATA":Form("MC%s", MCTYPE.c_str())), (isPbPb?"PbPb":"PP")))->plotOn(framezoom, Name("dOS"), DataError(RooAbsData::SumW2), XErrorSize(0), 
													     MarkerColor(kBlack), LineColor(kBlack), MarkerSize(1.2));
    myws.pdf(Form("pdfMASS_Tot_%s}", (isPbPb?"PbPb":"PP")))->plotOn(framezoom,Name("PDF"),
								    Normalization(myws.data(Form("dOS_%s_%s", (isData?"DATA":Form("MC%s", MCTYPE.c_str())), (isPbPb?"PbPb":"PP")))->sumEntries(),
										  RooAbsReal::NumEvent), 
								    LineColor(kBlack), LineStyle(1),Precision(1e-4));
  }			
  
  setTDRStyle();
  
  TCanvas* cFig = new TCanvas(Form("cMassFig_%s", (isPbPb?"PbPb":"PP")), "cMassFig",800,800);
  TPad *pad1 = new TPad(Form("pad1_%s", (isPbPb?"PbPb":"PP")),"",0,0.23,1,1);
  TPad *pad2 = new TPad(Form("pad2_%s", (isPbPb?"PbPb":"PP")),"",0,0,1,.228);
  TLine * pline = new TLine(cut.dMuon.M.Min, 0.0, cut.dMuon.M.Max, 0.0);
  
  TPad *pad4 = new TPad("pad4","This is pad4",0.55,0.46,0.97,0.87);
  pad4->SetFillStyle(0);
  pad4->SetLeftMargin(0.28);
  pad4->SetRightMargin(0.10);
  pad4->SetBottomMargin(0.21);
  pad4->SetTopMargin(0.072);

  float txtSize = opt.oniaMode==1 ? 0.032 : 0.028;
  float dx = opt.oniaMode==1 ? 0.63 : 0.61;
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

  // Find maximum and minimum points of Plot to rescale Y axis
  Double_t binW  = (cut.dMuon.M.Max-cut.dMuon.M.Min)/nbins;
  Double_t mNMax = myws.pdf(Form("pdfMASS_Tot_%s", (isPbPb?"PbPb":"PP")))->asTF(*myws.var("invMass"))->GetMaximumX(cut.dMuon.M.Min, cut.dMuon.M.Max); 
  Double_t mNMin = myws.pdf(Form("pdfMASS_Tot_%s", (isPbPb?"PbPb":"PP")))->asTF(*myws.var("invMass"))->GetMinimumX(cut.dMuon.M.Min, cut.dMuon.M.Max);
  myws.var("invMass")->setRange("YMaxBin",mNMax-binW*0.5,mNMax+binW*0.5);
  myws.var("invMass")->setRange("YMinBin",mNMin-binW*0.5,mNMin+binW*0.5);
  Double_t YMax = (myws.pdf(Form("pdfMASS_Tot_%s", (isPbPb?"PbPb":"PP")))->createIntegral(RooArgSet(*myws.var("invMass")), NormSet(RooArgSet(*myws.var("invMass"))), Range("YMaxBin"))->getVal()*
		   myws.data(Form("dOS_%s_%s", (isData?"DATA":Form("MC%s", MCTYPE.c_str())), (isPbPb?"PbPb":"PP")))->sumEntries());
  Double_t YMin = (myws.pdf(Form("pdfMASS_Tot_%s", (isPbPb?"PbPb":"PP")))->createIntegral(RooArgSet(*myws.var("invMass")), NormSet(RooArgSet(*myws.var("invMass"))), Range("YMinBin"))->getVal()*
		   myws.data(Form("dOS_%s_%s", (isData?"DATA":Form("MC%s", MCTYPE.c_str())), (isPbPb?"PbPb":"PP")))->sumEntries());
  if(setLogScale){ YMin=max(YMin,1.0); frame->GetYaxis()->SetRangeUser(max(YMin/pow((YMax/YMin), 0.05),1.0), YMax*pow((YMax/YMin), 0.4)); } 
  else{ frame->GetYaxis()->SetRangeUser(max(YMin-(YMax-YMin)*0.05,0.0), YMax+(YMax-YMin)*0.4); }
 
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

  if (!zoomPsi) { printParameters(myws, pad1, isPbPb); }
  pad1->SetLogy(setLogScale);

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
      t->DrawLatex(0.19, 0.86-dy, Form("<pt_{J/#psi}> = %.2f#pm%.2f GeV/c", myws.var(Form("ptJpsi%s", (isPbPb?"PbPb":"PP")))->getValV(), myws.var(Form("ptJpsi%s", (isPbPb?"PbPb":"PP")))->getError())); dy+=0.045;
      if (opt.inExcStat) {
	t->DrawLatex(0.19, 0.86-dy, Form("<pt_{#psi(2S)}> = %.2f#pm%.2f GeV/c", myws.var(Form("ptPsi2S%s", (isPbPb?"PbPb":"PP")))->getValV(), myws.var(Form("ptPsi2S%s", (isPbPb?"PbPb":"PP")))->getError())); dy+=0.045;
      }
      t->DrawLatex(0.19, 0.86-dy, Form("<pt_{bkg}> = %.2f#pm%.2f GeV/c", myws.var(Form("ptBkg%s", (isPbPb?"PbPb":"PP")))->getValV(), myws.var(Form("ptBkg%s", (isPbPb?"PbPb":"PP")))->getError())); dy+=0.045;
    }
  } 

  TLegend* leg = new TLegend(0.5175, 0.7802, 0.7180, 0.8809); leg->SetTextSize(0.03);
  leg->AddEntry(frame->findObject("dOS"), (incSS?"Opposite Charge":"Data"),"pe");
  if (incSS) { leg->AddEntry(frame->findObject("dSS"),"Same Charge","pe"); }
  leg->AddEntry(frame->findObject("PDF"),"Total fit","l");
  leg->AddEntry(frame->findObject("BKG"),"Background","fl");
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
   
  CMS_lumi(pad1, isPbPb ? 105 : 104, 33, label);
  gStyle->SetTitleFontSize(0.05);
  
  pad1->Update();
  cFig->cd(); 
     
  if (zoomPsi) {
    framezoom->SetName("zoom_frame_PbPb");
    framezoom->SetTitle("");
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

    cFig->cd(); 
  }

  //---plot pull
  pad2->Draw();
  pad2->cd();
    
  frame2->SetTitle("");
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
  frame2->GetYaxis()->SetRangeUser(-7.0, 7.0);

  frame2->Draw(); 
  
  // *** Print chi2/ndof 
  printChi2(myws, pad2, hpull, 
	    "invMass", 
	    Form("dOS_%s_%s", (isData?"DATA":Form("MC%s", MCTYPE.c_str())), (isPbPb?"PbPb":"PP")), 
	    Form("pdfMASS_Tot_%s", (isPbPb?"PbPb":"PP"))
	    );
  
  pline->Draw("same");
  pad2->Update();
  
  if(isData){
    gSystem->mkdir("./Plots/DATA/root/", kTRUE); 
    cFig->SaveAs(Form("./Plots/DATA/root/DATA_%s_%sPrompt_pt%.0f%.0f_rap%.0f%.0f_cent%d%d_%d_%d_%s.root", (opt.oniaMode==1?"Psi2SJpsi":"Upsilon"), (isPbPb?"PbPb":"PP"), (cut.dMuon.Pt.Min*10.0), (cut.dMuon.Pt.Max*10.0), (cut.dMuon.AbsRap.Min*10.0), (cut.dMuon.AbsRap.Max*10.0), cut.Centrality.Start, cut.Centrality.End, opt.PbPb.RunNb.Start, opt.PbPb.RunNb.End, addLabel));
    gSystem->mkdir("./Plots/DATA/png/", kTRUE);
    cFig->SaveAs(Form("./Plots/DATA/png/DATA_%s_%sPrompt_pt%.0f%.0f_rap%.0f%.0f_cent%d%d_%d_%d_%s.png", (opt.oniaMode==1?"Psi2SJpsi":"Upsilon"), (isPbPb?"PbPb":"PP"), (cut.dMuon.Pt.Min*10.0), (cut.dMuon.Pt.Max*10.0), (cut.dMuon.AbsRap.Min*10.0), (cut.dMuon.AbsRap.Max*10.0), cut.Centrality.Start, cut.Centrality.End, opt.PbPb.RunNb.Start, opt.PbPb.RunNb.End, addLabel));
    gSystem->mkdir("./Plots/DATA/pdf/", kTRUE);
    cFig->SaveAs(Form("./Plots/DATA/pdf/DATA_%s_%sPrompt_pt%.0f%.0f_rap%.0f%.0f_cent%d%d_%d_%d_%s.pdf", (opt.oniaMode==1?"Psi2SJpsi":"Upsilon"), (isPbPb?"PbPb":"PP"), (cut.dMuon.Pt.Min*10.0), (cut.dMuon.Pt.Max*10.0), (cut.dMuon.AbsRap.Min*10.0), (cut.dMuon.AbsRap.Max*10.0), cut.Centrality.Start, cut.Centrality.End, opt.PbPb.RunNb.Start, opt.PbPb.RunNb.End, addLabel));
    gSystem->mkdir("./Plots/DATA/pdf/", kTRUE);
    //cFig->SaveAs(Form("./Plots/DATA/pdf/DATA_%s_%sPrompt_pt%.0f%.0f_rap%.0f%.0f_cent%d%d_%d_%d.pdf", (opt.oniaMode==1?"Psi2SJpsi":"Upsilon"), (isPbPb?"PbPb":"PP"), (cut.dMuon.Pt.Min*10.0), (cut.dMuon.Pt.Max*10.0), (cut.dMuon.AbsRap.Min*10.0), (cut.dMuon.AbsRap.Max*10.0), cut.Centrality.Start, cut.Centrality.End, opt.PbPb.RunNb.Start, opt.PbPb.RunNb.End));
  } else {
    gSystem->mkdir("./Plots/MC/root/", kTRUE); 
    cFig->SaveAs(Form("./Plots/MC/root/MC%s_%s_%sPrompt_pt%.0f%.0f_rap%.0f%.0f_cent%d%d_%d_%d_%s.root", MCTYPE.c_str(), (opt.oniaMode==1?"Psi2SJpsi":"Upsilon"), (isPbPb?"PbPb":"PP"), (cut.dMuon.Pt.Min*10.0), (cut.dMuon.Pt.Max*10.0), (cut.dMuon.AbsRap.Min*10.0), (cut.dMuon.AbsRap.Max*10.0), cut.Centrality.Start, cut.Centrality.End, opt.PbPb.RunNb.Start, opt.PbPb.RunNb.End, addLabel));
    gSystem->mkdir("./Plots/MC/png/", kTRUE);
    cFig->SaveAs(Form("./Plots/MC/png/MC%s_%s_%sPrompt_pt%.0f%.0f_rap%.0f%.0f_cent%d%d_%d_%d_%s.png", MCTYPE.c_str(), (opt.oniaMode==1?"Psi2SJpsi":"Upsilon"), (isPbPb?"PbPb":"PP"), (cut.dMuon.Pt.Min*10.0), (cut.dMuon.Pt.Max*10.0), (cut.dMuon.AbsRap.Min*10.0), (cut.dMuon.AbsRap.Max*10.0), cut.Centrality.Start, cut.Centrality.End, opt.PbPb.RunNb.Start, opt.PbPb.RunNb.End, addLabel));
    gSystem->mkdir("./Plots/MC/pdf/", kTRUE);
    cFig->SaveAs(Form("./Plots/MC/pdf/MC%s_%s_%sPrompt_pt%.0f%.0f_rap%.0f%.0f_cent%d%d_%d_%d_%s.pdf", MCTYPE.c_str(), (opt.oniaMode==1?"Psi2SJpsi":"Upsilon"), (isPbPb?"PbPb":"PP"), (cut.dMuon.Pt.Min*10.0), (cut.dMuon.Pt.Max*10.0), (cut.dMuon.AbsRap.Min*10.0), (cut.dMuon.AbsRap.Max*10.0), cut.Centrality.Start, cut.Centrality.End, opt.PbPb.RunNb.Start, opt.PbPb.RunNb.End, addLabel));
    cFig->Close();
  }

}

#endif // #ifndef drawMassPlot_C


void printParameters(RooWorkspace myws, TPad* Pad, bool isPbPb)
{
  Pad->cd();
  TLatex *t = new TLatex(); t->SetNDC(); t->SetTextSize(0.026); float dy = 0.025; 
  RooArgSet* Parameters =  myws.pdf(Form("pdfMASS_Tot_%s", (isPbPb?"PbPb":"PP")))->getParameters(myws.data(Form("dOS_DATA_%s", (isPbPb?"PbPb":"PP"))));
  TIterator* parIt = Parameters->createIterator(); 
  for (RooRealVar* it = (RooRealVar*)parIt->Next(); it!=NULL; it = (RooRealVar*)parIt->Next() ) {
    stringstream ss(it->GetName()); string s1, s2, s3, label; 
    getline(ss, s1, '_'); getline(ss, s2, '_'); getline(ss, s3, '_');
    if(s1=="RFrac2Svs1S"){s1="R_{#psi(2S)/J/#psi}";} if(s2=="PbPbvsPP"){s2="PbPb/PP";}
    if(s1.find("sigma")!=std::string::npos || s1.find("lambda")!=std::string::npos || s1.find("alpha")!=std::string::npos){s1=Form("#%s",s1.c_str());}
    if(s2=="Jpsi"){s2="J/#psi";} else if(s2=="Psi2S"){s2="#psi(2S)";} else if(s2=="Bkg"){s2="bkg";}
    if(s3!=""){label=Form("%s_{%s}^{%s}", s1.c_str(), s2.c_str(), s3.c_str());} else {label=Form("%s^{%s}", s1.c_str(), s2.c_str());}
    if(s1=="invMass"){continue;} else if(s1=="MassRatio"){continue;}
    else if(s1=="N"){ t->DrawLatex(0.69, 0.75-dy, Form("%s = %.0f#pm%.0f ", label.c_str(), it->getValV(), it->getError())); dy+=0.045; }
    else if(s1.find("sigma")!=std::string::npos){ t->DrawLatex(0.69, 0.75-dy, Form("%s = %.2f#pm%.2f MeV/c^{2}", label.c_str(), it->getValV()*1000., it->getError()*1000.)); dy+=0.045; } 
    else if(s1.find("lambda")!=std::string::npos){ t->DrawLatex(0.69, 0.75-dy, Form("%s = %.4f#pm%.4f", label.c_str(), it->getValV(), it->getError())); dy+=0.045; }
    else if(s1.find("m")!=std::string::npos){ t->DrawLatex(0.69, 0.75-dy, Form("%s = %.5f#pm%.5f GeV/c^{2}", label.c_str(), it->getValV(), it->getError())); dy+=0.045; }
    else { t->DrawLatex(0.69, 0.75-dy, Form("%s = %.4f#pm%.4f", label.c_str(), it->getValV(), it->getError())); dy+=0.045; }
  }
}


void printChi2(RooWorkspace& myws, TPad* Pad, RooHist* hpull, string varLabel, string dataLabel, string pdfLabel) 
{
  double chi2=0; unsigned int ndof=0;
  Pad->cd();
  TLatex *t = new TLatex(); t->SetNDC(); t->SetTextSize(0.1); 
  unsigned int nbins = hpull->GetN();
  TH1 *hdatact = myws.data(dataLabel.c_str())->createHistogram("hdatact", *myws.var(varLabel.c_str()), Binning(nbins));
  unsigned int nFitPar = myws.pdf(pdfLabel.c_str())->getParameters(*myws.data(dataLabel.c_str()))->selectByAttrib("Constant",kFALSE)->getSize(); 
  double* ypulls = hpull->GetY();
  unsigned int nFullBins = 0;
  for (unsigned int i = 0; i < nbins; i++) {
    if (hdatact->GetBinContent(i+1) != 0) {
      chi2 += ypulls[i]*ypulls[i];
      nFullBins++;
    }
  }
  ndof = nFullBins - nFitPar;
  t->DrawLatex(0.7, 0.8, Form("#chi^{2}/ndof = %.0f / %d", chi2, ndof));
};
