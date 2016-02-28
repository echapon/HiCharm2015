#include "Utilities/initClasses.h"

void fitMeanPTJpsi2015(RooWorkspace& myws, struct InputOpt opt, struct KinCuts cut, bool isPbPb, int nbins)
{
  // SAVE CURRENT PARAMETER VALUES AND SET THEM CONSTANT 
  
  myws.saveSnapshot(Form("pdfFIT_%s", (isPbPb?"PbPb":"PP")), *myws.pdf(Form("pdf_%s", (isPbPb?"PbPb":"PP")))->getParameters(myws.data(Form("dataOS_%s", (isPbPb?"PbPb":"PP")))) , kTRUE); 
  myws.pdf(Form("pdf_%s", (isPbPb?"PbPb":"PP")))->getParameters(myws.data(Form("DHProf_%s", (isPbPb?"PbPb":"PP"))))->setAttribAll("Constant", kTRUE);

  // C r e a t e  Mean PT  m o d e l

  Double_t    WBin    = (((RooDataHist*)myws.data(Form("DHProf_%s", (isPbPb?"PbPb":"PP"))))->numEntries())/(cut.dMuon.M.Max - cut.dMuon.M.Min);
  RooRealVar* binW    = new RooRealVar("binWidth","binWidth", WBin, WBin, WBin);
  
  RooRealVar* ptPSI2S = NULL;	
  if(opt.inExcStat){ ptPSI2S = new RooRealVar(Form("<pt_{#psi(2S)}>^{%s}", (isPbPb?"PbPb":"PP")),"<pt_{#psi(2S)}>", 20.0, 0.0, 40.0); }
  RooRealVar* ptJPSI  = new RooRealVar(Form("<pt_{J/#psi}>^{%s}", (isPbPb?"PbPb":"PP")),"<pt_{J/#psi}>", 20.0, 0.0, 40.0);
  
  RooRealVar* ptBKG   = new RooRealVar(Form("<pt_{bkg}>^{%s}", (isPbPb?"PbPb":"PP")),"<pt_{bkg}>", 20.0, 0.0, 40.0);
  RooRealVar*  cBKG   = new RooRealVar(Form("cpt_{bkg}^{%s}", (isPbPb?"PbPb":"PP")),"cpt_{bkg}", 0.2, -1.0, 1.0);


  RooFormulaVar*  meanPTPSI2S = NULL;
  if(opt.inExcStat){ meanPTPSI2S = new RooFormulaVar(Form("meanPTPSI2S_%s", (isPbPb?"PbPb":"PP")), "meanPTPSI2S", "@0*@1", RooArgSet(*ptPSI2S, *binW)); }
  RooFormulaVar*  meanPTJPSI  = new RooFormulaVar(Form("meanPTJPSI_%s", (isPbPb?"PbPb":"PP")), "meanPTJPSI", "@0*@1", RooArgSet(*ptJPSI, *binW));
  RooFormulaVar*  meanPTBKG   = new RooFormulaVar(Form("meanPTBKG_%s", (isPbPb?"PbPb":"PP")), "meanPTBKG", "@0*@1*( 1.0 + @2*( @3 - 3.1 ) )", 
						  RooArgSet(*ptBKG, *binW, *cBKG, *myws.var("invariantMass")));
	
  RooFormulaVar*  totPDF = NULL;
  if(opt.inExcStat){
    if (opt.doSimulFit && isPbPb) {
      totPDF = new RooFormulaVar(Form("totPDF_%s", (isPbPb?"PbPb":"PP")),"totPDF", "(@0*@4 + @2*@1*@0*@5 + @3*@6)", 
				 RooArgList( *myws.var(Form("N_{J/#psi}^{%s}", (isPbPb?"PbPb":"PP"))),
					     *myws.var("R_{#frac{2S}{1S}}^{PP}"),
					     *myws.var("R_{#frac{2S}{1S}}^{PbPb/PP}"),
					     *myws.var(Form("N_{bkg}^{%s}", (isPbPb?"PbPb":"PP"))),
					     *myws.pdf(Form("sigJPSI_%s", (isPbPb?"PbPb":"PP"))),
					     *myws.pdf(Form("sigPSI2S_%s", (isPbPb?"PbPb":"PP"))),
					     *myws.pdf(Form("bkgPDF_%s", (isPbPb?"PbPb":"PP")))
					     ));
    } else {
      totPDF = new RooFormulaVar(Form("totPDF_%s", (isPbPb?"PbPb":"PP")),"totPDF", "(@0*@3 + @1*@0*@4 + @2*@5)", 
				 RooArgList( *myws.var(Form("N_{J/#psi}^{%s}", (isPbPb?"PbPb":"PP"))),
					     *myws.var(Form("R_{#frac{2S}{1S}}^{%s}", (isPbPb?"PbPb":"PP"))),
					     *myws.var(Form("N_{bkg}^{%s}", (isPbPb?"PbPb":"PP"))),
					     *myws.pdf(Form("sigJPSI_%s", (isPbPb?"PbPb":"PP"))),
					     *myws.pdf(Form("sigPSI2S_%s", (isPbPb?"PbPb":"PP"))),
					     *myws.pdf(Form("bkgPDF_%s", (isPbPb?"PbPb":"PP")))
					     ));	
    }
  } else {
    totPDF = new RooFormulaVar(Form("totPDF_%s", (isPbPb?"PbPb":"PP")),"totPDF", "(@0*@2 + @1*@3)", 
			       RooArgList( *myws.var(Form("N_{J/#psi}^{%s}", (isPbPb?"PbPb":"PP"))),
					   *myws.var(Form("N_{bkg}^{%s}", (isPbPb?"PbPb":"PP"))),
					   *myws.pdf(Form("sigJPSI_%s", (isPbPb?"PbPb":"PP"))),
					   *myws.pdf(Form("bkgPDF_%s", (isPbPb?"PbPb":"PP")))
					   ));
  }

  RooFormulaVar*  alphaPSI2S = NULL;
  if(opt.inExcStat){ 
    if (opt.doSimulFit && isPbPb) {
      alphaPSI2S = new RooFormulaVar(Form("alphaPSI2S_%s", (isPbPb?"PbPb":"PP")),"alphaPSI2S", "(@2*@1*@0*@3)/(@4)", 
				     RooArgList( *myws.var(Form("N_{J/#psi}^{%s}", (isPbPb?"PbPb":"PP"))),
						 *myws.var("R_{#frac{2S}{1S}}^{PP}"),
						 *myws.var("R_{#frac{2S}{1S}}^{PbPb/PP}"),
						 *myws.pdf(Form("sigPSI2S_%s", (isPbPb?"PbPb":"PP"))),
						 *totPDF
						 ));
    } else {
      alphaPSI2S = new RooFormulaVar(Form("alphaPSI2S_%s", (isPbPb?"PbPb":"PP")),"alphaPSI2S", "(@1*@0*@2)/(@3)", 
				     RooArgList( *myws.var(Form("N_{J/#psi}^{%s}", (isPbPb?"PbPb":"PP"))),
						 *myws.var(Form("R_{#frac{2S}{1S}}^{%s}", (isPbPb?"PbPb":"PP"))),
						 *myws.pdf(Form("sigPSI2S_%s", (isPbPb?"PbPb":"PP"))),
						 *totPDF
						 ));
    } 
  }
  RooFormulaVar*  alphaJPSI = new RooFormulaVar(Form("alphaJPSI_%s", (isPbPb?"PbPb":"PP")),"alphaJPSI", "(@0*@1)/(@2)", 
						RooArgList( *myws.var(Form("N_{J/#psi}^{%s}", (isPbPb?"PbPb":"PP"))),
							    *myws.pdf(Form("sigJPSI_%s", (isPbPb?"PbPb":"PP"))),
							    *totPDF
							    ));	
  RooFormulaVar*  alphaBKG = new RooFormulaVar(Form("alphaBKG_%s", (isPbPb?"PbPb":"PP")),"alphaBKG", "(@0*@1)/(@2)", 
					       RooArgList(*myws.var(Form("N_{bkg}^{%s}", (isPbPb?"PbPb":"PP"))),
							  *myws.pdf(Form("bkgPDF_%s", (isPbPb?"PbPb":"PP"))),
							  *totPDF
							  ));

  RooFormulaVar*  ptFUNPSI2S = NULL; 
  if(opt.inExcStat){ ptFUNPSI2S = new RooFormulaVar(Form("ptFUNPSI2S_%s", (isPbPb?"PbPb":"PP")),"ptFUNPSI2S", "@0*@1", RooArgList(*alphaPSI2S, *meanPTPSI2S)); }
  RooFormulaVar*  ptFUNJPSI = new RooFormulaVar(Form("ptFUNJPSI_%s", (isPbPb?"PbPb":"PP")),"ptFUNJPSI", "@0*@1", RooArgList(*alphaJPSI, *meanPTJPSI));
  RooFormulaVar*  ptFUNBKG  = new RooFormulaVar(Form("ptFUNBKG_%s", (isPbPb?"PbPb":"PP")),"ptFUNBKG", "@0*@1", RooArgList(*alphaBKG, *meanPTBKG));
  RooFormulaVar*  ptFUN     = NULL;
  if(opt.inExcStat){ 
    ptFUN = new RooFormulaVar(Form("ptFUN_%s", (isPbPb?"PbPb":"PP")), "ptFUN", "@0 + @1 + @2", RooArgList(*ptFUNJPSI, *ptFUNPSI2S, *ptFUNBKG));
  } else {
    ptFUN = new RooFormulaVar(Form("ptFUN_%s", (isPbPb?"PbPb":"PP")), "ptFUN", "@0 + @1", RooArgList(*ptFUNJPSI, *ptFUNBKG));
  }
  RooFitResult* fitPtProfile = ptFUN->chi2FitTo(*((RooDataHist*)myws.data(Form("DHProf_%s", (isPbPb?"PbPb":"PP")))), Save(), NumCPU(8), PrintEvalErrors(1));
  
  myws.import(*ptFUN);
  
  RooPlot*   framed     = myws.var("invariantMass")->frame(Bins(nbins), Range(cut.dMuon.M.Min, cut.dMuon.M.Max));
  ((RooDataHist*)myws.data(Form("DHProf_%s", (isPbPb?"PbPb":"PP"))))->plotOn(framed, Name(Form("DHProfFIT_%s", (isPbPb?"PbPb":"PP"))), DataError(RooAbsData::SumW2), XErrorSize(0), 
									     MarkerColor(kBlack), LineColor(kBlack), MarkerSize(1.2));
  ((RooFormulaVar*)myws.obj(Form("ptFUN_%s", (isPbPb?"PbPb":"PP"))))->plotOn(framed,Name(Form("ptFUNFIT_%s", (isPbPb?"PbPb":"PP"))), LineColor(kBlack), LineStyle(1), Normalization(1.0/WBin, RooAbsReal::NumEvent));
  
  setTDRStyle();
  TCanvas* cFig = new TCanvas(Form("cFigMeanPT_%s", (isPbPb?"PbPb":"PP")), "cFig",800,800);
  cFig->cd();
  framed->GetYaxis()->SetRangeUser(ptBKG->getValV()-2.0, ptBKG->getValV()+2.0);
  framed->Draw();
  //Drawing the title
  TString label;
  if (isPbPb) {
    label = Form("%s [%s]", "PbPb", "HIOniaL1DoubleMu0");
  } else {
    label = Form("%s [%s]", "PP", "DoubleMu0");
  }
  CMS_lumi(cFig, isPbPb ? 105 : 104, 33, label);
  cFig->Update();

  
  // RESTORE PREVIOUS PARAMETER VALUES 
  
  myws.loadSnapshot(Form("pdfFIT_%s", (isPbPb?"PbPb":"PP"))); 
  myws.Print();
 
}
