#include "Utilities/initClasses.h"

void buildModelJpsi2015(RooWorkspace& w, struct InputOpt opt,  int sigModel, int bkgModel, bool isPbPb)
{
  // C r e a t e   m o d e l  
  
  bool doRatio = true;
  int nt= 0; 
  if (isPbPb) { nt = 200000; }
  else        { nt = 300000; }

  // Signal Model
  RooRealVar *nSigJPSI = new RooRealVar(Form("N_{J/#psi}^{%s}", (isPbPb?"PbPb":"PP")),"N_{J/#psi}",0,nt);
  RooRealVar *nSigPSI2S = new RooRealVar(Form("N_{#psi(2S)}^{%s}", (isPbPb?"PbPb":"PP")),"N_{#psi(2S)}",0,nt);
  if(!opt.inExcStat){nSigPSI2S=NULL;}
  RooRealVar *mass = (RooRealVar*) w.var("invariantMass");

  RooRealVar *f2Svs1S   = NULL; RooRealVar *fDRat2Svs1S = NULL;
  if(doRatio && opt.inExcStat){
    if (opt.doSimulFit && isPbPb) {
      fDRat2Svs1S = new RooRealVar("R_{#frac{2S}{1S}}^{PbPb/PP}", "R_{#frac{2S}{1S}}^{PbPb/PP}", 0.26, 0.0, 3.0);
      RooFormulaVar *tmp0 = new RooFormulaVar("R_{#frac{2S}{1S}}^{PbPb}","@0*@1", RooArgList(*w.var("R_{#frac{2S}{1S}}^{PP}"),*fDRat2Svs1S));
      fDRat2Svs1S->setConstant(kFALSE);
      f2Svs1S = (RooRealVar*)tmp0;
    } else {
      f2Svs1S   = new RooRealVar(Form("R_{#frac{2S}{1S}}^{%s}", (isPbPb?"PbPb":"PP")),"R_{#frac{2S}{1S}}",0.26,0.0,1.0);
    }
    RooFormulaVar *tmp1 = new RooFormulaVar(Form("N_{#psi(2S)}^{%s}", (isPbPb?"PbPb":"PP")),"@0*@1", RooArgList(*nSigJPSI,*f2Svs1S));
    f2Svs1S->setConstant(kFALSE);
    nSigPSI2S = (RooRealVar*)tmp1;
  }

  RooRealVar *meanSigJPSI   = new RooRealVar(Form("m_{J/#psi}^{%s}", (isPbPb?"PbPb":"PP")),"m_{J/#psi}",Mass.JPsi,Mass.JPsi-0.1,Mass.JPsi+0.1);
  RooConstVar *rat = new RooConstVar(Form("ratio_%s", (isPbPb?"PbPb":"PP")), "m_{#psi(2S)}/m_{J/#psi}", Mass.Psi2S/Mass.JPsi);
  RooRealVar *sigmaSigJPSI  = new RooRealVar(Form("#sigma_{CB1}^{%s}", (isPbPb?"PbPb":"PP")),"#sigma_{CB1}",0.04, 0.01, 0.09); //0.1  work with 0.11
  RooRealVar *sigmaSigJPSI2  = new RooRealVar(Form("#sigma_{CB2}^{%s}", (isPbPb?"PbPb":"PP")),"#sigma_{CB2}",0.02, 0.01, 0.07); //0.03 
  RooRealVar *sigmaSigJPSI3  = new RooRealVar(Form("#sigma_{CB3}^{%s}", (isPbPb?"PbPb":"PP")),"#sigma_{CB3}",0.03, 0.01, 0.07); //0.03 

  RooRealVar *alpha = new RooRealVar(Form("#alpha_{CB}^{%s}", (isPbPb?"PbPb":"PP")),"#alpha_{CB}",2.21, 0.5, 30.0*100.0);//30.*100.    work with 1.01 , 5.0
  RooRealVar *n  = new RooRealVar(Form("n_{CB}^{%s}", (isPbPb?"PbPb":"PP")),"n_{CB}",5., 0.5, 30.*10.);  
  RooRealVar *nW = new RooRealVar(Form("nW_{CB}^{%s}", (isPbPb?"PbPb":"PP")),"nW_{CB}",2., 0.5, 30.*100.); 
  
  if (opt.doSimulFit && isPbPb) {
    meanSigJPSI   = w.var("m_{J/#psi}^{PP}");
    //sigmaSigJPSI  = w.var("#sigma_{CB1}^{PP}");
    //sigmaSigJPSI2 = w.var("#sigma_{CB2}^{PP}");
    alpha = w.var("#alpha_{CB}^{PP}");
    n     = w.var("n_{CB}^{PP}");
    nW    = w.var("nW_{CB}^{PP}");
  }

  RooFormulaVar *meanSigPSI2S = new RooFormulaVar(Form("m_{#psi(2S)}^{%s}", (isPbPb?"PbPb":"PP")),"@0*@1",RooArgList(*meanSigJPSI,*rat));
  RooFormulaVar *sigmaSigPSI2S  = new RooFormulaVar(Form("#sigma_{CB1PSI}^{%s}", (isPbPb?"PbPb":"PP")),"@0*@1",RooArgList(*sigmaSigJPSI,*rat)); // 
  RooFormulaVar *sigmaSigPSI2S2  = new RooFormulaVar(Form("#sigma_{CB2PSI}^{%s}", (isPbPb?"PbPb":"PP")),"@0*@1",RooArgList(*sigmaSigJPSI2,*rat)); //
  RooFormulaVar *sigmaSigPSI2S3  = new RooFormulaVar(Form("#sigma_{CB3PSI}^{%s}", (isPbPb?"PbPb":"PP")),"@0*@1",RooArgList(*sigmaSigJPSI3,*rat)); //

  //////// Candidates for signal
  // Normal gaussians
  RooGaussian* signalG1  = NULL;
  RooGaussian* signalG2  = NULL;
  RooGaussian* signalG3  = NULL;
  RooGaussian* signalG4  = NULL;
  RooCBShape*  signalCB1 = NULL;
  RooCBShape*  signalCB2 = NULL;
  RooCBShape*  signalCB1WN = NULL;
  RooCBShape*  signalCB2WN = NULL;
  RooCBShape*  signalCB3WN = NULL;
  RooCBShape*  signalCB4WN = NULL;

  RooAbsPdf* sigJPSI  = NULL;   
  RooAbsPdf* sigPSI2S = NULL; 
  RooRealVar *coeffGaus = new RooRealVar(Form("f_{G}^{%s}", (isPbPb?"PbPb":"PP")),"coeffGaus",0.1, 0.0,1.0);
  RooRealVar *coeffCB = new RooRealVar(Form("f_{CB}^{%s}", (isPbPb?"PbPb":"PP")),"coeffCB",0.3, 0.0,1.0);
  RooRealVar *coeffCB2 = new RooRealVar(Form("f_{CB2}^{%s}", (isPbPb?"PbPb":"PP")),"coeffCB2",0.3, 0.0,1.0);

  if (opt.doSimulFit && isPbPb) {
    // coeffCB = w.var("f_{CB}^{PP}");
  }

  switch(sigModel) 
    {    
    case 1: 
      // Gaussian for JPsi and Psi(2S)
      sigJPSI  = new RooGaussian(Form("sigJPSI_%s", (isPbPb?"PbPb":"PP")),"Gaussian Sig1", *mass, *meanSigJPSI, *sigmaSigJPSI); //signalG1
      sigPSI2S  = new RooGaussian(Form("sigPSI2S_%s", (isPbPb?"PbPb":"PP")),"Gaussian Sig2", *mass, *meanSigPSI2S, *sigmaSigPSI2S); //signalG2
      cout << "you're fitting 2 signal peaks with Gaussian functions"<< endl;
      break;  
    case 2: 
       // Gaussian for JPsi and Psi(2S)
      signalG1  = new RooGaussian(Form("signalG1_%s", (isPbPb?"PbPb":"PP")),"Gaussian Sig1", *mass, *meanSigJPSI, *sigmaSigJPSI); //signalG1
      signalCB1WN = new RooCBShape(Form("signalCB1WN_%s", (isPbPb?"PbPb":"PP")),"FSR cb 1s", *mass, *meanSigJPSI, *sigmaSigJPSI2, *alpha, *nW);
      sigJPSI    = new RooAddPdf (Form("sigJPSI_%s", (isPbPb?"PbPb":"PP")), "sigJPSI", RooArgList(*signalG1, *signalCB1WN), *coeffGaus);
      signalG2  = new RooGaussian(Form("signalG2_%s", (isPbPb?"PbPb":"PP")),"Gaussian Sig2", *mass, *meanSigPSI2S, *sigmaSigPSI2S); //signalG2
      signalCB2WN = new RooCBShape(Form("signalCB2WN_%s", (isPbPb?"PbPb":"PP")),"FSR cb 1s", *mass, *meanSigPSI2S, *sigmaSigPSI2S2, *alpha, *nW);
      sigPSI2S  = new RooAddPdf (Form("sigPSI2S_%s", (isPbPb?"PbPb":"PP")), "sigPSI2S", RooArgList(*signalG2, *signalCB2WN), *coeffGaus);
      cout << "you're fitting 2 signal peaks with Gaussian functions"<< endl;
      break;   
    case 3: 
      // Gaussian + Cristal Ball for JPsi
      signalG1  = new RooGaussian(Form("signalG1_%s", (isPbPb?"PbPb":"PP")),"Gaussian Sig1", *mass, *meanSigJPSI, *sigmaSigJPSI); //signalG1
      signalCB1WN = new RooCBShape(Form("signalCB1WN_%s", (isPbPb?"PbPb":"PP")),"FSR cb 1s", *mass, *meanSigJPSI, *sigmaSigJPSI2, *alpha, *nW);
      sigJPSI    = new RooAddPdf (Form("sigJPSI_%s", (isPbPb?"PbPb":"PP")), "sigJPSI", RooArgList(*signalG1, *signalCB1WN), *coeffGaus);
      cout << "you're fitting 1 signal peak with a Gaussian function"<< endl;
      break;  
    case 4: 
      // Gaussian + Gaussian for JPsi
      signalG1  = new RooGaussian(Form("signalG1_%s", (isPbPb?"PbPb":"PP")),"Gaussian Sig1", *mass, *meanSigJPSI, *sigmaSigJPSI); //signalG1
      signalG2 = new RooGaussian(Form("signalG2_%s", (isPbPb?"PbPb":"PP")),"Gaussian Sig2", *mass, *meanSigJPSI, *sigmaSigJPSI2);
      sigJPSI    = new RooAddPdf (Form("sigJPSI_%s", (isPbPb?"PbPb":"PP")), "sigJPSI", RooArgList(*signalG1, *signalG2), *coeffGaus);
      cout << "you're fitting 1 signal peak with a Gaussian function"<< endl;
      break;   
    case 5: 
      // Crystal Ball for JPsi
      //sigJPSI = new RooGaussian(Form("sigJPSI_%s", (isPbPb?"PbPb":"PP")),"Gaussian Sig1", *mass, *meanSigJPSI, *sigmaSigJPSI); //signalG1
      sigJPSI = new RooCBShape(Form("sigJPSI_%s", (isPbPb?"PbPb":"PP")),"Gaussian Sig1", *mass, *meanSigJPSI, *sigmaSigJPSI, *alpha, *nW); //signalG1
      //signalCB1WN = new RooCBShape(Form("signalCB1WN_%s", (isPbPb?"PbPb":"PP")),"FSR cb 1s", *mass, *meanSigJPSI, *sigmaSigJPSI, *alpha, *nW);
      ////signalG1    = new RooGaussian(Form("sigJPSI_%s", (isPbPb?"PbPb":"PP")),"Gaussian Sig1", *mass, *meanSigJPSI, *sigmaSigJPSI); //signalG1
      //signalCB2WN = new RooCBShape(Form("signalCB2WN_%s", (isPbPb?"PbPb":"PP")),"FSR cb 1s", *mass, *meanSigJPSI, *sigmaSigJPSI2, *alpha, *nW);
      //sigJPSI    = new RooAddPdf (Form("sigJPSI_%s", (isPbPb?"PbPb":"PP")), "sigJPSI", RooArgList(*signalCB1WN, *signalCB2WN), *coeffCB);
      ////sigJPSI    = new RooAddPdf (Form("sigJPSI_%s", (isPbPb?"PbPb":"PP")), "sigJPSI", RooArgList(*signalG1, *signalCB2WN), *coeffCB);
      // signalCB3WN = new RooCBShape(Form("signalCB3WN_%s", (isPbPb?"PbPb":"PP")),"FSR cb 1s", *mass, *meanSigJPSI, *sigmaSigJPSI3, *alpha, *nW);
      //sigJPSI    = new RooAddPdf (Form("sigJPSI_%s", (isPbPb?"PbPb":"PP")), "sigJPSI", RooArgList(*signalCB1WN, *signalCB2WN, *signalCB3WN), RooArgList(*coeffCB, *coeffCB2));
      cout << "you're fitting 1 signal peak with a Crystal Ball function"<< endl;
      break;
    case 6:
      // Currently used in JPsi analysis
      // Sum of gaussian 1 and crystall ball 2 with wide n  for JPsi
      signalG1    = new RooGaussian(Form("signalG1_%s", (isPbPb?"PbPb":"PP")),"Gaussian Sig1", *mass, *meanSigJPSI, *sigmaSigJPSI); //signalG1
      signalCB2WN = new RooCBShape(Form("signalCB2WN_%s", (isPbPb?"PbPb":"PP")),"FSR cb 1s", *mass, *meanSigJPSI, *sigmaSigJPSI2, *alpha, *nW);
      sigJPSI    = new RooAddPdf (Form("sigJPSI_%s", (isPbPb?"PbPb":"PP")), "sigPDF", RooArgList(*signalG1, *signalCB2WN), *coeffGaus);
      cout << "you're fitting 1 signal peak with a sum of a Gaussian and a Crystal Ball function"<< endl;
      break;
    case 7:
      // Sum of gaussian 1 and a crystall ball for JPsi
      signalG1  = new RooGaussian(Form("signalG1_%s", (isPbPb?"PbPb":"PP")),"Gaussian Sig1", *mass, *meanSigJPSI, *sigmaSigJPSI);
      signalCB1 = new RooCBShape(Form("signalCB1_%s", (isPbPb?"PbPb":"PP")), "FSR cb 1s", *mass, *meanSigJPSI, *sigmaSigJPSI, *alpha, *n);
      sigJPSI    = new RooAddPdf (Form("sigJPSI_%s", (isPbPb?"PbPb":"PP")), "sigPDF", RooArgList(*signalG1, *signalCB1), *coeffGaus);
      cout << "you're fitting 1 signal peak with a sum of a Gaussian and a Crystal Ball function"<< endl;
      break;
    case 8:
      // Sum of gaussian 1 and a crystall ball for JPsi
      signalG1  = new RooGaussian(Form("signalG1_%s", (isPbPb?"PbPb":"PP")),"Gaussian Sig1", *mass, *meanSigJPSI, *sigmaSigJPSI);
      signalCB2 = new RooCBShape(Form("signalCB2_%s", (isPbPb?"PbPb":"PP")), "FSR cb 1s", *mass, *meanSigJPSI, *sigmaSigJPSI2, *alpha, *n);
      sigJPSI    = new RooAddPdf (Form("sigJPSI_%s", (isPbPb?"PbPb":"PP")), "sigPDF", RooArgList(*signalG1, *signalCB2), *coeffGaus);
      cout << "you're fitting 1 signal peak with a sum of a Gaussian and a Crystal Ball function"<< endl;
      break;
    case 9:
      // Sum of gaussian 1 and crystall ball with wide n for JPsi
      signalG1    = new RooGaussian(Form("signalG1_%s", (isPbPb?"PbPb":"PP")),"Gaussian Sig1", *mass, *meanSigJPSI, *sigmaSigJPSI);
      signalCB1WN = new RooCBShape(Form("signalCB1WN_%s", (isPbPb?"PbPb":"PP")),"FSR cb 1s", *mass, *meanSigJPSI, *sigmaSigJPSI, *alpha, *nW);
      sigJPSI    = new RooAddPdf (Form("sigJPSI_%s", (isPbPb?"PbPb":"PP")), "sigPDF", RooArgList(*signalG1, *signalCB1WN), *coeffGaus);
      cout << "you're fitting 1 signal peak with a sum of a Gaussian and a Crystal Ball function"<< endl;
      break;
    case 10: 
       // Gaussian for JPsi and Psi(2S)
      signalG1  = new RooGaussian(Form("signalG1_%s", (isPbPb?"PbPb":"PP")),"Gaussian Sig1", *mass, *meanSigJPSI, *sigmaSigJPSI); //signalG1
      signalG2 = new RooGaussian(Form("signalCB1WN_%s", (isPbPb?"PbPb":"PP")),"FSR cb 1s", *mass, *meanSigJPSI, *sigmaSigJPSI2);
      sigJPSI    = new RooAddPdf (Form("sigJPSI_%s", (isPbPb?"PbPb":"PP")), "sigJPSI", RooArgList(*signalG1, *signalG2), *coeffGaus);
      signalG3  = new RooGaussian(Form("signalG2_%s", (isPbPb?"PbPb":"PP")),"Gaussian Sig2", *mass, *meanSigPSI2S, *sigmaSigPSI2S); //signalG2
      signalG4 = new RooGaussian(Form("signalCB2WN_%s", (isPbPb?"PbPb":"PP")),"FSR cb 1s", *mass, *meanSigPSI2S, *sigmaSigPSI2S2);
      sigPSI2S  = new RooAddPdf (Form("sigPSI2S_%s", (isPbPb?"PbPb":"PP")), "sigPSI2S", RooArgList(*signalG3, *signalG4), *coeffGaus);
      cout << "you're fitting 2 signal peaks with Gaussian functions"<< endl;
      cout << "you're fitting 2 signal peaks with Gaussian functions"<< endl;
      break;   
    case 11: 
       // Gaussian for JPsi and Psi(2S)
      /*
      sigJPSI = new RooCBShape(Form("sigJPSI_%s", (isPbPb?"PbPb":"PP")),"FSR cb 1s", *mass, *meanSigJPSI, *sigmaSigJPSI, *alpha, *nW);
      sigPSI2S = new RooCBShape(Form("sigPSI2S_%s", (isPbPb?"PbPb":"PP")),"FSR cb 1s", *mass, *meanSigPSI2S, *sigmaSigPSI2S, *alpha, *nW);
      */
      
      signalCB1WN = new RooCBShape(Form("signalCB1WN_%s", (isPbPb?"PbPb":"PP")),"FSR cb 1s", *mass, *meanSigJPSI, *sigmaSigJPSI, *alpha, *nW);
      signalCB2WN = new RooCBShape(Form("signalCB2WN_%s", (isPbPb?"PbPb":"PP")),"FSR cb 1s", *mass, *meanSigJPSI, *sigmaSigJPSI2, *alpha, *nW);
      sigJPSI    = new RooAddPdf (Form("sigJPSI_%s", (isPbPb?"PbPb":"PP")), "sigJPSI", RooArgList(*signalCB1WN, *signalCB2WN), *coeffCB);
      signalCB3WN = new RooCBShape(Form("signalCB3WN_%s", (isPbPb?"PbPb":"PP")),"FSR cb 1s", *mass, *meanSigPSI2S, *sigmaSigPSI2S, *alpha, *nW);
      signalCB4WN = new RooCBShape(Form("signalCB4WN_%s", (isPbPb?"PbPb":"PP")),"FSR cb 1s", *mass, *meanSigPSI2S, *sigmaSigPSI2S2, *alpha, *nW);
      sigPSI2S  = new RooAddPdf (Form("sigPSI2S_%s", (isPbPb?"PbPb":"PP")), "sigPSI2S", RooArgList(*signalCB3WN, *signalCB4WN), *coeffCB);
      
      cout << "you're fitting 2 signal peaks with Gaussian functions"<< endl;
      break;
    case 12: 
      sigJPSI = new RooCBShape(Form("sigJPSI_%s", (isPbPb?"PbPb":"PP")),"FSR cb 1s", *mass, *meanSigJPSI, *sigmaSigJPSI, *alpha, *nW);
      sigPSI2S = new RooCBShape(Form("sigPSI2S_%s", (isPbPb?"PbPb":"PP")),"FSR cb 1s", *mass, *meanSigPSI2S, *sigmaSigPSI2S, *alpha, *nW);
      
      cout << "you're fitting 2 signal peaks with Gaussian functions"<< endl;
      break;
    case 13: 
      // Gaussian + Cristal Ball for JPsi
      signalG1  = new RooGaussian(Form("signalG1_%s", (isPbPb?"PbPb":"PP")),"Gaussian Sig1", *mass, *meanSigPSI2S, *sigmaSigPSI2S); //signalG1
      signalCB1WN = new RooCBShape(Form("signalCB1WN_%s", (isPbPb?"PbPb":"PP")),"FSR cb 1s", *mass, *meanSigPSI2S, *sigmaSigPSI2S2, *alpha, *nW);
      sigPSI2S    = new RooAddPdf (Form("sigJPSI_%s", (isPbPb?"PbPb":"PP")), "sigJPSI", RooArgList(*signalG1, *signalCB1WN), *coeffGaus);
      cout << "you're fitting 1 signal peak with a Gaussian function"<< endl;
      break;  
    case 14: 
      // Gaussian + Gaussian for JPsi
      signalG1  = new RooGaussian(Form("signalG1_%s", (isPbPb?"PbPb":"PP")),"Gaussian Sig1", *mass, *meanSigPSI2S, *sigmaSigPSI2S); //signalG1
      signalG2 = new RooGaussian(Form("signalG2_%s", (isPbPb?"PbPb":"PP")),"Gaussian Sig2", *mass, *meanSigPSI2S, *sigmaSigPSI2S2);
      sigPSI2S    = new RooAddPdf (Form("sigJPSI_%s", (isPbPb?"PbPb":"PP")), "sigJPSI", RooArgList(*signalG1, *signalG2), *coeffGaus);
      cout << "you're fitting 1 signal peak with a Gaussian function"<< endl;
      break;      
    default :
      cout<<"Donno what you are talking about! Pick another fit option for signal!"<<endl;
      break;
  }
      
  // Background Model  
  RooRealVar *nBkg   = new RooRealVar(Form("N_{bkg}^{%s}", (isPbPb?"PbPb":"PP")),"nbkgd",0,nt);
  RooRealVar *coeffPoll1 = new RooRealVar(Form("#lambda_{p1}^{%s}", (isPbPb?"PbPb":"PP")),"coeffPoll",0.05, -5.*100., 5.*100.);
  RooRealVar *coeffPoll2 = new RooRealVar(Form("#lambda_{p2}^{%s}", (isPbPb?"PbPb":"PP")),"coeffPoll1",0.05, -5.*100., 5.*100.);
  RooRealVar *coeffPoll3 = new RooRealVar(Form("#lambda_{p3}^{%s}", (isPbPb?"PbPb":"PP")),"coeffPoll2",0.05, -5.*100., 5.*100.);
  RooRealVar *coeffCheb1 = new RooRealVar(Form("#lambda_{CB1}^{%s}", (isPbPb?"PbPb":"PP")),"coeffCheb1",0.05,-1.5,1.5);
  RooRealVar *coeffCheb2 = new RooRealVar(Form("#lambda_{CB2}^{%s}", (isPbPb?"PbPb":"PP")),"coeffCheb2",0.05,-1.0,1.0);
  RooRealVar *coeffCheb3 = new RooRealVar(Form("#lambda_{CB3}^{%s}", (isPbPb?"PbPb":"PP")),"coeffCheb3",0.05,-1.0,1.0);
  RooRealVar *coeffExp  = new RooRealVar(Form("#lambda_{exp}^{%s}", (isPbPb?"PbPb":"PP")),"coeffExp",0.1, -3., 1.*100.);
  RooRealVar turnOn(Form("turnOn_%s", (isPbPb?"PbPb":"PP")),"turnOn", 0,10);
  RooRealVar width(Form("width_%s", (isPbPb?"PbPb":"PP")),"width", 0.1,10);// MB 2.63
  RooRealVar decay(Form("decay_%s", (isPbPb?"PbPb":"PP")),"decay",0,10);// MB: 3.39
  RooAbsPdf* bkgPDF = NULL;

  switch(bkgModel) {  
  case 1: 
    // 1st order polynomial
    bkgPDF = new RooPolynomial(Form("bkgPDF_%s", (isPbPb?"PbPb":"PP")), "bkgPoll1", *mass, *coeffPoll1);
    cout << "you're fitting the background with a 1st order polynomial"<< endl;
    break;
  case 2:
    // Exponential
    bkgPDF = new RooExponential(Form("bkgPDF_%s", (isPbPb?"PbPb":"PP")), "bkgExp", *mass, *coeffExp);
    cout << "you're fitting the background with an exponential"<< endl;
    break;
  case 3:
    // Chebychev 1st
    bkgPDF = new RooChebychev(Form("bkgPDF_%s", (isPbPb?"PbPb":"PP")), "bkgCheb1", *mass, *coeffCheb1);
    cout << "you're fitting the background with an Chebychev 1st"<< endl;
    break;
  case 4:
    // Chebychev 2nd
    bkgPDF = new RooChebychev(Form("bkgPDF_%s", (isPbPb?"PbPb":"PP")), "bkgCheb2", *mass, {*coeffCheb1, *coeffCheb2});
    cout << "you're fitting the background with an Chebychev 2nd"<< endl;
    break;
  case 5:
    // Chebychev 3rd
    bkgPDF = new RooChebychev(Form("bkgPDF_%s", (isPbPb?"PbPb":"PP")), "bkgCheb3", *mass, {*coeffCheb1, *coeffCheb2,*coeffCheb3});
    cout << "you're fitting the background with an Chebychev 3rd"<< endl;
    break;  
  case 6: 
    // 2st order polynomial
    bkgPDF = new RooPolynomial(Form("bkgPDF_%s", (isPbPb?"PbPb":"PP")), "bkgPoll2", *mass, {*coeffPoll1, *coeffPoll2});
    cout << "you're fitting the background with a 1st order polynomial"<< endl;
    break;  
  case 7: 
    // 3st order polynomial
    bkgPDF = new RooPolynomial(Form("bkgPDF_%s", (isPbPb?"PbPb":"PP")), "bkgPoll3", *mass, {*coeffPoll1, *coeffPoll2,*coeffPoll3});
    cout << "you're fitting the background with a 1st order polynomial"<< endl;
    break;
  case 8:
    // ErrFn
    width.setConstant(false);
    decay.setConstant(false);
    turnOn.setConstant(false);
    bkgPDF = new  RooGenericPdf(Form("bkgPDF_%s", (isPbPb?"PbPb":"PP")),"ErrFn", "exp(-@0/decay)*(TMath::Erf((@0-turnOn)/width)+1)", 
				RooArgList(*mass,turnOn,width,decay));
    cout << "you're fitting the background with an Errfn"<< endl;  
    default :
      cout<<"Donno what you are talking about! Pick another fit option for background!"<<endl;
      break;
  }

  // Total PDF = Signal + Background
  RooAbsPdf  *pdf  = new RooAddPdf (Form("pdf_%s", (isPbPb?"PbPb":"PP")),"total p.d.f.", RooArgList(*sigJPSI, *sigPSI2S, *bkgPDF), RooArgList(*nSigJPSI, *nSigPSI2S, *nBkg));

  w.import(*pdf);
  w.Print();

}
