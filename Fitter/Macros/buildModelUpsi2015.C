#include "allFunctions.h"

using namespace RooFit;

void buildModelUpsi2015(RooWorkspace& w, int sigModel, int bkgModel,bool fitExcited){
   // C r e a t e   m o d e l  
   int nt=100000;
   // cout << "you're building a model for the quarkonium resonance of mass = "<< M1S <<" GeV/c^{2},"endl;
   RooRealVar *nsig1f   = new RooRealVar("N_{#varUpsilon(1S)}","nsig1S",0,nt*10);
   RooRealVar* mass       = (RooRealVar*) w.var("invariantMass");

  

   RooRealVar *nsig2f = new RooRealVar("N_{#varUpsilon(2S)}","nsig2S",0,nt*10);
   RooRealVar *nsig3f = new RooRealVar("N_{#varUpsilon(3S)}","nsig3S",0,nt*10); 
   RooRealVar *f2Svs1S   = NULL;
   RooRealVar *f3Svs1S   = NULL;

   bool doRatio = true;

   if(doRatio && fitExcited){
   f2Svs1S   = new RooRealVar("R_{#frac{2S}{1S}}","f2Svs1S",0.26,-0.1,1.0);
   f3Svs1S   = new RooRealVar("R_{#frac{3S}{1S}}","f3Svs1S",0.13,-0.1,1.0);
   RooFormulaVar *tmp1 = new RooFormulaVar("N_{ #varUpsilon(2S)}","@0*@1", RooArgList(*nsig1f,*f2Svs1S));
   RooFormulaVar *tmp2 = new RooFormulaVar("N_{ #varUpsilon(3S)}","@0*@1", RooArgList(*nsig1f,*f3Svs1S));
   f2Svs1S->setConstant(kFALSE);
   f3Svs1S->setConstant(kFALSE);
   nsig2f = (RooRealVar*)tmp1;
   nsig3f = (RooRealVar*)tmp2;
   }



   if(!fitExcited){
   nsig2f = NULL;
   nsig3f = NULL;  
   }
   
   RooRealVar  *mean = new RooRealVar("m_{ #varUpsilon(1S)}","#Upsilon mean", Mass.Y1S, Mass.Y1S-0.2, Mass.Y1S+0.2);
   RooConstVar *rat2 = new RooConstVar("rat2", "rat2", Mass.Y2S/Mass.Y1S);
   RooConstVar *rat3 = new RooConstVar("rat3", "rat3", Mass.Y3S/Mass.Y1S);
   // scale mean and resolution by mass ratio
   RooFormulaVar *mean1S = new RooFormulaVar("mean1S","@0",RooArgList(*mean));
   RooFormulaVar *mean2S = new RooFormulaVar("mean2S","@0*@1", RooArgList(*mean,*rat2));
   RooFormulaVar *mean3S = new RooFormulaVar("mean3S","@0*@1", RooArgList(*mean,*rat3));

   // //detector resolution ?? where is this coming from?
   RooRealVar    *sigma1  = new RooRealVar("#sigma_{CB1}","#sigma_{CB1}",0.01,0.1); // 
   RooFormulaVar *sigma1S = new RooFormulaVar("sigma1S","@0"   ,RooArgList(*sigma1));
   RooFormulaVar *sigma2S = new RooFormulaVar("sigma2S","@0*@1",RooArgList(*sigma1,*rat2));
   RooFormulaVar *sigma3S = new RooFormulaVar("sigma3S","@0*@1",RooArgList(*sigma1,*rat3));
   RooRealVar *alpha  = new RooRealVar("#alpha_{CB}","tail shift",0.1,10);    // MC 5tev 1S pol2 
   RooRealVar *npow   = new RooRealVar("n_{CB}","power order",0.1,10);    // MC 5tev 1S pol2 
   RooRealVar *sigmaFraction = new RooRealVar("sigmaFraction","Sigma Fraction",0.,1.);
   // scale the sigmaGaus with sigma1S*scale=sigmaGaus now.
   RooRealVar    *scaleWidth = new RooRealVar("#sigma_{CB2}/#sigma_{CB1}","scaleWidth",1.,2.5);
   RooFormulaVar *sigmaGaus = new RooFormulaVar("sigmaGaus","@0*@1", RooArgList(*sigma1,*scaleWidth));
   RooFormulaVar *sigmaGaus2 = new RooFormulaVar("sigmaGaus","@0*@1*@2", RooArgList(*sigma1,*scaleWidth,*rat2));
   RooFormulaVar *sigmaGaus3 = new RooFormulaVar("sigmaGaus","@0*@1*@2", RooArgList(*sigma1,*scaleWidth,*rat3));
   RooGaussian* gauss1 = new RooGaussian("gaus1s","gaus1s",
         *nsig1f,
         *mass,    //mean
         *sigmaGaus); //sigma
   RooGaussian* gauss1b = new RooGaussian("gaus1sb","gaus1sb",
         *nsig1f,
         *mean,    //mean
         *sigma1); //sigma

   // declarations
   RooAbsPdf  *cb1S_1   = NULL;
   RooAbsPdf  *cb1S_2   = NULL;
   RooAbsPdf  *cb2S_1   = NULL;
   RooAbsPdf  *cb2S_2   = NULL;
   RooAbsPdf  *cb3S_1   = NULL;
   RooAbsPdf  *cb3S_2   = NULL;
   RooAbsPdf  *sig1S   = NULL;
   RooAbsPdf  *sig2S   = NULL;
   RooAbsPdf  *sig3S   = NULL;

   switch(sigModel){   
         
      case 1: //Crystal Balls
         sig1S   = new RooCBShape("cb1S_1", "FSR cb 1s",
               *mass,*mean1S,*sigma1,*alpha,*npow);

         sig2S   = new RooCBShape ("cb2S_1", "FSR cb 1s",
               *mass,*mean2S,*sigma2S,*alpha,*npow);
         sig3S   = new RooCBShape ("cb3S_1", "FSR cb 1s",
               *mass,*mean3S,*sigma3S,*alpha,*npow);
         cout << "you're fitting each signal peak with a Crystal Ball function"<< endl;
         break;
      case 2: //Two Crystal Balls

	 //----Upsilon 1S
         cb1S_1    = new RooCBShape ("cb1S_1", "FSR cb 1s",
               *mass,*mean1S,*sigma1,*alpha,*npow);

         cb1S_2    = new RooCBShape ("cb1S_2", "FSR cb 1s",
               *mass,*mean1S,*sigmaGaus,*alpha,*npow);
         sig1S  = new RooAddPdf  ("cbcb","1S mass pdf",
               RooArgList(*cb1S_1,*cb1S_2),*sigmaFraction);
         //----Upsilon 2S
         cb2S_1    = new RooCBShape ("cb2S_1", "FSR cb 2s", 
               *mass,*mean2S,*sigma2S,*alpha,*npow); 
         cb2S_2    = new RooCBShape ("cb2S_2", "FSR cb 2s", 
               *mass,*mean2S,*sigmaGaus2,*alpha,*npow); 
         sig2S  = new RooAddPdf  ("sig2S","2S mass pdf",
               RooArgList(*cb2S_1,*cb2S_2),*sigmaFraction);

         //----Upsilon 3S
         cb3S_1    = new RooCBShape ("cb3S_1", "FSR cb 3s", 
               *mass,*mean3S,*sigma3S,*alpha,*npow); 
         cb3S_2    = new RooCBShape ("cb3S_2", "FSR cb 3s", 
               *mass,*mean3S,*sigmaGaus3,*alpha,*npow); 
         sig3S  = new RooAddPdf  ("sig3S","3S mass pdf",
               RooArgList(*cb3S_1,*cb3S_2),*sigmaFraction); // = cb3S1*sigmaFrac + cb3S2*(1-sigmaFrac)
         cout << "you're fitting each signal peak with a Double Crystal Ball function"<< endl;
         break;

      //---Cases below are fore single peak fitting
      case 3: //Gaussian
         sig1S  = new RooGaussian ("g1", "gaus 1s",
               *mass,*mean1S,*sigma1);
         cout << "you're fitting 1 signal peak with a Gaussian function"<< endl;
         break;

      case 4: //Crystal Ball
         cb1S_1    = new RooCBShape ("cb1S_1", "FSR cb 1s",
               *mass,*mean1S,*sigma1,*alpha,*npow);
         sig1S  = new RooAddPdf ("cbg", "cbgaus 1s",
               RooArgList(*cb1S_1),*sigma1);
         cout << "you're fitting 1 signal peak with a sum of  a Crystal Ball function"<< endl;
         break;

      case 5: //Gaussian + Crystal Ball
         cb1S_1    = new RooCBShape ("cb1S_1", "FSR cb 1s",
               *mass,*mean1S,*sigma1,*alpha,*npow);
         sig1S  = new RooAddPdf ("cbg", "cbgaus 1s",
               RooArgList(*gauss1,*cb1S_1),*sigmaFraction);
         cout << "you're fitting 1 signal peak with a sum of a Gaussian and a Crystal Ball function"<< endl;
         break;

      case 6: //Two Gaussians 
         sig1S  = new RooAddPdf ("cb1S_1", "cbgaus 1s",
               RooArgList(*gauss1,*gauss1b),*sigmaFraction);
         cout << "you're fitting each signal peak with a Double Gaussian function"<< endl;
         break;
   }
   // bkg Chebychev
   RooRealVar *nbkgd   = new RooRealVar("n_{Bkgd}","nbkgd",0,nt);
   RooRealVar *bkg_a1  = new RooRealVar("a1_bkg", "bkg_{a1}", 0, -5, 5);
   RooRealVar *bkg_a2  = new RooRealVar("a2_Bkg", "bkg_{a2}", 0, -2, 2);
   RooRealVar *bkg_a3  = new RooRealVar("a3_Bkg", "bkg_{a3}", 0, -0.9, 2);

   //  likesign
   RooRealVar *nLikesignbkgd = new RooRealVar("NLikesignBkg","nlikesignbkgd",nt*0.75,0,10*nt);
   // *************************************************** bkgModel

   RooRealVar turnOn("turnOn","turnOn", 0,10);
   RooRealVar width("width","width", 0.1,10);// MB 2.63
   RooRealVar decay("decay","decay",0,10);// MB: 3.39


   width.setConstant(false);
   decay.setConstant(false);
   turnOn.setConstant(false);


   alpha->setConstant(false);
   npow->setConstant(false);
   sigma1->setConstant(false);
   scaleWidth->setConstant(false);
   sigmaFraction->setConstant(false);

   //thisPdf: form of the bkg pdf
   //pdf_combinedbkgd; // total bkg pdf. usually form*normalization  (so that you can do extended ML fits)

   RooAbsPdf *pdf_combinedbkgd = NULL;;
   RooRealVar *fPol   = new RooRealVar("F_{pol}","fraction of polynomial distribution",0.0,1);
   RooAbsPdf *ChebPdf = NULL;
   RooAbsPdf *ErrPdf = NULL;
   RooAbsPdf* ExpPdf = NULL;

   switch (bkgModel) 
   {
      case 1 : //use polynomial
	 ChebPdf = new RooChebychev("ChebPdf","ChebPdf",
		*mass, RooArgList(*bkg_a1,*bkg_a2));
	 pdf_combinedbkgd = new RooAddPdf ("bkgPdf","total combined pdf",
		RooArgList(*ChebPdf));
        break;

      case 2 : //use error function to fit the OS directly
         bkg_a3->setConstant(true);
         pdf_combinedbkgd            = new  RooGenericPdf("bkgPdf","bkgPdf",
               "exp(-@0/decay)*(TMath::Erf((@0-turnOn)/width)+1)",
               RooArgList(*mass,turnOn,width,decay));
         break;

      case 3 : //use pol 2+ErfExp to fit the OS directly

         ChebPdf  = new RooChebychev("ChebPdf","ChebPdf",
               *mass, RooArgList(*bkg_a1,*bkg_a2));
         ErrPdf     = new  RooGenericPdf("ErrPdf","ErrPdf",
               "exp(-@0/decay)*(TMath::Erf((@0-turnOn)/width)+1)",
               RooArgList(*mass,turnOn,width,decay));
         pdf_combinedbkgd   = new RooAddPdf ("bkgPdf","total combined background pdf",
               RooArgList(*ChebPdf,*ErrPdf),
               RooArgList(*fPol));

         break;
      case 4 : //use ( error function + polynomial 1) to fit the OS directly

         bkg_a3->setConstant(true);
         bkg_a2->setConstant(true);
         ChebPdf  = new RooChebychev("ChebPdf","ChebPdf",
               *mass, RooArgList(*bkg_a1,*bkg_a2,*bkg_a3));
         ErrPdf     = new  RooGenericPdf("ErrPdf","ErrPdf",
               "exp(-@0/decay)*(TMath::Erf((@0-turnOn)/width)+1)",
               RooArgList(*mass,turnOn,width,decay));
         pdf_combinedbkgd   = new RooAddPdf ("bkgPdf","total combined background pdf",
               RooArgList(*ChebPdf,*ErrPdf),
               RooArgList(*fPol));
         break;
      case 5: // NOT WORKING
         ChebPdf  = new RooChebychev("ChebPdf","ChebPdf",
               *mass, RooArgList(*bkg_a1,*bkg_a2));
         ExpPdf     = new  RooGenericPdf("ExpPdf","ExpPdf",
               "exp(-@0/decay)",
               RooArgList(*mass,decay));
         pdf_combinedbkgd  =  new RooAddPdf ("bkgPdf","total combined background pdf",
               RooArgList(*ChebPdf,*ExpPdf),
               RooArgList(*fPol));
         break;
      default :
         cout<<"Donno what you are talking about! Pick another fit option!"<<endl;
         break;
   }

   //###### the nominal fit with default pdf 

   // can remove the double crystal ball in pbpb: just commenting out and copying an appropriate version
   RooAbsPdf  *pdf             = new RooAddPdf ("pdf","total p.d.f.",
         RooArgList(*sig1S,*sig2S,*sig3S,*pdf_combinedbkgd),
         RooArgList(*nsig1f,*nsig2f,*nsig3f,*nbkgd));
   //  nsig3f->setVal(0); nsig3f->setConstant();

   w.import(*pdf);
   w.Print();
}
