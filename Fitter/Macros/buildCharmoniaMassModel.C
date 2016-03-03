#include "Utilities/initClasses.h"

void fixPsi2StoJpsi(map<string, string>& parIni, bool isPbPb);
void fixPbPbtoPP(map<string, string>& parIni);
void setDefaultParameters(map<string, string> &parIni, bool isPbPb);
bool addSignalMassModel(RooWorkspace& ws, string object, MassModel model, map<string,string> parIni, bool isPbPb); 
bool addBackgroundMassModel(RooWorkspace& ws, string object, MassModel model, map<string,string> parIni, bool isPbPb);
bool defineCtauResolModel(RooWorkspace& ws, CtauModel model, map<string,string> parIni, bool isPbPb); 
bool addSignalCtauModel(RooWorkspace& ws, string object, CtauModel model, map<string,string> parIni, bool isPbPb); 
bool addBackgroundCtauModel(RooWorkspace& ws, string object, CtauModel model, map<string,string> parIni, bool isPbPb);


bool buildCharmoniaMassModel(RooWorkspace& ws, struct InputOpt opt, struct CharmModel model, map<string, string>  parIni, bool isPbPb)
{

  // If the initial parameters are empty, set defaul parameter values
  if (parIni.size()==0) { setDefaultParameters(parIni, isPbPb); }
  // Fix all psi2S parameters to jpsi
  fixPsi2StoJpsi(parIni, isPbPb);

  if (opt.inExcStat) {
    if (opt.doSimulFit && isPbPb) {
      // Fix mean, alpha and n parameters in PbPb to PP values 
      fixPbPbtoPP(parIni);
      // Asign NPsi2S(PbPb) = DoubleRatio(PbPb/PP) * SingleRatioPP(Psi2S/Jpsi) * NJpsi(PbPb) 
      ws.factory(Form("RooFormulaVar::%s('@0*@1*@2',{%s,%s,%s})", "N_Psi2S_PbPb", 
		      parIni["RFrac2Svs1S_PbPbvsPP"].c_str(), 
		      "RFrac2Svs1S_PP", 
		      parIni[Form("N_Jpsi_%s", (isPbPb?"PbPb":"PP"))].c_str() 
		      )); 
      ws.var("RFrac2Svs1S_PbPbvsPP")->setConstant(kFALSE);
    } else {

      // Asign N(Psi2S) = SingleRatio(Psi2S/Jpsi) * N(Jpsi)
      ws.factory(Form("RooFormulaVar::%s('@0*@1',{%s,%s})", Form("N_Psi2S_%s", (isPbPb?"PbPb":"PP")), 
		      parIni[Form("RFrac2Svs1S_%s", (isPbPb?"PbPb":"PP"))].c_str(), 
		      parIni[Form("N_Jpsi_%s", (isPbPb?"PbPb":"PP"))].c_str() 
		      )); 
      ws.var(Form("RFrac2Svs1S_%s", (isPbPb?"PbPb":"PP")))->setConstant(kFALSE);
    }
    parIni[Form("N_Jpsi_%s", (isPbPb?"PbPb":"PP"))] = Form("N_Jpsi_%s", (isPbPb?"PbPb":"PP"));
    parIni[Form("N_Psi2S_%s", (isPbPb?"PbPb":"PP"))] = Form("N_Psi2S_%s", (isPbPb?"PbPb":"PP")); 
  }

  // C r e a t e   m o d e l  

  if(!addSignalMassModel(ws, "Jpsi", model.Jpsi.Mass, parIni, isPbPb)) { cout << "[ERROR] Adding Jpsi Mass Model failed" << endl; return false; }
  if (opt.inExcStat) { 
    if (!addSignalMassModel(ws, "Psi2S", model.Psi2S.Mass, parIni, isPbPb)) { cout << "[ERROR] Adding Psi(2S) Mass Model failed" << endl; return false; }  
  }
  if(!addBackgroundMassModel(ws, "Bkg", model.Bkg.Mass, parIni, isPbPb)) { cout << "[ERROR] Adding Background Mass Model failed" << endl; return false; }
    
  // Total PDF = Signal + Background
  if (opt.inExcStat) {
    ws.factory(Form("SUM::%s(%s*%s, %s*%s, %s*%s)", Form("pdfMASS_Tot_%s", (isPbPb?"PbPb":"PP")),
		    parIni[Form("N_Jpsi_%s", (isPbPb?"PbPb":"PP"))].c_str(),
		    Form("pdfMASS_Jpsi_%s", (isPbPb?"PbPb":"PP")),
		    parIni[Form("N_Psi2S_%s", (isPbPb?"PbPb":"PP"))].c_str(),
		    Form("pdfMASS_Psi2S_%s", (isPbPb?"PbPb":"PP")),
		    parIni[Form("N_Bkg_%s", (isPbPb?"PbPb":"PP"))].c_str(),
		    Form("pdfMASS_Bkg_%s", (isPbPb?"PbPb":"PP"))
		    ));
  } else {
    ws.factory(Form("SUM::%s(%s*%s, %s*%s)", Form("pdfMASS_Tot_%s", (isPbPb?"PbPb":"PP")),
		    parIni[Form("N_Jpsi_%s", (isPbPb?"PbPb":"PP"))].c_str(),
		    Form("pdfMASS_Jpsi_%s", (isPbPb?"PbPb":"PP")),
		    parIni[Form("N_Bkg_%s", (isPbPb?"PbPb":"PP"))].c_str(),
		    Form("pdfMASS_Bkg_%s", (isPbPb?"PbPb":"PP"))
		    ));
  }
  ws.pdf(Form("pdfMASS_Tot_%s", (isPbPb?"PbPb":"PP")))->setNormRange("MassWindow");

  // save the initial values of the model we've just created
  RooAbsPdf *themodel = ws.pdf(Form("pdfMASS_Tot_%s", (isPbPb?"PbPb":"PP")));
  RooRealVar *x = ws.var("invMass");
  RooArgSet* params = (RooArgSet*) themodel->getParameters(*x) ;
  ws.saveSnapshot(Form("pdfMASS_Tot_%s_parIni", (isPbPb?"PbPb":"PP")),*params,kTRUE) ;
  
  //ws.Print();
  return true;
};

bool addBackgroundMassModel(RooWorkspace& ws, string object, MassModel model, map<string,string> parIni, bool isPbPb) 
{
  cout << Form("[INFO] Implementing %s Background Mass Model", object.c_str()) << endl;
  
  switch(model) 
    {  
    case (MassModel::FirstOrderPolynomial): 
      if (!( parIni.count(Form("lambda1_%s_%s", object.c_str(), (isPbPb?"PbPb":"PP"))) )) { 
	cout << Form("[ERROR] Initial parameters where not found for %s Background First Order Polynomial in %s", object.c_str(), (isPbPb?"PbPb":"PP")) << endl; return false;
      } 
      ws.factory(Form("Polynomial::%s(%s, {%s})", Form("pdfMASS_%s_%s", object.c_str(), (isPbPb?"PbPb":"PP")), "invMass", 
		      parIni[Form("lambda1_%s_%s", object.c_str(), (isPbPb?"PbPb":"PP"))].c_str()
		      ));
      cout << Form("[INFO] %s Background 1st Order Polynomial PDF in %s included", object.c_str(), (isPbPb?"PbPb":"PP")) << endl; break;
 
    case (MassModel::SecondOrderPolynomial): 
      if (!( parIni.count(Form("lambda1_%s_%s", object.c_str(), (isPbPb?"PbPb":"PP"))) && parIni.count(Form("lambda2_%s_%s", object.c_str(), (isPbPb?"PbPb":"PP"))) )) { 
	cout << Form("[ERROR] Initial parameters where not found for %s Background Second Order Polynomial in %s", object.c_str(), (isPbPb?"PbPb":"PP")) << endl; return false;
      } 
      ws.factory(Form("Polynomial::%s(%s, {%s, %s})", Form("pdfMASS_%s_%s", object.c_str(), (isPbPb?"PbPb":"PP")), "invMass", 
		      parIni[Form("lambda1_%s_%s", object.c_str(), (isPbPb?"PbPb":"PP"))].c_str(), 
		      parIni[Form("lambda2_%s_%s", object.c_str(), (isPbPb?"PbPb":"PP"))].c_str()
		      ));
      cout << Form("[INFO] %s Background 2nd Order Polynomial PDF in %s included", object.c_str(), (isPbPb?"PbPb":"PP")) << endl; break; 

    case (MassModel::ThirdOrderPolynomial): 
      if (!( parIni.count(Form("lambda1_%s_%s", object.c_str(), (isPbPb?"PbPb":"PP"))) && parIni.count(Form("lambda2_%s_%s", object.c_str(), (isPbPb?"PbPb":"PP")))
	     && parIni.count(Form("lambda3_%s_%s", object.c_str(), (isPbPb?"PbPb":"PP"))) )) { 
	cout << Form("[ERROR] Initial parameters where not found for %s Background Third Order Polynomial in %s", object.c_str(), (isPbPb?"PbPb":"PP")) << endl; return false;
      } 
      ws.factory(Form("Polynomial::%s(%s, {%s, %s, %s})", Form("pdfMASS_%s_%s", object.c_str(), (isPbPb?"PbPb":"PP")), "invMass", 
		      parIni[Form("lambda1_%s_%s", object.c_str(), (isPbPb?"PbPb":"PP"))].c_str(), 
		      parIni[Form("lambda2_%s_%s", object.c_str(), (isPbPb?"PbPb":"PP"))].c_str(), 
		      parIni[Form("lambda3_%s_%s", object.c_str(), (isPbPb?"PbPb":"PP"))].c_str()
		      ));
      cout << Form("[INFO] %s Background 3rd Order Polynomial PDF in %s included", object.c_str(), (isPbPb?"PbPb":"PP")) << endl; break; 

    case (MassModel::FourthOrderPolynomial): 
      if (!( parIni.count(Form("lambda1_%s_%s", object.c_str(), (isPbPb?"PbPb":"PP"))) && parIni.count(Form("lambda2_%s_%s", object.c_str(), (isPbPb?"PbPb":"PP"))) 
	     && parIni.count(Form("lambda3_%s_%s", object.c_str(), (isPbPb?"PbPb":"PP"))) && parIni.count(Form("lambda4_%s_%s", object.c_str(), (isPbPb?"PbPb":"PP"))) )) { 
	cout << Form("[ERROR] Initial parameters where not found for %s Background Fourth Order Polynomial in %s", object.c_str(), (isPbPb?"PbPb":"PP")) << endl; return false;
      } 
      ws.factory(Form("Polynomial::%s(%s, {%s, %s, %s, %s})", Form("pdfMASS_%s_%s", object.c_str(), (isPbPb?"PbPb":"PP")), "invMass", 
		      parIni[Form("lambda1_%s_%s", object.c_str(), (isPbPb?"PbPb":"PP"))].c_str(), 
		      parIni[Form("lambda2_%s_%s", object.c_str(), (isPbPb?"PbPb":"PP"))].c_str(), 
		      parIni[Form("lambda3_%s_%s", object.c_str(), (isPbPb?"PbPb":"PP"))].c_str(), 
		      parIni[Form("lambda4_%s_%s", object.c_str(), (isPbPb?"PbPb":"PP"))].c_str()
		      ));
      cout << Form("[INFO] %s Background 4th Order Polynomial PDF in %s included", object.c_str(), (isPbPb?"PbPb":"PP")) << endl; break; 

    case (MassModel::FirstOrderChebychev): 
      if (!( parIni.count(Form("lambda1_%s_%s", object.c_str(), (isPbPb?"PbPb":"PP"))) )) { 
	cout << Form("[ERROR] Initial parameters where not found for %s Background First Order Chebychev in %s", object.c_str(), (isPbPb?"PbPb":"PP")) << endl; return false;
      } 
      ws.factory(Form("Chebychev::%s(%s, {%s})", Form("pdfMASS_%s_%s", object.c_str(), (isPbPb?"PbPb":"PP")), "invMass", 
		      parIni[Form("lambda1_%s_%s", object.c_str(), (isPbPb?"PbPb":"PP"))].c_str()
		      ));
      cout << Form("[INFO] %s Background 1st Order Chebychev PDF in %s included", object.c_str(), (isPbPb?"PbPb":"PP")) << endl; break;
 
    case (MassModel::SecondOrderChebychev): 
      if (!( parIni.count(Form("lambda1_%s_%s", object.c_str(), (isPbPb?"PbPb":"PP"))) && parIni.count(Form("lambda2_%s_%s", object.c_str(), (isPbPb?"PbPb":"PP"))) )) { 
	cout << Form("[ERROR] Initial parameters where not found for %s Background Second Order Chebychev in %s", object.c_str(), (isPbPb?"PbPb":"PP")) << endl; return false;
      } 
      ws.factory(Form("Chebychev::%s(%s, {%s, %s})", Form("pdfMASS_%s_%s", object.c_str(), (isPbPb?"PbPb":"PP")), "invMass", 
		      parIni[Form("lambda1_%s_%s", object.c_str(), (isPbPb?"PbPb":"PP"))].c_str(), 
		      parIni[Form("lambda2_%s_%s", object.c_str(), (isPbPb?"PbPb":"PP"))].c_str()
		      ));
      cout << Form("[INFO] %s Background 2nd Order Chebychev PDF in %s included", object.c_str(), (isPbPb?"PbPb":"PP")) << endl; break; 

    case (MassModel::ThirdOrderChebychev): 
      if (!( parIni.count(Form("lambda1_%s_%s", object.c_str(), (isPbPb?"PbPb":"PP"))) && parIni.count(Form("lambda2_%s_%s", object.c_str(), (isPbPb?"PbPb":"PP")))
	     && parIni.count(Form("lambda3_%s_%s", object.c_str(), (isPbPb?"PbPb":"PP"))) )) { 
	cout << Form("[ERROR] Initial parameters where not found for %s Background Third Order Chebychev in %s", object.c_str(), (isPbPb?"PbPb":"PP")) << endl; return false;
      } 
      ws.factory(Form("Chebychev::%s(%s, {%s, %s, %s})", Form("pdfMASS_%s_%s", object.c_str(), (isPbPb?"PbPb":"PP")), "invMass", 
		      parIni[Form("lambda1_%s_%s", object.c_str(), (isPbPb?"PbPb":"PP"))].c_str(), 
		      parIni[Form("lambda2_%s_%s", object.c_str(), (isPbPb?"PbPb":"PP"))].c_str(), 
		      parIni[Form("lambda3_%s_%s", object.c_str(), (isPbPb?"PbPb":"PP"))].c_str()
		      ));
      cout << Form("[INFO] %s Background 3rd Order Polynomial PDF in %s included", object.c_str(), (isPbPb?"PbPb":"PP")) << endl; break; 

    case (MassModel::FourthOrderChebychev): 
      if (!( parIni.count(Form("lambda1_%s_%s", object.c_str(), (isPbPb?"PbPb":"PP"))) && parIni.count(Form("lambda2_%s_%s", object.c_str(), (isPbPb?"PbPb":"PP"))) 
	     && parIni.count(Form("lambda3_%s_%s", object.c_str(), (isPbPb?"PbPb":"PP"))) && parIni.count(Form("lambda4_%s_%s", object.c_str(), (isPbPb?"PbPb":"PP"))) )) { 
	cout << Form("[ERROR] Initial parameters where not found for %s Background Fourth Order Chebychev in %s", object.c_str(), (isPbPb?"PbPb":"PP")) << endl; return false;
      } 
      ws.factory(Form("Polynomial::%s(%s, {%s, %s, %s, %s})", Form("pdfMASS_%s_%s", object.c_str(), (isPbPb?"PbPb":"PP")), "invMass", 
		      parIni[Form("lambda1_%s_%s", object.c_str(), (isPbPb?"PbPb":"PP"))].c_str(), 
		      parIni[Form("lambda2_%s_%s", object.c_str(), (isPbPb?"PbPb":"PP"))].c_str(), 
		      parIni[Form("lambda3_%s_%s", object.c_str(), (isPbPb?"PbPb":"PP"))].c_str(), 
		      parIni[Form("lambda4_%s_%s", object.c_str(), (isPbPb?"PbPb":"PP"))].c_str()
		      ));
      cout << Form("[INFO] %s Background 4th Order Chebychev PDF in %s included", object.c_str(), (isPbPb?"PbPb":"PP")) << endl; break; 

    case (MassModel::Exponential): 
      if (!( parIni.count(Form("lambda1_%s_%s", object.c_str(), (isPbPb?"PbPb":"PP"))) )) { 
	cout << Form("[ERROR] Initial parameters where not found for %s Background Exponential in %s", object.c_str(), (isPbPb?"PbPb":"PP")) << endl; return false;
      } 
      ws.factory(Form("Exponential::%s(%s, %s)", Form("pdfMASS_%s_%s", object.c_str(), (isPbPb?"PbPb":"PP")), "invMass", 
		      parIni[Form("lambda1_%s_%s", object.c_str(), (isPbPb?"PbPb":"PP"))].c_str()
		      ));
      cout << Form("[INFO] %s Background Exponential PDF in %s included", object.c_str(), (isPbPb?"PbPb":"PP")) << endl; break;
      
    default :
      cout<< Form("[ERROR] Selected Background Mass Model for %s has not been implemented", object.c_str()) << endl; return false;
    }
  
  return true;
};


bool addSignalMassModel(RooWorkspace& ws, string object, MassModel model, map<string,string> parIni, bool isPbPb) 
{
  cout << Form("[INFO] Implementing %s Mass Model", object.c_str()) << endl;
  
  switch(model) 
    {    
    case (MassModel::SingleGaussian): 
      if (!( parIni.count(Form("m_%s_%s", object.c_str(), (isPbPb?"PbPb":"PP"))) && parIni.count(Form("sigma1_%s_%s", object.c_str(), (isPbPb?"PbPb":"PP"))) )) { 
	cout << Form("[ERROR] Initial parameters where not found for %s Single Gaussian Model in %s", object.c_str(), (isPbPb?"PbPb":"PP")) << endl; return false;
      }
      ws.factory(Form("Gaussian::%s(%s, %s, %s)", Form("pdfMASS_%s_%s", object.c_str(), (isPbPb?"PbPb":"PP")), "invMass", 
		      parIni[Form("m_%s_%s", object.c_str(), (isPbPb?"PbPb":"PP"))].c_str(), 
		      parIni[Form("sigma1_%s_%s", object.c_str(), (isPbPb?"PbPb":"PP"))].c_str()
		      ));
    
      cout << Form("[INFO] %s Single Gaussian PDF in %s included", object.c_str(), (isPbPb?"PbPb":"PP")) << endl; break;  
      
    case (MassModel::DoubleGaussian): 
      if (!( parIni.count(Form("m_%s_%s", object.c_str(), (isPbPb?"PbPb":"PP"))) && parIni.count(Form("sigma1_%s_%s", object.c_str(), (isPbPb?"PbPb":"PP"))) 
	     && parIni.count(Form("sigma2_%s_%s", object.c_str(), (isPbPb?"PbPb":"PP"))) && parIni.count(Form("f_%s_%s", object.c_str(), (isPbPb?"PbPb":"PP"))) )) { 
	cout << Form("[ERROR] Initial parameters where not found for %s Double Gaussian Model in %s", object.c_str(), (isPbPb?"PbPb":"PP")) << endl; return false; 
      }
      ws.factory(Form("SUM::%s(%s*%s, %s)", Form("pdfMASS_Jpsi_%s", (isPbPb?"PbPb":"PP")),
		      parIni[Form("f_Jpsi_%s", (isPbPb?"PbPb":"PP"))].c_str(),
		      Form("Gaussian::%s(%s, %s, %s)", Form("pdfMASSG1_%s_%s", object.c_str(), (isPbPb?"PbPb":"PP")), "invMass", 
			   parIni[Form("m_%s_%s", object.c_str(), (isPbPb?"PbPb":"PP"))].c_str(), 
			   parIni[Form("sigma1_%s_%s", object.c_str(), (isPbPb?"PbPb":"PP"))].c_str()
			   ),
		      Form("Gaussian::%s(%s, %s, %s)", Form("pdfMASSG2_%s_%s", object.c_str(), (isPbPb?"PbPb":"PP")), "invMass", 
			   Form("m_%s_%s", object.c_str(), (isPbPb?"PbPb":"PP")), 
			   parIni[Form("sigma2_%s_%s", object.c_str(), (isPbPb?"PbPb":"PP"))].c_str()
			   )
		      ));
      cout << Form("[INFO] %s Double Gaussian PDF in %s included", object.c_str(), (isPbPb?"PbPb":"PP")) << endl; break; 

    case (MassModel::SingleCrystalBall):  
      if (!( parIni.count(Form("m_%s_%s", object.c_str(), (isPbPb?"PbPb":"PP"))) && parIni.count(Form("sigma1_%s_%s", object.c_str(), (isPbPb?"PbPb":"PP")))
	     && parIni.count(Form("alpha_%s_%s", object.c_str(), (isPbPb?"PbPb":"PP"))) && parIni.count(Form("n_%s_%s", object.c_str(), (isPbPb?"PbPb":"PP"))) )) {
	cout << Form("[ERROR] Initial parameters where not found for %s Single Crystal Ball Model in %s", object.c_str(), (isPbPb?"PbPb":"PP")) << endl; return false; 
      }      
      ws.factory(Form("CBShape::%s(%s, %s, %s, %s, %s)", Form("pdfMASS_%s_%s", object.c_str(), (isPbPb?"PbPb":"PP")), "invMass", 
		      parIni[Form("m_%s_%s", object.c_str(), (isPbPb?"PbPb":"PP"))].c_str(), 
		      parIni[Form("sigma1_%s_%s", object.c_str(), (isPbPb?"PbPb":"PP"))].c_str(),
		      parIni[Form("alpha_%s_%s", object.c_str(), (isPbPb?"PbPb":"PP"))].c_str(),
		      parIni[Form("n_%s_%s", object.c_str(), (isPbPb?"PbPb":"PP"))].c_str()
		      ));
      cout << Form("[INFO] %s Single Crystal Ball PDF in %s included", object.c_str(), (isPbPb?"PbPb":"PP")) << endl; break;
      
    case (MassModel::DoubleCrystalBall): 
      if (!( parIni.count(Form("m_%s_%s", object.c_str(), (isPbPb?"PbPb":"PP"))) && parIni.count(Form("sigma1_%s_%s", object.c_str(), (isPbPb?"PbPb":"PP")))
	     && parIni.count(Form("sigma2_%s_%s", object.c_str(), (isPbPb?"PbPb":"PP"))) && parIni.count(Form("alpha_%s_%s", object.c_str(), (isPbPb?"PbPb":"PP"))) 
	     && parIni.count(Form("n_%s_%s", object.c_str(), (isPbPb?"PbPb":"PP"))) && parIni.count(Form("f_%s_%s", object.c_str(), (isPbPb?"PbPb":"PP"))) )) { 
	cout << Form("[ERROR] Initial parameters where not found for %s Double Crystal Ball Model in %s", object.c_str(), (isPbPb?"PbPb":"PP")) << endl; return false; 
      }  
      ws.factory(Form("SUM::%s(%s*%s, %s)", Form("pdfMASS_%s_%s", object.c_str(), (isPbPb?"PbPb":"PP")),
		      parIni[Form("f_%s_%s", object.c_str(), (isPbPb?"PbPb":"PP"))].c_str(),
		      Form("CBShape::%s(%s, %s, %s, %s, %s)", Form("pdfMASSCB1_%s_%s", object.c_str(), (isPbPb?"PbPb":"PP")), "invMass", 
			   parIni[Form("m_%s_%s", object.c_str(), (isPbPb?"PbPb":"PP"))].c_str(), 
			   parIni[Form("sigma1_%s_%s", object.c_str(), (isPbPb?"PbPb":"PP"))].c_str(),
			   parIni[Form("alpha_%s_%s", object.c_str(), (isPbPb?"PbPb":"PP"))].c_str(),
			   parIni[Form("n_%s_%s", object.c_str(), (isPbPb?"PbPb":"PP"))].c_str()
			   ),
		      Form("CBShape::%s(%s, %s, %s, %s, %s)", Form("pdfMASSCB2_%s_%s", object.c_str(), (isPbPb?"PbPb":"PP")), "invMass", 
			   Form("m_%s_%s", object.c_str(), (isPbPb?"PbPb":"PP")), 
			   parIni[Form("sigma2_%s_%s", object.c_str(), (isPbPb?"PbPb":"PP"))].c_str(),
			   Form("alpha_%s_%s", object.c_str(), (isPbPb?"PbPb":"PP")),
			   Form("n_%s_%s", object.c_str(), (isPbPb?"PbPb":"PP"))
			   )
		      ));   
      cout << Form("[INFO] %s Double Crystal Ball PDF in %s included", object.c_str(), (isPbPb?"PbPb":"PP")) << endl; break;
      
    case (MassModel::GaussianAndCrystalBall):
      if (!( parIni.count(Form("m_%s_%s", object.c_str(), (isPbPb?"PbPb":"PP"))) && parIni.count(Form("sigma1_%s_%s", object.c_str(), (isPbPb?"PbPb":"PP")))
	     && parIni.count(Form("sigma2_%s_%s", object.c_str(), (isPbPb?"PbPb":"PP"))) && parIni.count(Form("alpha_%s_%s", object.c_str(), (isPbPb?"PbPb":"PP"))) 
	     && parIni.count(Form("n_%s_%s", object.c_str(), (isPbPb?"PbPb":"PP"))) && parIni.count(Form("f_%s_%s", object.c_str(), (isPbPb?"PbPb":"PP")))  )) { 
	cout << Form("[ERROR] Initial parameters where not found for %s Gaussian and Crystal Ball Model in %s", object.c_str(), (isPbPb?"PbPb":"PP")) << endl; return false;
      }  
      ws.factory(Form("SUM::%s(%s*%s, %s)", Form("pdfMASS_%s_%s", object.c_str(), (isPbPb?"PbPb":"PP")),
		      parIni[Form("f%s%s", object.c_str(), (isPbPb?"PbPb":"PP"))].c_str(),
		      Form("Gaussian::%s(%s, %s, %s)", Form("pdfMASSG1_%s_%s", object.c_str(), (isPbPb?"PbPb":"PP")), "invMass", 
			   parIni[Form("m_%s_%s", object.c_str(), (isPbPb?"PbPb":"PP"))].c_str(), 
			   parIni[Form("sigma1_%s_%s", object.c_str(), (isPbPb?"PbPb":"PP"))].c_str()
			   ),
		      Form("CBShape::%s(%s, %s, %s, %s, %s)", Form("pdfMASSCB1_%s_%s", object.c_str(), (isPbPb?"PbPb":"PP")), "invMass", 
			   Form("m_%s_%s", object.c_str(), (isPbPb?"PbPb":"PP")), 
			   parIni[Form("sigma2_%s_%s", object.c_str(), (isPbPb?"PbPb":"PP"))].c_str(),
			   Form("alpha_%s_%s", object.c_str(), (isPbPb?"PbPb":"PP")),
			   Form("n_%s_%s", object.c_str(), (isPbPb?"PbPb":"PP"))
			   )
		      ));   
      cout << Form("[INFO] %s Gaussian and Crystal Ball PDF in %s included", object.c_str(), (isPbPb?"PbPb":"PP")) << endl; break;
	
    default :
      cout<< "[ERROR] Selected Signal Mass Model has not been implemented"<< endl; return false;
    }
  
  return true;
};

void fixPbPbtoPP(map<string, string>& parIni)
{
  parIni["m_Jpsi_PbPb"]  = Form("RooFormulaVar::%s('@0',{%s})", "m_Jpsi_PbPb", "m_Jpsi_PP");
  parIni["m_Psi2S_PbPb"] = Form("RooFormulaVar::%s('@0',{%s})", "m_Psi2S_PbPb", "m_Psi2S_PP");
  //parIni["sigma1_Jpsi_PbPb"]  = Form("RooFormulaVar::%s('@0',{%s})", "sigma1_Jpsi_PbPb", "sigma1_JpsiPP");
  //parIni["sigma1_Psi2S_PbPb"] = Form("RooFormulaVar::%s('@0',{%s})", "sigma1_Psi2S_PbPb", "sigma1_Psi2SPP");
  //parIni["sigma2_Jpsi_PbPb"]  = Form("RooFormulaVar::%s('@0',{%s})", "sigma2_Jpsi_PbPb", "sigma2_JpsiPP");
  //parIni["sigma2_Psi2S_PbPb"] = Form("RooFormulaVar::%s('@0',{%s})", "sigma2_Jpsi_PbPb", "sigma2_Psi2SPP");
  parIni["alpha_Jpsi_PbPb"]  = Form("RooFormulaVar::%s('@0',{%s})", "alpha_Jpsi_PbPb", "alpha_Jpsi_PP");
  parIni["alpha_Psi2S_PbPb"] = Form("RooFormulaVar::%s('@0',{%s})", "alpha_Psi2S_PbPb", "alpha_Jpsi_PP");
  parIni["n_Jpsi_PbPb"] = Form("RooFormulaVar::%s('@0',{%s})", "n_Jpsi_PbPb", "n_Jpsi_PP");
  parIni["n_Psi2S_PbPb"] = Form("RooFormulaVar::%s('@0',{%s})", "n_Psi2S_PbPb", "n_Jpsi_PP");
  //parIni["f_Jpsi_PbPb"] = Form("RooFormulaVar::%s('@0',{%s})", "f_Jpsi_PbPb", "f_Jpsi_PP");
  //parIni["f_Psi2S_PbPb"] = Form("RooFormulaVar::%s('@0',{%s})", "f_Psi2S_PbPb", "f_Jpsi_PP");  
};

void fixPsi2StoJpsi(map<string, string>& parIni, bool isPbPb)
{
  cout << "[INFO] Constraining Psi(2S) parameters to Jpsi using PDF Mass Ratio" << endl;
  Double_t MassRatio = (Mass.Psi2S/Mass.JPsi);
  parIni[Form("m_Psi2S_%s", (isPbPb?"PbPb":"PP"))]      = Form("RooFormulaVar::%s('@0*@1',{MassRatio[%.4f],%s})", Form("m_Psi2S_%s", (isPbPb?"PbPb":"PP")), MassRatio, Form("m_Jpsi_%s", (isPbPb?"PbPb":"PP") ));
  parIni[Form("sigma1_Psi2S_%s", (isPbPb?"PbPb":"PP"))] = Form("RooFormulaVar::%s('@0*@1',{MassRatio,%s})", Form("sigma1_Psi2S_%s", (isPbPb?"PbPb":"PP")), Form("sigma1_Jpsi_%s", (isPbPb?"PbPb":"PP") ));
  parIni[Form("sigma2_Psi2S_%s", (isPbPb?"PbPb":"PP"))] = Form("RooFormulaVar::%s('@0*@1',{MassRatio,%s})", Form("sigma2_Psi2S_%s", (isPbPb?"PbPb":"PP")), Form("sigma2_Jpsi_%s", (isPbPb?"PbPb":"PP") ));
  parIni[Form("alpha_Psi2S_%s", (isPbPb?"PbPb":"PP"))]  = Form("RooFormulaVar::%s('@0',{%s})", Form("alpha_Psi2S_%s", (isPbPb?"PbPb":"PP")), Form("alpha_Jpsi_%s", (isPbPb?"PbPb":"PP")));
  parIni[Form("n_Psi2S_%s", (isPbPb?"PbPb":"PP"))]      = Form("RooFormulaVar::%s('@0',{%s})", Form("n_Psi2S_%s", (isPbPb?"PbPb":"PP")), Form("n_Jpsi_%s", (isPbPb?"PbPb":"PP")));
  parIni[Form("f_Psi2S_%s", (isPbPb?"PbPb":"PP"))]      = Form("RooFormulaVar::%s('@0',{%s})", Form("f_Psi2S_%s", (isPbPb?"PbPb":"PP")), Form("f_Jpsi_%s", (isPbPb?"PbPb":"PP")));
};

void setDefaultParameters(map<string, string> &parIni, bool isPbPb)
{
  
  // MASS FIT PARAMETERS, TO BE REMOVED SOON
  parIni[Form("N_Jpsi_%s", (isPbPb?"PbPb":"PP"))]  = Form("%s[%d,%d]", Form("N_Jpsi_%s", (isPbPb?"PbPb":"PP")), 0, 350000);
  parIni[Form("N_Psi2S_%s", (isPbPb?"PbPb":"PP"))] = Form("%s[%d,%d]", Form("N_Psi2S_%s", (isPbPb?"PbPb":"PP")), 0, 350000);
  parIni[Form("N_Bkg_%s", (isPbPb?"PbPb":"PP"))]   = Form("%s[%d,%d]", Form("N_Bkg_%s", (isPbPb?"PbPb":"PP")), 0, 350000);

  parIni[Form("m_Jpsi_%s", (isPbPb?"PbPb":"PP"))]      = Form("%s[%.4f,%.4f,%.4f]", Form("m_Jpsi_%s", (isPbPb?"PbPb":"PP")), Mass.JPsi, Mass.JPsi-0.1, Mass.JPsi+0.1);
  parIni[Form("sigma1_Jpsi_%s", (isPbPb?"PbPb":"PP"))] = Form("%s[%.4f,%.4f,%.4f]", Form("sigma1_Jpsi_%s", (isPbPb?"PbPb":"PP")), 0.04, 0.01, 0.09);
  parIni[Form("sigma2_Jpsi_%s", (isPbPb?"PbPb":"PP"))] = Form("%s[%.4f,%.4f,%.4f]", Form("sigma2_Jpsi_%s", (isPbPb?"PbPb":"PP")), 0.02, 0.01, 0.07);
  parIni[Form("alpha_Jpsi_%s", (isPbPb?"PbPb":"PP"))]  = Form("%s[%.4f,%.4f,%.4f]", Form("alpha_Jpsi_%s", (isPbPb?"PbPb":"PP")), 2.21, 0.5, 30.0*100.0);
  parIni[Form("n_Jpsi_%s", (isPbPb?"PbPb":"PP"))]      = Form("%s[%.4f,%.4f,%.4f]", Form("n_Jpsi_%s", (isPbPb?"PbPb":"PP")), 2., 0.5, 30.0*100.0);
  parIni[Form("f_Jpsi_%s", (isPbPb?"PbPb":"PP"))]      = Form("%s[%.4f,%.4f,%.4f]", Form("f_Jpsi_%s", (isPbPb?"PbPb":"PP")), 0.3, 0.0, 1.0);

  fixPsi2StoJpsi(parIni, isPbPb);

  parIni["RFrac2Svs1S_PbPbvsPP"] = Form("%s[%.4f,%.4f,%.4f]", "RFrac2Svs1S_PbPbvsPP", 0.26, 0.0, 3.0);
  parIni[Form("RFrac2Svs1S_%s", (isPbPb?"PbPb":"PP"))] = Form("%s[%.4f,%.4f,%.4f]", Form("RFrac2Svs1S_%s", (isPbPb?"PbPb":"PP")), 0.26, 0.0, 1.0);

  parIni[Form("lambda1_Bkg_%s", (isPbPb?"PbPb":"PP"))] = Form("%s[%.4f,%.4f,%.4f]", Form("lambda1_Bkg_%s", (isPbPb?"PbPb":"PP")), 0.05, -1.0, 1.0);
  parIni[Form("lambda2_Bkg_%s", (isPbPb?"PbPb":"PP"))] = Form("%s[%.4f,%.4f,%.4f]", Form("lambda2_Bkg_%s", (isPbPb?"PbPb":"PP")), 0.05, -1.0, 1.0);
  parIni[Form("lambda3_Bkg_%s", (isPbPb?"PbPb":"PP"))] = Form("%s[%.4f,%.4f,%.4f]", Form("lambda3_Bkg_%s", (isPbPb?"PbPb":"PP")), 0.05, -1.0, 1.0);
  parIni[Form("lambda4_Bkg_%s", (isPbPb?"PbPb":"PP"))] = Form("%s[%.4f,%.4f,%.4f]", Form("lambda4_Bkg_%s", (isPbPb?"PbPb":"PP")), 0.05, -1.0, 1.0);
 
};
