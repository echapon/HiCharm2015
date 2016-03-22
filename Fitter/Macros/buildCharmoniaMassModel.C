#include "Utilities/initClasses.h"

void fixPsi2StoJpsi(map<string, string>& parIni, bool isPbPb);
void fixPbPbtoPP(map<string, string>& parIni);
void setDefaultParameters(map<string, string> &parIni, bool isPbPb, double numEntries);
bool addSignalMassModel(RooWorkspace& ws, string object, MassModel model, map<string,string> parIni, bool isPbPb); 
bool addBackgroundMassModel(RooWorkspace& ws, string object, MassModel model, map<string,string> parIni, bool isPbPb);


bool buildCharmoniaMassModel(RooWorkspace& ws, struct CharmModel model, map<string, string>  parIni, 
                             bool isPbPb,                 // Determine if we are working with PbPb (True) or PP (False)
                             bool doSimulFit,             // Do simultaneous fit
                             bool incBkg,                 // Include background model
                             bool incJpsi,                // Include Jpsi model
                             bool incPsi2S,               // Include Psi(2S) model
                             string label,                // pdf label
                             double  numEntries = 300000. // Number of entries in the dataset
                             )
{

  // If the initial parameters are empty, set defaul parameter values
  setDefaultParameters(parIni, isPbPb, numEntries);

  // Fix all psi2S parameters to jpsi
  if (incJpsi && incPsi2S) {
    fixPsi2StoJpsi(parIni, isPbPb);
  }

  // Let's define the single and double ratio variables
  if (incPsi2S && incJpsi) {
    if (doSimulFit && isPbPb) {

      // Fix mean, alpha and n parameters in PbPb to PP values 
      fixPbPbtoPP(parIni);

      // create parameters related to the double ratio
      ws.factory( parIni["RFrac2Svs1S_PbPbvsPP"].c_str() );                     // Double Ratio
      ws.factory( parIni[Form("N_Jpsi_%s", (isPbPb?"PbPb":"PP"))].c_str() );    // Number of Jpsi in PbPb

      // Asign NPsi2S(PbPb) = DoubleRatio(PbPb/PP) * SingleRatioPP(Psi2S/Jpsi) * NJpsi(PbPb) 
      ws.factory(Form("RooFormulaVar::%s('@0*@1*@2',{%s,%s,%s})", "N_Psi2S_PbPb", 
		      "RFrac2Svs1S_PbPbvsPP", 
		      "RFrac2Svs1S_PP", 
		      Form("N_Jpsi_%s", (isPbPb?"PbPb":"PP"))
		      )); 
      ws.var("RFrac2Svs1S_PbPbvsPP")->setConstant(kFALSE);

    } else {

      // create parameters related to the double ratio
      ws.factory( parIni[Form("RFrac2Svs1S_%s", (isPbPb?"PbPb":"PP"))].c_str() );
      ws.factory( parIni[Form("N_Jpsi_%s", (isPbPb?"PbPb":"PP"))].c_str() );

      // Asign N(Psi2S) = SingleRatio(Psi2S/Jpsi) * N(Jpsi)
      ws.factory(Form("RooFormulaVar::%s('@0*@1',{%s,%s})", Form("N_Psi2S_%s", (isPbPb?"PbPb":"PP")), 
		      Form("RFrac2Svs1S_%s", (isPbPb?"PbPb":"PP")),
		      Form("N_Jpsi_%s", (isPbPb?"PbPb":"PP")) 
		      )); 
      ws.var(Form("RFrac2Svs1S_%s", (isPbPb?"PbPb":"PP")))->setConstant(kFALSE);

    }
    
    // As we have already declared the Number of Psi2S and Jpsi before, just use their names in their PDFs
    parIni[Form("N_Jpsi_%s", (isPbPb?"PbPb":"PP"))] = Form("N_Jpsi_%s", (isPbPb?"PbPb":"PP"));
    parIni[Form("N_Psi2S_%s", (isPbPb?"PbPb":"PP"))] = Form("N_Psi2S_%s", (isPbPb?"PbPb":"PP")); 
  }

  // C r e a t e   m o d e l  

  if (incJpsi) {
    if(!addSignalMassModel(ws, "Jpsi", model.Jpsi.Mass, parIni, isPbPb)) { cout << "[ERROR] Adding Jpsi Mass Model failed" << endl; return false; }
  }
  if (incPsi2S) { 
    if (!addSignalMassModel(ws, "Psi2S", model.Psi2S.Mass, parIni, isPbPb)) { cout << "[ERROR] Adding Psi(2S) Mass Model failed" << endl; return false; }  
  }
  if (incBkg) {
    if(!addBackgroundMassModel(ws, "Bkg", model.Bkg.Mass, parIni, isPbPb)) { cout << "[ERROR] Adding Background Mass Model failed" << endl; return false; }
  }
  // Total PDF
  string pdfName = Form("pdfMASS_Tot_%s", (isPbPb?"PbPb":"PP"));
  if (!label.empty())pdfName+=Form("_%s",label.c_str());
  if (incJpsi && incPsi2S && incBkg) {
    ws.factory(Form("SUM::%s(%s*%s, %s*%s, %s*%s)", pdfName.c_str(),
		    parIni[Form("N_Jpsi_%s", (isPbPb?"PbPb":"PP"))].c_str(),
		    Form("pdfMASS_Jpsi_%s", (isPbPb?"PbPb":"PP")),
		    parIni[Form("N_Psi2S_%s", (isPbPb?"PbPb":"PP"))].c_str(),
		    Form("pdfMASS_Psi2S_%s", (isPbPb?"PbPb":"PP")),
		    parIni[Form("N_Bkg_%s", (isPbPb?"PbPb":"PP"))].c_str(),
		    Form("pdfMASS_Bkg_%s", (isPbPb?"PbPb":"PP"))
		    ));
  }
  if (incJpsi && incPsi2S && !incBkg) {
    ws.factory(Form("SUM::%s(%s*%s, %s*%s)", pdfName.c_str(),
		    parIni[Form("N_Jpsi_%s", (isPbPb?"PbPb":"PP"))].c_str(),
		    Form("pdfMASS_Jpsi_%s", (isPbPb?"PbPb":"PP")),
		    parIni[Form("N_Psi2S_%s", (isPbPb?"PbPb":"PP"))].c_str(),
		    Form("pdfMASS_Psi2S_%s", (isPbPb?"PbPb":"PP"))
		    ));
  }
  if (incJpsi && !incPsi2S && incBkg) {
    ws.factory(Form("SUM::%s(%s*%s, %s*%s)", pdfName.c_str(),
		    parIni[Form("N_Jpsi_%s", (isPbPb?"PbPb":"PP"))].c_str(),
		    Form("pdfMASS_Jpsi_%s", (isPbPb?"PbPb":"PP")),
		    parIni[Form("N_Bkg_%s", (isPbPb?"PbPb":"PP"))].c_str(),
		    Form("pdfMASS_Bkg_%s", (isPbPb?"PbPb":"PP"))
		    ));
  }
  if (!incJpsi && incPsi2S && incBkg) {
    ws.factory(Form("SUM::%s(%s*%s, %s*%s)", pdfName.c_str(),
		    parIni[Form("N_Psi2S_%s", (isPbPb?"PbPb":"PP"))].c_str(),
		    Form("pdfMASS_Psi2S_%s", (isPbPb?"PbPb":"PP")),
		    parIni[Form("N_Bkg_%s", (isPbPb?"PbPb":"PP"))].c_str(),
		    Form("pdfMASS_Bkg_%s", (isPbPb?"PbPb":"PP"))
		    ));
  }
  if (incJpsi && !incPsi2S && !incBkg) {
    ws.factory(Form("SUM::%s(%s*%s)", pdfName.c_str(),
		    parIni[Form("N_Jpsi_%s", (isPbPb?"PbPb":"PP"))].c_str(),
		    Form("pdfMASS_Jpsi_%s", (isPbPb?"PbPb":"PP"))
                    ));
  }
  if (!incJpsi && incPsi2S && !incBkg) {
    ws.factory(Form("SUM::%s(%s*%s)", pdfName.c_str(),
		    parIni[Form("N_Psi2S_%s", (isPbPb?"PbPb":"PP"))].c_str(),
		    Form("pdfMASS_Psi2S_%s", (isPbPb?"PbPb":"PP"))
		    ));
  }
  if (!incJpsi && !incPsi2S && incBkg) {
    ws.factory(Form("SUM::%s(%s*%s)", pdfName.c_str(),
		    parIni[Form("N_Bkg_%s", (isPbPb?"PbPb":"PP"))].c_str(),
		    Form("pdfMASS_Bkg_%s", (isPbPb?"PbPb":"PP"))
		    ));
  }
  if (!incJpsi && !incPsi2S && !incBkg) {
    cout << "[ERROR] User did not include any model, please fix your input settings!" << endl; return false;
  }
  ws.pdf(pdfName.c_str())->setNormRange("MassWindow");

  // save the initial values of the model we've just created
  RooAbsPdf *themodel = ws.pdf(pdfName.c_str());
  RooRealVar *x = ws.var("invMass");
  RooArgSet* params = (RooArgSet*) themodel->getParameters(*x) ;
  pdfName+="_parIni";
  ws.saveSnapshot(pdfName.c_str(),*params,kTRUE) ;
  
  //ws.Print();
  return true;
};

bool addBackgroundMassModel(RooWorkspace& ws, string object, MassModel model, map<string,string> parIni, bool isPbPb) 
{
  cout << Form("[INFO] Implementing %s Background Mass Model", object.c_str()) << endl;
  
  switch(model) 
    {  
    case (MassModel::FirstOrderChebychev): 

      // check that all input parameters are defined 
      if (!( 
            parIni.count(Form("lambda1_%s_%s", object.c_str(), (isPbPb?"PbPb":"PP"))) 
             )) { 
	cout << Form("[ERROR] Initial parameters where not found for %s Background First Order Chebychev in %s", object.c_str(), (isPbPb?"PbPb":"PP")) << endl; return false;
      } 

      // create the variables for this model
      ws.factory( parIni[Form("lambda1_%s_%s", object.c_str(), (isPbPb?"PbPb":"PP"))].c_str() );

      // create the PDF           
      ws.factory(Form("Chebychev::%s(%s, {%s})", Form("pdfMASS_%s_%s", object.c_str(), (isPbPb?"PbPb":"PP")), "invMass", 
		      Form("lambda1_%s_%s", object.c_str(), (isPbPb?"PbPb":"PP"))
		      ));

      cout << Form("[INFO] %s Background 1st Order Chebychev PDF in %s included", object.c_str(), (isPbPb?"PbPb":"PP")) << endl; break;
 
    case (MassModel::SecondOrderChebychev): 

      // check that all input parameters are defined 
      if (!( 
            parIni.count(Form("lambda1_%s_%s", object.c_str(), (isPbPb?"PbPb":"PP"))) && 
            parIni.count(Form("lambda2_%s_%s", object.c_str(), (isPbPb?"PbPb":"PP"))) 
             )) { 
	cout << Form("[ERROR] Initial parameters where not found for %s Background Second Order Chebychev in %s", object.c_str(), (isPbPb?"PbPb":"PP")) << endl; return false;
      } 

      // create the variables for this model
      ws.factory( parIni[Form("lambda1_%s_%s", object.c_str(), (isPbPb?"PbPb":"PP"))].c_str() );
      ws.factory( parIni[Form("lambda2_%s_%s", object.c_str(), (isPbPb?"PbPb":"PP"))].c_str() );

      // create the PDF           
      ws.factory(Form("Chebychev::%s(%s, {%s, %s})", Form("pdfMASS_%s_%s", object.c_str(), (isPbPb?"PbPb":"PP")), "invMass", 
		      Form("lambda1_%s_%s", object.c_str(), (isPbPb?"PbPb":"PP")), 
		      Form("lambda2_%s_%s", object.c_str(), (isPbPb?"PbPb":"PP"))
		      ));

      cout << Form("[INFO] %s Background 2nd Order Chebychev PDF in %s included", object.c_str(), (isPbPb?"PbPb":"PP")) << endl; break; 

    case (MassModel::ThirdOrderChebychev): 

      // check that all input parameters are defined 
      if (!( 
            parIni.count(Form("lambda1_%s_%s", object.c_str(), (isPbPb?"PbPb":"PP"))) && 
            parIni.count(Form("lambda2_%s_%s", object.c_str(), (isPbPb?"PbPb":"PP"))) && 
            parIni.count(Form("lambda3_%s_%s", object.c_str(), (isPbPb?"PbPb":"PP"))) 
             )) { 
	cout << Form("[ERROR] Initial parameters where not found for %s Background Third Order Chebychev in %s", object.c_str(), (isPbPb?"PbPb":"PP")) << endl; return false;
      }

      // create the variables for this model 
      ws.factory( parIni[Form("lambda1_%s_%s", object.c_str(), (isPbPb?"PbPb":"PP"))].c_str() );
      ws.factory( parIni[Form("lambda2_%s_%s", object.c_str(), (isPbPb?"PbPb":"PP"))].c_str() );
      ws.factory( parIni[Form("lambda3_%s_%s", object.c_str(), (isPbPb?"PbPb":"PP"))].c_str() );

      // create the PDF                 
      ws.factory(Form("Chebychev::%s(%s, {%s, %s, %s})", Form("pdfMASS_%s_%s", object.c_str(), (isPbPb?"PbPb":"PP")), "invMass", 
		      Form("lambda1_%s_%s", object.c_str(), (isPbPb?"PbPb":"PP")), 
		      Form("lambda2_%s_%s", object.c_str(), (isPbPb?"PbPb":"PP")), 
		      Form("lambda3_%s_%s", object.c_str(), (isPbPb?"PbPb":"PP"))
		      ));

      cout << Form("[INFO] %s Background 3rd Order Polynomial PDF in %s included", object.c_str(), (isPbPb?"PbPb":"PP")) << endl; break; 

    case (MassModel::FourthOrderChebychev): 

      // check that all input parameters are defined 
      if (!( 
            parIni.count(Form("lambda1_%s_%s", object.c_str(), (isPbPb?"PbPb":"PP"))) && 
            parIni.count(Form("lambda2_%s_%s", object.c_str(), (isPbPb?"PbPb":"PP"))) && 
            parIni.count(Form("lambda3_%s_%s", object.c_str(), (isPbPb?"PbPb":"PP"))) && 
            parIni.count(Form("lambda4_%s_%s", object.c_str(), (isPbPb?"PbPb":"PP"))) 
             )) { 
	cout << Form("[ERROR] Initial parameters where not found for %s Background Fourth Order Chebychev in %s", object.c_str(), (isPbPb?"PbPb":"PP")) << endl; return false;
      }

      // create the variables for this model        
      ws.factory( parIni[Form("lambda1_%s_%s", object.c_str(), (isPbPb?"PbPb":"PP"))].c_str() );
      ws.factory( parIni[Form("lambda2_%s_%s", object.c_str(), (isPbPb?"PbPb":"PP"))].c_str() );
      ws.factory( parIni[Form("lambda3_%s_%s", object.c_str(), (isPbPb?"PbPb":"PP"))].c_str() );
      ws.factory( parIni[Form("lambda4_%s_%s", object.c_str(), (isPbPb?"PbPb":"PP"))].c_str() );

      // create the PDF                 
      ws.factory(Form("Polynomial::%s(%s, {%s, %s, %s, %s})", Form("pdfMASS_%s_%s", object.c_str(), (isPbPb?"PbPb":"PP")), "invMass", 
		      Form("lambda1_%s_%s", object.c_str(), (isPbPb?"PbPb":"PP")), 
		      Form("lambda2_%s_%s", object.c_str(), (isPbPb?"PbPb":"PP")), 
		      Form("lambda3_%s_%s", object.c_str(), (isPbPb?"PbPb":"PP")), 
		      Form("lambda4_%s_%s", object.c_str(), (isPbPb?"PbPb":"PP"))
		      ));

      cout << Form("[INFO] %s Background 4th Order Chebychev PDF in %s included", object.c_str(), (isPbPb?"PbPb":"PP")) << endl; break; 

    case (MassModel::FifthOrderChebychev): 

      // check that all input parameters are defined 
      if (!( 
            parIni.count(Form("lambda1_%s_%s", object.c_str(), (isPbPb?"PbPb":"PP"))) && 
            parIni.count(Form("lambda2_%s_%s", object.c_str(), (isPbPb?"PbPb":"PP"))) && 
            parIni.count(Form("lambda3_%s_%s", object.c_str(), (isPbPb?"PbPb":"PP"))) && 
            parIni.count(Form("lambda4_%s_%s", object.c_str(), (isPbPb?"PbPb":"PP"))) && 
            parIni.count(Form("lambda5_%s_%s", object.c_str(), (isPbPb?"PbPb":"PP"))) 
             )) { 
	cout << Form("[ERROR] Initial parameters where not found for %s Background Fifth Order Chebychev in %s", object.c_str(), (isPbPb?"PbPb":"PP")) << endl; return false;
      }

      // create the variables for this model        
      ws.factory( parIni[Form("lambda1_%s_%s", object.c_str(), (isPbPb?"PbPb":"PP"))].c_str() );
      ws.factory( parIni[Form("lambda2_%s_%s", object.c_str(), (isPbPb?"PbPb":"PP"))].c_str() );
      ws.factory( parIni[Form("lambda3_%s_%s", object.c_str(), (isPbPb?"PbPb":"PP"))].c_str() );
      ws.factory( parIni[Form("lambda4_%s_%s", object.c_str(), (isPbPb?"PbPb":"PP"))].c_str() );
      ws.factory( parIni[Form("lambda5_%s_%s", object.c_str(), (isPbPb?"PbPb":"PP"))].c_str() );

      // create the PDF                 
      ws.factory(Form("Polynomial::%s(%s, {%s, %s, %s, %s, %s})", Form("pdfMASS_%s_%s", object.c_str(), (isPbPb?"PbPb":"PP")), "invMass", 
		      Form("lambda1_%s_%s", object.c_str(), (isPbPb?"PbPb":"PP")), 
		      Form("lambda2_%s_%s", object.c_str(), (isPbPb?"PbPb":"PP")), 
		      Form("lambda3_%s_%s", object.c_str(), (isPbPb?"PbPb":"PP")), 
		      Form("lambda4_%s_%s", object.c_str(), (isPbPb?"PbPb":"PP")), 
		      Form("lambda5_%s_%s", object.c_str(), (isPbPb?"PbPb":"PP"))
		      ));

      cout << Form("[INFO] %s Background 5th Order Chebychev PDF in %s included", object.c_str(), (isPbPb?"PbPb":"PP")) << endl; break; 

    case (MassModel::SixthOrderChebychev): 

      // check that all input parameters are defined 
      if (!( 
            parIni.count(Form("lambda1_%s_%s", object.c_str(), (isPbPb?"PbPb":"PP"))) && 
            parIni.count(Form("lambda2_%s_%s", object.c_str(), (isPbPb?"PbPb":"PP"))) && 
            parIni.count(Form("lambda3_%s_%s", object.c_str(), (isPbPb?"PbPb":"PP"))) && 
            parIni.count(Form("lambda4_%s_%s", object.c_str(), (isPbPb?"PbPb":"PP"))) && 
            parIni.count(Form("lambda5_%s_%s", object.c_str(), (isPbPb?"PbPb":"PP"))) && 
            parIni.count(Form("lambda6_%s_%s", object.c_str(), (isPbPb?"PbPb":"PP"))) 
             )) { 
	cout << Form("[ERROR] Initial parameters where not found for %s Background Sixth Order Chebychev in %s", object.c_str(), (isPbPb?"PbPb":"PP")) << endl; return false;
      }

      // create the variables for this model        
      ws.factory( parIni[Form("lambda1_%s_%s", object.c_str(), (isPbPb?"PbPb":"PP"))].c_str() );
      ws.factory( parIni[Form("lambda2_%s_%s", object.c_str(), (isPbPb?"PbPb":"PP"))].c_str() );
      ws.factory( parIni[Form("lambda3_%s_%s", object.c_str(), (isPbPb?"PbPb":"PP"))].c_str() );
      ws.factory( parIni[Form("lambda4_%s_%s", object.c_str(), (isPbPb?"PbPb":"PP"))].c_str() );
      ws.factory( parIni[Form("lambda5_%s_%s", object.c_str(), (isPbPb?"PbPb":"PP"))].c_str() );
      ws.factory( parIni[Form("lambda6_%s_%s", object.c_str(), (isPbPb?"PbPb":"PP"))].c_str() );

      // create the PDF                 
      ws.factory(Form("Polynomial::%s(%s, {%s, %s, %s, %s, %s, %s})", Form("pdfMASS_%s_%s", object.c_str(), (isPbPb?"PbPb":"PP")), "invMass", 
		      Form("lambda1_%s_%s", object.c_str(), (isPbPb?"PbPb":"PP")), 
		      Form("lambda2_%s_%s", object.c_str(), (isPbPb?"PbPb":"PP")), 
		      Form("lambda3_%s_%s", object.c_str(), (isPbPb?"PbPb":"PP")), 
		      Form("lambda4_%s_%s", object.c_str(), (isPbPb?"PbPb":"PP")), 
		      Form("lambda5_%s_%s", object.c_str(), (isPbPb?"PbPb":"PP")), 
		      Form("lambda6_%s_%s", object.c_str(), (isPbPb?"PbPb":"PP"))
		      ));

      cout << Form("[INFO] %s Background 6th Order Chebychev PDF in %s included", object.c_str(), (isPbPb?"PbPb":"PP")) << endl; break; 

    case (MassModel::Exponential): 

      // check that all input parameters are defined 
      if (!( 
            parIni.count(Form("lambda1_%s_%s", object.c_str(), (isPbPb?"PbPb":"PP"))) 
             )) { 
	cout << Form("[ERROR] Initial parameters where not found for %s Background Exponential in %s", object.c_str(), (isPbPb?"PbPb":"PP")) << endl; return false;
      }

      // create the variables for this model        
      ws.factory( parIni[Form("lambda1_%s_%s", object.c_str(), (isPbPb?"PbPb":"PP"))].c_str() );

      // create the PDF                 
      ws.factory(Form("Exponential::%s(%s, %s)", Form("pdfMASS_%s_%s", object.c_str(), (isPbPb?"PbPb":"PP")), "invMass", 
		      Form("lambda1_%s_%s", object.c_str(), (isPbPb?"PbPb":"PP"))
		      ));

      cout << Form("[INFO] %s Background Exponential PDF in %s included", object.c_str(), (isPbPb?"PbPb":"PP")) << endl; break;
      
    default :
      
      cout << "[ERROR] Selected Background Mass Model: " << parIni[Form("Model_%s_%s", object.c_str(), (isPbPb?"PbPb":"PP"))] << " has not been implemented" << endl; return false;
    
    }
  
  return true;
};


bool addSignalMassModel(RooWorkspace& ws, string object, MassModel model, map<string,string> parIni, bool isPbPb) 
{
  cout << Form("[INFO] Implementing %s Mass Model", object.c_str()) << endl;
  
  switch(model) 
    {    
    case (MassModel::SingleGaussian): 
      
      // check that all input parameters are defined
      if (!( 
            parIni.count(Form("m_%s_%s", object.c_str(), (isPbPb?"PbPb":"PP"))) && 
            parIni.count(Form("sigma1_%s_%s", object.c_str(), (isPbPb?"PbPb":"PP"))) 
             )) { 
	cout << Form("[ERROR] Initial parameters where not found for %s Single Gaussian Model in %s", object.c_str(), (isPbPb?"PbPb":"PP")) << endl; return false;
      }

      // create the variables for this model        
      ws.factory( parIni[Form("m_%s_%s", object.c_str(), (isPbPb?"PbPb":"PP"))].c_str() );
      ws.factory( parIni[Form("sigma1_%s_%s", object.c_str(), (isPbPb?"PbPb":"PP"))].c_str() );

      // create the PDF                       
      ws.factory(Form("Gaussian::%s(%s, %s, %s)", Form("pdfMASS_%s_%s", object.c_str(), (isPbPb?"PbPb":"PP")), "invMass", 
		      Form("m_%s_%s", object.c_str(), (isPbPb?"PbPb":"PP")), 
		      Form("sigma1_%s_%s", object.c_str(), (isPbPb?"PbPb":"PP"))
		      ));
    
      cout << Form("[INFO] %s Single Gaussian PDF in %s included", object.c_str(), (isPbPb?"PbPb":"PP")) << endl; break;  
      
    case (MassModel::DoubleGaussian): 

      // check that all input parameters are defined
      if (!( 
            parIni.count(Form("m_%s_%s", object.c_str(), (isPbPb?"PbPb":"PP"))) && 
            parIni.count(Form("sigma1_%s_%s", object.c_str(), (isPbPb?"PbPb":"PP"))) && 
            parIni.count(Form("sigma2_%s_%s", object.c_str(), (isPbPb?"PbPb":"PP"))) && 
            parIni.count(Form("f_%s_%s", object.c_str(), (isPbPb?"PbPb":"PP"))) 
             )) { 
	cout << Form("[ERROR] Initial parameters where not found for %s Double Gaussian Model in %s", object.c_str(), (isPbPb?"PbPb":"PP")) << endl; return false; 
      }

      // create the variables for this model              
      ws.factory( parIni[Form("f_%s_%s", object.c_str(), (isPbPb?"PbPb":"PP"))].c_str() );
      ws.factory( parIni[Form("m_%s_%s", object.c_str(), (isPbPb?"PbPb":"PP"))].c_str() ); 
      ws.factory( parIni[Form("sigma1_%s_%s", object.c_str(), (isPbPb?"PbPb":"PP"))].c_str() );
      ws.factory( parIni[Form("sigma2_%s_%s", object.c_str(), (isPbPb?"PbPb":"PP"))].c_str() );

      // create the two PDFs             
      ws.factory(Form("Gaussian::%s(%s, %s, %s)", Form("pdfMASSG1_%s_%s", object.c_str(), (isPbPb?"PbPb":"PP")), "invMass", 
                      Form("m_%s_%s", object.c_str(), (isPbPb?"PbPb":"PP")), 
                      Form("sigma1_%s_%s", object.c_str(), (isPbPb?"PbPb":"PP"))
                      ));    
      ws.factory(Form("Gaussian::%s(%s, %s, %s)", Form("pdfMASSG2_%s_%s", object.c_str(), (isPbPb?"PbPb":"PP")), "invMass", 
                      Form("m_%s_%s", object.c_str(), (isPbPb?"PbPb":"PP")), 
                      Form("sigma2_%s_%s", object.c_str(), (isPbPb?"PbPb":"PP"))
                      ));               

      // Sum the PDFs to get the signal PDF
      ws.factory(Form("SUM::%s(%s*%s, %s)", Form("pdfMASS_%s_%s", object.c_str(), (isPbPb?"PbPb":"PP")),
		      Form("f_%s_%s", object.c_str(), (isPbPb?"PbPb":"PP")),
		      Form("pdfMASSG1_%s_%s", object.c_str(), (isPbPb?"PbPb":"PP")),
		      Form("pdfMASSG2_%s_%s", object.c_str(), (isPbPb?"PbPb":"PP"))
		      ));

      cout << Form("[INFO] %s Double Gaussian PDF in %s included", object.c_str(), (isPbPb?"PbPb":"PP")) << endl; break; 

    case (MassModel::SingleCrystalBall):  

      // check that all input parameters are defined
      if (!( 
            parIni.count(Form("m_%s_%s", object.c_str(), (isPbPb?"PbPb":"PP"))) && 
            parIni.count(Form("sigma1_%s_%s", object.c_str(), (isPbPb?"PbPb":"PP"))) && 
            parIni.count(Form("alpha_%s_%s", object.c_str(), (isPbPb?"PbPb":"PP"))) && 
            parIni.count(Form("n_%s_%s", object.c_str(), (isPbPb?"PbPb":"PP"))) 
             )) {
	cout << Form("[ERROR] Initial parameters where not found for %s Single Crystal Ball Model in %s", object.c_str(), (isPbPb?"PbPb":"PP")) << endl; return false; 
      }

      // create the variables for this model             
      ws.factory( parIni[Form("m_%s_%s", object.c_str(), (isPbPb?"PbPb":"PP"))].c_str() );
      ws.factory( parIni[Form("sigma1_%s_%s", object.c_str(), (isPbPb?"PbPb":"PP"))].c_str() );
      ws.factory( parIni[Form("alpha_%s_%s", object.c_str(), (isPbPb?"PbPb":"PP"))].c_str() );
      ws.factory( parIni[Form("n_%s_%s", object.c_str(), (isPbPb?"PbPb":"PP"))].c_str() );

      // create the PDF              
      ws.factory(Form("CBShape::%s(%s, %s, %s, %s, %s)", Form("pdfMASS_%s_%s", object.c_str(), (isPbPb?"PbPb":"PP")), "invMass", 
		      Form("m_%s_%s", object.c_str(), (isPbPb?"PbPb":"PP")), 
		      Form("sigma1_%s_%s", object.c_str(), (isPbPb?"PbPb":"PP")),
		      Form("alpha_%s_%s", object.c_str(), (isPbPb?"PbPb":"PP")),
		      Form("n_%s_%s", object.c_str(), (isPbPb?"PbPb":"PP"))
		      ));

      cout << Form("[INFO] %s Single Crystal Ball PDF in %s included", object.c_str(), (isPbPb?"PbPb":"PP")) << endl; break;
      
    case (MassModel::DoubleCrystalBall): 
      
      // check that all input parameters are defined
        if (!(
              parIni.count(Form("m_%s_%s", object.c_str(), (isPbPb?"PbPb":"PP"))) &&
              parIni.count(Form("sigma1_%s_%s", object.c_str(), (isPbPb?"PbPb":"PP"))) &&
              parIni.count(Form("sigma2_%s_%s", object.c_str(), (isPbPb?"PbPb":"PP"))) &&
              parIni.count(Form("alpha_%s_%s", object.c_str(), (isPbPb?"PbPb":"PP"))) &&
              parIni.count(Form("n_%s_%s", object.c_str(), (isPbPb?"PbPb":"PP"))) &&
              parIni.count(Form("f_%s_%s", object.c_str(), (isPbPb?"PbPb":"PP")))
              )) {
          cout << Form("[ERROR] Initial parameters where not found for %s Double Crystal Ball Model in %s", object.c_str(), (isPbPb?"PbPb":"PP")) << endl; return false;
        }
        
      // create the variables for this model             
      ws.factory( parIni[Form("f_%s_%s", object.c_str(), (isPbPb?"PbPb":"PP"))].c_str() );
      ws.factory( parIni[Form("m_%s_%s", object.c_str(), (isPbPb?"PbPb":"PP"))].c_str() );
      ws.factory( parIni[Form("sigma1_%s_%s", object.c_str(), (isPbPb?"PbPb":"PP"))].c_str() );
      ws.factory( parIni[Form("alpha_%s_%s", object.c_str(), (isPbPb?"PbPb":"PP"))].c_str() );
      ws.factory( parIni[Form("n_%s_%s", object.c_str(), (isPbPb?"PbPb":"PP"))].c_str() );
      ws.factory( parIni[Form("sigma2_%s_%s", object.c_str(), (isPbPb?"PbPb":"PP"))].c_str() );

      // create the two PDFs
      ws.factory(Form("CBShape::%s(%s, %s, %s, %s, %s)", Form("pdfMASSCB1_%s_%s", object.c_str(), (isPbPb?"PbPb":"PP")), "invMass", 
                      Form("m_%s_%s", object.c_str(), (isPbPb?"PbPb":"PP")), 
                      Form("sigma1_%s_%s", object.c_str(), (isPbPb?"PbPb":"PP")),
                      Form("alpha_%s_%s", object.c_str(), (isPbPb?"PbPb":"PP")),
                      Form("n_%s_%s", object.c_str(), (isPbPb?"PbPb":"PP"))
                      ));   
      ws.factory(Form("CBShape::%s(%s, %s, %s, %s, %s)", Form("pdfMASSCB2_%s_%s", object.c_str(), (isPbPb?"PbPb":"PP")), "invMass", 
                      Form("m_%s_%s", object.c_str(), (isPbPb?"PbPb":"PP")), 
                      Form("sigma2_%s_%s", object.c_str(), (isPbPb?"PbPb":"PP")),
                      Form("alpha_%s_%s", object.c_str(), (isPbPb?"PbPb":"PP")),
                      Form("n_%s_%s", object.c_str(), (isPbPb?"PbPb":"PP"))
                      ));   

      // Sum the PDFs to get the signal PDF
      ws.factory(Form("SUM::%s(%s*%s, %s)", Form("pdfMASS_%s_%s", object.c_str(), (isPbPb?"PbPb":"PP")),
		      Form("f_%s_%s", object.c_str(), (isPbPb?"PbPb":"PP")),
		      Form("pdfMASSCB1_%s_%s", object.c_str(), (isPbPb?"PbPb":"PP")),
		      Form("pdfMASSCB2_%s_%s", object.c_str(), (isPbPb?"PbPb":"PP"))
		      ));   

      cout << Form("[INFO] %s Double Crystal Ball PDF in %s included", object.c_str(), (isPbPb?"PbPb":"PP")) << endl; break;
      
    case (MassModel::GaussianAndCrystalBall):

      // check that all input parameters are defined
      if (!( 
            parIni.count(Form("m_%s_%s", object.c_str(), (isPbPb?"PbPb":"PP"))) && 
            parIni.count(Form("sigma1_%s_%s", object.c_str(), (isPbPb?"PbPb":"PP"))) && 
            parIni.count(Form("sigma2_%s_%s", object.c_str(), (isPbPb?"PbPb":"PP"))) && 
            parIni.count(Form("alpha_%s_%s", object.c_str(), (isPbPb?"PbPb":"PP"))) && 
            parIni.count(Form("n_%s_%s", object.c_str(), (isPbPb?"PbPb":"PP"))) && 
            parIni.count(Form("f_%s_%s", object.c_str(), (isPbPb?"PbPb":"PP")))  
             )) { 
	cout << Form("[ERROR] Initial parameters where not found for %s Gaussian and Crystal Ball Model in %s", object.c_str(), (isPbPb?"PbPb":"PP")) << endl; return false;
      } 

      // create the variables for this model             
      ws.factory( parIni[Form("f_%s_%s", object.c_str(), (isPbPb?"PbPb":"PP"))].c_str() );
      ws.factory( parIni[Form("m_%s_%s", object.c_str(), (isPbPb?"PbPb":"PP"))].c_str() );
      ws.factory( parIni[Form("sigma1_%s_%s", object.c_str(), (isPbPb?"PbPb":"PP"))].c_str() );
      ws.factory( parIni[Form("sigma2_%s_%s", object.c_str(), (isPbPb?"PbPb":"PP"))].c_str() );
      ws.factory( parIni[Form("alpha_%s_%s", object.c_str(), (isPbPb?"PbPb":"PP"))].c_str() );
      ws.factory( parIni[Form("n_%s_%s", object.c_str(), (isPbPb?"PbPb":"PP"))].c_str() );
      
      // create the two PDFs 
      ws.factory(Form("Gaussian::%s(%s, %s, %s)", Form("pdfMASSG1_%s_%s", object.c_str(), (isPbPb?"PbPb":"PP")), "invMass", 
                      Form("m_%s_%s", object.c_str(), (isPbPb?"PbPb":"PP")), 
                      Form("sigma1_%s_%s", object.c_str(), (isPbPb?"PbPb":"PP"))
                      ));   

      ws.factory(Form("CBShape::%s(%s, %s, %s, %s, %s)", Form("pdfMASSCB1_%s_%s", object.c_str(), (isPbPb?"PbPb":"PP")), "invMass", 
                      Form("m_%s_%s", object.c_str(), (isPbPb?"PbPb":"PP")), 
                      Form("sigma2_%s_%s", object.c_str(), (isPbPb?"PbPb":"PP")),
                      Form("alpha_%s_%s", object.c_str(), (isPbPb?"PbPb":"PP")),
                      Form("n_%s_%s", object.c_str(), (isPbPb?"PbPb":"PP"))
                      ));   
 
      // Sum the PDFs to get the signal PDF 
      ws.factory(Form("SUM::%s(%s*%s, %s)", Form("pdfMASS_%s_%s", object.c_str(), (isPbPb?"PbPb":"PP")),
		      Form("f_%s_%s", object.c_str(), (isPbPb?"PbPb":"PP")),
		      Form("pdfMASSG1_%s_%s", object.c_str(), (isPbPb?"PbPb":"PP")),
		      Form("pdfMASSCB1_%s_%s", object.c_str(), (isPbPb?"PbPb":"PP"))
		      ));   

      cout << Form("[INFO] %s Gaussian and Crystal Ball PDF in %s included", object.c_str(), (isPbPb?"PbPb":"PP")) << endl; break;
	
    default :

      cout << "[ERROR] Selected Signal Mass Model: " << parIni[Form("Model_%s_%s", object.c_str(), (isPbPb?"PbPb":"PP"))] << " has not been implemented" << endl; return false;

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
  parIni[Form("m_Psi2S_%s", (isPbPb?"PbPb":"PP"))]      = Form("RooFormulaVar::%s('@0*@1',{MassRatio[%.6f],%s})", Form("m_Psi2S_%s", (isPbPb?"PbPb":"PP")), MassRatio, Form("m_Jpsi_%s", (isPbPb?"PbPb":"PP") ));
  parIni[Form("sigma1_Psi2S_%s", (isPbPb?"PbPb":"PP"))] = Form("RooFormulaVar::%s('@0*@1',{MassRatio,%s})", Form("sigma1_Psi2S_%s", (isPbPb?"PbPb":"PP")), Form("sigma1_Jpsi_%s", (isPbPb?"PbPb":"PP") ));
  parIni[Form("sigma2_Psi2S_%s", (isPbPb?"PbPb":"PP"))] = Form("RooFormulaVar::%s('@0*@1',{MassRatio,%s})", Form("sigma2_Psi2S_%s", (isPbPb?"PbPb":"PP")), Form("sigma2_Jpsi_%s", (isPbPb?"PbPb":"PP") ));
  parIni[Form("alpha_Psi2S_%s", (isPbPb?"PbPb":"PP"))]  = Form("RooFormulaVar::%s('@0',{%s})", Form("alpha_Psi2S_%s", (isPbPb?"PbPb":"PP")), Form("alpha_Jpsi_%s", (isPbPb?"PbPb":"PP")));
  parIni[Form("n_Psi2S_%s", (isPbPb?"PbPb":"PP"))]      = Form("RooFormulaVar::%s('@0',{%s})", Form("n_Psi2S_%s", (isPbPb?"PbPb":"PP")), Form("n_Jpsi_%s", (isPbPb?"PbPb":"PP")));
  parIni[Form("f_Psi2S_%s", (isPbPb?"PbPb":"PP"))]      = Form("RooFormulaVar::%s('@0',{%s})", Form("f_Psi2S_%s", (isPbPb?"PbPb":"PP")), Form("f_Jpsi_%s", (isPbPb?"PbPb":"PP")));
};

void setDefaultParameters(map<string, string> &parIni, bool isPbPb, double numEntries)
{

  cout << "[INFO] Setting user undefined initial parameters to their default values" << endl;

  // DEFAULT SINGLE AND DOUBLE RATIO PARAMETERS
  if (parIni.count("RFrac2Svs1S_PbPbvsPP")==0 || parIni["RFrac2Svs1S_PbPbvsPP"]=="") { 
    parIni["RFrac2Svs1S_PbPbvsPP"] = Form("%s[%.4f,%.4f,%.4f]", "RFrac2Svs1S_PbPbvsPP", 0.26, 0.0, 3.0);
  }
  if (parIni.count(Form("RFrac2Svs1S_%s", (isPbPb?"PbPb":"PP")))==0 || parIni[Form("RFrac2Svs1S_%s", (isPbPb?"PbPb":"PP"))]=="") {
    parIni[Form("RFrac2Svs1S_%s", (isPbPb?"PbPb":"PP"))] = Form("%s[%.4f,%.4f,%.4f]", Form("RFrac2Svs1S_%s", (isPbPb?"PbPb":"PP")), 0.26, 0.0, 1.0);
  }

  // DEFAULT RANGE OF NUMBER OF EVENTS
  if (parIni.count(Form("N_Jpsi_%s", (isPbPb?"PbPb":"PP")))==0 || parIni[Form("N_Jpsi_%s", (isPbPb?"PbPb":"PP"))]=="") { 
    parIni[Form("N_Jpsi_%s", (isPbPb?"PbPb":"PP"))]  = Form("%s[%.10f,%.10f]", Form("N_Jpsi_%s", (isPbPb?"PbPb":"PP")), 0.0, numEntries*10.0);
  }
  if (parIni.count(Form("N_Psi2S_%s", (isPbPb?"PbPb":"PP")))==0 || parIni[Form("N_Psi2S_%s", (isPbPb?"PbPb":"PP"))]=="") { 
    parIni[Form("N_Psi2S_%s", (isPbPb?"PbPb":"PP"))]  = Form("%s[%.10f,%.10f]", Form("N_Psi2S_%s", (isPbPb?"PbPb":"PP")), 0.0, numEntries*10.0);
  }
  if (parIni.count(Form("N_Bkg_%s", (isPbPb?"PbPb":"PP")))==0 || parIni[Form("N_Bkg_%s", (isPbPb?"PbPb":"PP"))]=="") { 
    parIni[Form("N_Bkg_%s", (isPbPb?"PbPb":"PP"))]  = Form("%s[%.10f,%.10f]", Form("N_Bkg_%s", (isPbPb?"PbPb":"PP")), 0.0, numEntries);
  }

  // DEFAULT SIGNAL MASS MODEL PARAMETERS 
  if (parIni.count(Form("m_Jpsi_%s", (isPbPb?"PbPb":"PP")))==0 || parIni[Form("m_Jpsi_%s", (isPbPb?"PbPb":"PP"))]=="") {
    parIni[Form("m_Jpsi_%s", (isPbPb?"PbPb":"PP"))] = Form("%s[%.6f,%.6f,%.6f]", Form("m_Jpsi_%s", (isPbPb?"PbPb":"PP")), Mass.JPsi, Mass.JPsi-0.1, Mass.JPsi+0.1);
  }
  if (parIni.count(Form("m_Psi2S_%s", (isPbPb?"PbPb":"PP")))==0 || parIni[Form("m_Psi2S_%s", (isPbPb?"PbPb":"PP"))]=="") {
    parIni[Form("m_Psi2S_%s", (isPbPb?"PbPb":"PP"))] = Form("%s[%.6f,%.6f,%.6f]", Form("m_Psi2S_%s", (isPbPb?"PbPb":"PP")), Mass.Psi2S, Mass.Psi2S-0.1, Mass.Psi2S+0.1);
  }
  if (parIni.count(Form("sigma1_Jpsi_%s", (isPbPb?"PbPb":"PP")))==0 || parIni[Form("sigma1_Jpsi_%s", (isPbPb?"PbPb":"PP"))]=="") {
    parIni[Form("sigma1_Jpsi_%s", (isPbPb?"PbPb":"PP"))] = Form("%s[%.4f,%.4f,%.4f]", Form("sigma1_Jpsi_%s", (isPbPb?"PbPb":"PP")), 0.02, 0.005, 0.07);
  }
  if (parIni.count(Form("rSigma21_Jpsi_%s", (isPbPb?"PbPb":"PP")))==0) {
    if (parIni.count(Form("sigma2_Jpsi_%s", (isPbPb?"PbPb":"PP")))==0 || parIni[Form("sigma2_Jpsi_%s", (isPbPb?"PbPb":"PP"))]=="") {
      parIni[Form("sigma2_Jpsi_%s", (isPbPb?"PbPb":"PP"))] = Form("%s[%.4f,%.4f,%.4f]", Form("sigma2_Jpsi_%s", (isPbPb?"PbPb":"PP")), 0.04, 0.01, 0.10);
    }
  } else {
    if (parIni[Form("rSigma21_Jpsi_%s", (isPbPb?"PbPb":"PP"))]=="") {
      parIni[Form("rSigma21_Jpsi_%s", (isPbPb?"PbPb":"PP"))] = Form("%s[%.4f,%.4f,%.4f]", Form("rSigma21_Jpsi_%s", (isPbPb?"PbPb":"PP")), 2.0, 1.0, 4.0);
    }
    parIni[Form("sigma2_Jpsi_%s", (isPbPb?"PbPb":"PP"))] = Form("RooFormulaVar::%s('@0*@1',{%s,%s})", Form("sigma2_Jpsi_%s", (isPbPb?"PbPb":"PP")),
                                                                parIni[Form("rSigma21_Jpsi_%s", (isPbPb?"PbPb":"PP"))].c_str(), Form("sigma1_Jpsi_%s", (isPbPb?"PbPb":"PP") ));
  }      
  if (parIni.count(Form("sigma1_Psi2S_%s", (isPbPb?"PbPb":"PP")))==0 || parIni[Form("sigma1_Psi2S_%s", (isPbPb?"PbPb":"PP"))]=="") {
    parIni[Form("sigma1_Psi2S_%s", (isPbPb?"PbPb":"PP"))] = Form("%s[%.4f,%.4f,%.4f]", Form("sigma1_Psi2S_%s", (isPbPb?"PbPb":"PP")), 0.02, 0.005, 0.07);
  }
  if (parIni.count(Form("rSigma21_Psi2S_%s", (isPbPb?"PbPb":"PP")))==0) {
    if (parIni.count(Form("sigma2_Psi2S_%s", (isPbPb?"PbPb":"PP")))==0 || parIni[Form("sigma2_Psi2S_%s", (isPbPb?"PbPb":"PP"))]=="") {
      parIni[Form("sigma2_Psi2S_%s", (isPbPb?"PbPb":"PP"))] = Form("%s[%.4f,%.4f,%.4f]", Form("sigma2_Psi2S_%s", (isPbPb?"PbPb":"PP")), 0.04, 0.01, 0.10);
    }
  } else {
    if (parIni[Form("rSigma21_Psi2S_%s", (isPbPb?"PbPb":"PP"))]=="") {
      parIni[Form("rSigma21_Psi2S_%s", (isPbPb?"PbPb":"PP"))] = Form("%s[%.4f,%.4f,%.4f]", Form("rSigma21_Psi2S_%s", (isPbPb?"PbPb":"PP")), 2.0, 1.0, 4.0);
    }
    parIni[Form("sigma2_Psi2S_%s", (isPbPb?"PbPb":"PP"))] = Form("RooFormulaVar::%s('@0*@1',{%s,%s})", Form("sigma2_Psi2S_%s", (isPbPb?"PbPb":"PP")), 
                                                                 parIni[Form("rSigma21_Psi2S_%s", (isPbPb?"PbPb":"PP"))].c_str(), Form("sigma1_Psi2S_%s", (isPbPb?"PbPb":"PP") ));
  }
  if (parIni.count(Form("alpha_Jpsi_%s", (isPbPb?"PbPb":"PP")))==0 || parIni[Form("alpha_Jpsi_%s", (isPbPb?"PbPb":"PP"))]=="") {
    parIni[Form("alpha_Jpsi_%s", (isPbPb?"PbPb":"PP"))] = Form("%s[%.4f,%.4f,%.4f]", Form("alpha_Jpsi_%s", (isPbPb?"PbPb":"PP")), 2.21, 0.5, 30.0);
  }
  if (parIni.count(Form("alpha_Psi2S_%s", (isPbPb?"PbPb":"PP")))==0 || parIni[Form("alpha_Psi2S_%s", (isPbPb?"PbPb":"PP"))]=="") {
    parIni[Form("alpha_Psi2S_%s", (isPbPb?"PbPb":"PP"))] = Form("%s[%.4f,%.4f,%.4f]", Form("alpha_Psi2S_%s", (isPbPb?"PbPb":"PP")), 2.21, 0.5, 30.0);
  }
  if (parIni.count(Form("n_Jpsi_%s", (isPbPb?"PbPb":"PP")))==0 || parIni[Form("n_Jpsi_%s", (isPbPb?"PbPb":"PP"))]=="") {
    parIni[Form("n_Jpsi_%s", (isPbPb?"PbPb":"PP"))] = Form("%s[%.4f,%.4f,%.4f]", Form("n_Jpsi_%s", (isPbPb?"PbPb":"PP")), 2., 0.5, 30.0);
  }
  if (parIni.count(Form("n_Psi2S_%s", (isPbPb?"PbPb":"PP")))==0 || parIni[Form("n_Psi2S_%s", (isPbPb?"PbPb":"PP"))]=="") {
    parIni[Form("n_Psi2S_%s", (isPbPb?"PbPb":"PP"))] = Form("%s[%.4f,%.4f,%.4f]", Form("n_Psi2S_%s", (isPbPb?"PbPb":"PP")), 2., 0.5, 30.0);
  }
  if (parIni.count(Form("f_Jpsi_%s", (isPbPb?"PbPb":"PP")))==0 || parIni[Form("f_Jpsi_%s", (isPbPb?"PbPb":"PP"))]=="") {
    parIni[Form("f_Jpsi_%s", (isPbPb?"PbPb":"PP"))] = Form("%s[%.4f,%.4f,%.4f]", Form("f_Jpsi_%s", (isPbPb?"PbPb":"PP")), 0.3, 0.0, 1.0);
  }
  if (parIni.count(Form("f_Psi2S_%s", (isPbPb?"PbPb":"PP")))==0 || parIni[Form("f_Psi2S_%s", (isPbPb?"PbPb":"PP"))]=="") {
    parIni[Form("f_Psi2S_%s", (isPbPb?"PbPb":"PP"))] = Form("%s[%.4f,%.4f,%.4f]", Form("f_Psi2S_%s", (isPbPb?"PbPb":"PP")), 0.3, 0.0, 1.0);
  }

  // DEFAULT BACKGROUND MASS MODEL PARAMETERS
  if (parIni[Form("Model_Bkg_%s",(isPbPb?"PbPb":"PP"))].find("Chebychev")!=std::string::npos) {
    if (parIni.count(Form("lambda1_Bkg_%s", (isPbPb?"PbPb":"PP")))==0 || parIni[Form("lambda1_Bkg_%s", (isPbPb?"PbPb":"PP"))]=="") { 
      parIni[Form("lambda1_Bkg_%s", (isPbPb?"PbPb":"PP"))] = Form("%s[%.4f,%.4f,%.4f]", Form("lambda1_Bkg_%s", (isPbPb?"PbPb":"PP")), 0.0, -1.0, 1.0);
    }
    if (parIni.count(Form("lambda2_Bkg_%s", (isPbPb?"PbPb":"PP")))==0 || parIni[Form("lambda2_Bkg_%s", (isPbPb?"PbPb":"PP"))]=="") { 
      parIni[Form("lambda2_Bkg_%s", (isPbPb?"PbPb":"PP"))] = Form("%s[%.4f,%.4f,%.4f]", Form("lambda2_Bkg_%s", (isPbPb?"PbPb":"PP")), 0.0, -1.0, 1.0);
    }
    if (parIni.count(Form("lambda3_Bkg_%s", (isPbPb?"PbPb":"PP")))==0 || parIni[Form("lambda3_Bkg_%s", (isPbPb?"PbPb":"PP"))]=="") { 
      parIni[Form("lambda3_Bkg_%s", (isPbPb?"PbPb":"PP"))] = Form("%s[%.4f,%.4f,%.4f]", Form("lambda3_Bkg_%s", (isPbPb?"PbPb":"PP")), 0.0, -1.0, 1.0);
    }
    if (parIni.count(Form("lambda4_Bkg_%s", (isPbPb?"PbPb":"PP")))==0 || parIni[Form("lambda4_Bkg_%s", (isPbPb?"PbPb":"PP"))]=="") { 
      parIni[Form("lambda4_Bkg_%s", (isPbPb?"PbPb":"PP"))] = Form("%s[%.4f,%.4f,%.4f]", Form("lambda4_Bkg_%s", (isPbPb?"PbPb":"PP")), 0.0, -1.0, 1.0);
    }
    if (parIni.count(Form("lambda5_Bkg_%s", (isPbPb?"PbPb":"PP")))==0 || parIni[Form("lambda5_Bkg_%s", (isPbPb?"PbPb":"PP"))]=="") { 
      parIni[Form("lambda5_Bkg_%s", (isPbPb?"PbPb":"PP"))] = Form("%s[%.4f,%.4f,%.4f]", Form("lambda5_Bkg_%s", (isPbPb?"PbPb":"PP")), 0.0, -1.0, 1.0);
    }
    if (parIni.count(Form("lambda6_Bkg_%s", (isPbPb?"PbPb":"PP")))==0 || parIni[Form("lambda6_Bkg_%s", (isPbPb?"PbPb":"PP"))]=="") { 
      parIni[Form("lambda6_Bkg_%s", (isPbPb?"PbPb":"PP"))] = Form("%s[%.4f,%.4f,%.4f]", Form("lambda6_Bkg_%s", (isPbPb?"PbPb":"PP")), 0.0, -1.0, 1.0);
    }
  } 
  else if (parIni[Form("Model_Bkg_%s",(isPbPb?"PbPb":"PP"))].find("Exponential")!=std::string::npos) {
    if (parIni.count(Form("lambda1_Bkg_%s", (isPbPb?"PbPb":"PP")))==0 || parIni[Form("lambda1_Bkg_%s", (isPbPb?"PbPb":"PP"))]=="") { 
      parIni[Form("lambda1_Bkg_%s", (isPbPb?"PbPb":"PP"))] = Form("%s[%.4f,%.4f,%.4f]", Form("lambda1_Bkg_%s", (isPbPb?"PbPb":"PP")), 0.05, -100.0, 100.0);
    }
  }
 
};
