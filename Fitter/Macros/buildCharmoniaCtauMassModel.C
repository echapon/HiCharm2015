#include "Utilities/initClasses.h"


void setDefaultParameters(map<string, string> &parIni, int nt, bool isPbPb);
bool addSignalMassModel(RooWorkspace& ws, string object, MassModel model, map<string,string> parIni, bool isPbPb); 
bool addBackgroundMassModel(RooWorkspace& ws, string object, MassModel model, map<string,string> parIni, bool isPbPb);
bool defineCtauResolModel(RooWorkspace& ws, CtauModel model, map<string,string> parIni, bool isPbPb); 
bool addSignalCtauModel(RooWorkspace& ws, string object, CtauModel model, map<string,string> parIni, bool isPbPb); 
bool addBackgroundCtauModel(RooWorkspace& ws, string object, CtauModel model, map<string,string> parIni, bool isPbPb);


bool buildCharmoniaCtauMassModel(RooWorkspace& ws, struct InputOpt opt, struct CharmModel model, bool isPbPb, bool do2DFit)
{
  
  int nt= 0; 
  if (isPbPb) { nt = 250000; }
  else        { nt = 300000; }

  map<string, string> parIni;
  if (parIni.size()==0) { setDefaultParameters(parIni, nt, isPbPb); } 

  if (opt.inExcStat) {
    if (opt.doSimulFit && isPbPb) {
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
	
      ws.factory(Form("RooFormulaVar::%s('@0*@1*@2',{%s,%s,%s})", "N_Psi2S_PbPb", 
		      parIni["RFrac2Svs1S_PbPbvsPP"].c_str(), 
		      "RFrac2Svs1S_PP", 
		      parIni[Form("N_Jpsi_%s", (isPbPb?"PbPb":"PP"))].c_str() 
		      )); 
      ws.var("RFrac2Svs1S_PbPbvsPP")->setConstant(kFALSE);
      parIni[Form("N_Jpsi_%s", (isPbPb?"PbPb":"PP"))] = Form("N_Jpsi_%s", (isPbPb?"PbPb":"PP"));
      parIni[Form("N_Psi2S_%s", (isPbPb?"PbPb":"PP"))] = Form("N_Psi2S_%s", (isPbPb?"PbPb":"PP")); 
    } else {
      ws.factory(Form("RooFormulaVar::%s('@0*@1',{%s,%s})", Form("N_Psi2S_%s", (isPbPb?"PbPb":"PP")), 
		      parIni[Form("RFrac2Svs1S_%s", (isPbPb?"PbPb":"PP"))].c_str(), 
		      parIni[Form("N_Jpsi_%s", (isPbPb?"PbPb":"PP"))].c_str() 
		      )); 
      ws.var(Form("RFrac2Svs1S_%s", (isPbPb?"PbPb":"PP")))->setConstant(kFALSE);
    }
    parIni[Form("N_Jpsi_%s", (isPbPb?"PbPb":"PP"))] = Form("N_Jpsi_%s", (isPbPb?"PbPb":"PP"));
    parIni[Form("N_Psi2S_%s", (isPbPb?"PbPb":"PP"))] = Form("N_Psi2S_%s", (isPbPb?"PbPb":"PP")); 
  }

  if (do2DFit==false) {

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

  } else {

    // THIS IS THE 2D FIT PART, IT IS VERY PRELIMINARY, NEEDS TO BE TUNED 
    if(!defineCtauResolModel(ws, model.CtauRes, parIni, isPbPb)) { cout << "[ERROR] Defining the Ctau Resolution Model failed" << endl; return false; }
    
    if(!addSignalMassModel(ws, "Jpsi", model.Jpsi.Mass, parIni, isPbPb)) { cout << "[ERROR] Adding Prompt and NonPrompt Jpsi Mass Model failed" << endl; return false; }
    if(!addSignalCtauModel(ws, "JpsiP", model.Jpsi.Ctau.Prompt, parIni, isPbPb)) { cout << "[ERROR] Adding Prompt Jpsi Ctau Model failed" << endl; return false; }
    if(!addSignalCtauModel(ws, "JpsiNP", model.Jpsi.Ctau.NonPrompt, parIni, isPbPb)) { cout << "[ERROR] Adding NonPrompt Jpsi Ctau Model failed" << endl; return false; } 
    ws.factory(Form("SUM::%s(%s*%s, %s)", Form("pdfCTAU_Jpsi_%s", (isPbPb?"PbPb":"PP")),
		    parIni[Form("f_JpsiNP_%s", (isPbPb?"PbPb":"PP"))].c_str(),
		    Form("pdfCTAU_JpsiNP_%s", (isPbPb?"PbPb":"PP")),
		    Form("pdfCTAU_JpsiP_%s", (isPbPb?"PbPb":"PP"))
		    ));
    ws.factory(Form("PROD::%s(%s, %s)", Form("pdfCTAUMASS_Jpsi_%s", (isPbPb?"PbPb":"PP")),
		    Form("pdfCTAU_Jpsi_%s", (isPbPb?"PbPb":"PP")),
		    Form("pdfMASS_Jpsi_%s", (isPbPb?"PbPb":"PP"))
		    ));
    if (opt.inExcStat) { 
      if (!addSignalMassModel(ws, "Psi2S", model.Psi2S.Mass, parIni, isPbPb)) { cout << "[ERROR] Adding Psi(2S) Mass Model failed" << endl; return false; }
      if (!addSignalCtauModel(ws, "Psi2SP", model.Psi2S.Ctau.Prompt, parIni, isPbPb)) { cout << "[ERROR] Adding Prompt Psi(2S) Mass Model failed" << endl; return false; } 
      if (!addSignalCtauModel(ws, "Psi2SNP", model.Psi2S.Ctau.NonPrompt, parIni, isPbPb)) { cout << "[ERROR] Adding NonPrompt Psi(2S) Mass Model failed" << endl; return false; }   
      ws.factory(Form("SUM::%s(%s*%s, %s)", Form("pdfCTAU_Psi2S_%s", (isPbPb?"PbPb":"PP")),
		      parIni[Form("f_Psi2SNP_%s", (isPbPb?"PbPb":"PP"))].c_str(),
		      Form("pdfCTAU_Psi2SNP_%s", (isPbPb?"PbPb":"PP")),
		      Form("pdfCTAU_Psi2SP_%s", (isPbPb?"PbPb":"PP"))
		      ));
      ws.factory(Form("PROD::%s(%s, %s)", Form("pdfCTAUMASS_Psi2S_%s", (isPbPb?"PbPb":"PP")),
		      Form("pdfCTAU_Bkg_%s", (isPbPb?"PbPb":"PP")),
		      Form("pdfMASS_Bkg_%s", (isPbPb?"PbPb":"PP"))
		      ));
    }
    if(!addBackgroundMassModel(ws, "Bkg", model.Bkg.Mass, parIni, isPbPb)) { cout << "[ERROR] Adding Background Mass Model failed" << endl; return false; }
    if(!addBackgroundCtauModel(ws, "BkgP", model.Bkg.Ctau.Prompt, parIni, isPbPb)) { cout << "[ERROR] Adding Prompt Background Ctau Model failed" << endl; return false; }
    if(!addBackgroundCtauModel(ws, "BkgNP", model.Bkg.Ctau.NonPrompt, parIni, isPbPb)) { cout << "[ERROR] Adding NonPrompt Background Ctau Model failed" << endl; return false; }
    ws.factory(Form("SUM::%s(%s*%s, %s)", Form("pdfCTAU_Bkg_%s", (isPbPb?"PbPb":"PP")),
		    parIni[Form("f_BkgNP_%s", (isPbPb?"PbPb":"PP"))].c_str(),
		    Form("pdfCTAU_BkgNP_%s", (isPbPb?"PbPb":"PP")),
		    Form("pdfCTAU_BkgP_%s", (isPbPb?"PbPb":"PP"))
		    ));
    ws.factory(Form("PROD::%s(%s, %s)", Form("pdfCTAUMASS_Bkg_%s", (isPbPb?"PbPb":"PP")),
		    Form("pdfCTAU_Bkg_%s", (isPbPb?"PbPb":"PP")),
		    Form("pdfMASS_Bkg_%s", (isPbPb?"PbPb":"PP"))
		    ));
    
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
      ws.factory(Form("SUM::%s(%s*%s, %s*%s, %s*%s)", Form("pdfCTAUMASS_Tot_%s", (isPbPb?"PbPb":"PP")),
		      parIni[Form("N_Jpsi_%s", (isPbPb?"PbPb":"PP"))].c_str(),
		      Form("pdfCTAUMASS_Jpsi_%s", (isPbPb?"PbPb":"PP")),
		      parIni[Form("N_Psi2S_%s", (isPbPb?"PbPb":"PP"))].c_str(),
		      Form("pdfCTAUMASS_Psi2S_%s", (isPbPb?"PbPb":"PP")),
		      parIni[Form("N_Bkg_%s", (isPbPb?"PbPb":"PP"))].c_str(),
		      Form("pdfCTAUMASS_Bkg_%s", (isPbPb?"PbPb":"PP"))
		      ));
    } else {
      ws.factory(Form("SUM::%s(%s*%s, %s*%s)", Form("pdfMASSTot%s", (isPbPb?"PbPb":"PP")),
		      parIni[Form("N_Jpsi_%s", (isPbPb?"PbPb":"PP"))].c_str(),
		      Form("pdfMASS_Jpsi_%s", (isPbPb?"PbPb":"PP")),
		      parIni[Form("N_Bkg_%s", (isPbPb?"PbPb":"PP"))].c_str(),
		      Form("pdfMASS_Bkg_%s", (isPbPb?"PbPb":"PP"))
		      ));
      ws.factory(Form("SUM::%s(%s*%s, %s*%s)", Form("pdfCTAUMASS_Tot_%s", (isPbPb?"PbPb":"PP")),
		      parIni[Form("N_Jpsi_%s", (isPbPb?"PbPb":"PP"))].c_str(),
		      Form("pdfCTAUMASS_Jpsi_%s", (isPbPb?"PbPb":"PP")),
		      parIni[Form("N_Bkg_%s", (isPbPb?"PbPb":"PP"))].c_str(),
		      Form("pdfCTAUMASS_Bkg_%s", (isPbPb?"PbPb":"PP"))
		      ));
    }
    ws.pdf(Form("pdfMASS_Tot_%s", (isPbPb?"PbPb":"PP")))->setNormRange("MassWindow");
    ws.pdf(Form("pdfCTAUMASS_Tot_%s", (isPbPb?"PbPb":"PP")))->setNormRange("MassWindow");

  }

  ws.Print();
  return true;
};


void setDefaultParameters(map<string, string> &parIni, int nt, bool isPbPb)
{
  
  // MASS FIT PARAMETERS, TO BE REMOVED SOON
  parIni[Form("N_Jpsi_%s", (isPbPb?"PbPb":"PP"))]  = Form("%s[%d,%d]", Form("N_Jpsi_%s", (isPbPb?"PbPb":"PP")), 0, nt);
  parIni[Form("N_Psi2S_%s", (isPbPb?"PbPb":"PP"))] = Form("%s[%d,%d]", Form("N_Psi2S_%s", (isPbPb?"PbPb":"PP")), 0, nt);
  parIni[Form("N_Bkg_%s", (isPbPb?"PbPb":"PP"))]   = Form("%s[%d,%d]", Form("N_Bkg_%s", (isPbPb?"PbPb":"PP")), 0, nt);

  parIni[Form("m_Jpsi_%s", (isPbPb?"PbPb":"PP"))]      = Form("%s[%.4f,%.4f,%.4f]", Form("m_Jpsi_%s", (isPbPb?"PbPb":"PP")), Mass.JPsi, Mass.JPsi-0.1, Mass.JPsi+0.1);
  parIni[Form("sigma1_Jpsi_%s", (isPbPb?"PbPb":"PP"))] = Form("%s[%.4f,%.4f,%.4f]", Form("sigma1_Jpsi_%s", (isPbPb?"PbPb":"PP")), 0.04, 0.01, 0.09);
  parIni[Form("sigma2_Jpsi_%s", (isPbPb?"PbPb":"PP"))] = Form("%s[%.4f,%.4f,%.4f]", Form("sigma2_Jpsi_%s", (isPbPb?"PbPb":"PP")), 0.02, 0.01, 0.07);
  parIni[Form("alpha_Jpsi_%s", (isPbPb?"PbPb":"PP"))]  = Form("%s[%.4f,%.4f,%.4f]", Form("alpha_Jpsi_%s", (isPbPb?"PbPb":"PP")), 2.21, 0.5, 30.0*100.0);
  parIni[Form("n_Jpsi_%s", (isPbPb?"PbPb":"PP"))]      = Form("%s[%.4f,%.4f,%.4f]", Form("n_Jpsi_%s", (isPbPb?"PbPb":"PP")), 2., 0.5, 30.0*100.0);
  parIni[Form("f_Jpsi_%s", (isPbPb?"PbPb":"PP"))]      = Form("%s[%.4f,%.4f,%.4f]", Form("f_Jpsi_%s", (isPbPb?"PbPb":"PP")), 0.3, 0.0, 1.0);

  cout << "[INFO] Constraining Psi(2S) parameters to Jpsi using PDF Mass Ratio" << endl;
  Double_t MassRatio = (Mass.Psi2S/Mass.JPsi);
  parIni[Form("m_Psi2S_%s", (isPbPb?"PbPb":"PP"))]      = Form("RooFormulaVar::%s('@0*@1',{MassRatio[%.4f],%s})", Form("m_Psi2S_%s", (isPbPb?"PbPb":"PP")), 
							       MassRatio, Form("m_Jpsi_%s", (isPbPb?"PbPb":"PP") ));
  parIni[Form("sigma1_Psi2S_%s", (isPbPb?"PbPb":"PP"))] = Form("RooFormulaVar::%s('@0*@1',{MassRatio,%s})", Form("sigma1_Psi2S_%s", (isPbPb?"PbPb":"PP")), Form("sigma1_Jpsi_%s", (isPbPb?"PbPb":"PP") ));
  parIni[Form("sigma2_Psi2S_%s", (isPbPb?"PbPb":"PP"))] = Form("RooFormulaVar::%s('@0*@1',{MassRatio,%s})", Form("sigma2_Psi2S_%s", (isPbPb?"PbPb":"PP")), Form("sigma2_Jpsi_%s", (isPbPb?"PbPb":"PP") ));
  parIni[Form("alpha_Psi2S_%s", (isPbPb?"PbPb":"PP"))]  = Form("RooFormulaVar::%s('@0',{%s})", Form("alpha_Psi2S_%s", (isPbPb?"PbPb":"PP")), Form("alpha_Jpsi_%s", (isPbPb?"PbPb":"PP")));
  parIni[Form("n_Psi2S_%s", (isPbPb?"PbPb":"PP"))]      = Form("RooFormulaVar::%s('@0',{%s})", Form("n_Psi2S_%s", (isPbPb?"PbPb":"PP")), Form("n_Jpsi_%s", (isPbPb?"PbPb":"PP")));
  parIni[Form("f_Psi2S_%s", (isPbPb?"PbPb":"PP"))]      = Form("RooFormulaVar::%s('@0',{%s})", Form("f_Psi2S_%s", (isPbPb?"PbPb":"PP")), Form("f_Jpsi_%s", (isPbPb?"PbPb":"PP")));

  parIni["RFrac2Svs1S_PbPbvsPP"] = Form("%s[%.4f,%.4f,%.4f]", "RFrac2Svs1S_PbPbvsPP", 0.26, 0.0, 3.0);
  parIni[Form("RFrac2Svs1S_%s", (isPbPb?"PbPb":"PP"))] = Form("%s[%.4f,%.4f,%.4f]", Form("RFrac2Svs1S_%s", (isPbPb?"PbPb":"PP")), 0.26, 0.0, 1.0);

  parIni[Form("lambda1_Bkg_%s", (isPbPb?"PbPb":"PP"))] = Form("%s[%.4f,%.4f,%.4f]", Form("lambda1_Bkg_%s", (isPbPb?"PbPb":"PP")), 0.05, -1.0, 1.0);
  parIni[Form("lambda2_Bkg_%s", (isPbPb?"PbPb":"PP"))] = Form("%s[%.4f,%.4f,%.4f]", Form("lambda2_Bkg_%s", (isPbPb?"PbPb":"PP")), 0.05, -1.0, 1.0);
  parIni[Form("lambda3_Bkg_%s", (isPbPb?"PbPb":"PP"))] = Form("%s[%.4f,%.4f,%.4f]", Form("lambda3_Bkg_%s", (isPbPb?"PbPb":"PP")), 0.05, -1.0, 1.0);
  parIni[Form("lambda4_Bkg_%s", (isPbPb?"PbPb":"PP"))] = Form("%s[%.4f,%.4f,%.4f]", Form("lambda4_Bkg_%s", (isPbPb?"PbPb":"PP")), 0.05, -1.0, 1.0);

  // CTAU FIT PARAMETERS, TO BE REMOVED SOON
  parIni[Form("ctau1_ctauRes_%s", (isPbPb?"PbPb":"PP"))] = Form("%s[%.4f,%.4f,%.4f]", Form("ctau1_Res_%s", (isPbPb?"PbPb":"PP")), 0., -0.01, 0.01);
  parIni[Form("sigma1_ctauRes_%s", (isPbPb?"PbPb":"PP"))] = Form("%s[%.4f,%.4f,%.4f]", Form("sigma1_ctauRes_%s", (isPbPb?"PbPb":"PP")), 0.8, 0.6, 2.0);
  parIni[Form("ctau2_ctauRes_%s", (isPbPb?"PbPb":"PP"))] = Form("%s[%.4f,%.4f,%.4f]", Form("ctau2_ctauRes_%s", (isPbPb?"PbPb":"PP")), 0., -0.01, 0.01);
  parIni[Form("sigma2_ctauRes_%s", (isPbPb?"PbPb":"PP"))] = Form("%s[%.4f,%.4f,%.4f]", Form("sigma2_ctauRes_%s", (isPbPb?"PbPb":"PP")), 2.3, 1.1, 5.5);
  parIni[Form("f_ctauRes_%s", (isPbPb?"PbPb":"PP"))] = Form("%s[%.4f,%.4f,%.4f]", Form("f_ctauRes_%s", (isPbPb?"PbPb":"PP")), 0.05, 0., 1.);
  parIni[Form("f_JpsiNP_%s", (isPbPb?"PbPb":"PP"))] = Form("%s[%.4f,%.4f,%.4f]", Form("f_JpsiNP_%s", (isPbPb?"PbPb":"PP")), 0.2, 0., 1.);
  parIni[Form("f_Psi2SNP_%s", (isPbPb?"PbPb":"PP"))] = Form("%s[%.4f,%.4f,%.4f]", Form("f_Psi2SNP_%s", (isPbPb?"PbPb":"PP")), 0.2, 0., 1.);
  parIni[Form("f_BkgNP_%s", (isPbPb?"PbPb":"PP"))] = Form("%s[%.4f,%.4f,%.4f]", Form("f_BkgNP_%s", (isPbPb?"PbPb":"PP")), 0.3, 0., 1.);
  parIni[Form("lambdaDSS_BkgNP_%s", (isPbPb?"PbPb":"PP"))] = Form("%s[%.4f,%.4f,%.4f]", Form("lambdaDSS_BkgNP_%s", (isPbPb?"PbPb":"PP")), 0.42, 0.05, 1.5);
  parIni[Form("lambdaDF_BkgNP_%s", (isPbPb?"PbPb":"PP"))] = Form("%s[%.4f,%.4f,%.4f]", Form("lambdaDF_BkgNP_%s", (isPbPb?"PbPb":"PP")), 0.79, 0.001, 1.5);
  parIni[Form("lambdaDDS_BkgNP_%s", (isPbPb?"PbPb":"PP"))] = Form("%s[%.4f,%.4f,%.4f]", Form("lambdaDDS_BkgNP_%s", (isPbPb?"PbPb":"PP")), 0.69, 0.02, 5.0);
  parIni[Form("fDFSS_BkgNP_%s", (isPbPb?"PbPb":"PP"))] = Form("%s[%.4f,%.4f,%.4f]", Form("fDFSS_BkgNP_%s", (isPbPb?"PbPb":"PP")), 0.9, 0., 1.);
  parIni[Form("fDLIV_BkgNP_%s", (isPbPb?"PbPb":"PP"))] = Form("%s[%.4f,%.4f,%.4f]", Form("fDLIV_BkgNP_%s", (isPbPb?"PbPb":"PP")), 0.9, 0., 1.);
  parIni[Form("lambdaDSS_JpsiNP_%s", (isPbPb?"PbPb":"PP"))] = Form("%s[%.4f,%.4f,%.4f]", Form("lambdaDSS_JpsiNP_%s", (isPbPb?"PbPb":"PP")), 0.8, 0.01, 2.0);
  parIni[Form("lambdaDSS_Psi2SNP_%s", (isPbPb?"PbPb":"PP"))] = Form("%s[%.4f,%.4f,%.4f]", Form("lambdaDSS_Psi2SNP_%s", (isPbPb?"PbPb":"PP")), 0.8, 0.01, 2.0);
  
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


bool defineCtauResolModel(RooWorkspace& ws, CtauModel model, map<string,string> parIni, bool isPbPb) 
{ 
  cout << "[INFO] Implementing Ctau Resolution Model" << endl;
  
  switch(model) 
    {  
    case (CtauModel::SingleGaussianResolution):  
      if (!( parIni.count(Form("ctau1_ctauRes_%s", (isPbPb?"PbPb":"PP"))) && parIni.count(Form("sigma1_ctauRes_%s", (isPbPb?"PbPb":"PP"))) )) { 
	cout << Form("[ERROR] Initial parameters where not found for Single Gaussian Ctau Resolution Model in %s", (isPbPb?"PbPb":"PP")) << endl; return false; 
      }
      ws.factory(Form("GaussModel::%s(%s, %s, %s, one[1.0], %s)", Form("pdfCTAU_ctauRes_%s", (isPbPb?"PbPb":"PP")), "ctau", 
		      parIni[Form("ctau1_ctauRes_%s", (isPbPb?"PbPb":"PP"))].c_str(), 
		      parIni[Form("sigma1_ctauRes_%s", (isPbPb?"PbPb":"PP"))].c_str(),
		      "ctauErr"
		      ));
      cout << Form("[INFO] Single Gaussian Ctau Resolution PDF in %s included", (isPbPb?"PbPb":"PP")) << endl; break;
      
    case (CtauModel::DoubleGaussianResolution):  
      if (!( parIni.count(Form("ctau1_ctauRes_%s", (isPbPb?"PbPb":"PP"))) && parIni.count(Form("sigma1_ctauRes_%s", (isPbPb?"PbPb":"PP")))
	     && parIni.count(Form("ctau2_ctauRes_%s", (isPbPb?"PbPb":"PP"))) && parIni.count(Form("sigma2_ctauRes_%s", (isPbPb?"PbPb":"PP"))) 
	     && parIni.count(Form("f_ctauRes_%s", (isPbPb?"PbPb":"PP"))) )) { 
	cout << Form("[ERROR] Initial parameters where not found for Double Gaussian Ctau Resolution Model in %s", (isPbPb?"PbPb":"PP")) << endl; return false; 
      }
      ws.factory(Form("GaussModel::%s(%s, %s, %s, one[1.0], %s)", Form("pdfCTAUG1_ctauRes_%s", (isPbPb?"PbPb":"PP")), "ctau", 
		      parIni[Form("ctau1_ctauRes_%s", (isPbPb?"PbPb":"PP"))].c_str(), 
		      parIni[Form("sigma1_ctauRes_%s", (isPbPb?"PbPb":"PP"))].c_str(),
		      "ctauErr"
		      ));
      ws.factory(Form("GaussModel::%s(%s, %s, %s, one[1.0], %s)", Form("pdfCTAUG2_ctauRes_%s", (isPbPb?"PbPb":"PP")), "ctau", 
		      parIni[Form("ctau2_ctauRes_%s", (isPbPb?"PbPb":"PP"))].c_str(), 
		      parIni[Form("sigma2_ctauRes_%s", (isPbPb?"PbPb":"PP"))].c_str(),
		      "ctauErr"
		      ));
      ws.factory(Form("AddModel::%s({%s, %s}, {%s})", Form("pdfCTAU_ctauRes_%s", (isPbPb?"PbPb":"PP")), 
		      Form("pdfCTAUG1_ctauRes_%s", (isPbPb?"PbPb":"PP")), 
		      Form("pdfCTAUG2_ctauRes_%s", (isPbPb?"PbPb":"PP")),  
		      parIni[Form("f_ctauRes_%s", (isPbPb?"PbPb":"PP"))].c_str()
		      ));
      cout << Form("[INFO] Double Gaussian Ctau Resolution PDF in %s included", (isPbPb?"PbPb":"PP")) << endl; break;

    default :
      cout<< "[ERROR] Selected Ctau Resolution Model has not been implemented"<< endl; return false;
    }
  
  return true;
};

bool addBackgroundCtauModel(RooWorkspace& ws, string object, CtauModel model, map<string,string> parIni, bool isPbPb) 
{
  cout << Form("[INFO] Implementing %s Background Ctau Model", object.c_str()) << endl;
  
  switch(model) 
    {  
    case (CtauModel::TripleDecay): 
      if (!( parIni.count(Form("lambdaDSS_%s_%s", object.c_str(), (isPbPb?"PbPb":"PP"))) && parIni.count(Form("lambdaDF_%s_%s", object.c_str(), (isPbPb?"PbPb":"PP"))) 
	     && parIni.count(Form("lambdaDDS_%s_%s", object.c_str(), (isPbPb?"PbPb":"PP"))) && parIni.count(Form("fDFSS_%s_%s", object.c_str(), (isPbPb?"PbPb":"PP")))
	     && parIni.count(Form("fDLIV_%s_%s", object.c_str(), (isPbPb?"PbPb":"PP"))) )) { 
	cout << Form("[ERROR] Initial parameters where not found for %s Background Triple Decay Ctau Model in %s", object.c_str(), (isPbPb?"PbPb":"PP")) << endl; return false; 
      }
      ws.factory(Form("Decay::%s(%s, %s, %s, RooDecay::SingleSided)", Form("pdfCTAUDSS_%s_%s", object.c_str(), (isPbPb?"PbPb":"PP")), "ctau", 
		      parIni[Form("lambdaDSS_%s_%s", object.c_str(), (isPbPb?"PbPb":"PP"))].c_str(),
		      Form("pdfCTAU_ctauRes_%s", (isPbPb?"PbPb":"PP"))
		      ));
      ws.factory(Form("Decay::%s(%s, %s, %s, RooDecay::Flipped)", Form("pdfCTAUDF_%s_%s", object.c_str(), (isPbPb?"PbPb":"PP")), "ctau", 
		      parIni[Form("lambdaDF_%s_%s", object.c_str(), (isPbPb?"PbPb":"PP"))].c_str(),
		      Form("pdfCTAU_ctauRes_%s", (isPbPb?"PbPb":"PP"))
		      ));
      ws.factory(Form("Decay::%s(%s, %s, %s, RooDecay::DoubleSided)", Form("pdfCTAUDDS_%s_%s", object.c_str(), (isPbPb?"PbPb":"PP")), "ctau", 
		      parIni[Form("lambdaDDS_%s_%s", object.c_str(), (isPbPb?"PbPb":"PP"))].c_str(),
		      Form("pdfCTAU_ctauRes_%s", (isPbPb?"PbPb":"PP"))
		      ));
      ws.factory(Form("SUM::%s(%s*%s, %s)", Form("pdfCTAU1_%s_%s", object.c_str(), (isPbPb?"PbPb":"PP")), 
		      parIni[Form("fDFSS_%s_%s", object.c_str(), (isPbPb?"PbPb":"PP"))].c_str(),
		      Form("pdfCTAUDSS_%s_%s", object.c_str(), (isPbPb?"PbPb":"PP")),
		      Form("pdfCTAUDF_%s_%s", object.c_str(), (isPbPb?"PbPb":"PP"))
		      ));
      ws.factory(Form("SUM::%s(%s*%s, %s)", Form("pdfCTAU_%s_%s", object.c_str(), (isPbPb?"PbPb":"PP")), 
		      parIni[Form("fDLIV_%s_%s", object.c_str(), (isPbPb?"PbPb":"PP"))].c_str(),
		      Form("pdfCTAU1_%s_%s", object.c_str(), (isPbPb?"PbPb":"PP")),
		      Form("pdfCTAUDDS_%s_%s", object.c_str(), (isPbPb?"PbPb":"PP"))
		      ));

      cout << Form("[INFO] %s Background Triple Decay Ctau PDF in %s included", object.c_str(), (isPbPb?"PbPb":"PP")) << endl; break; 

    case (CtauModel::Delta): 
      ws.factory(Form("SUM::%s(%s)", Form("pdfCTAU_%s_%s", object.c_str(), (isPbPb?"PbPb":"PP")), 
		      Form("pdfCTAU_ctauRes_%s", (isPbPb?"PbPb":"PP"))
		      ));
		      
      cout << Form("[INFO] %s Background Delta Ctau PDF in %s included", object.c_str(), (isPbPb?"PbPb":"PP")) << endl; break; 

    default :
      cout<< "[ERROR] Selected Background Ctau Model has not been implemented"<< endl; return false;
    }
  
  return true;
};

bool addSignalCtauModel(RooWorkspace& ws, string object, CtauModel model, map<string,string> parIni, bool isPbPb) 
{
  cout << Form("[INFO] Implementing %s Signal Ctau Model", object.c_str()) << endl;
  
  switch(model) 
    {  
    case (CtauModel::SingleSidedDecay): 
      if (!( parIni.count(Form("lambdaDSS_%s_%s", object.c_str(), (isPbPb?"PbPb":"PP"))) )) { 
	cout << Form("[ERROR] Initial parameters where not found for %s Signal Single Sided Decay Ctau Model in %s", object.c_str(), (isPbPb?"PbPb":"PP")) << endl; return false; 
      }
      ws.factory(Form("Decay::%s(%s, %s, %s, RooDecay::SingleSided)", Form("pdfCTAU_%s_%s", object.c_str(), (isPbPb?"PbPb":"PP")), "ctau", 
		      parIni[Form("lambdaDSS_%s_%s", object.c_str(), (isPbPb?"PbPb":"PP"))].c_str(),
		      Form("pdfCTAU_ctauRes_%s", (isPbPb?"PbPb":"PP"))
		      ));

      cout << Form("[INFO] %s Signal Single Sided Decay Ctau PDF in %s included", object.c_str(), (isPbPb?"PbPb":"PP")) << endl; break; 

    case (CtauModel::Delta): 
      ws.factory(Form("SUM::%s(%s)", Form("pdfCTAU_%s_%s", object.c_str(), (isPbPb?"PbPb":"PP")), 
		      Form("pdfCTAU_ctauRes_%s", (isPbPb?"PbPb":"PP"))
		      ));
		      
      cout << Form("[INFO] %s Signal Delta Ctau PDF in %s included", object.c_str(), (isPbPb?"PbPb":"PP")) << endl; break; 

    default :
      cout<< "[ERROR] Selected Signal Ctau Model has not been implemented"<< endl; return false;
    }
  
  return true;
};
