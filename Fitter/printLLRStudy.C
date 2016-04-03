// a macro that simply extracts the ratio of Psi2S/Jpsi from a list of files containing fit results

#include "TFile.h"
#include "TMath.h"
#include "RooWorkspace.h"
#include "RooRealVar.h"
#include "RooAbsPdf.h"
#include <iostream>
#include <fstream>
#include <dirent.h>
#include <vector>
#include <map>
#include <set>

using namespace std;

// let's define a few useful small class
struct model_t {
   string fileName;
   string binName;
   string modelName;
   int npar;
   double nll;
   int cnt;

   bool operator<(model_t other) const { // needed for std::set
      if (npar<other.npar) return true;
      if (npar>other.npar) return false;
      return fileName<other.fileName;
   }
};
typedef vector<model_t> vecModels_t;
typedef set<model_t> setModels_t;



vector<string> printNLL(map< string, setModels_t > content, string dirPath, string type, string dirLabel) ;
void setLines(vector<string>& strLin, vector<string> lin); 
void printLines(vector<string> strLin, ofstream& fout); 
bool findFiles(string dirPath, vector<string>& fileNames); 
void splitString(string stringOriginal, const string Key, string& stringWithoutKey, string& stringWithKey); 
bool readFiles(string dirPath, vector<string> fileNames, map<string, setModels_t>& content, string type);
bool extractNLL(string fileName, model_t& value) ;


void printLLRStudy(
                   string dirLabel="BkgShapeStudy",  // Name of the working directory (currently hardcoded to work with DATA)
                   string type = "Bkg"               // Type of the LLR test, available options are: "Bkg" , "Jpsi" and "Psi2S"
                   ) 
{
  
  vector<string> fileNames;
  string dirPath = Form("./Output/%s/result/DATA/",dirLabel.c_str());
  if (!findFiles(dirPath, fileNames)) { return; } 

  cout << "[INFO] Creating " << ((type=="Bkg")?"Background":"Signal") << " Study summary!" << endl;
  
  // Group the files based on their background model
  map<string, setModels_t> content;
  if (!readFiles(dirPath, fileNames, content, type)) { return; }

  // Loop over each kinematic bin and compute the LLR/AIC tests
  vector<string> bestModelFiles = printNLL(content, dirPath, type, dirLabel); 
  
  cout << "[INFO] " << ((type=="Bkg")?"Background":"Signal") << " Study summary file done!" << endl; 
    
  cout << "The files for the best models are: " << endl;
  for (vector<string>::iterator it = bestModelFiles.begin(); it != bestModelFiles.end(); it++) {
     cout << *it << endl;
  }

};


vector<string> printNLL(map< string, setModels_t > content, string dirPath, string type, string dirLabel) 
{ 
  vector<string> ans;

  string outputFile = Form("%s/../LLRTest_%s.txt",dirPath.c_str(), type.c_str());
  ofstream fout( outputFile );
  string outputFileTexTable = Form("%s/../LLRTest_%s_TexTables.txt",dirPath.c_str(), type.c_str());
  ofstream foutTexTable( outputFileTexTable );
  map< string, setModels_t>::iterator contIt;

  for ( contIt = content.begin(); contIt != content.end(); contIt++) {
      
    string      binName = contIt->first;
    setModels_t binCont = contIt->second;
    
    cout << "Analzing Kinematic Bin: " << binName << endl;
    cout << " " << endl;
    fout << "Analzing Kinematic Bin: " << binName << endl;
    fout << " " << endl;

    vecModels_t modelNLLB;
    setModels_t::iterator modelRow;

    for (modelRow = binCont.begin(); modelRow != binCont.end(); modelRow++) {  

      vector<string> strLin;
      unsigned int iB=0;

      setModels_t::iterator modelCol;
      for (modelCol = binCont.begin(); modelCol != binCont.end(); modelCol++) {

        string  modelNameA = modelCol->modelName;
        int     nParA      = modelCol->npar;
        double  NLLA       = modelCol->nll;
        double  AICA       = 2*(nParA + NLLA);

        if (modelRow==binCont.begin()) {
          modelNLLB.push_back(*modelCol);
        }

        vector<string> lin;

        if (nParA>=modelRow->npar) {
          string  modelNameB = modelNLLB[iB].modelName;
          int     nParB      = modelNLLB[iB].npar;
          double  NLLB       = modelNLLB[iB].nll;
          double  AICB       = 2*(nParB + NLLB);
          if (modelNameA==modelNameB) {
            lin.push_back("|------------------------");
            lin.push_back("| "+modelNameA);
            lin.push_back(Form("|    NLL: %.2f  ", NLLA));
            lin.push_back(Form("|    AIC: %.2f  ", AICA));
            lin.push_back("|------------------------");
          } else if (nParA >= nParB) {
            double  diffNLL    = -2.0*(NLLA - NLLB);
            double  diffNPar   =  2.0*(nParA-nParB);
            double  probChi2   = 100.*TMath::Prob(diffNLL, diffNPar);
            if (diffNLL<0) probChi2 = 100.;
            if (probChi2>5. && (nParA-nParB)<=2) modelNLLB[iB].cnt++;
            lin.push_back("| "+modelNameA);
            lin.push_back(Form("|    NLL: %.2f  ", NLLA));
            lin.push_back(Form("|    Diff: %.2f  ", diffNLL));
            lin.push_back(Form("|    Prob: %.1f%s   ", probChi2, "%"));
            lin.push_back(Form("|    AIC: %.2f  ", -(AICA-AICB)));
            lin.push_back("|------------------------");
          } 
          setLines(strLin, lin); iB++;
        } 
      }
      for (unsigned int j=0; j<strLin.size(); j++) { strLin[j] = strLin[j] + "|"; }
      printLines(strLin, fout);
    }      

    // which is the best model for this bin?
    string bestModelFile="NOTFOUND"; int minok=999; int maxpar=0;
    for (vecModels_t::iterator it=modelNLLB.begin(); it!=modelNLLB.end();it++) {
       if (it->npar > maxpar) maxpar = it->npar;
       if (it->cnt>=2 && it->npar<minok) {
          bestModelFile=it->fileName;
          minok = it->npar;
       }
    }
    if (minok==999) { // sometimes the best model is one of the two highest orders...
       for (vecModels_t::iterator it=modelNLLB.begin(); it!=modelNLLB.end();it++) {
          int npar = it->npar;
          if (it->cnt>=maxpar-npar && npar<minok) {
             bestModelFile=it->fileName+" WARNING, HIGH ORDER";
             minok = it->npar;
          }
       }
    }

    cout << endl << " And the winner is... " << bestModelFile << endl << endl << endl;
    fout << endl << " And the winner is... " << bestModelFile << endl << endl << endl;
    ans.push_back(bestModelFile);

    vector<string> TexTable;
    TexTable.push_back("\\begin{table}");
    TexTable.push_back("\\centering");
    string iniLinLatex = "\\begin{tabular}{ c  c";
    string header = "N & NLL ";
    for (unsigned int iM=1; iM<=binCont.size(); iM++) { 
      if ((binCont.size()-iM)>=2) {
        iniLinLatex = iniLinLatex + " c ";
        header = header + "& " + Form("p(H0: N = %d)",iM) + " ";
      }
    }
    iniLinLatex = iniLinLatex + "}";
    header = header + "\\\\";
    TexTable.push_back(iniLinLatex);
    TexTable.push_back(header);
    TexTable.push_back("\\hline");
    
    vector<string> strLinLatex;
    unsigned int iR=0;
    int nParMax = 0;
    setModels_t::iterator model;
    for (model = binCont.begin(); model != binCont.end(); model++) { 
      if (model->npar > nParMax) { nParMax = model->npar; }
    }
    for (modelRow = binCont.begin(); modelRow != binCont.end(); modelRow++) {  

      unsigned int iA=0;

      setModels_t::iterator modelCol;
      for (modelCol = binCont.begin(); modelCol != binCont.end(); modelCol++) {

        string  modelNameA = modelCol->modelName;
        int     nParA      = modelCol->npar;
        double  NLLA       = modelCol->nll;
        double  AICA       = 2*(nParA + NLLA);

        if (modelRow==binCont.begin()) {
          modelNLLB.push_back(*modelCol);
          strLinLatex.push_back( Form("%d & %.2f", ++iA, NLLA) );
        }

        if (nParA<modelRow->npar) {
          int     nParRow = modelRow->npar;
          double  NLLRow  = modelRow->nll;
          if ((nParMax - nParA)>=2) {
            double  diffNLL    = -2.0*(NLLRow - NLLA);
            double  diffNPar   =  2.0*(nParRow-nParA);
            double  probChi2   = 100.*TMath::Prob(diffNLL, diffNPar);
            if (diffNLL<0) probChi2 = 100.;
            if (bestModelFile==modelCol->fileName && (nParRow-nParA)<=2) {
              strLinLatex[iR] = strLinLatex[iR] + " & " + Form("\\textbf{%.1f%s}", probChi2, "$\\%$");
            } else {
              strLinLatex[iR] = strLinLatex[iR] + " & " + Form("%.1f%s", probChi2, "$\\%$");
            } 
          }
        } else {
          if ((nParMax - nParA)>=2) {
            strLinLatex[iR] = strLinLatex[iR] + " & ";
          }
        } 
      }
      iR++;
    }      
    for (unsigned int j=0; j<strLinLatex.size(); j++) {  
      strLinLatex[j] = strLinLatex[j] + "\\\\";
      TexTable.push_back(strLinLatex[j]);
    }
    TexTable.push_back("\\end{tabular}");
    TexTable.push_back(Form("\\label{tab:LLRTEST_%s_%s}", dirLabel.c_str(), binName.c_str()));
    string rapStr, centStr, ptStr, modelStr, colStr;
    if (bestModelFile.find("ExpChebychev")!=std::string::npos){ modelStr = "exponential chebychev polynomials"; }
    else if (bestModelFile.find("Chebychev")!=std::string::npos){ modelStr = "chebychev polynomials"; }
    else { modelStr = "I DON'T KNOW"; }
    if (binName.find("rap016")!=std::string::npos){ rapStr = "$\\abs{y} <$ 1.6"; }
    else if (binName.find("rap1624")!=std::string::npos){ rapStr = "1.6 $\\leq \\abs{y} <$ 2.4"; }
    if (binName.find("pt30300")!=std::string::npos){ ptStr = "3.0 $\\leq \\PT <$ 30.0 $\\GeVc$"; }
    else if (binName.find("pt3065")!=std::string::npos){ ptStr = "3.0 $\\leq \\PT <$ 6.5 $\\GeVc$"; }
    else if (binName.find("pt6590")!=std::string::npos){ ptStr = "6.5 $\\leq \\PT <$ 9.0 $\\GeVc$"; }
    else if (binName.find("pt65120")!=std::string::npos){ ptStr = "6.5 $\\leq \\PT <$ 12.0 $\\GeVc$"; }
    else if (binName.find("pt65300")!=std::string::npos){ ptStr = "6.5 $\\leq \\PT <$ 30.0 $\\GeVc$"; }
    else if (binName.find("pt90120")!=std::string::npos){ ptStr = "9.0 $\\leq \\PT <$ 12.0 $\\GeVc$"; }
    else if (binName.find("pt120150")!=std::string::npos){ ptStr = "12.0 $\\leq \\PT <$ 15.0 $\\GeVc$"; }
    else if (binName.find("pt120300")!=std::string::npos){ ptStr = "12.0 $\\leq \\PT <$ 30.0 $\\GeVc$"; }
    else if (binName.find("pt150200")!=std::string::npos){ ptStr = "15.0 $\\leq \\PT <$ 20.0 $\\GeVc$"; }
    else if (binName.find("pt200300")!=std::string::npos){ ptStr = "20.0 $\\leq \\PT <$ 30.0 $\\GeVc$"; }
    if (binName.find("cent020")!=std::string::npos){ centStr = "centratility bin 0-10$\\%$"; }
    else if (binName.find("cent2040")!=std::string::npos){ centStr = "centratility bin 10-20$\\%$"; }
    else if (binName.find("cent4060")!=std::string::npos){ centStr = "centratility bin 20-30$\\%$"; }
    else if (binName.find("cent6080")!=std::string::npos){ centStr = "centratility bin 30-40$\\%$"; }
    else if (binName.find("cent80100")!=std::string::npos){ centStr = "centratility bin 40-50$\\%$"; }
    else if (binName.find("cent100200")!=std::string::npos){ centStr = "centratility bin 50-1000$\\%$"; }
    else if (binName.find("cent040")!=std::string::npos){ centStr = "centratility bin 0-20$\\%$"; }
    else if (binName.find("cent4080")!=std::string::npos){ centStr = "centratility bin 20-40$\\%$"; }
    else if (binName.find("cent80200")!=std::string::npos){ centStr = "centratility bin 40-100$\\%$"; }
    else if (binName.find("cent0200")!=std::string::npos){ centStr = "centratility bin 0-100$\\%$"; }
    if (binName.find("PP")!=std::string::npos){ colStr = "pp"; }
    else if (binName.find("PbPb")!=std::string::npos){ colStr = "PbPb"; }
    TexTable.push_back(Form("\\caption{Negative loglikelihoods for fits with %s of orders 1-5 of %s data in %s %s. In addition the p-values of the LLR-test for the null-hypothesis are listed. Tests of which the null-hypothesis cannot be rejected for two consecutive orders are highlighted in bold, together with the corresponding order.}", modelStr.c_str(), colStr.c_str(), rapStr.c_str(), (colStr=="pp" ? Form("and %s", ptStr.c_str()) : Form(", %s and %s", ptStr.c_str(), centStr.c_str()))));
    TexTable.push_back("\\end{table}");
    printLines(TexTable, foutTexTable);
    foutTexTable << endl; foutTexTable << endl;

  } // bin loop

  return ans;
};


void setLines(vector<string>& strLin, vector<string> lin) 
{
  const  string  empty  = "                          ";
  if (strLin.size() < lin.size()) {
    strLin = lin;  
  } else {
    for (unsigned int i = 0; i < lin.size(); i++) {
      strLin.at(i) = strLin.at(i) + lin.at(i);
    }
  }
  for (unsigned int i = 0; i < lin.size(); i++) {
    strLin.at(i).append(empty, 0, (25-lin.at(i).length()));
  }
    
  return;
};


void printLines(vector<string> strLin, ofstream& fout) 
{
  for (vector<string>::iterator line = strLin.begin(); line < strLin.end(); line++) {
      cout << *line << endl;
      fout << *line << endl;
  }
};


bool extractNLL(string fileName, model_t& value) 
{
  value.npar = -1;
  value.nll = 1e99;
  TFile *f = new TFile( fileName.c_str() );
  if (!f) {
    cout << "[Error] " << fileName << " not found" << endl; return false;
  }
  RooWorkspace *ws = (RooWorkspace*) f->Get("workspace");
  if (!ws) {
    f->Close(); delete f;
    cout << "[ERROR] Workspace not found in " << fileName << endl; return false;
  }    
  RooAbsReal *nll = NULL;
  double NLL = 0;
  int npar = 0;
  if (fileName.find("PP")!=std::string::npos) {
    nll = ws->pdf("pdfMASS_Tot_PP")->createNLL(*ws->data("dOS_DATA_PP"));
    npar = ws->pdf("pdfMASS_Bkg_PP")->getParameters(*ws->data("dOS_DATA_PP"))->getSize();
  }
  if (fileName.find("PbPb")!=std::string::npos) {
    nll = ws->pdf("pdfMASS_Tot_PbPb")->createNLL(*ws->data("dOS_DATA_PbPb"));
    npar = ws->pdf("pdfMASS_Bkg_PbPb")->getParameters(*ws->data("dOS_DATA_PbPb"))->getSize();
  }
  
  if (!nll) {
    delete ws;
    f->Close(); delete f;
    cout << "[ERROR] NLL was not found" << endl; return false;
  } else {
    NLL = nll->getVal();
  }

  value.npar = npar;
  value.nll = NLL;

  delete nll; delete ws;
  f->Close(); delete f;
        
  return true;
}



bool readFiles(string dirPath, vector<string> fileNames, map<string, setModels_t>& content, string type)
{
  for (vector<string>::iterator it = fileNames.begin(); it < fileNames.end(); it++) {
    string fileName = *it;

    string binName;   
    string modelName; 
    splitString(fileName, type, binName, modelName); 
       
    model_t modelNLL;
    modelNLL.binName = binName;
    modelNLL.modelName = modelName;
    modelNLL.fileName = fileName;
    modelNLL.cnt=0;
    if (extractNLL(dirPath+fileName, modelNLL) && modelName!="") {
      if (content.find(binName) == content.end()) content[binName] = setModels_t();
      content[binName].insert(modelNLL);
    } 
  }
  if (content.size()==0) {
    cout << "[ERROR] No NLL values were found in the input files" << endl; return false;
  }
  return true;
};

 
void splitString(string stringOriginal, const string Key, string& stringWithoutKey, string& stringWithKey) 
{
  string tmp  = stringOriginal;   

  if (tmp.find(Key+"_")==std::string::npos) {
    stringWithoutKey = "";
    stringWithKey = stringOriginal;
    return;
  }

  tmp.erase(0, tmp.find(Key+"_")+(Key+"_").length()); stringWithKey = tmp;
  stringWithKey.erase(stringWithKey.find("_"), stringWithKey.size());
      
  tmp.erase(0, tmp.find("_"));
  stringWithoutKey = stringOriginal;
  stringWithoutKey = stringWithoutKey.erase(stringWithoutKey.find(Key+"_")-1, stringWithoutKey.size())+tmp.erase(tmp.find(".root"), tmp.size());
};


bool findFiles(string dirPath, vector<string>& fileNames) 
{

  DIR *dpdf;
  struct dirent *epdf;

  // Open the working directory
  dpdf = opendir(dirPath.c_str());
  // Search for all the files inside this directory
  if (dpdf != NULL){
    while ( (epdf=readdir(dpdf)) ){
      if (strcmp(epdf->d_name,".")!=0 && strcmp(epdf->d_name,"..")!=0 ) {
        std::cout << "[INFO] Adding file: " << epdf->d_name << std::endl;
        string filePath = epdf->d_name; 
        fileNames.push_back(filePath);
      }
    }
  } else {
    cout << "[ERROR] Working directory was not found!" << endl; return false;
  }
  return true;

};
