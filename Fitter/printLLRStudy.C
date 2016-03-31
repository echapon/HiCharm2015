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



vector<string> printNLL(map< string, setModels_t > content, string outputFile) ;
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
  string outputFile = Form("%s/../LLRTest_%s.txt",dirPath.c_str(), type.c_str());
  vector<string> bestModelFiles = printNLL(content, outputFile); 
  
  cout << "[INFO] " << ((type=="Bkg")?"Background":"Signal") << " Study summary file done!" << endl; 
    
  cout << "The files for the best models are: " << endl;
  for (vector<string>::iterator it = bestModelFiles.begin(); it != bestModelFiles.end(); it++) {
     cout << *it << endl;
  }

};


vector<string> printNLL(map< string, setModels_t > content, string outputFile) 
{ 
  vector<string> ans;

  ofstream fout( outputFile );
  map< string, setModels_t>::iterator contIt;

  for ( contIt = content.begin(); contIt != content.end(); contIt++) {
      
    string      binName = contIt->first;
    setModels_t binCont = contIt->second;
    
    cout << "Analzing Kinematic Bin: " << binName << endl;
    cout << " " << endl;
    fout << "Analzing Kinematic Bin: " << binName << endl;
    fout << " " << endl;

    vecModels_t modelNLLB;
    setModels_t::iterator it;
    for (it = binCont.begin(); it != binCont.end(); it++) {  

      model_t modelNLL = *it;
      vector<string> strLin; int i = 0;
      setModels_t::iterator modelNLLA;
      for (modelNLLA = it; modelNLLA != binCont.end(); modelNLLA++) {

        string  modelNameA = modelNLLA->modelName;
        int     nParA      = modelNLLA->npar;
        double  NLLA       = modelNLLA->nll;
        double  AICA       = 2*(nParA + NLLA);

        if (it==binCont.begin()) modelNLLB.push_back(*modelNLLA);

        string  modelNameB = modelNLLB[i].modelName;
        int     nParB      = modelNLLB[i].npar;
        double  NLLB       = modelNLLB[i].nll;
        double  AICB       = 2*(nParB + NLLB);

        vector<string> lin;
        if (modelNameA==modelNameB) {
          lin.push_back("|------------------------");
          lin.push_back("| "+modelNameA);
          lin.push_back(Form("|    NLL: %.2f  ", NLLA));
          lin.push_back(Form("|    AIC: %.2f  ", AICA));
          lin.push_back("|------------------------");

        } else if (nParA >= nParB) {
          double diffNLL  = -2.0*(NLLA - NLLB);
          double diffNPar =  2.0*(nParA-nParB);
          double probChi2 = 100.*TMath::Prob(diffNLL, diffNPar);
          if (diffNLL<0) probChi2 = 100.;
          if (probChi2>5.) modelNLLB[i].cnt++;
          else modelNLLB[i].cnt--; // put a penalty when the higher order model is better
        
          lin.push_back("| "+modelNameA);
          lin.push_back(Form("|    NLL: %.2f  ", NLLA));
          lin.push_back(Form("|    Diff: %.2f  ", diffNLL));
          lin.push_back(Form("|    Prob: %.1f%s   ", probChi2, "%"));
          lin.push_back(Form("|    AIC: %.2f  ", -(AICA-AICB)));
          lin.push_back("|------------------------");
        }
        setLines(strLin, lin);
        i = i + 1;
 
      }
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
    cout << *line << "|" << endl;
    fout << *line << "|" << endl;
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
