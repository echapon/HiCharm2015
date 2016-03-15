// a macro that simply extracts the ratio of Psi2S/Jpsi from a list of files containing fit results

#include "TFile.h"
#include "TMath.h"
#include "RooWorkspace.h"
#include "RooRealVar.h"
#include "RooAbsPdf.h"
#include <iostream>
#include <fstream>
#include <dirent.h>

using namespace std;


bool printNLL(map< string, map<pair<int, string>,  double> > content, string outputFile) ;
void setLines(vector<string>& strLin, vector<string> lin); 
void printLines(vector<string> strLin, ofstream& fout); 
bool findFiles(string dirPath, vector<string>& fileNames); 
void splitString(string stringOriginal, const string Key, string& stringWithoutKey, string& stringWithKey); 
bool readFiles(string dirPath, vector<string> fileNames, map< string, map<pair<int, string>,  double> >& content, string type);
bool extractNLL(string fileName, pair<int,double>& value) ;


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
  map< string, map<pair<int, string>,  double> > content;
  if (!readFiles(dirPath, fileNames, content, type)) { return; }

  // Loop over each kinematic bin and compute the LLR/AIC tests
  string outputFile = Form("%s/../LLRTest_%s.txt",dirPath.c_str(), type.c_str());
  printNLL(content, outputFile); 
  
  cout << "[INFO] " << ((type=="Bkg")?"Background":"Signal") << " Study summary file done!" << endl; 
    

};


bool printNLL(map< string, map<pair<int, string>,  double> > content, string outputFile) 
{ 

  ofstream fout( outputFile );
  map< string, map<pair<int, string>,  double> >::iterator contIt;

  for ( contIt = content.begin(); contIt != content.end(); contIt++) {
      
    string                          binName = contIt->first;
    map<pair<int, string>,  double> binCont = contIt->second;
    
    cout << "Analzing Kinematic Bin: " << binName << endl;
    cout << " " << endl;
    fout << "Analzing Kinematic Bin: " << binName << endl;
    fout << " " << endl;
      
    vector<pair<pair<int, string>, double>> modelNLLB;
    map<pair<int, string>,  double>::iterator modelNLL;
    for (modelNLL = binCont.begin(); modelNLL != binCont.end(); modelNLL++) {  

      vector<string> strLin; int i = 0;
      map<pair<int, string>,  double>::iterator modelNLLA;
      for (modelNLLA = modelNLL; modelNLLA != binCont.end(); modelNLLA++) {
                   
        pair<int, string> keyA = modelNLLA->first;
        string  modelNameA = keyA.second;
        int     nParA      = keyA.first;
        double  NLLA       = modelNLLA->second;
        double  AIC        = 2*(nParA + NLLA);
          
        if (modelNLL==binCont.begin()) {
          modelNLLB.push_back(make_pair(modelNLLA->first, modelNLLA->second));

        }

        pair<int, string> keyB = modelNLLB[i].first;
        string  modelNameB = keyB.second;
        int     nParB      = keyB.first;
        double  NLLB       = modelNLLB[i].second;

        vector<string> lin;
        if (modelNameA==modelNameB) {
          lin.push_back("|------------------------");
          lin.push_back("| "+modelNameA);
          lin.push_back(Form("|    NLL: %.4f", NLLA));
          lin.push_back(Form("|    AIC: %.4f", AIC));
          lin.push_back("|------------------------");

        } else if (nParA >= nParB) {
          double diffNLL  = -2.0*(NLLA - NLLB);
          double diffNPar =  2.0*(nParA-nParB);
          double probChi2 = TMath::Prob(diffNLL, diffNPar);
        
          lin.push_back("| "+modelNameA);
          lin.push_back(Form("|    NLL: %.4f", NLLA));
          lin.push_back(Form("|    Diff: %.4f", diffNLL));
          lin.push_back(Form("|    Prob: %.6f", probChi2));
          lin.push_back(Form("|    AIC: %.4f", AIC));
          lin.push_back("|------------------------");
        }
        setLines(strLin, lin);
        i = i + 1;
 
      }
      printLines(strLin, fout); 
    }
    cout << " " << endl; cout << " " << endl; cout << " " << endl;
    fout << " " << endl; fout << " " << endl; fout << " " << endl;
  }

  return true;
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


bool extractNLL(string fileName, pair<int,double>& value) 
{
  value = make_pair( -1 , -9999999.9 );
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

  value = make_pair( npar , NLL );

  delete nll; delete ws;
  f->Close(); delete f;
        
  return true;
}



bool readFiles(string dirPath, vector<string> fileNames, map< string, map<pair<int, string>,  double> >& content, string type)
{
  for (vector<string>::iterator it = fileNames.begin(); it < fileNames.end(); it++) {
    string fileName = *it;

    string binName;   
    string modelName; 
    splitString(fileName, type, binName, modelName); 
       
    pair< int , double > modelNLL;
    if (extractNLL(dirPath+fileName, modelNLL) && modelName!="") {
      int     nPar = modelNLL.first;
      double  NLL  = modelNLL.second;
      content[binName][make_pair(nPar, modelName)] = NLL;
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
