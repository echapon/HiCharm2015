// a macro that simply extracts the ratio of Psi2S/Jpsi from a list of files containing fit results

#include "TFile.h"
#include "TMath.h"
#include "RooWorkspace.h"
#include "RooRealVar.h"
#include <iostream>
#include <fstream>
#include <dirent.h>

using namespace std;


bool extractParameter(string fileName, const char* parName, pair<double,double>& value);


void printPsi2SvsJpsiRatioSummary(string dirPATH="./Output/BinStudy/result/DATA/",string saveLabel="TestBinning") {
  
  vector<string> filenames;
  DIR *dpdf;
  struct dirent *epdf;

  // Open the working directory
  dpdf = opendir(dirPATH.c_str());
  // Search for all the files inside this directory
  if (dpdf != NULL){
    while ( (epdf=readdir(dpdf)) ){
      if (strcmp(epdf->d_name,".")!=0 && strcmp(epdf->d_name,"..")!=0 ) {
        std::cout << "[INFO] Adding file: " << epdf->d_name << std::endl;
        filenames.push_back(epdf->d_name);
      }
   }
  } else {
    cout << "[ERROR] Working directory was not found!" << endl; return;
  }

  cout << "[INFO] Creating Psi2S/Jpsi summary!" << endl;

  // Group the files based on PP and PbPb
  map< string, map< string , string > > listFilesNames;
  for (vector<string>::iterator it = filenames.begin(); it < filenames.end(); it++) {
    if (it->find("PbPb")!=std::string::npos) {
      string fileTAG = *it; 
      fileTAG.erase(fileTAG.find("_PbPb"), string("_PbPb").length());
      listFilesNames[fileTAG]["PbPb"] = (dirPATH + *it);
    } 
    else if (it->find("PP")!=std::string::npos) {
      string fileTAG = *it; 
      fileTAG.erase(fileTAG.find("_PP"), string("_PP").length());
      listFilesNames[fileTAG]["PP"] = (dirPATH + *it);
    } 
    else {
      cout << "[ERROR] The input file name does not contain the word PbPb or PP!" << endl; return;
    }
  }

  vector<string> lines;
  lines.push_back( "BIN ; (Psi(2S)/Jpsi)_PP ; (Psi(2S)/Jpsi)_PbPb; (Psi(2S)/Jpsi)_(PbPb/PP)"); 

  // Loop over the map of files 
  for ( map<string, map<string,string> >::iterator it = listFilesNames.begin(); it != listFilesNames.end(); it++) {

    pair<double,double> RFrac2Svs1S_PbPbvsPP;
    pair<double,double> RFrac2Svs1S_PbPb;
    pair<double,double> RFrac2Svs1S_PP;

    map<string,string> file = it->second;
    
    if (file.count("PbPb")>0) {
      string fileName = file["PbPb"]; 
      if (!extractParameter(fileName, "RFrac2Svs1S_PbPbvsPP", RFrac2Svs1S_PbPbvsPP)) { return; }
      if (!extractParameter(fileName, "RFrac2Svs1S_PbPb", RFrac2Svs1S_PbPb)) { return; }
    }
    if (file.count("PP")>0) {
      string fileName = file["PP"]; 
      if (!extractParameter(fileName, "RFrac2Svs1S_PP", RFrac2Svs1S_PP)) { return; }
    }

    // If there is a missing parameter, compute it using the other two parameters
    if ( RFrac2Svs1S_PbPbvsPP.first==-999.9 && (RFrac2Svs1S_PbPb.first!=-999.9 && RFrac2Svs1S_PP.first!=-999.9) ){
      
      double SR_PbPb_Value = RFrac2Svs1S_PbPb.first;
      double SR_PbPb_Error = RFrac2Svs1S_PbPb.second;
      double SR_PP_Value   = RFrac2Svs1S_PP.first;
      double SR_PP_Error   = RFrac2Svs1S_PP.second;

      double DR_Value = ( SR_PbPb_Value / SR_PP_Value );
      // We assume no correlation between the single ratios in PbPb and PP
      double DR_Error = TMath::Sqrt( 
                                    TMath::Power( ( SR_PbPb_Error / SR_PP_Value ) , 2.0 ) + 
                                    TMath::Power( ( (SR_PbPb_Value*SR_PP_Error)  / (SR_PP_Value*SR_PP_Value) ) , 2.0 ) 
                                     );
      
      RFrac2Svs1S_PbPbvsPP = make_pair( DR_Value , DR_Error );
    }
    else if ( RFrac2Svs1S_PbPb.first==-999.9 && (RFrac2Svs1S_PbPbvsPP.first!=-999.9 && RFrac2Svs1S_PP.first!=-999.9) ){
      
      double DR_Value      = RFrac2Svs1S_PbPbvsPP.first;
      double DR_Error      = RFrac2Svs1S_PbPbvsPP.second; 
      double SR_PP_Value   = RFrac2Svs1S_PP.first;
      double SR_PP_Error   = RFrac2Svs1S_PP.second;

      double SR_PbPb_Value = ( DR_Value * SR_PP_Value ); 
      // We assume no correlation between the Double Ratio and the Single Ratio in PP
      double SR_PbPb_Error = TMath::Sqrt( 
                                         TMath::Power( ( DR_Error * SR_PP_Value ) , 2.0 ) + 
                                         TMath::Power( ( DR_Value * SR_PP_Error ) , 2.0 ) 
                                          );
      
      RFrac2Svs1S_PbPb = make_pair( SR_PbPb_Value , SR_PbPb_Error );
    } 
    else if ( RFrac2Svs1S_PbPb.first!=-999.9 && RFrac2Svs1S_PbPbvsPP.first!=-999.9 && RFrac2Svs1S_PP.first!=-999.9 ){
      cout << "[INFO] All parameters where found! " << endl;  
    }
    else {
      cout << "[ERROR] No parameters where found!" << endl; return;
    }

    // Clean the name of your file ( maybe later we can parse it! )
    string cleanFileName = it->first;
    cleanFileName.erase( cleanFileName.find("FIT_"), string("FIT_").length());
    cleanFileName.erase( cleanFileName.find("Psi2SJpsi_"), string("Psi2SJpsi_").length());
    cleanFileName.erase( cleanFileName.find(".root"), string(".root").length());

    // Save the output line for this bin
    lines.push_back( Form("%s ; %.5f±%.5f(%.2f %%) ; %.5f±%.5f(%.2f %%) ; %.5f±%.5f(%.2f %%) ", cleanFileName.c_str(), RFrac2Svs1S_PP.first , RFrac2Svs1S_PP.second , (RFrac2Svs1S_PP.second/RFrac2Svs1S_PP.first)*100. , RFrac2Svs1S_PbPb.first , RFrac2Svs1S_PbPb.second , (RFrac2Svs1S_PbPb.second/RFrac2Svs1S_PbPb.first)*100. , RFrac2Svs1S_PbPbvsPP.first , RFrac2Svs1S_PbPbvsPP.second , (RFrac2Svs1S_PbPbvsPP.second/RFrac2Svs1S_PbPbvsPP.first)*100.) );

  }

  // Write all the lines over the output file
  string saveName= Form("../RatiosSummary_%s.csv",saveLabel.c_str());
  ofstream fout( (dirPATH+saveName) );
  for (vector<string>::iterator it = lines.begin(); it < lines.end(); it++) {
    fout << it->c_str() << endl;
  } 

  cout << "[INFO] Psi2S/Jpsi summary file done!" << endl; 
    

}


bool extractParameter(string fileName, const char* parName, pair<double,double>& value) 
{
  TFile *f = new TFile( fileName.c_str() );
  if (!f) {
    cout << "[Error] " << fileName << " not found" << endl; return false;
  }
  RooWorkspace *ws = (RooWorkspace*) f->Get("workspace");
  if (!ws) {
    cout << "[ERROR] Workspace not found in " << fileName << endl; return false;
  }    
  RooRealVar *var = ws->var(parName);
  if (!var) {
    value = make_pair( -999.9 , -999.9 ); return true;
  }
  
  value = make_pair( var->getValV() , var->getError() );
        
  return true;
}
