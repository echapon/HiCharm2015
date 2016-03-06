#include "Macros/Utilities/initClasses.h"
#include "Macros/tree2DataSet.C"
#include "Macros/fitCharmonia.C"


bool parseFile(string FileName, vector< map<string, string> >& data);
bool parseString(string input, string delimiter, vector<double>& output);

bool iniWorkEnv(map<string,string>& DIR, const string workDirName);
bool existDir(string dir);
bool readFile(string FileName, vector< vector<string> >& content, const int nCol=-1);
bool getInputFileNames(const string InputTrees, map<string, vector<string> >& InputFileCollection);
bool setParameters(map<string, string> row,  string mass, string ctau, struct KinCuts& cut, map<string, string>& parIni);
bool addParamters(string InputFile,  vector< struct KinCuts > cutVector, vector< map<string, string> >&  parIniVector);


void fitter(
            const string workDirName="Test", // Working directory
            bool isPbPb      = false,        // isPbPb = false for pp, true for PbPb
            bool doSimulFit  = true,         // Do simultaneous fit
            bool inExcStat   = false,        // if inExcStat is true, then the excited states are fitted
            bool zoomPsi     = false,        // Zoom Psi(2S) peak on extra pad
            bool setLogScale = true,         // Draw plot with log scale
            bool incSS       = false,        // Include Same Sign data
            int  nbins       = 74,           // Number of bins used for fitting
            int  numCores    = 2             // Number of cores used for fitting
            ) 
{
  // -------------------------------------------------------------------------------
  // STEP 0: INITIALIZE THE FITTER WORK ENVIROMENT
  // The work enviroment is divided as follows:
  /*
    main |-> Macros: Contain all the macros
         |-> Input   |-> <WorkDir> : Contain Input File, Bin and Parameter List for a given work directory (e.g. 20160201)
	 |-> Output  |-> <WorkDir> : Contain Output Plots and Results for a given work directory (e.g. 20160201)
	 |-> DataSet : Contain all the datasets (MC and Data)
  */

  map<string,string> DIR;
  if(!iniWorkEnv(DIR, workDirName)){ return; } 
 
  // -------------------------------------------------------------------------------
  // STEP 1: CREATE/LOAD THE ROODATASETS
  /*
    Input : List of TTrees with format:  TAG <tab> FILE_NAME
    Output: Collection of RooDataSets splitted by tag name, including OS and SS dimuons.
  */

  const string InputTrees = DIR["input"] + "InputTrees.txt";
  map<string, vector<string> > InputFileCollection;
  if(!getInputFileNames(InputTrees, InputFileCollection)){ return; }
  map<string, RooWorkspace> Workspace;
  for(map<string, vector<string> >::iterator FileCollection=InputFileCollection.begin(); FileCollection!=InputFileCollection.end(); ++FileCollection) {
    string         TAG            = FileCollection->first; 
    vector<string> InputFileNames = FileCollection->second; 
    string         OutputFileName;
    OutputFileName= DIR["dataset"] + "DATASET_" + TAG + ".root"; 
    if(TAG.find("DATA")!=std::string::npos) { 
      if(!tree2DataSet(Workspace["DATA"], InputFileNames, TAG, OutputFileName)){ return; }
    }
    if(TAG.find("MC")!=std::string::npos) { 
      if(!tree2DataSet(Workspace["MC"],   InputFileNames, TAG, OutputFileName)){ return; } 
    }
  }

  // -------------------------------------------------------------------------------
  // STEP 2: LOAD THE INITIAL PARAMETERS
  /*
    Input : List of initial parameters with format PT <tab> RAP <tab> CEN <tab> iniPar ... 
    Output: two vectors with one entry per kinematic bin filled with the cuts and initial parameters
  */

  string InputFile;
  vector< struct KinCuts >       cutVector;
  vector< map<string, string> >  parIniVector;
  vector< map<string, string> > data;

  // Add initial parameters for PbPb background models
  InputFile = (DIR["input"] + "InitialParam_MASS_BKG_PbPb.csv");
  data.clear(); 
  if(!parseFile(InputFile, data)) { return; }
  for(vector< map<string, string> >::iterator row=data.begin(); row!=data.end(); ++row) {
    struct KinCuts cut; map<string, string> parIni;
    if(!setParameters(*row, "2.2-4.5", "-100.0->100.0", cut, parIni)) { return; }
    cutVector.push_back(cut);  parIniVector.push_back(parIni);
  }

  // Add initial parameters for PP background models
  InputFile = (DIR["input"] + "InitialParam_MASS_BKG_PP.csv");
  if (!addParamters(InputFile, cutVector, parIniVector)) { return; }

  // Add initial parameters for PbPb jpsi models
  InputFile = (DIR["input"] + "InitialParam_MASS_JPSI_PbPb.csv");
  if (!addParamters(InputFile, cutVector, parIniVector)) { return; }

  // Add initial parameters for PP jpsi models
  InputFile = (DIR["input"] + "InitialParam_MASS_JPSI_PP.csv");
  if (!addParamters(InputFile, cutVector, parIniVector)) { return; }

  // Add initial parameters for PbPb psi(2S) models
  InputFile = (DIR["input"] + "InitialParam_MASS_PSI2S_PbPb.csv");
  if (!addParamters(InputFile, cutVector, parIniVector)) { return; }

  // Add initial parameters for PP jpsi models
  InputFile = (DIR["input"] + "InitialParam_MASS_PSI2S_PP.csv");
  if (!addParamters(InputFile, cutVector, parIniVector)) { return; }
  
  unsigned int nBins = data.size();

  // -------------------------------------------------------------------------------  
  // STEP 3: FIT THE DATASETS
  /*
    Input : 
              -> The cuts and initial parameters per kinematic bin
	      -> The workspace with the full datasets included.
    Output: 
              -> Plots (png, pdf and root format) of each fit.
	      -> The local workspace used for each fit.
  */
  for (unsigned int i=0; i<nBins; i++) {
    string outputDir = DIR["output"];
    // DO SOMETHING WITH DATA
    if (Workspace.count("DATA")>0) { 
      // DATA datasets were loaded
      if (doSimulFit) {
        // If do simultaneous fits, then just fits once
        if (!fitCharmonia( Workspace["DATA"], cutVector.at(i), parIniVector.at(i), outputDir, "DATA",
                           isPbPb,           // isPbPb = false for pp, true for PbPb
                           zoomPsi,          // Zoom Psi(2S) peak on extra pad
                           setLogScale,      // Draw plot with log scale
                           incSS,            // Include Same Sign data
                           false,            // Compute the mean PT (NEED TO FIX)
                           inExcStat,        // if inExcStat is true, then the excited states are fitted
                           true,             // Do simultaneous fit
                           nbins,            // number of bins
                           numCores 
                           )
            ) { return; } 
      }
      else {
        // If don't want simultaneous fits, then fit first PbPb and then PP separately
        if (!fitCharmonia( Workspace["DATA"], cutVector.at(i), parIniVector.at(i), outputDir, "DATA",
                           true,             // isPbPb = false for pp, true for PbPb
                           zoomPsi,          // Zoom Psi(2S) peak on extra pad
                           setLogScale,      // Draw plot with log scale
                           incSS,            // Include Same Sign data
                           false,            // Compute the mean PT (NEED TO FIX)
                           inExcStat,        // if inExcStat is true, then the excited states are fitted
                           false,            // Do simultaneous fit
                           nbins,            // number of bins
                           numCores 
                           )
            ) { return; } 
        if (!fitCharmonia( Workspace["DATA"], cutVector.at(i), parIniVector.at(i), outputDir, "DATA",
                           false,            // isPbPb = false for pp, true for PbPb
                           zoomPsi,          // Zoom Psi(2S) peak on extra pad
                           setLogScale,      // Draw plot with log scale
                           incSS,            // Include Same Sign data
                           false,            // Compute the mean PT (NEED TO FIX)
                           inExcStat,        // if inExcStat is true, then the excited states are fitted
                           false,            // Do simultaneous fit
                           nbins,            // number of bins
                           numCores 
                           )
            ) { return; } 
      }
    }
  }
};

bool addParamters(string InputFile,  vector< struct KinCuts > cutVector, vector< map<string, string> >&  parIniVector)
{
  vector< map<string, string> > data;
  if(!parseFile(InputFile, data)) { return false; }
  if (data.size()!=cutVector.size()) { cout << "[ERROR] The initial parameters in file " << InputFile << " are not consistent with previous files!" << endl; return false; }
  for (unsigned int i=0; i<data.size(); i++) {
    struct KinCuts cut;
    if (!setParameters(data.at(i), "2.2-4.5", "-100.0->100.0", cut, parIniVector.at(i))) { return false; };
    if (!isEqualKinCuts(cut, cutVector.at(i))) { cout << "[ERROR] The bins in file " << InputFile << " are not consistent with previous files!" << endl; return false; }
  }
  return true;
};

bool setParameters(map<string, string> row,  string mass, string ctau, struct KinCuts& cut, map<string, string>& parIni)
{

  // set initial paramters
  cut.sMuon.Pt.Min  =  0.0;    
  cut.sMuon.Pt.Max  = 100000.0;
  cut.sMuon.Eta.Min = -2.4;   
  cut.sMuon.Eta.Max = 2.4;
  cut.dMuon.ctauErr.Min = 0.0; 
  cut.dMuon.ctauErr.Max = 10.0;
  cut.dMuon.ctau.Min = -3.0;   
  cut.dMuon.ctau.Max = 5.0;    
  cut.dMuon.M.Min = 2.0; 
  cut.dMuon.M.Max = 5.0;  
  cut.dMuon.AbsRap.Min = 0.0;
  cut.dMuon.AbsRap.Max = 2.4;
  cut.dMuon.Pt.Min  =  0.0;
  cut.dMuon.Pt.Max  =  1000.0;
  cut.Centrality.Start = 0;
  cut.Centrality.End = 200;

  // set paramters from file
  row["mass"] = mass;
  row["ctau"] = ctau;
  for(map<string, string>::iterator col=row.begin(); col!=row.end(); ++col) {
    string label = col->first;
    if (label=="rap") {
      std::vector<double> v; 
      if(!parseString(col->second, "-", v)) { return false; }
      cut.dMuon.AbsRap.Min = v.at(0); 
      cut.dMuon.AbsRap.Max = v.at(1);
    } 
    else if (label=="pt"){
      std::vector<double> v; 
      if(!parseString(col->second, "-", v)) { return false; }
      cut.dMuon.Pt.Min = v.at(0); 
      cut.dMuon.Pt.Max = v.at(1);
    } 
    else if (label=="cent"){
      std::vector<double> v; 
      if(!parseString(col->second, "-", v)) { return false; }
      cut.Centrality.Start = (int) (2.0*v.at(0)); 
      cut.Centrality.End = (int) (2.0*v.at(1));
    } 
    else if (label=="mass"){
      std::vector<double> v; 
      if(!parseString(col->second, "-", v)) { return false; }
      cut.dMuon.M.Min = v.at(0); 
      cut.dMuon.M.Max = v.at(1);
    } 
    else if (label=="ctau"){
      std::vector<double> v; 
      if(!parseString(col->second, "->", v)) { return false; }
      cut.dMuon.ctau.Min = v.at(0); 
      cut.dMuon.ctau.Max = v.at(1);
    } 
    else if (label.find("Model")!=std::string::npos){ 
      parIni[col->first] = col->second;
    } 
    else {
      if (col->second != "") {
	string value = col->second;
	// check that initial parameters format is correct: [ num, num num ]
	if ((value.find("[")==std::string::npos)||(value.find("]")==std::string::npos)) {
	  // Special cases like parameter constrains could be set here but for now, let's keep it simple

	  cout << "[ERROR] Either ']' or '[' are missing in the initial parameter values!" << endl; return false;
	} else {
	  value.erase(value.find("["), string("[").length());
	  value.erase(value.find("]"), string("]").length());
	}
	std::vector<double> v; 
	if(!parseString(value, ",", v)) { return false; }
	// everything seems alright, then proceed to save the values
	if (v.size()==1) { 
	  // if only one value is given i.e. [ num ], consider it a constant value
	  if (col->first.find("N_")!=std::string::npos){
	    parIni[col->first] = Form("%s[ %.0f, %.0f, %.0f]", col->first.c_str(), v.at(0), v.at(0), v.at(0));
	  } else {
	    parIni[col->first] = Form("%s[ %.6f, %.6f, %.6f]", col->first.c_str(), v.at(0), v.at(0), v.at(0));
	  }
	} else {
	  parIni[col->first] = col->first + col->second;
	}
      }
    }
  }

  return true;
};

bool parseString(string input, string delimiter, vector<double>& output)
{
  // remove spaces from input string 
  input.erase(std::remove(input.begin(), input.end(), ' '), input.end());
  // proceed to parse input string
  char *end;
  while(input!="") {
    double d = strtod(input.c_str(), &end);
    if (end != input) {
      output.push_back(d);
    } else {
      cout << "[ERROR] The conversion from string to double failed!"; return false;
    }
    input = end; 
    if(input.find(delimiter.c_str())!= std::string::npos){ input.erase(input.find(delimiter.c_str()), delimiter.length()); }
  }
  return true;
}

bool parseFile(string FileName, vector< map<string, string> >& data)
{
  vector< vector<string> > content; 
  if(!readFile(FileName, content)){ return false; }
  vector<string> header = content.at(0);
  if (header.size()==0) { cout << "[ERROR] The header is null!" << endl; return false; }
  for(vector<string>::iterator rHeader=header.begin(); rHeader!=header.end(); ++rHeader) {
    if (*rHeader=="") { cout << "[ERROR] A column has no label!" << endl; return false; }
  }

  for(vector< vector<string> >::iterator row=content.begin()+1; row!=content.end(); ++row) {
    map<string, string> col;
    for (unsigned int i=0; i<header.size(); i++) {
      if (i<row->size()) {
	col[header.at(i)] = row->at(i);
      } else {
	col[header.at(i)] = "";
      }
    }
    data.push_back(col);
  }

  return true;
}
					


bool getInputFileNames(const string InputTrees, map<string, vector<string> >& InputFileCollection)
{
  vector< vector<string> > content; 
  if(!readFile(InputTrees, content, 2)){ return false; }
  for(vector< vector<string> >::iterator row=content.begin(); row!=content.end(); ++row) {
    // remove spaces
    row->at(0).erase(std::remove(row->at(0).begin(), row->at(0).end(), ' '), row->at(0).end());
    row->at(1).erase(std::remove(row->at(1).begin(), row->at(1).end(), ' '), row->at(1).end());
    // store the filenames mapped by the tag
    InputFileCollection[row->at(0)].push_back(row->at(1));
  }
  return true;
};

bool readFile(string FileName, vector< vector<string> >& content, const int nCol)
{
  cout << "[INFO] Reading file: " << FileName << endl; 
  ifstream myfile(FileName.c_str());
  if (myfile.is_open()){ 
    string line;
    while ( getline(myfile, line) ){
      stringstream row(line);
      vector<string> cols; int i=0;
      while (true){
	string col; getline(row, col, ';');
	if ( (nCol>=0) ? (i>=nCol) : (col=="") ){ break; }
	cols.push_back(col);
	i++;
      }
      content.push_back(cols);
    }
  } else {
    cout << "[ERROR] File: " << FileName << " was not found!" << endl; return false;
  }
  return true;
};

bool iniWorkEnv(map<string,string>& DIR, const string workDirName) 
{
  cout << "[INFO] Initializing the work enviroment" << endl;
  DIR["main"]   = gSystem->ExpandPathName(gSystem->pwd());
  DIR["macros"] = DIR["main"] + "/Macros/";
  if (existDir(DIR["macros"].c_str())==false){ 
    cout << "[ERROR] Input directory: " << DIR["macros"] << " does not exist!" << endl; 
    return false; 
  }     
  DIR["input"]  = DIR["main"] + "/Input/" + workDirName + "/";
if (existDir(DIR["input"])==false){ 
    cout << "[ERROR] Input directory: " << DIR["input"] << " does not exist!" << endl; 
    return false; 
  }     
  DIR["output"] = DIR["main"] + "/Output/" + workDirName + "/";
  if (existDir(DIR["output"].c_str())==false){ 
    cout << "[INFO] Output directory: " << DIR["output"] << " does not exist, will create it!" << endl;  
    gSystem->mkdir(DIR["output"].c_str(), kTRUE);
  }
  DIR["dataset"] = DIR["main"] + "/DataSet/";
  if (existDir(DIR["dataset"])==false){ 
    cout << "[INFO] DataSet directory: " << DIR["dataset"] << " does not exist, will create it!" << endl;  
    gSystem->mkdir(DIR["dataset"].c_str(), kTRUE);
  }
  return true; 
};

bool existDir(string dir)
{
  void * dirp = gSystem->OpenDirectory(dir.c_str());
  if (dirp){
    gSystem->FreeDirectory(dirp);
    return true;
  }
  return false;
};

