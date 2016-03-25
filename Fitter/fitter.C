#include "Macros/Utilities/initClasses.h"
#include "Macros/tree2DataSet.C"
#include "Macros/fitCharmonia.C"


bool checkSettings(bool fitData, bool fitPbPb, bool fitPP, bool incJpsi, bool incPsi2S, bool incBkg, bool cutCtau,  
                   bool doSimulFit, bool wantPureSMC, bool setLogScale, bool zoomPsi, bool incSS, int numCores, int nBins);


bool parseFile(string FileName, vector< map<string, string> >& data);
bool parseString(string input, string delimiter, vector<double>& output);

bool iniWorkEnv(map<string,string>& DIR, const string workDirName);
bool existDir(string dir);
bool readFile(string FileName, vector< vector<string> >& content, const int nCol=-1, int nRow=-1);
bool getInputFileNames(const string InputTrees, map<string, vector<string> >& InputFileCollection);
bool setParameters(map<string, string> row, struct KinCuts& cut, map<string, string>& parIni);
bool addParameters(string InputFile,  vector< struct KinCuts >& cutVector, vector< map<string, string> >&  parIniVector);


void fitter(
            const string workDirName="Test", // Working directory
            // Select the type of datasets to fit
            bool fitData     = true,         // Fits Data if true, otherwise fits MC
            bool fitPbPb     = true,         // Fits PbPb datasets
            bool fitPP       = true,         // Fits PP datasets
            // Select the type of object to fit
            bool incJpsi     = true,         // Includes Jpsi model
            bool incPsi2S    = true,         // Includes Psi(2S) model
            bool incBkg      = true,         // Includes Background model
            // Select the fitting options
            bool cutCtau     = false,        // Apply prompt ctau cuts
            bool doSimulFit  = false,        // Do simultaneous fit
            bool wantPureSMC = false,        // Flag to indicate if we want to fit pure signal MC
            int  numCores    = 2,            // Number of cores used for fitting
            // Select the drawing options
            bool setLogScale = true,         // Draw plot with log scale
            bool incSS       = false,        // Include Same Sign data
            bool zoomPsi     = false,        // Zoom Psi(2S) peak on extra pad
            int  nBins       = 54            // Number of bins used for plotting
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

  if (!checkSettings(fitData, fitPbPb, fitPP, incJpsi, incPsi2S, incBkg, cutCtau, doSimulFit, wantPureSMC, setLogScale, zoomPsi, incSS, numCores, nBins)) { return; }

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
  
  TObjArray* aDSTAG = new TObjArray(); // Array to store the different tags in the list of trees
  aDSTAG->SetOwner(true);
  map<string, RooWorkspace> Workspace;

  for(map<string, vector<string> >::iterator FileCollection=InputFileCollection.begin(); FileCollection!=InputFileCollection.end(); ++FileCollection) {
    // Get the file tag which has the following format: DSTAG_COLL , i.e. DATA_PP 
    string FILETAG = FileCollection->first;  
    string DSTAG   = FILETAG;
    if (FILETAG.size()) {
      DSTAG.erase(DSTAG.find("_"));
    } else {
      cout << "[ERROR] FILETAG is empty!" << endl;
    }
    // Extract the filenames
    vector<string> InputFileNames = FileCollection->second; 
    string         OutputFileName;
    // If we have data, check if the user wants to fit data
    if ( (FILETAG.find("DATA")!=std::string::npos) && fitData==true ) {
      if ( (FILETAG.find("PP")!=std::string::npos)   && !fitPP   ) continue; // If we find PP, check if the user wants PP
      if ( (FILETAG.find("PbPb")!=std::string::npos) && !fitPbPb ) continue; // If we find PbPb, check if the user wants PbPb
      OutputFileName = DIR["dataset"] + "DATASET_" + FILETAG + ".root";
      if(!tree2DataSet(Workspace[DSTAG], InputFileNames, FILETAG, OutputFileName)){ return; }
      if (!aDSTAG->FindObject(DSTAG.c_str())) aDSTAG->Add(new TObjString(DSTAG.c_str()));
    }
    // If we find MC, check if the user wants to fit MC
    if ( (FILETAG.find("MC")!=std::string::npos) && fitData==false ) {
      if ( (FILETAG.find("PP")!=std::string::npos)    && !fitPP    ) continue; // If we find PP, check if the user wants PP
      if ( (FILETAG.find("PbPb")!=std::string::npos)  && !fitPbPb  ) continue; // If we find PbPb, check if the user wants PbPb
      if ( (FILETAG.find("JPSI")!=std::string::npos)  && !incJpsi  ) continue; // If we find Jpsi MC, check if the user wants to include Jpsi
      if ( (FILETAG.find("PSI2S")!=std::string::npos) && !incPsi2S ) continue; // If we find Psi2S MC, check if the user wants to include Psi2S
      OutputFileName = DIR["dataset"] + "DATASET_" + FILETAG + ".root";
      if(!tree2DataSet(Workspace[DSTAG], InputFileNames, FILETAG, OutputFileName)){ return; }
      if (!aDSTAG->FindObject(DSTAG.c_str())) aDSTAG->Add(new TObjString(DSTAG.c_str()));
      if (wantPureSMC)
      {
        OutputFileName = DIR["dataset"] + "DATASET_" + FILETAG + "_PureS" + ".root";
        if(!tree2DataSet(Workspace[Form("%s_PureS",DSTAG.c_str())], InputFileNames, FILETAG, OutputFileName)){ return; }
      }
    } 
  }
  if (Workspace.size()==0) {
    cout << "[ERROR] No onia tree files were found matching the user's input settings!" << endl; return;
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
 
  if (fitPbPb && incBkg) {
    // Add initial parameters for PbPb background models
    InputFile = (DIR["input"] + "InitialParam_MASS_BKG_PbPb.csv");
    if (!addParameters(InputFile, cutVector, parIniVector)) { return; }
  } 
  if (fitPbPb && incJpsi) {
    // Add initial parameters for PbPb jpsi models
    InputFile = (DIR["input"] + "InitialParam_MASS_JPSI_PbPb.csv");
    if (!addParameters(InputFile, cutVector, parIniVector)) { return; }
  } 
  if (fitPbPb && incPsi2S) {
    // Add initial parameters for PbPb psi(2S) models
    InputFile = (DIR["input"] + "InitialParam_MASS_PSI2S_PbPb.csv");
    if (!addParameters(InputFile, cutVector, parIniVector)) { return; }
  } 
  if (fitPP && incBkg) {
    // Add initial parameters for PP background models
    InputFile = (DIR["input"] + "InitialParam_MASS_BKG_PP.csv");
    if (!addParameters(InputFile, cutVector, parIniVector)) { return; }
  } 
  if (fitPP && incJpsi) {
    // Add initial parameters for PP jpsi models
    InputFile = (DIR["input"] + "InitialParam_MASS_JPSI_PP.csv");
    if (!addParameters(InputFile, cutVector, parIniVector)) { return; }
  } 
  if (fitPP && incPsi2S) {
    // Add initial parameters for PP psi(2S) models
    InputFile = (DIR["input"] + "InitialParam_MASS_PSI2S_PP.csv");
    if (!addParameters(InputFile, cutVector, parIniVector)) { return; }
  }

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
  
  TIter nextDSTAG(aDSTAG);
  string outputDir = DIR["output"];
  for (unsigned int i=0; i<cutVector.size(); i++) {

    nextDSTAG.Reset();
    TObjString* soDSTAG(0x0);
    while ( (soDSTAG = static_cast<TObjString*>(nextDSTAG.Next())) )
    {
      TString DSTAG   = static_cast<TString>(soDSTAG->GetString());
      
      if (Workspace.count(DSTAG.Data())>0) {
        // DATA/MC datasets were loaded
        if (doSimulFit) {
          // If do simultaneous fits, then just fits once
          if (!fitCharmonia( Workspace[DSTAG.Data()], cutVector.at(i), parIniVector.at(i), outputDir,
                             // Select the type of datasets to fit
                             DSTAG.Data(),
                             false,           // dummy flag when fitting simultaneously since both PP and PbPb are used
                             // Select the type of object to fit
                             incJpsi,         // Includes Jpsi model
                             incPsi2S,        // Includes Psi(2S) model
                             incBkg,          // Includes Background model
                             // Select the fitting options
                             cutCtau,         // Apply prompt ctau cuts
                             true,            // Do simultaneous fit
                             wantPureSMC,     // Flag to indicate if we want to fit pure signal MC
                             numCores,        // Number of cores used for fitting
                             // Select the drawing options
                             setLogScale,     // Draw plot with log scale
                             incSS,           // Include Same Sign data
                             zoomPsi,         // Zoom Psi(2S) peak on extra pad
                             nBins,           // Number of bins used for plotting
                             false            // Compute the mean PT (NEED TO FIX)
                             )
              ) { return; }
        } else {
          // If don't want simultaneous fits, then fit PbPb or PP separately
          if ( DSTAG.Contains("MCJPSI")  ) { incJpsi = true;  incPsi2S = false; }
          if ( DSTAG.Contains("MCPSI2S") ) { incJpsi = false; incPsi2S = true;  }
            
          if (fitPbPb) {
            if (!fitCharmonia( Workspace[DSTAG.Data()], cutVector.at(i), parIniVector.at(i), outputDir,
                               // Select the type of datasets to fit
                               DSTAG.Data(),
                               true,            // In this case we are fitting PbPb
                               // Select the type of object to fit
                               incJpsi,         // Includes Jpsi model
                               incPsi2S,        // Includes Psi(2S) model
                               incBkg,          // Includes Background model
                               // Select the fitting options
                               cutCtau,         // Apply prompt ctau cuts
                               false,           // Do simultaneous fit
                               false,     // Flag to indicate if we want to fit pure signal MC
                               numCores,        // Number of cores used for fitting
                               // Select the drawing options
                               setLogScale,     // Draw plot with log scale
                               incSS,           // Include Same Sign data
                               zoomPsi,         // Zoom Psi(2S) peak on extra pad
                               nBins,           // Number of bins used for plotting
                               false            // Compute the mean PT (NEED TO FIX)
                               )
                ) { return; }
            if (DSTAG.Contains("MC") && wantPureSMC)
            {
              if (!fitCharmonia( Workspace[Form("%s_PureS",DSTAG.Data())], cutVector.at(i), parIniVector.at(i), outputDir,
                                // Select the type of datasets to fit
                                DSTAG.Data(),
                                true,            // In this case we are fitting PbPb
                                // Select the type of object to fit
                                incJpsi,         // Includes Jpsi model
                                incPsi2S,        // Includes Psi(2S) model
                                incBkg,          // Includes Background model
                                // Select the fitting options
                                cutCtau,         // Apply prompt ctau cuts
                                false,           // Do simultaneous fit
                                true,            // Flag to indicate if we want to fit pure signal MC
                                numCores,        // Number of cores used for fitting
                                // Select the drawing options
                                setLogScale,     // Draw plot with log scale
                                incSS,           // Include Same Sign data
                                zoomPsi,         // Zoom Psi(2S) peak on extra pad
                                nBins,           // Number of bins used for plotting
                                false            // Compute the mean PT (NEED TO FIX)
                                )
                  ) { return; }
            }
          }
          if (fitPP) {
            if (!fitCharmonia( Workspace[DSTAG.Data()], cutVector.at(i), parIniVector.at(i), outputDir,
                               // Select the type of datasets to fit
                               DSTAG.Data(),
                               false,           // In this case we are fitting PP
                               // Select the type of object to fit
                               incJpsi,         // Includes Jpsi model
                               incPsi2S,        // Includes Psi(2S) model
                               incBkg,          // Includes Background model
                               // Select the fitting options
                               cutCtau,         // Apply prompt ctau cuts
                               false,           // Do simultaneous fit
                               false,           // Flag to indicate if we want to fit pure signal MC
                               numCores,        // Number of cores used for fitting
                               // Select the drawing options
                               setLogScale,     // Draw plot with log scale
                               incSS,           // Include Same Sign data
                               zoomPsi,         // Zoom Psi(2S) peak on extra pad
                               nBins,           // Number of bins used for plotting
                               false            // Compute the mean PT (NEED TO FIX)
                               )
                ) { return; }
            if (DSTAG.Contains("MC") && wantPureSMC)
            {
              if (!fitCharmonia( Workspace[Form("%s_PureS",DSTAG.Data())], cutVector.at(i), parIniVector.at(i), outputDir,
                                // Select the type of datasets to fit
                                DSTAG.Data(),
                                false,           // In this case we are fitting PP
                                // Select the type of object to fit
                                incJpsi,         // Includes Jpsi model
                                incPsi2S,        // Includes Psi(2S) model
                                incBkg,          // Includes Background model
                                // Select the fitting options
                                cutCtau,         // Apply prompt ctau cuts
                                false,           // Do simultaneous fit
                                true,            // Flag to indicate if we want to fit pure signal MC
                                numCores,        // Number of cores used for fitting
                                // Select the drawing options
                                setLogScale,     // Draw plot with log scale
                                incSS,           // Include Same Sign data
                                zoomPsi,         // Zoom Psi(2S) peak on extra pad
                                nBins,           // Number of bins used for plotting
                                false            // Compute the mean PT (NEED TO FIX)
                                )
                  ) { return; }
            }
          }
        }
      } else {
        cout << "[ERROR] The workspace for " << DSTAG.Data() << " was not found!" << endl; return;
      }  
    }
  }
  
  delete aDSTAG;
};


bool addParameters(string InputFile,  vector< struct KinCuts >& cutVector, vector< map<string, string> >&  parIniVector)
{
  vector< map<string, string> >  data;
  if(!parseFile(InputFile, data)) { return false; }
  
  if (cutVector.size()==0) {
    for(vector< map<string, string> >::iterator row=data.begin(); row!=data.end(); ++row) {
      struct KinCuts cut; map<string, string> parIni;
      if(!setParameters(*row, cut, parIni)) { return false; }
      cutVector.push_back(cut);  parIniVector.push_back(parIni);
    }
  }
  else {
    if (data.size()!=cutVector.size()) { cout << "[ERROR] The initial parameters in file " << InputFile << " are not consistent with previous files!" << endl; return false; }
    for (unsigned int i=0; i<data.size(); i++) {
      struct KinCuts cut;
      if (!setParameters(data.at(i), cut, parIniVector.at(i))) { return false; };
      if (!isEqualKinCuts(cut, cutVector.at(i))) { cout << "[ERROR] The bins in file " << InputFile << " are not consistent with previous files!" << endl; return false; }
    }
  }

  return true;
};


bool setParameters(map<string, string> row, struct KinCuts& cut, map<string, string>& parIni)
{

  // set initial parameters
  cut.sMuon.Pt.Min  =  0.0;    
  cut.sMuon.Pt.Max  = 100000.0;
  cut.sMuon.Eta.Min = -2.4;   
  cut.sMuon.Eta.Max = 2.4;
  cut.dMuon.ctauErr.Min = -100.0; 
  cut.dMuon.ctauErr.Max = 100.0;
  cut.dMuon.ctau.Min = -100.0;   
  cut.dMuon.ctau.Max = 100.0;    
  cut.dMuon.M.Min = 2.0; 
  cut.dMuon.M.Max = 5.0;  
  cut.dMuon.AbsRap.Min = 0.0;
  cut.dMuon.AbsRap.Max = 2.4;
  cut.dMuon.Pt.Min  =  0.0;
  cut.dMuon.Pt.Max  =  1000.0;
  cut.Centrality.Start = 0;
  cut.Centrality.End = 200;

  // set parameters from file
  for(map<string, string>::iterator col=row.begin(); col!=row.end(); ++col) {
    string label = col->first;
    if (label=="rap") {
      if (col->second=="" || col->second.find("-")==std::string::npos) {
        cout << "[ERROR] Input column 'rap' has invalid value: " << col->second << endl; return false;
      }  
      std::vector<double> v; 
      if(!parseString(col->second, "-", v)) { return false; }
      if (v.size()!=2) {
        cout << "[ERROR] Input column 'rap' has incorrect number of values, it should have 2 values but has: " << v.size() << endl; return false;
      }  
      cut.dMuon.AbsRap.Min = v.at(0); 
      cut.dMuon.AbsRap.Max = v.at(1);
    } 
    else if (label=="pt"){
      if (col->second=="" || col->second.find("-")==std::string::npos) {
        cout << "[ERROR] Input column 'pt' has invalid value: " << col->second << endl; return false;
      }
      std::vector<double> v; 
      if(!parseString(col->second, "-", v)) { return false; }
      if (v.size()!=2) {
        cout << "[ERROR] Input column 'pt' has incorrect number of values, it should have 2 values but has: " << v.size() << endl; return false;
      }  
      cut.dMuon.Pt.Min = v.at(0); 
      cut.dMuon.Pt.Max = v.at(1);
    } 
    else if (label=="cent"){
      if (col->second=="" || col->second.find("-")==std::string::npos) {
        cout << "[ERROR] Input column 'cent' has invalid value: " << col->second << endl; return false;
      }
      std::vector<double> v; 
      if(!parseString(col->second, "-", v)) { return false; }
      if (v.size()!=2) {
        cout << "[ERROR] Input column 'cent' has incorrect number of values, it should have 2 values but has: " << v.size() << endl; return false;
      }  
      cut.Centrality.Start = (int) (2.0*v.at(0)); 
      cut.Centrality.End = (int) (2.0*v.at(1));
    } 
    else if (label=="mass"){
      if (col->second=="" || col->second.find("-")==std::string::npos) {
        cout << "[ERROR] Input column 'mass' has invalid value: " << col->second << endl; return false;
      }
      std::vector<double> v; 
      if(!parseString(col->second, "-", v)) { return false; }
      if (v.size()!=2) {
        cout << "[ERROR] Input column 'mass' has incorrect number of values, it should have 2 values but has: " << v.size() << endl; return false;
      }  
      cut.dMuon.M.Min = v.at(0); 
      cut.dMuon.M.Max = v.at(1);
    } 
    else if (label=="ctau"){
      if (col->second=="" || col->second.find("->")!=std::string::npos) {
        cout << "[ERROR] Input column 'ctau' has invalid value: " << col->second << endl; return false;
      }
      std::vector<double> v; 
      if(!parseString(col->second, "->", v)) { return false; }
      if (v.size()!=2) {
        cout << "[ERROR] Input column 'ctau' has incorrect number of values, it should have 2 values but has: " << v.size() << endl; return false;
      }  
      cut.dMuon.ctau.Min = v.at(0); 
      cut.dMuon.ctau.Max = v.at(1);
    } 
    else if (label.find("Model")!=std::string::npos){
      if (col->second=="") {
        cout << "[ERROR] Input column 'Model' has empty value" << endl; return false;
      }
      parIni[col->first] = col->second;
    } 
    else {
      if (col->second != "") {
	string value = col->second;
	// check that initial parameters format is correct: [ num, num, num ]
	if ((value.find("[")==std::string::npos)||(value.find("]")==std::string::npos)) {
	  // Special cases like parameter constrains could be set here but for now, let's keep it simple

	  cout << "[ERROR] Either ']' or '[' are missing in the initial parameter values!" << endl; return false;
	} else {
	  value.erase(value.find("["), string("[").length());
	  value.erase(value.find("]"), string("]").length());
	}
	std::vector<double> v; 
	if(!parseString(value, ",", v)) { return false; }
        if (v.size()>3 || v.size()<1) {
          cout << "[ERROR] Initial parameter " << col->first << " has incorrect number of values, it has: " << v.size() << endl; return false;
        } 
	// everything seems alright, then proceed to save the values
	if (v.size()==1) { 
	  // if only one value is given i.e. [ num ], consider it a constant value
	  parIni[col->first] = Form("%s[ %.6f, %.6f, %.6f]", col->first.c_str(), v.at(0), v.at(0), v.at(0));
	} else {
	  parIni[col->first] = col->first + col->second;
	}
      } else {
        parIni[col->first] = "";
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
};


bool parseFile(string FileName, vector< map<string, string> >& data)
{
  vector< vector<string> > content, tmp; 
  if(!readFile(FileName, tmp, -1, 1)){ return false; }
  vector<string> header = tmp.at(0);
  if (header.size()==0) { cout << "[ERROR] The header is null!" << endl; return false; }
  if(!readFile(FileName, content, header.size())){ return false; }
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
};					


bool getInputFileNames(const string InputTrees, map<string, vector<string> >& InputFileCollection)
{
  vector< vector<string> > content; 
  if(!readFile(InputTrees, content, 2)){ return false; }
  for(vector< vector<string> >::iterator row=content.begin(); row!=content.end(); ++row) {
    // remove spaces
    row->at(0).erase(std::remove(row->at(0).begin(), row->at(0).end(), ' '), row->at(0).end());
    row->at(1).erase(std::remove(row->at(1).begin(), row->at(1).end(), ' '), row->at(1).end());
    // remove tabs
    row->at(0).erase(std::remove(row->at(0).begin(), row->at(0).end(), '\t'), row->at(0).end());
    row->at(1).erase(std::remove(row->at(1).begin(), row->at(1).end(), '\t'), row->at(1).end());
    if (row->at(0)!="" && row->at(1)=="") { cout << "[ERROR] There is an empty file name in your InputTrees.txt, please fix it" << endl; return false; }
    if (row->at(0)=="" && row->at(1)!="") { cout << "[ERROR] There is an empty file tag in your InputTrees.txt, please fix it" << endl; return false; }
    if (row->at(0)!="" && row->at(1)!="") {
      // store the filenames mapped by the tag
      InputFileCollection[row->at(0)].push_back(row->at(1));
    }
  }
  return true;
};


bool readFile(string FileName, vector< vector<string> >& content, const int nCol, int nRow)
{
  if (nCol==0 || nRow==0) { 
    cout << "[WARNING] Ignoring content of File: " << FileName << endl; return true; 
  }
  if (nRow!=1) { cout << "[INFO] Reading file: " << FileName << endl; }
  ifstream myfile(FileName.c_str());
  if (myfile.is_open()){ 
    string line;
    while ( getline(myfile, line) ){
      if (nRow==0) break; else {nRow=nRow-1;}
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
  bool exist = false;
  void * dirp = gSystem->OpenDirectory(dir.c_str());
  if (dirp){
    gSystem->FreeDirectory(dirp);
    exist = true;
  }
  return exist;
};


bool checkSettings(bool fitData, bool fitPbPb, bool fitPP, bool incJpsi, bool incPsi2S, bool incBkg, bool cutCtau, 
                   bool doSimulFit, bool wantPureSMC, bool setLogScale, bool zoomPsi, bool incSS, int numCores, int nBins)
{ 
  cout << "[INFO] Checking user settings " << endl;

  if (!fitPbPb && !fitPP) {
    cout << "[ERROR] No datasets were selected, please select either PbPb or PP!" << endl; return false;
  }
  if (!incJpsi && !incPsi2S && !incBkg) {
    cout << "[ERROR] No models were included for fitting, please include either Jpsi, Psi2S or Bkg" << endl; return false;
  }
  if (!fitData && (!incJpsi && !incPsi2S)) {
    cout << "[ERROR] We can't fit only background on MC, please include also either Jpsi or Psi2S!" << endl; return false;
  }
  if (fitData && wantPureSMC) {
    cout << "[WARNING] Pure signal MC will not be fitted since you have selected DATA!" << endl;
  }
  if (wantPureSMC && (!incPsi2S && !incJpsi)) {
    cout << "[ERROR] Pure signal MC can only be fitted if at least one signal model is included, please fix your input settings!" << endl;
  }
  if (numCores<=0 || numCores>32) {
    cout << "[ERROR] Invalid number of cores: " << numCores << " , please select an interger number between 1 and 32." << endl; return false;
  }
  if (nBins<=0) {
    cout << "[ERROR] Invalid number of bins: " << nBins << " , please select an interger number > 0 ." << endl; return false;
  }
  if (doSimulFit && (!fitPbPb || !fitPP)) {
    cout << "[ERROR] Can't perform simultaneous fit without both PP and PbPb, please fix your fitter settings!" << endl; return false;
  }
  if (doSimulFit && (!incPsi2S || !incJpsi)) {
    cout << "[ERROR] Simultaneous fits require both Psi2S and Jpsi, please fix your fitter settings!" << endl; return false;
  }
  if (doSimulFit && !fitData) {
    cout << "[ERROR] Can't perform simultaneous fits on MC, please fix your fitter settings!" << endl; return false;
  }
  if (setLogScale && zoomPsi) {
    cout << "[ERROR] Zooming on Psi(2S) is currently not supported with logarithmic scale." << endl;
  }

  cout << "[INFO] All user setting are correct " << endl;
  return true;
};

