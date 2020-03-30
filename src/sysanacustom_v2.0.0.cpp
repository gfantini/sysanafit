/*
 *     This is intended to load a file .root processed by root2dat_(v06)
 *     The important thing (required!) is that there exist histos
 *     h1000/ h700 / h2000 / h2200 / h5000 / h3000 / hrege  [histo name]
 *     DAT  / MC   / MC    / DAT   / DAT   / DAT   / DAT    [binning]
 *
 *     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 *     NEW HISTOGRAMS CREATED (ON THE FLY BY LAMBDA FUNCTION)
 *     h700norm   normalized smearing matrix such that sum_reco h700norm(true,reco) = 1
 *     hTheoryMCsmeared
 *     hTheoryDATsmeared
 *
 *     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 *     MEANING OF THE ARGUMENTS
 *
 *     folder     working folder: all files will be read / written there
 *     rootfile   name of .root containing input histo (should be inside "folder")
 *
 *
 *     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 *     MEANING OF RETURNED VALUE
 *     0          terminated with error
 *     1          everything went well
 *
 *
 *     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 *     MEANING OF VERSIONS
 *     1.0.0      release of the 1st working fitting procedure for zetaSL
 *     1.1.0      added implementation of gamma function (fit converges!)
 *     1.2.0      first implementation of the decoherencelib.hh -> BUG! the fit for SL does not converge
 *     1.2.1      bug fixed: was applying regeneration to non-rebinned histo
 *                BUG FOUND!! the smearing matrix used in decoherencelib is NOT NORMALIZED! As it was in v1.0.0
 *                Moved the old decoherencelib.hh -> decoherencelib.bug.hh and fixed the bug.
 *     2.0.0      First implementation of all the models in the GFitinterf constructor class -> fit gives reasonable output
 *                Some discrepancy wrt FITINTERF. -> this can depend on the different cut config (pz) or numerical approximations
 *                Added ::Save(<outputRootFile>,[index=0]) method to save fit result in <outputRootFile>/ModelName/... FitResult(_%d) %d = index or omitted
 *                   It is safe to run ::Save many times on the same root file. It will properly overwrite previous results.
 *     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 *     WORK TO BE DONE
 *     0. implement basic functionalities in stand-alone functions (?)
 *     1. implement all the models and fit the default file
 *     2. set up a ROOT procedure to extract ALL the relevant information from the fit
 *        like: likelihood contour / profiled likelihood / par. values / errors / fit status code / ...
 *        and store ALL OF THESE THINGS in the rootfile with appropriate names
 *
 *    ****** USER MANUAL
 *    1. define GFitinterf object providing
 *       - iInputFile = .root file with all the input histograms
 *       - model code = SL/00/gamma/omega (will be used to set the fitting function and output parameters)
 *                      This defines the "model" set of variables
 *       - [iOutputFile] eventuale .root diverso dall'input per salvare risultati (ad ora non salva nulla)
 *    2. object.Fit()
 *       - this defines the "Result" set of variables, comprising the result of the fit and histograms of data and fit distributions
 *         Those histograms can be retrieved with GetFitResultHdata GetFitResultHtheory
 *         The fit result root object can be read through GetFitResult
 */

#include "../inc/decoherencelib.hh"

// definition of verbosity flags
const int V_DEBUG = 0;
const int V_WARNING = 1;
const int V_ERROR = 2;
const int V_FATAL = 9001; // over 9000

class GFitinterf
{
	/* "Input" set of variables */
	bool Input_IsSet = false;
	TString InputFile;      // The InputFile should contain
	TString OutputFile;     // h1000/ h700 / h2000 / h2200 / h5000 / h3000 / hrege  [histo name]
	TFile *pInputFile;      // DAT  / MC   / MC    / DAT   / DAT   / DAT   / DAT    [binning]
	TFile *pOutputFile;
	
	/* "Model" set of variables */
	bool Model_IsSet = false;
	TString ModelName;
	Int_t Npar;        // number of fit parameters
	Double_t *Par;     // Npar objects, store here the initial value (can be modified..)
	Double_t *ParStep; // --> ParStep == 0. means SET CONSTANT
	TString *ParName;  // Npar objects
	Int_t FitRange;    // Fit is performed in |Dt|/tauS from 0 ... FitRange
	
	/* "Result" set of variables */
	bool Result_IsSet = false;
	bool bFitOk;       // did the fit converge? t/f
	TH1D* hData;        // data histogram after bkg subtraction, with errors comprising STAT (+) DATA/MC (+) EFFI (init in constructor)
	TH1D* hFit;         // fit  histogram (disregard errors!!) (init in constructor)
	ROOT::Fit::FitResult Result; // here the result of the fit is stored
	ROOT::Fit::FitResult ResultTmp; // here the result when the profile likelihood is being computed is temporarily stored
	
	/* "Profile" stuff */
	bool Profile_IsSet = false;
	TGraph Profile;
public:
  GFitinterf(TString iInputFile,TString Model,TString iOutputFile,const Int_t iFitRange); // allowed models are SL / 00 / gamma / omega
  ~GFitinterf();
  bool Fit(bool bOverwriteResult = true,const int verbosity = V_ERROR,const bool silent = true); // bOverwriteResult = true/false --> save fit result in Result / ResultTmp
  Double_t GetParameter(TString iParName); // get best fit parameter
  Double_t GetSymmetricError(TString iParName); // get best fit error (symmetric)
  /*
    Double_t GetSymmetricError(TString iParName); // tbi
    Double_t GetAsymmetricErrorLeft(TString iParName); // tbi
    Double_t GetAsymmetricErrorRight(TString iParName); // tbi
  */
  Double_t GetMinFCN(); // get value of FCN for best fit
  ROOT::Fit::FitResult GetFitResult(); // return the Result variable
  TH1D* GetFitResultHdata(bool print = false); // return hData
  TH1D* GetFitResultHtheory(bool print = false); // return hFit
  TGraph GetProfiledLikelihood(TString iParName,Double_t iParMin,Double_t iParMax,Double_t iParStep); /* return graph of FCN as a function of iParName, 
												       * minimized wrt all the other parameters */
  void Save(const char* iRootFileName,const int index = 0); // save file/ModelName/ FitResult_%d hData_%d hFit_%d, omit _%d if index = 0
  void PrintInfo();
private:
  Int_t GetParameterNumber(TString iParName); // returns -1 on failure
};

GFitinterf::GFitinterf(TString iInputFile,TString Model,TString iOutputFile = "UPDATE*INPUT",const Int_t iFitRange = 12)
{
	// *** setting up all the INPUT-block variables
	if(iOutputFile == "UPDATE*INPUT")
	{
		InputFile = iInputFile;
		OutputFile = iInputFile;
		pInputFile = new TFile(InputFile,"UPDATE");
		pOutputFile = pInputFile;
	}else{
		InputFile = iInputFile;
		OutputFile = iOutputFile;
		pInputFile = new TFile(InputFile,"READ");
		pOutputFile = new TFile(OutputFile,"UPDATE");
	}
	if( !pInputFile->IsZombie() && !pOutputFile->IsZombie() ) Input_IsSet = true;
	if(!Input_IsSet) cerr << "[GFinterf] ERROR: Input is not set" << endl;
	// *** end

	// *** setting up all the MODEL-block variables
	// this means all the settings about the fit
	ModelName = Model;
	if(iFitRange > 0 && iFitRange <= 12){
		FitRange = iFitRange;
	}else{
		cerr << "[GFinterf] FATAL: iFitRange should be in 0..12 range." << endl;
		exit(0);
	}
	if(ModelName == "SL")
	{

		Npar = 3;
		Par = new Double_t[Npar];
		ParStep = new Double_t[Npar];
		ParName = new TString[Npar];
		
		Par[0] = 0;
		ParStep[0] = 0.;
		ParName[0] = "ModelCode"; // 0=SL 1=00 2=gamma 3=omega
		
		Par[1] = 3.8e4;
		ParStep[1] = 1e1;
		ParName[1] = "Norm";
			
		Par[2] = 0.;
		ParStep[2] = 1e-3;
		ParName[2] = "ZetaSL";
			
		Model_IsSet = true;

	}
	else if (ModelName == "00")
	{
		Npar = 3;
		Par = new Double_t[Npar];
		ParStep = new Double_t[Npar];
		ParName = new TString[Npar];
		
		Par[0] = 1;
		ParStep[0] = 0.;
		ParName[0] = "ModelCode"; // 0=SL 1=00 2=gamma 3=omega
		
		Par[1] = 3.8e4;
		ParStep[1] = 1e1;
		ParName[1] = "Norm";
		
		Par[2] = 0.;
		ParStep[2] = 1e-7;
		ParName[2] = "Zeta00";
		
		Model_IsSet = true;
	}
	else if (ModelName == "gamma")
	{
		Npar = 3;
		Par = new Double_t[Npar];
		ParStep = new Double_t[Npar];
		ParName = new TString[Npar];
		
		Par[0] = 2; // 0=SL 1=00 2=gamma 3=omega
		ParStep[0] = 0.;
		ParName[0] = "ModelCode";
		
		Par[1] = 3.8e4;
		ParStep[1] = 1e1;
		ParName[1] = "Norm";
		
		Par[2] = 0.;
		ParStep[2] = 1e-21;
		ParName[2] = "gamma";
		
		Model_IsSet = true;
	}
	else if (ModelName == "omega")
	{
		Npar = 4;
		Par = new Double_t[Npar];
		ParStep = new Double_t[Npar];
		ParName = new TString[Npar];
		
		Par[0] = 3; // 0=SL 1=00 2=gamma 3=omega
		ParStep[0] = 0.;
		ParName[0] = "ModelCode";
		
		Par[1] = 3.8e4;
		ParStep[1] = 1e1;
		ParName[1] = "Norm";
		
		Par[2] = 0.;
		ParStep[2] = 1e-4;
		ParName[2] = "ReOmega";
		
		Par[3] = 0.;
		ParStep[3] = 1e-4;
		ParName[3] = "ImOmega";
		
		Model_IsSet = true;
	}
	else
	{
		cerr << "[GFinterf] PANIC: Model is not set. Invalid model " << ModelName <<endl;
		exit(0); // crash the software
	}
	// *** end

	// *** setting the RESULT histograms
	hData = new TH1D("hData","data after background subtraction (events/#tau_S)",FitRange,0,FitRange);
	hFit = new TH1D("hFit","best fit of data (events/#tau_S)",FitRange,0,FitRange);
}
void GFitinterf::PrintInfo()
{
	if(Model_IsSet)
	{
		cout << "--------------------------------------------------------" << endl;
		cout << "[GFinterf::PrintInfo] Model is set." << endl;
		cout << "[GFinterf::PrintInfo] ModelName \t = " << ModelName << endl;
		cout << "[GFinterf::PrintInfo] Npar \t\t = " << Npar << endl;
		for(auto i=0;i<Npar;i++)
		{
			cout << "[GFinterf::PrintInfo] Par["<< i << "] \t\t = " << Par[i] << endl;
			cout << "[GFinterf::PrintInfo] ParStep["<< i << "] \t = " << ParStep[i] << endl;
			cout << "[GFinterf::PrintInfo] ParName["<< i << "] \t = " << ParName[i] << endl;
		}
	}else{
		cout << "--------------------------------------------------------" << endl;
		cout << "[GFinterf::PrintInfo] Model is NOT set." << endl;
	}
	
	if(Input_IsSet)
	{
		cout << "--------------------------------------------------------" << endl;
		cout << "[GFinterf::PrintInfo] InputFile is " << InputFile << endl;
		cout << "[GFinterf::PrintInfo] OutputFile is " << OutputFile << endl;
	}
}

GFitinterf::~GFitinterf()
{
	pInputFile->Close();

	if(!pOutputFile->IsZombie()){
		pOutputFile->Close();
	}
	//pInputFile->~TFile();
	//pOutputFile->~TFile();
	//delete [] Par;
	//delete [] ParStep;
	//delete [] ParName;
}
Int_t GFitinterf::GetParameterNumber(TString iParName)
{
	if(Model_IsSet){
		for(int kk = 0; kk < Npar; kk++){
			if(ParName[kk] == iParName)return kk;
		}
		cerr << "[GFitinterf::GetParameterNumber] ERROR: Parameter " << iParName << " not found." << endl;
	}else{
		cerr << "[GFitinterf::GetParameterNumber] ERROR: Model_IsSet = false" << endl;
	}
	return -1;
}
Double_t GFitinterf::GetParameter(TString iParName)
{
	Int_t num = this->GetParameterNumber(iParName);
	return Result.GetParams()[num];
}

Double_t GFitinterf::GetSymmetricError(TString iParName)
{
	Int_t num = this->GetParameterNumber(iParName);
	return Result.GetErrors()[num];
}
Double_t GFitinterf::GetMinFCN()
{
	return Result.MinFcnValue(); // return value of whatever is minimized @ minimum
}

TGraph GFitinterf::GetProfiledLikelihood(TString iParName,Double_t iParMin,Double_t iParMax,Double_t iParStep)
{
	TGraph g;
	
	// check that fit was already performed and went ok
	if(!bFitOk)cerr << "[GetProfiledLikelihood] bFitOk = FALSE" << endl;
	
	ROOT::Fit::FitResult * pFitResult;
	
	// save here old settings for the profile parameter to be restored at the end
	Int_t myParNum = this->GetParameterNumber(iParName);
	Double_t myOldParStep = ParStep[myParNum];
	Double_t myOldPar = Par[myParNum];
	
	Double_t myParameterNow = iParMin;
	int i = 0;
	while(myParameterNow <= iParMax)
	{
		// compute profile FCN
		
		// set "iParName" == constant to ParameterNow
		ParStep[myParNum] = 0.;
		Par[myParNum] = myParameterNow;
		
		// perform fit and save to ResultTmp
		this->Fit(false,V_ERROR);
		
		// get result and add to TGraph
		g.SetPoint(i,myParameterNow,ResultTmp.MinFcnValue());
		
		// update myParameterNow
		myParameterNow += iParStep;
		i++;
	}
	
	// restore old parameter settings
	ParStep[myParNum] = myOldParStep;
	Par[myParNum] = myOldPar;
	
	// set TGraph names on axis
	g.GetXaxis()->SetTitle(iParName);
	g.GetYaxis()->SetTitle("profiled FCN value");
	
	return g;
}
ROOT::Fit::FitResult GFitinterf::GetFitResult()
{
  if(Result_IsSet){
    return Result;
  }else{
	  Error("GFitinterf::GetFitResult()","Result_IsSet = false! Panic!"); 
    exit(0);
  }
}
TH1D* GFitinterf::GetFitResultHdata(bool print)
{
	if(Result_IsSet == false){
		Error("GFitinterf::GetFitResultHdata()","Result_IsSet = false! Panic!");
		exit(0);
	}
	if(print)for(auto i = 1; i<=hData->GetNbinsX(); i++)cout << "hData[" << i << "] / error = " << hData->GetBinContent(i) << " / " << hData->GetBinError(i)<< endl;
	return hData;
}

TH1D* GFitinterf::GetFitResultHtheory(bool print)
{
	if(Result_IsSet == false){
		Error("GFitinterf::GetFitResultHdata()","Result_IsSet = false! Panic!");
		exit(0);
	}
	if(print)for(auto i = 1; i<=hFit->GetNbinsX(); i++)cout << "hFit[" << i << "] = " << hFit->GetBinContent(i)<< endl;
	return hFit;
}

void GFitinterf::Save(const char* iRootFileName,const int index){
	if(Model_IsSet == false || Result_IsSet == false){
		cerr << "[GFitinterf::Save] FATAL. Model or Result is not set." << endl;
		exit(0);
	}
	TFile* fp = new TFile(iRootFileName,"UPDATE");
	if(fp->IsZombie()){
		cerr << "[GFitinterf]::Save FATAL: could not open " << iRootFileName << endl;
		exit(0);
	}
	// --> invent some way to keep track of different model fits... (folders?)
	// if not already existing, make folder with the ModelName = SL/00/gamma/omega
	if(fp->GetDirectory(ModelName.Data()) == 0){
		fp->mkdir(ModelName.Data());
	}
	fp->cd(ModelName.Data());
	if(index == 0)
	{
		// save FilterResult hData hFit
		gDirectory->WriteObject(&Result,"FitResult","SingleKey | WriteDelete");
		gDirectory->WriteObject(hData,"hData","SingleKey | WriteDelete");
		gDirectory->WriteObject(hFit,"hFit","SingleKey | WriteDelete");
	}else{
		gDirectory->WriteObject(&Result,Form("FitResult_%d",index),"SingleKey | WriteDelete");
		gDirectory->WriteObject(hData,Form("hData_%d",index),"SingleKey | WriteDelete");
		gDirectory->WriteObject(hFit,Form("hFit_%d",index),"SingleKey | WriteDelete");
	}
	fp->Close();
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
bool GFitinterf::Fit(bool bOverwriteResult,const int verbosity,const bool silent){
	// 1.     load histos to RAM
	TH1D *h1000 = (TH1D*)pInputFile->Get("h1000"); // data distribution
	TH2D *h700  = (TH2D*)pInputFile->Get("h700");  // smearing matrix (non-normalized!)
	TH1D *h2000 = (TH1D*)pInputFile->Get("h2000"); // efficiency [MC-bin]
	TH1D *h2200 = (TH1D*)pInputFile->Get("h2200"); // efficiency [DAT-bin]
	TH1D *h5000 = (TH1D*)pInputFile->Get("h5000"); // 4pi-background
	TH1D *h3000 = (TH1D*)pInputFile->Get("h3000"); // DATA-MC correction
	TH1D *hrege = (TH1D*)pInputFile->Get("hrege"); // regeneration correction
	
	
	// 2.     build likelihood function to fit data
	auto likelihood = [&](const Double_t *par)
	{// temporary version of function to minimize, just to check that it works
		const bool debug = false;
		if(debug)cerr << "[Likelihood] par[0]:"<< par[0] << " par[1]:" << par[1] <<" par[2] = " << par[2] << endl;
		TH1D hTheory = GetTheoryHistogram(par,280,0.,70.); // store here the integral of theoretical function (par[0] defines the model)
		if(debug)for(int jj = 1; jj <= 280; jj++)cerr << "THEORY[" << jj << "] = " << hTheory.GetBinContent(jj) << endl;
		TH1D hTheoryEffi = MultiplyHistogram(hTheory,*h2000,280,0.,70.); // apply efficiency correction bin-by-bin
		if(debug)for(int jj = 1; jj <= 280; jj++)cerr << "THEORY+EFFI[" << jj << "] = " << hTheoryEffi.GetBinContent(jj) << endl;
		TH1D hTheorySmeared = GetSmearedHistogram(hTheoryEffi,*h700); // apply smearing (yay!)
		if(debug)for(int jj = 1; jj <= 280; jj++)cerr << "THEORY+EFFI+SMEARING[" << jj << "] = " << hTheorySmeared.GetBinContent(jj) << endl;
		TH1D hTheorySmearedRebin = *(TH1D*)hTheorySmeared.Clone();
		hTheorySmearedRebin.Rebin(4); // rebin from 0.25 --> 1.
		if(debug)for(int jj = 1; jj <= 12; jj++)cerr << "THEORY+EFFI+SMEARING+REBIN[" << jj << "] = " << hTheorySmearedRebin.GetBinContent(jj) << endl;
		TH1D hTheoryRegen = MultiplyHistogram(hTheorySmearedRebin,*hrege,12,0.,12.); // apply regeneration correction
		if(debug)for(int jj = 1; jj <= 12; jj++)cerr << "THEORY+..+REGEN[" << jj << "] = " << hTheoryRegen.GetBinContent(jj) << endl;
		TH1D hTheoryDatMc = MultiplyHistogram(hTheoryRegen,*h3000,12,0.,12.); // apply DATA/MC correction
		if(debug)for(int jj = 1; jj <= 12; jj++)cerr << "THEORY+..+REGEN+DAT/MC[" << jj << "] = " << hTheoryDatMc.GetBinContent(jj) << endl;

		Double_t chisquare = 0.;
		Double_t bin_value_data = 0.;
		Double_t bin_value_theory = 0.;
		Double_t bin_error = 0.;
		
		if(debug)cout << "- - - - - - - - - - - - - - - - - - - " << endl;
		if(debug)cout << "[Likelihood] CHISQUARE COMPUTATION" << endl;
		if(debug)cout << "- - - - - - - - - - - - - - - - - - - " << endl;

		for(int bin_reco = 1; bin_reco <= FitRange; bin_reco++)
		{
			bin_value_data   = h1000->GetBinContent(bin_reco) - h5000->GetBinContent(bin_reco); // background subtraction
			bin_value_theory = par[1]*hTheoryDatMc.GetBinContent(bin_reco);
			bin_error = bin_value_data*sqrt(  (h1000->GetBinContent(bin_reco) + h5000->GetBinContent(bin_reco))/pow(bin_value_data,2) + pow(h3000->GetBinError(bin_reco)/h3000->GetBinContent(bin_reco) ,2) + pow(h2200->GetBinError(bin_reco)/h2200->GetBinContent(bin_reco) ,2) ); // 3 error terms come from STAT (+) DATA/MC (+) EFFI
			if(bOverwriteResult){ // overwrite the histograms in the Result part It was initialized by the constructor.
				hData->SetBinContent(bin_reco,bin_value_data);
				hData->SetBinError(bin_reco,bin_error);
				hFit->SetBinContent(bin_reco,bin_value_theory);
			}
			chisquare += pow((bin_value_data-bin_value_theory)/bin_error,2);
			if(debug)cout << std::setprecision(8) << "BIN: "<< bin_reco << "\tDAT: " << bin_value_data << "\tTHEORY: "<< bin_value_theory << "\tERR: " << bin_error << endl;
		}
		return chisquare;
	};
	
	
	// wrap chi2 funciton in a function object for the fit
	// Npar is the number of fit parameters (size of array par)
	ROOT::Math::Functor fcn(likelihood,Npar);
	
	
	// 3.     fit data
	ROOT::Fit::Fitter  fitter;
	fitter.SetFCN(fcn,Par);

	// configure fitter according to class construction
	for(int jj = 0; jj < Npar; jj++)
	{
		if(ParStep[jj] == 0.){ // parameter should be set to constant
			fitter.Config().ParSettings(jj).Set(ParName[jj].Data(),Par[jj]);
		}else{
			fitter.Config().ParSettings(jj).SetName(ParName[jj].Data());
			fitter.Config().ParSettings(jj).SetStepSize(ParStep[jj]);
			// fitter.Config().ParSettings(jj).SetLowerLimit(/**/); // EVENTUALLY SET LOWER LIMIT...
		}
	}
	
	/*
	 // SET OF CONFIGURATIONS FOR ZETA SL
	fitter.Config().ParSettings(0).Set("ModelCode",0.);  // set as fixed parameter
	fitter.Config().ParSettings(1).SetName("Norm");
	fitter.Config().ParSettings(2).SetName("ZetaSL");
	fitter.Config().ParSettings(3).Set("disabled",-1);   // set as fixed
	//	fitter.Config().ParSettings(1).SetLowerLimit(0.);
	fitter.Config().ParSettings(1).SetStepSize(10);
	fitter.Config().ParSettings(2).SetStepSize(0.001);
     */
	
	/*
	// SET OF CONFIGURATIONS FOR GAMMA
	fitter.Config().ParSettings(0).Set("ModelCode",2.);   // set as fixed parameter 2. <-> Gamma model
	fitter.Config().ParSettings(1).SetName("Norm");       // parameter 1 <-> always the normalization
	fitter.Config().ParSettings(2).SetName("Gamma [GeV]");// parameter 2 <-> 1st physical parameter
	fitter.Config().ParSettings(3).Set("disabled",-1);    // parameter 3 <-> 2nd physical parameter or disabled (in this case)
	fitter.Config().ParSettings(1).SetStepSize(10);
	fitter.Config().ParSettings(2).SetStepSize(1e-21);
	 */
	
	bool my_bFitOk = fitter.FitFCN();
	if( !my_bFitOk ) Error("sysanacustom","Fit did not converge.");
	if( false == fitter.CalculateMinosErrors())cerr << "ERROR: impossible to compute MINOS uncertainty." << endl;
	
	// 4.     retrieve results
	const ROOT::Fit::FitResult & my_result = fitter.Result(); // Rene Brun la sa lunga
	if(!silent)my_result.Print(std::cout); // print results on screen
	
	// 5.     store results
	if(bOverwriteResult == true)
	{
		// nominal fit, save result normally
		Result_IsSet = true;
		bFitOk = my_bFitOk;
		Result = fitter.Result();
	}else{
		// it means it is used by some other method (like profiled likelihood) and should be exported as temporary result
		ResultTmp = fitter.Result();
	}
	
	// 6.     return result (fit converged? t/f)
	return my_bFitOk;
}

// DEBUG STUFF - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

void NormalizationTest(double zetasl){
	double normalization = ComputeIntegralSL(zetasl,0.,12.);
	cout << "Integral in 0. <-> 12. is \t" << normalization << endl;
	double tmp = 0.;
	double tot = 0.;
	for(int ii = 0.; ii < 280; ii++)
	{
		tmp = ComputeIntegralSL(zetasl,0. + ii*0.25,0. + (ii+1)*0.25)/normalization;
		cout << "VAL " << 0.+0.25*ii << " <-> " << 0.+0.25*(ii+1) << "\t" << tmp << endl;
		if(ii < 12*4)tot += tmp;
	}
	cout << "tot = " << tot << endl;
}
