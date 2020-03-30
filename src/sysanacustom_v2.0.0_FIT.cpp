/* 
 * Questa macro fitta (con impostazioni di default) tutti e quattro i parametri
 * e li stampa a schermo
 *
 * PARAMETRI IN INPUT
 * inputFileName = nome del file root che contiene tutti gli istogrammi (cfr. sysanacustom)
 */


#include "sysanacustom_v2.0.0.cpp"
// fit one file with one model (not considering PhysConstVariation)
void FitOneFile_Model(const char* inputFileName,const char* model,bool print = false)
{
	if(print)cout << "FitOneFile: " << inputFileName << " " << model << endl;
	// init fitter
	GFitinterf fitter(inputFileName,model);
	fitter.Fit(); // fit silent
	fitter.Save(inputFileName); // save result
	if(print)fitter.GetFitResult().Print(std::cout);
	
	// check if there is TTree of physical constants inside
	TFile* pFile = new TFile(inputFileName,"READ");
	TTree* pTree = (TTree*)pFile->Get("PhysicalConstants");
	Double_t TauS,TauL,DeltaM,modeta,phieta;
	if(pTree != NULL){ // Tree found --> gotta loop over entries
		// IMPLEMENT METHOD OF *** to load directly from TTree? Nah.
		pTree->SetBranchAddress("TauS",&TauS);
		pTree->SetBranchAddress("TauL",&TauL);
		pTree->SetBranchAddress("DeltaM",&DeltaM);
		pTree->SetBranchAddress("modeta",&modeta);
		pTree->SetBranchAddress("phieta",&phieta);
		for(int i=0;i<pTree->GetEntries();i++){
			pTree->GetEntry(i);
			ResetGlobalPhysicalConstants(TauS,TauL,DeltaM,modeta,phieta);
			fitter.Fit();
			fitter.Save(inputFileName,i);
			if(i%100 == 0)cout << "Fit " << i << " / " <<pTree->GetEntries() << endl;
		}
	}
	if(print)cout << "*****************" << endl;
}

// fit one file with all models
void FitOneFile(const char* inputFileName)
{
	FitOneFile_Model(inputFileName,"SL");
	FitOneFile_Model(inputFileName,"00");
	FitOneFile_Model(inputFileName,"gamma");
	FitOneFile_Model(inputFileName,"omega");
}
