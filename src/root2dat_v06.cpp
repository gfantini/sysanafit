/*
 v06 Guido Fantini 20.07.2017
 + load new histograms into OUTPUT rootfile (all the necessary histo to perform fit)
   as far as I remember this one still does NOT produce .dat output
    // to be implemented
    + h2000 = (TH1D) h200/h100: 280 bin in 0. <-> 70. = efficienza con bin MC
      ++ [OK] h_MCeffi | "h2000"
    + h2200 = (TH1D) h2000 rebinnato: 50 bin in 0. <-> 50. = efficienza con bin DAT
      ->>>>  Mi rifiuto di scrivere un istogramma con binning custom: sarà 0. <-> 70. con 70 bin (DAT binning)
      ++ [OK] h_DATeffi | "h2200"
    + h5000 = (TH1D) fondo 4pi: capire come viene generato
      ++ [OK] h_DAT4pi | "h5000"
    + h3000 = (TH1D) data/mc correction: (cost) da 3000v.dat / 3000e.dat in fit/dat
              path is NOT hard coded, can be changed at runtime as argument to the macro
      ++ [OK] h_DATdatamontecarlo | "h3000"
    + hrege = regeneration: 12 bin in 0. <-> 12. from fit/dat/regen_default.dat
              regen_default.dat is produced from root2dat_v0x from (const) ./inc/regen.dat
              The same content will be added to rootfile.
              Now the path to load from can be passed as an argument
	 ++ [OK] h_DATregen || "hrege"
 + in order to modify the 4pi subtraction when the mpp cut is changed, or during 4pi cut variation tests, or the regeneration correction
   new input were added
 ----------------------------------------------------------
 Old histos that were there: (filled by hmaker_v02)
 h1000   = data histo (TH1F) with INTEGER bin content
 h700    = smearing matrix (TH2F)                        <--- low numerical accuracy!
 h100    = MC (all)  (TH1F) with INTEGER bin content
 h200    = MC (sel)  (TH1F) with INTEGER bin content
 ==========================================================
 + OPEN BUGS / PROBLEMS
      + write documentation about this in the logbook and/or everywhere...
      + 4pi bkg correction as a function of the window implemented HERE (reading the input coefficients for cuts & co..)
      + regeneration correction variation implemented HERE -> use 0.12 / 0.24 values for Kregeneration constant
 
 + CLOSED ISSUES
      - if I run twice this macro on the same .root file and add different versions
        of an histogram, when I retrieve it, which histo do I get?
        I don't care, as long as I use histo.Write("",TObject::kOverwrite)
        I don't produce multiple versions of the same object
 
 v05 Guido Fantini
 + new prototype, loads configfile
 + makes in outputfolder regendefault.dat and 5000v.dat based on configfile
 [OK] tested!!
 + BUG fix on errDAT.dat and errMC.dat files. Now they contain
 errDAT = (N-B)sqrt(1/N + 1/B)
 errMC = (N-B)sqrt( [σ(dat-MC) / dat-MC]**2 + [σ(ε) / ε]**2 )
 and need to be combined in quadrature 'cause the final result is a product (N-B)*ε*(dat-MC)
 
 v04 Guido Fantini 02-06-2016 [looks OK]
 Generates (in the given outputfolder) ./dat/<file>.dat as the older versions
 + makes ./errDAT.dat
 + makes ./errMC.dat
 Those 2 files, combined in quadrature, give the error on data-bkg
 
 v03 Guido Fantini 30-05-2016 TESTED OK (inserted in sysana.sh)
 Same as v02 (generates ./dat/<file>.dat)
 + generates ./1000toterr.dat <-- in which there is total error for data histo
 -> needs ./fit/dat/3000v.dat ./fit/dat/3000e.dat to load DATA/MC correction value AND error
 + file2vector function
 
 v02 Guido Fantini 20-05-2016
 Non-interactive interface to generate 1000v 1000e 2000v 2200e 700v .dat
 Needs input filename (.root) <-- TH1F 1000 100 200 TH2F 700 should be inside
 Needs output folder: write path including final "/"
 
 v01 Guido Fantini
 Interactive interface to generate .dat
 
 */
#include "../inc/phys_const.inc"
#include "../inc/gflib_v08.h"

#define MATRIX_NBIN 280
struct matrix
{
	float el[MATRIX_NBIN][MATRIX_NBIN];
};

int load_TH1F(const char* rootfile,const char* histoname,int nbin,int* output)
{
	TFile* myfile = new TFile(rootfile,"READ");
	TH1F* myhisto = (TH1F*)myfile->Get(histoname);
	int j;
	for(j=1;j<=nbin;j++)output[j-1]=myhisto->GetBinContent(j);
	cout << rootfile << " underflow "<< myhisto->GetBinContent(0) << endl;
	cout << rootfile << " overflow " << myhisto->GetBinContent(nbin+1) << endl;
	myfile->Close();
	delete myfile;
	
	return 1;
}

int load_TH1F(const char* rootfile,const char* histoname,int nbin,double* output)
{
	TFile* myfile = new TFile(rootfile,"READ");
	TH1F* myhisto = (TH1F*)myfile->Get(histoname);
	int j;
	for(j=1;j<=nbin;j++)output[j-1]=myhisto->GetBinContent(j);
	cout << rootfile << " underflow "<< myhisto->GetBinContent(0) << endl;
	cout << rootfile << " overflow " << myhisto->GetBinContent(nbin+1) << endl;
	myfile->Close();
	delete myfile;
	
	return 1;
}

int load_TH2F(const char* rootfile,const char* matrixname,int nbin,matrix* out)
{
	TFile* myfile = new TFile(rootfile,"READ");
	TH2F* mymatrix = (TH2F*)myfile->Get(matrixname);
	Int_t bin;
	
	int i,j;
	for(i=1;i<=nbin;i++)
	{
		for(j=1;j<=nbin;j++)
		{
	  // procedure to follow, dunno why
	  bin=mymatrix->GetBin(i,j);
	  out->el[i-1][j-1]=mymatrix->GetBinContent(bin);
		}
	}
	myfile->Close();
	delete myfile;
	return 1;
}
/*
 int vector2file(const char* filename,int nbin,int* input)
 {
 ofstream myfile;
 int i;
 myfile.open(filename);
 for(i=0;i<nbin;i++)myfile << input[i] << endl;
 myfile.close();
 return 1;
 }
 int vector2file(const char* filename,int nbin,double* input)
 {
 ofstream myfile;
 int i;
 myfile.open(filename);
 for(i=0;i<nbin;i++)myfile << setprecision(20) << input[i] << endl;
 myfile.close();
 return 1;
 }
 int file2vector(const char* filename,int nbin,double* output)
 {
 ifstream myfile;
 int i;
 myfile.open(filename,ios::in);
 for(i=0;i<nbin;i++)myfile >> output[i];
 myfile.close();
 return 1;
 }
 */
int matrix2file(const char* filename,int nbin,matrix input)
{
	ofstream myfile;
	int i,j;
	myfile.open(filename);
	for(i=0;i<nbin;i++)
	{
		for(j=0;j<nbin;j++)
		{
	  myfile << input.el[i][j] << endl;
		}
	}
	myfile.close();
	return 1;
}
//--------------------------------------------------------------
//       END OF GENERAL PURPOSE FUNCTIONS
//--------------------------------------------------------------


#define NBIN_DAT 50
#define NBIN_MC 280
#define NBIN_FIT 12

// verbosity level definition
const int ERROR_V = 0;
const int WARNING_V = 1;
const int INFO_V = 2;
const int DEBUG_V = 3;

void root2dat_v06()
{
	cout << "HELP for root2dat_v06" << endl;
	cout << "---------------------" << endl;
	cout << "No, you ain't gettin help :]      ...kidding" << endl;
	
	cout << "You call: root2dat_v06(CONFIGFILE,FOLDER,ROOTFILE,[VERBOSITY],[REGEN_FILENAME],[DMCV_FILENAME],[DMCE_FILENAME],[see hmaker_h700_v03])" << endl;
	cout << "CONFIGFILE:\t\t\t\t ASCII column of integer telling what to do with the ordered cuts described in gflib" << endl;
	cout << "FOLDER:\t\t\t\t\t path to folder where there is I/O root file" << endl;
	cout << "ROOTFILE:\t\t\t\t just the rootfile name like something.root, with no slashes. The output of the macro will be stored there." << endl;
	cout << "VERBOSITY THRESHOLDS: "<< endl;
	cout << "\t error:\t " << ERROR_V << endl;
	cout << "\t warn:\t " << WARNING_V << endl;
	cout << "\t info:\t " << INFO_V << endl;
	cout << "\t debug:\t " << DEBUG_V << endl;
	cout << "REGEN_FILENAME:\t\t\t\t path wrt FOLDER where the regeneration file is" << endl;
	cout << "\t\t\t\t default is regen_default.dat" << endl;
	cout << "DATA-MC CORRECTION FILENAME:\t\t path wrt folder where the DAT/MC correction values are stored (.dat)" << endl;
	cout << "\t\t\t\t default is h3000v.dat" << endl;
	cout << "DATA-MC CORRECTION ERROR FILENAME:\t path wrt folder where the error on the DAT/MC correction values are stored (.dat)" << endl;
	cout << "\t\t\t\t default is h3000e.dat" << endl;
	cout << "Kregeneration: fractional variation (coherent) of regeneration correction. Default is 0.00, normal is +-0.12 +-0.24" << endl;
	//cout << "OUTPUTFOLDER: folder where the output files will be placed." << endl;
	cout << "- - - - - - - - - - - - " << endl;
	cout << "This will add to .root a set of histograms (overwriting each time):" << endl;
	cout << "+ h2000 \t efficienza, h200/h100 with binomial error [MC-binning]" << endl;
	cout << "+ h2200 \t efficienza, h200/h100 with binomial error [DAT-binning]" << endl;
	cout << "+ h5000 \t 4pi background <-- HARD CODED!!" << endl;
	cout << "+ h3000 \t DATA-MC correction with error loaded from DMCV_FILENAME / DMCE_FILENAME" << endl;
	cout << "+ hrege \t regeneration correction loaded from REGEN_FILENAME (12 bins in 0. <-> 12.)" << endl;
	
}

// root [0] .L rootdat_v06.cpp
// root [1] root2dat_v06("default.config","./default/","input.root");
void root2dat_v06(const char* configfile,const char* folder,const char* rootfile,const int verbosity = 0,const char* regen_filename = "../fit/dat/regen_default.dat",const char* datamc_correction_filename = "../fit/dat/3000v.dat",const char* datamc_correction_error_filename = "../fit/dat/3000e.dat",const int cut_to_vary=0,const int nsigma=0,const int nSigma4pi = 0,const double Kregeneration=0/* +-0.12,+-0.24*/,const double Kresolution = 0.)
{// always put final "/" to outputfolder ex.: "./dat/" is ok "./dat" is NOT ok
	int i,j,k;
	// setting the output rootfile name = to input_filename but in outputfolder
	char rootfilename[100];
	sprintf(rootfilename,"%s%s",folder,rootfile);
	if(verbosity >= INFO_V){
		cout << "root2dat_v06 [INFO]: I/O folder " << folder << endl;
		cout << "root2dat_v06 [INFO]: I/O file name " << rootfile << endl;
		cout << "root2dat_v06 [INFO]: configfile path " << configfile << endl;
		//cout << "root2dat_v06 [INFO]: output folder " << outputfolder << endl;
		//cout << "root2dat_v06 [INFO]: output file name set to " << output_rootfilename << endl;
	}
	
	// opening rootfile
	TFile *fin = new TFile(rootfilename,"UPDATE");
	if(fin->IsZombie() && verbosity >= ERROR_V)cerr <<"root2dat_v06 [ERROR]: Failed to open I/O rootfile " << rootfilename << endl;
	
	// read TH1F of data (h1000) from input
	TH1F* h_dati = (TH1F*)fin->Get("h1000");
	// read TH1F of efficiency from input
	TH1F* h_MCsel = (TH1F*)fin->Get("h200");
	TH1F* h_MCtot = (TH1F*)fin->Get("h100");
	// read Smearing Matrix
	TH2F* h_MCSmMat = (TH2F*)fin->Get("h700");
	// get DAT binning
	const int DATnbin = h_dati->GetNbinsX();
	TAxis *DATaxis = h_dati->GetXaxis();
	const double DATxmin = DATaxis->GetXmin();
	const double DATxmax = DATaxis->GetXmax();
	// get MC binning
	const int MCnbin = h_MCSmMat->GetNbinsX();
	if(MCnbin != h_MCSmMat->GetNbinsY() && verbosity >= ERROR_V)cout << "root2dat_v06 [ERROR] SmMat bin number mismatch!" << endl;
	if(MCnbin != h_MCsel->GetNbinsX() && verbosity >= ERROR_V)cout << "root2dat_v06 [ERROR] MC-sel: bin number mismatch!" << endl;
	if(MCnbin != h_MCtot->GetNbinsX() && verbosity >= ERROR_V)cout << "root2dat_v06 [ERROR] MC-tot: bin number mismatch!" << endl;
	TAxis * MCaxis = h_MCsel->GetXaxis();
	const double MCxmin = MCaxis->GetXmin();
	if(MCxmin != h_MCSmMat->GetXaxis()->GetXmin() && verbosity >= ERROR_V)cout << "root2dat_v06 [ERROR] SmMat Xmin mismatch!" << endl;
	if(MCxmin != h_MCSmMat->GetYaxis()->GetXmin() && verbosity >= ERROR_V)cout << "root2dat_v06 [ERROR] SmMat Ymin mismatch!" << endl;
	if(MCxmin != h_MCtot->GetXaxis()->GetXmin() && verbosity >= ERROR_V)cout << "root2dat_v06 [ERROR] MC-tot: Xmin mismatch!" << endl;
	const double MCxmax = MCaxis->GetXmax();
	if(MCxmax != h_MCSmMat->GetXaxis()->GetXmax() && verbosity >= ERROR_V)cout << "root2dat_v06 [ERROR] SmMat Xmax mismatch!" << endl;
	if(MCxmax != h_MCSmMat->GetYaxis()->GetXmax() && verbosity >= ERROR_V)cout << "root2dat_v06 [ERROR] SmMat Ymax mismatch!" << endl;
	if(MCxmax != h_MCtot->GetXaxis()->GetXmax() && verbosity >= ERROR_V)cout << "root2dat_v06 [ERROR] MC-tot: Xmax mismatch!" << endl;
	// print binning info
	if(verbosity >= DEBUG_V){
		cout << "root2dat_v06 [DEBUG]: DAT binning is " << DATxmin << "<-->" << DATxmax << " with " << DATnbin << " bins." << endl;
		cout << "root2dat_v06 [DEBUG]: MC binning is << " << MCxmin << "<-->" << MCxmax << " with " << MCnbin << " bins." << endl;
	}
	
	// compute efficiency + efficiency uncertainty [h2000 / h2200]
	TH1D* h_MCeffi = new TH1D("h2000","Efficiency MC-sel/MC-tot",MCnbin,MCxmin,MCxmax);
	TH1D* h_DATeffi = new TH1D("h2200","Efficiency MC-sel/MC-tot (DAT bin)",DATnbin,DATxmin,DATxmax);
	double Nsel_tmp = 0.;
	double Ntot_tmp = 0.;
	double Effi_tmp = 0.;
	for(i=1;i<=MCnbin;i++)
	{
		Nsel_tmp = (double)h_MCsel->GetBinContent(i);
		Ntot_tmp = (double)h_MCtot->GetBinContent(i);
		Effi_tmp = Nsel_tmp / Ntot_tmp;
		if(verbosity >= DEBUG_V)
		{
			cerr << "[h2000] #N Effi / Err on bin_" << i << ":\t" << Effi_tmp << "/" << sqrt(Effi_tmp*(1.-Effi_tmp)/Ntot_tmp) << endl;
		}
		h_MCeffi->SetBinContent(i,Effi_tmp);
		h_MCeffi->SetBinError(i,sqrt(Effi_tmp*(1.-Effi_tmp)/Ntot_tmp));
	}
	h_MCeffi->Write("",TObject::kOverwrite);
	// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
	
	// now if I just rebin the original histos... then I can repeat the procedure
	h_MCsel->Rebin(4);
	h_MCtot->Rebin(4);
	if(MCnbin%4!=0 && verbosity >= DEBUG_V)cerr << "[root2dat_v06] ERROR: MCnbin not a multiple of 4!" << endl;
	for(i=1;i<=MCnbin/4;i++)
	{
		Nsel_tmp = (double)h_MCsel->GetBinContent(i);
		Ntot_tmp = (double)h_MCtot->GetBinContent(i);
		Effi_tmp = Nsel_tmp / Ntot_tmp;
		if(verbosity >= DEBUG_V)
		{
			cerr << "[h2200] #N Effi / Err on bin_" << i << ":\t" << Effi_tmp << "/" << sqrt(Effi_tmp*(1.-Effi_tmp)/Ntot_tmp) << endl;
		}
		h_DATeffi->SetBinContent(i,Effi_tmp);
		h_DATeffi->SetBinError(i,sqrt(Effi_tmp*(1.-Effi_tmp)/Ntot_tmp));
	}
	h_DATeffi->Write("",TObject::kOverwrite);
	// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
	// h_DATeffi->Draw("E");
	
	// 4-PI BACKGROUND HISTOGRAM [HARD-CODED]
	TH1D* h_DAT4pi = new TH1D("h5000","4pi background (DAT bin)",DATnbin,DATxmin,DATxmax);
	double quattropioni_val[DATnbin]; for(int kkk=0; kkk<DATnbin; kkk++)quattropioni_val[kkk]=0.; // this is initialized to 0.
	double quattropioni_err[DATnbin]; for(int kkk=0; kkk<DATnbin; kkk++)quattropioni_err[kkk]=0.; // this also
	double cut_mpipi = TRKMASSCUT; // catching the mass cut (if != from default -> correct 4pi bkg)
	if(cut_to_vary == 0)cut_mpipi+=nsigma*TRKMASSRES;
	quattropioni_val[0] = (N4PI+nSigma4pi*RES4PI)*cut_mpipi/(double)TRKMASSCUT;
	quattropioni_err[0] = RES4PI; // error on nominal bkg (do not think to much about it..)
	cerr << "[root2dat_v06] INFO: TRKMASSCUT = " << TRKMASSCUT << endl;
	cerr << "[root2dat_v06] INFO: TRKMASSRES = " << TRKMASSRES << endl;
	cerr << "[root2dat_v06] INFO: cut_mpipi = " << cut_mpipi << endl;
	cerr << "[root2dat_v06] INFO: nSigma4pi = " << nSigma4pi << endl;
	if(nSigma4pi != 0 && cut_to_vary == 0){
		cerr << "[root2dat_v06] FATAL: variation of mpp cut AND 4pi, this is impossible! Check input parameters." << endl;
		exit(0);
	}
	cerr << "[root2dat_v06] WARNING: 4pi value = N4PI 4pi error = RES4PI" << endl;
	for(int kkk = 1; kkk <= DATnbin; kkk++)
	{
		h_DAT4pi->SetBinContent(kkk,quattropioni_val[kkk-1]);
		h_DAT4pi->SetBinError(kkk,quattropioni_err[kkk-1]);
	}
	h_DAT4pi->Write("",TObject::kOverwrite);
	// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
	
	// define output file stream to use later to load from .dat files
	ifstream IFS;
	
	// DATA-MC CORRECTION HISTOGRAM
	TH1D* h_DATdatamontecarlo = new TH1D("h3000","DATA-MC correction (DAT bin)",DATnbin,DATxmin,DATxmax);
	IFS.open( Form("%s/%s",folder,datamc_correction_filename) );
	if(!IFS.is_open())cerr << "[root2dat_v06] ERROR: Could not open file " << Form("%s/%s",folder,datamc_correction_filename) << endl;
	if(verbosity > INFO_V)cerr << "[root2dat_v06] INFO: DATA-MC (values) FILE is " << Form("%s/%s",folder,datamc_correction_filename) << endl;
	if(verbosity > INFO_V)cerr << "[root2dat_v06] INFO: DATA-MC (errors) FILE is " << Form("%s/%s",folder,datamc_correction_error_filename) << endl;
	double datamc_correction_factor[DATnbin]; for(int kkk=0; kkk<DATnbin; kkk++)datamc_correction_factor[kkk]=0.;
	double datamc_correction_error[DATnbin];  for(int kkk=0; kkk<DATnbin; kkk++)datamc_correction_error[kkk]=0.; // init all to 0.
	for(int kkk = 0; kkk < 50; kkk++)IFS >> datamc_correction_factor[kkk]; // load DAT-MC correction factor
	IFS.close();
	IFS.open( Form("%s/%s",folder,datamc_correction_error_filename) );
	if(!IFS.is_open())cerr << "[root2dat_v06] ERROR: Could not open file " << Form("%s/%s",folder,datamc_correction_error_filename) << endl;
	for(int kkk = 0; kkk < 50; kkk++)IFS >> datamc_correction_error[kkk]; // load DAT-MC correction error
	IFS.close();
	for(int kkk = 0; kkk < DATnbin; kkk++)
	{    // apply loaded corrections into TH1D
		h_DATdatamontecarlo->SetBinContent(kkk+1,datamc_correction_factor[kkk]);
		h_DATdatamontecarlo->SetBinError(kkk+1,datamc_correction_error[kkk]);
	}
	h_DATdatamontecarlo->Write("",TObject::kOverwrite);
	// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
	
	// REGENERATION CORRECTION HISTOGRAM
	TH1D* h_DATregen = new TH1D("hrege","REGENERATION correction (DAT bin)",DATnbin,DATxmin,DATxmax);
	double regen_default[DATnbin]; for(int kkk=0; kkk<DATnbin; kkk++)regen_default[kkk]=0.; // init all to 0.
	IFS.open(Form("%s/%s",folder,regen_filename)); // open file to load from it
	if(!IFS.is_open()){
		cerr << "[root2dat_v06] FATAL: Could not open file " << Form("%s/%s",folder,regen_filename) << endl;
		exit(0);
	}
	if(verbosity > INFO_V)cerr << "[root2dat_v06] INFO: REGENERATION FILE is " << Form("%s/%s",folder,regen_filename) << endl;
	int num = 0;
	for(int kkk=0;kkk<NBIN_FIT;kkk++)
	{
		IFS >> num >> regen_default[num]; // yes, this syntax works
		if(num != kkk)cerr << "[root2dat_v06] ERROR! REGENERATION FILE MISALIGNMENT!" << endl;
	}
	IFS.close();
	double myoldregen[12]={1.00034,1.00031,1.00052,1.00112,1.00277,1.00742,1.02057,1.05169,1.05605,1.04286,1.03072,1.02111};
	bool regen_is_different = false;
	for(int mybin=0;mybin<NBIN_FIT;mybin++)
		if(regen_default[mybin]!=myoldregen[mybin]){
			cerr << "FATAL: I am quite sure there is some mismatch btw regeneration coefficients loaded and correct ones.. please check!" << endl;
			for(int iii=0;iii<NBIN_FIT;iii++)cerr << "[0] " << regen_default[iii] << " <-regen_default | myoldregen -> " << myoldregen[iii] << endl;
			exit(0);
		}
	// apply modification if required by systematics analysis....
	for(int mybin=1;mybin<=NBIN_FIT; mybin++)regen_default[mybin-1]*=(1.00 + Kregeneration);
	cerr << "[root2dat_v06] INFO: Applied regen_default * " << 1.00 + Kregeneration << endl;
	for(int mybin=1; mybin<=NBIN_FIT; mybin++) // paste result into TH1D
	{
		h_DATregen->SetBinContent(mybin,regen_default[mybin-1]);
	}
	h_DATregen->Write("",TObject::kOverwrite);
	// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
	
	if(verbosity >= DEBUG_V)
	{
		cout << "root2dat_v06 [DEBUG]: -------" << endl;
		cout << "root2dat_v06 [DEBUG]: From now on, start doing what v05 did..." << endl;
	}
	
	fin->Close("R"); // with R it is more safe to close a file (recursively closes related things and frees memory
}
