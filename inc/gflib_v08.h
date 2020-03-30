/*
v08 04-02-2017 -> 03.05.2017 --> 04.09.2017
    [test] Modify class fitresult to become compliant with multiple parameters
    CUTARRAY[8] --> BUGFIX inserito RES4PI invece di sqrt(N4PI)
    CUTARRAY[8] --> BUGFIX inserito valore del fit sidebands mpipi N4pi = 78 +- 9
    ---> BUGFIX 2 --> inserito valore del fit * 1/2 perché nelle sideband mpipi entrano il doppio degli eventi: uno per vertice
 
v07 18-11-2016 (12-12-2016 ok)
    [ok] Insert new class to compare parametric 2d plots

v06 14-07-2016
    [ok] is_selected, GetCutResult, GetCutUnits, GetCutResolutions updated to REMOVE pz cut
         NCUT updated to 10 -> 2 free slot for 4pi and regen bkg variation
	 --> Someday this non-proper use of NCUT has to be fixed somehow

v05 18-06-2016
    [ok] New class GFloop to load events automatically looping through files
    09-07-2016
    Included cut on pull to comply with v04
    [NO] DeltaT12reco (without kin fit) (some comments in the code)
    [ok] GetCutResult -> writes an array cutresult[NCUT] with 0 if FAILED 1 if PASSED for each cut

v04 (sviluppo in corso...)
    30-05-2016
    13-06-2016 update
    [ok] Function to smear MC reco p(trk) by f*gauss(m=0,sigma=1)
         Inserted random generator as global variable randgen2
    [ok] Array of cut names for plots available with FillCutNames
    [ok] Inserted file2vector and vector2file functions from root2dat_v03
    [ok] Inserted constant FSMEARMC = 2 x 10**-3
    [ok] Included DeltaT12 (data and MC)
    09-07-2016 -> update cut pull
    
v03 (definitiva)
    19-05-2016
    [ok] function to  SetBranchAddresses to global variables
    [ok] function to apply cuts (vector dependent)
    [ok] function to SetCutArray from config file
    [ok] function to SetDefaultCutArray

v02
  09-05-2016 is_selected(event *e) <--- inserted pointer to struct to load my version of lorentz variables my<LORENTZ_VAR_NAME>
  This function loads all my variables.

*/


#ifndef GFLIB
#define GFLIB

#define NCUT 10// number of cuts
#define TRKMASSCUT 5 // 5 MeV cut on |m(pipi)-mK0|
#define EPMISSCUT 10 // 10 MeV cut on EPmiss
#define MISSMASSCUT_LOW -50 // MeV**2
#define MISSMASSCUT_HIGH 10 // MeV**2
#define PSTARCUT 10 // MeV
#define GLOBALFITCHISQUARE 30 //-logL > 30 discard
#define GLOBALFITPULL 10 // reject |in-out|/sigma(out) > 10
#define COSPIONCUT -0.975 // reject too open pions
//#define PZCUT 2 // regect K (dal pz-con-segno più alto) if pz<2 MeV
#define N4PI 78 // good old 54 || old 78 NOTICE! this number should be 1/2 of that in the sideband fit (78) (there are 2 vertices)
#define REGENDEFAULT 1.

#define TRKMASSRES 0.9
#define EPMISSRES 1.3
#define MISSMASSRES_LOW 4.0
#define MISSMASSRES_HIGH 4.0
#define PSTARRES 1.3
#define GLOBALFITCHISQUARERES 2.0
#define GLOBALFITPULLRES 1.0 // <------------------- this cut was not varied
#define COSPIONRES 0.001
//#define PZRES 2.0
#define RES4PI 4. //7.348 // modified! = sqrt(54) // modified!! = 0.5*8.89 (from sideband fit)
#define REGENEXTRAPERCENTAGE 12

#define FSMEARMC 0.002


// set global variables to be fetched event by event
UShort_t NRun;         //number of run
UInt_t NEv;            //number of event (given the run)
Int_t NPar;            //number of parameters of global fit
UShort_t tagMask;      //?.. MC classification
Float_t pKS[2][3];     //KS track momentum
Float_t pKL[2][3];     //KL track momentum
Float_t KSLMass[2];    //m(KS_L) - MK0 from tracks in pion hp
Float_t KSLEPMiss[2];  //KS (KL) sqrt(emiss**2+pmiss**2) defined in k2id-18
Float_t KSLMissMass[2];//KS (KL) emiss**2-pmiss**2
Float_t KSLEstar[2];   //p*S/L - p0*
Float_t FIn[3][2];     //Input parameters of global fit (ADS)
Float_t FOut[3][2];    //Results of fit to the vertex (ADS) [#par][0/1] 0: central value [1]: error #par= 0-> lKS 1->lKL
Float_t FChi2;         //Chi2 from global fit  
Float_t pKStot[3];     // needed for Phi kinematics
Float_t pKLtag[3];     //   "     "   "      "

// MC variables
Float_t pKSMC[2][3];   //
Float_t pKLMC[3][3];   // <--- problema, perché è di [3][3]? è così anche nel .ntu
Float_t pKStotMC[3];
Float_t pKLtotMC[3];
UInt_t KSwordMC;
UInt_t KLwordMC;
Float_t xPhiMC[3];
Float_t xKSMC[3];
Float_t xKLMC[3];

//Random generator
TRandom2 *randgen2 = new TRandom2();


class event
{
 public:
  TVector3 p3KS[2],p3KL[2];         // pions 3-momenta
  TLorentzVector p4KS,p4KL;         // kaons 4-momenta from tracks pion hp
  TLorentzVector p4KSmiss,p4KLmiss; //
  TLorentzVector p4Phi;             // phi 4-momentum (event-by-event)
  Double_t BetaGammaS,BetaGammaL;
  Float_t DistKS,DistKL;
  Float_t DeltaT; // [ns]
  Double_t DeltaT12; // t1 - t2 [ns] pz(K1) > pz(K2) (after kin fit)
  //Double_t DeltaT12reco; // t1 - t2 [ns] pz(K1) > pz(K2) (NO kin fit)
  
  // true variables end by MC
  TVector3 p3KSMC[2],p3KLMC[2];         // pions 3-momenta (MC)
  TLorentzVector p4KSMC,p4KLMC;         // kaons 4-momenta (MC)
  TLorentzVector p4KSmissMC,p4KLmissMC; //
  TLorentzVector p4PhiMC;               // phi 4-momentum (MC)
  
  Double_t BetaGammaSMC,BetaGammaLMC;   
  Float_t DistKSMC,DistKLMC;
  Float_t DeltaTMC; // absolute value
  Double_t DeltaT12MC; // t1-t2 (could be <0)
  
  double myKSLMass[2];
  double myKSLEPMiss[2];
  double myKSLEstar[2];
  double myKSLMissMass[2];
  
 public:
  void fill_data(); // fill reco variables
  void fill_mc();   // fill true variables
};

class GFloop
{
  char MCmode;
  char DATmode;
  int nfile; // number of file 1-419
  int nev; // number of event 0 - (eventTOT-1)
  int eventTOT; // total number of events in <filename>
  TFile* myfile;
  TTree* mytree;
  char filename[100];
public:
  GFloop(char mode); // constructor mode = M/D for MC/DATA
  int LoadNext(int verbose); // fill global variables with ntuple variables
  //                              returns 0 if END of files
  //                              verbose = 0 (no output) | 1 (some output)
};

//----------------------------------------------------------
// Version 07

class plot2D
{
  int N;
  int ToggleErr; // 1 if ex ey filled or not
 public:
  
  Double_t x[1000];
  Double_t ex[1000];
  Double_t y[1000];
  Double_t ey[1000];
  Double_t p0;
  Double_t p1;
  TGraph *plot;
  void Load(const char* filename,int nrow,int ToggleErr);
  void LoadVector(double *input_x,double *input_y,int num);
  void DumpN();
  void Dump();
  void RescaleXY(double scalex,double scaley);
  void ShiftX(double shift);
  Double_t GetXmin();
  Double_t GetXmax();
  // integrate with trapezoidal method from xmin to xmax
  Double_t TrapezoidIntegrate(double xmin,double xmax); 
  Double_t BinSum(double xmin,double xmax,double binwidth);
  void Replot(); // updates the plot variable
};
void plot2D::Load(const char* filename,int nrow,int errortoggle){
  if(errortoggle!=0 && errortoggle!=1){
    cout << "Plot2D warning: invalid errortoggle. Setting to 0." << endl;
    ToggleErr = 0.;
  }else{
    ToggleErr = errortoggle;
  }
  ifstream myfile;
  int i;
  myfile.open(filename,ios::in);
  if(!myfile.is_open()){
    cout << "plot2d::Load ERROR! Unable to open file " << filename << endl;
  }
  if(errortoggle==0){ // no errorbar
    for(i=0;i<nrow;i++)myfile >> x[i] >> y[i];
    plot = new TGraph(nrow,x,y);
  }else if(errortoggle==1){// with errorbar
    for(i=0;i<nrow;i++)myfile >> x[i] >> ex[i] >> y[i] >> ey[i];
    }
  myfile.close();
  N = nrow;
}
void plot2D::LoadVector(double *input_x,double *input_y,int num)
{
  int i;
  N = num;
  ToggleErr = 0.;
  for(i=0;i<N;i++){
    x[i] = input_x[i];
    y[i] = input_y[i];
  }
  plot = new TGraph(N,x,y);
}
void plot2D::DumpN()
{
  cout << "plot2d::DumpN [DEBUG]" << endl;
  cout << "N = " << N << endl;
}
void plot2D::Dump()
{
  int i;
  for(i=0;i<N;i++){
    cout << "["<< i << "] " << endl;
    cout << "   [x] = " << x[i] << endl;
    cout << "   [y] = " << y[i] << endl;
  }
}
void plot2D::RescaleXY(double scalex,double scaley)
{
  int i;
  for(i=0;i<N;i++){
    x[i]*=scalex;
    y[i]*=scaley;
  }
  plot->~TGraph();
  plot = new TGraph(N,x,y);

}
void plot2D::ShiftX(double shift){
  int i;
  for(i=0;i<N;i++){
    x[i]+=shift;
  }
  plot->~TGraph();
  plot = new TGraph(N,x,y);
}
Double_t plot2D::GetXmin()
{
  Double_t tmp = x[0];
  int i;
  for(i=1;i<N;i++){
    if(x[i]<tmp)tmp=x[i];
  }
  return tmp;
}
Double_t plot2D::GetXmax(){
  Double_t tmp2 = x[0];
  int i;
  for(i=1;i<N;i++){
    if(x[i]>tmp2)tmp2=x[i];
  }
  return tmp2;
}
Double_t plot2D::TrapezoidIntegrate(double xmin,double xmax)
{
  Double_t min = GetXmin();
  Double_t max = GetXmax();
  // security checks
  if(xmin < min || xmax > max || xmin > xmax){
    cout << "plot2D::TrapezoidIntegrate() FATAL ERROR! Invalid argument!" << endl;
    cout << "Should be (xmin,xmax) = "<< min << " , " << max << endl;
    return 0.;
  }
  int i;
  for(i=1;i<N;i++){
    if(x[i]<x[i-1]){
      cout << "plot2D::Integrate ERROR! x[i] values not ordered. Input again." << endl;
      return 0.;
    }
  }
  // find index of array corresponding to integration area
  int imin,imax; // imin = just after xmin / imax = just before xmax
  int flagmin=0;
  int flagmax=0;
  for(i=0;i<N;i++){
    if(x[i]>xmin && flagmin==0){
      imin=i;
      flagmin++;
    }
    if(x[i]>xmax && x[i-1]<=xmax && flagmax==0){
      imax=i;
      flagmax++;
    }
  }
  // consistency check
  if(imin > imax){
    cout << "plot2D::TrapezoidIntegrate() FATAL CONSISTENCY VIOLATION! BUG FOUND." << endl;
    return 0.;
  }
  // compute integral
  double b,B,h;
  Double_t integral = 0.;
  for(i=imin;i<imax;i++){
    h = x[i+1]-x[i];
    b = y[i];
    B = y[i+1];
    integral+=(b+B)*h*0.5;
  }
  // extrema
  // low
  b = y[imin-1]+(xmin-x[imin-1])*(y[imin]-y[imin-1])/(x[imin]-x[imin-1]);
  B = y[imin];
  h = x[imin]-xmin;
  integral+=(b+B)*h*0.5;
  // up
  b = y[imax];
  B = y[imax]+(xmax-x[imax])*(y[imax+1]-y[imax])/(x[imax+1]-x[imax]);
  h = xmax-x[imax];
  integral+=(b+B)*h*0.5;
  return integral;
}
Double_t plot2D::BinSum(double xmin,double xmax,double binwidth)
{
  // cut & paste from Integrate------------------------------------
  Double_t min = GetXmin(); // get 1st bin central position
  Double_t max = GetXmax(); // get last bin central position
  // security checks
  if(xmin < min - binwidth/2. || xmax > max + binwidth/2.|| xmin > xmax){
    cout << "plot2D::BinSum() FATAL ERROR! Invalid argument!" << endl;
    cout << "Should be (xmin,xmax) = "<< min - binwidth/2. << " , " << max + binwidth/2. << endl;
    return 0.;
  }
  int i;
  for(i=1;i<N;i++){
    if(x[i]<x[i-1]){
      cout << "plot2D::BinSum ERROR! x[i] values not ordered. Input again." << endl;
      return 0.;
    }
  }
  // find index of array corresponding to integration area
  int imin,imax; // imin(imax) = index of bin where xmin(xmax) belong
  // In case of ambiguity, the bin on the left.
  int flagmin=0;
  int flagmax=0;

  int debugflag[2];
  for(i=0;i<N;i++){
    // DEBUG
    //if(x[i]-binwidth/2.< xmax){debugflag[0]=1;}else{debugflag[0]=0;}
    //if(xmax <= x[i]+binwidth/2.){debugflag[1]=1;}else{debugflag[1]=0;}
    //cout << "DEBUG: " << x[i]-binwidth/2. << " -- " << x[i]+binwidth/2. << " FLAGS: " << debugflag[0] << debugflag[1] << endl;

    if(x[i]-binwidth/2.< xmin && xmin <= x[i]+binwidth/2.){
      imin=i;
      flagmin++;
    }
    if(x[i]-binwidth/2.<xmax && xmax <= x[i]+binwidth/2.){
      imax=i;
      flagmax++;
    }
  }

  // consistency check
  if(imin > imax || flagmin!=1 || flagmax != 1){
    cout << "plot2D::BinSum() FATAL CONSISTENCY VIOLATION! BUG FOUND." << endl;
    cout << "-- imin "<< imin << endl;
    cout << "-- imax "<< imax << endl;
    cout << "-- flagmin " << flagmin << endl;
    cout << "-- flagmax " << flagmax << endl;
    cout << "-- xmin / xmax " << xmin << " / " << xmax << endl;
    cout << "-- x[0] " << x[0] << endl;
    cout << "-- x[N-1] " << x[N-1] << endl;
    return 0.;
  }
  if(imin < 0){cout << "plot2D::BinSum() Error: imin < 0" << endl;}
  if(imax > N){cout << "plot2D::BinSum() Error: imax > MAX = " << N << endl;}

  // last consistency check on binwidth
  double w;
  for(i=imin;i<imax;i++){
    w = x[i+1]-x[i];
    if( pow(w/binwidth - 1.,2) > 1E-6 ){
      cout << "plot2D::BinSum() FATAL CONSISTENCY CHECK VIOLATION! Bin width non compatible with data." << endl;
      cout << "Setting integral to 0." << endl;
      return 0.;
    }
  }
  
  Double_t integral = 0.;
  double delta;
  // compute left border integral
  delta = x[imin]+binwidth/2.- xmin;
  // DEBUG
  //cout << "DELTA(left)/binwidth = " << delta/binwidth << endl;
  //cout << "imin = " << imin << endl;
  if(delta < 0.){ // consistency check
    cout << "plot2D::BinSum() FATAL ERROR! delta(SX) < 0" << endl;
    return 0.;
  }
  integral+=y[imin]*delta/binwidth;

  // compute right border integral
  delta = xmax - (x[imax]-binwidth/2.);
  // DEBUG
  //cout << "DELTA(right)/binwidth = " << delta/binwidth << endl;
  //cout << "imax = " << imax << endl;
  if(delta < 0.){
    cout << "plot2D::BinSum() FATAL ERROR! delta(DX) < 0" << endl;
    return 0.;
  }
  integral+=y[imax]*delta/binwidth;
  //cout << "DEBUG: "<< x[imin] << " ... " << x[imax] << endl;
  //cout << "DEBUG: integral = " << integral << endl;
  
  // compute the rest
  for(i=imin+1;i<imax;i++){
    integral+=y[i];
  }
  return integral;
}
void plot2D::Replot(){
  plot->~TGraph();
  plot = new TGraph(N,x,y);
}

Double_t I_zetaSL(Double_t *x,Double_t *par)
{
  // x = variable vector
  // par = parameter vector
  Double_t xx = x[0];
  Double_t zetaSL = par[0];
  // DEBUG ------ INCOMPLETO
  return exp(-1.1*xx);

}

// Function for systematics <--- not useful anymore
// fills cutarray with: cut + nsigma resolutions
int UnpackCutword(int ncut,int cutword,double* nsigma)
{
  int i,mask,tmp,res;
  if(cutword > pow(10,ncut))
    {
      cout << "ERROR! Cutword and ncut not compatible." << endl;
      cout << "Setting all nsigma[i] = 0" << endl;
      for(i=0;i<ncut;i++)nsigma[i]=0.;
      return 0;
    }
  tmp = cutword;
  for(i=0;i<ncut;i++)
    {
      mask = pow(10,ncut-1-i);
      res = tmp/mask;
      tmp = tmp - res*mask;
      nsigma[i] = res;
    }
  return 1;
}

// based on global variables defined on top of this lib fills cutarray
// with the values of cuts shifted by nsigma resolutions
// loaded from configfile (.dat): one line = one element of nsigma
void SetCutArray(const char* configfile,double* cutarray)
{
  double nsigma[NCUT];
  ifstream config;
  int i;
  config.open(configfile,ios::in);
  if(config.is_open())
    {
      for(i=0;i<NCUT;i++){
	config>>nsigma[i];
	if(!config.good() && !(i-NCUT == -1 && config.eof()))
	  {
	    cout << "SetCutArray: some problem reading line " << i << " occurred." << endl;
	    nsigma[i] = 0.;
	  }
	cout << "SetCutArray: nsigma[" << i << "] = " << nsigma[i]<< endl;
      }
    }
  else
    {
      cout << "SetCutArray: unable to open configuration file " << configfile << endl;
      cout << "SetCutArray: setting all nsigma[*] = 0" << endl;
      for(i=0;i<NCUT;i++)nsigma[i]=0.;
    }
  config.close();
  cutarray[0] = TRKMASSCUT + nsigma[0]*TRKMASSRES;
  cutarray[1] = EPMISSCUT + nsigma[1]*EPMISSRES;
  cutarray[2] = MISSMASSCUT_LOW + nsigma[2]*MISSMASSRES_LOW;
  cutarray[3] = MISSMASSCUT_HIGH + nsigma[3]*MISSMASSRES_HIGH;
  cutarray[4] = PSTARCUT + nsigma[4]*PSTARRES;
  cutarray[5] = GLOBALFITCHISQUARE + nsigma[5]*GLOBALFITCHISQUARERES;
  cutarray[6] = GLOBALFITPULL + nsigma[6]*GLOBALFITPULLRES;
  cutarray[7] = COSPIONCUT + nsigma[7]*COSPIONRES;
  // BACKGROUND VARIATION [> v06]
  cutarray[8] = N4PI + nsigma[8]*RES4PI;
  cutarray[9] = REGENDEFAULT + nsigma[9]*REGENEXTRAPERCENTAGE/100.;
}

void SetDefaultCutArray(double* cutarray)
{
  cutarray[0] = TRKMASSCUT;
  cutarray[1] = EPMISSCUT;
  cutarray[2] = MISSMASSCUT_LOW;
  cutarray[3] = MISSMASSCUT_HIGH;
  cutarray[4] = PSTARCUT;
  cutarray[5] = GLOBALFITCHISQUARE;
  cutarray[6] = GLOBALFITPULL;
  cutarray[7] = COSPIONCUT;
  cutarray[8] = N4PI;
  cutarray[9] = REGENDEFAULT;
}

void GetCutUnits(char string[NCUT][10])
{
  sprintf(string[0],"MeV");
  sprintf(string[1],"MeV");
  sprintf(string[2],"MeV^{2}");
  sprintf(string[3],"MeV^{2}");
  sprintf(string[4],"MeV");
  sprintf(string[5]," ");
  sprintf(string[6]," ");
  sprintf(string[7]," ");
  sprintf(string[8],"events");  
  sprintf(string[9],"%%");
}

void GetCutResolutions(double* resarray)
{
  resarray[0] = TRKMASSRES;
  resarray[1] = EPMISSRES;
  resarray[2] = MISSMASSRES_LOW;
  resarray[3] = MISSMASSRES_HIGH;
  resarray[4] = PSTARRES;
  resarray[5] = GLOBALFITCHISQUARERES;
  resarray[6] = GLOBALFITPULLRES;
  resarray[7] = COSPIONRES;
  resarray[8] = RES4PI;
  resarray[9] = REGENEXTRAPERCENTAGE;
}

void FillCutNames(int ncut,char *string)
{
  switch(ncut){
  case 0:
    sprintf(string,"m(#pi+#pi-)");
    break;
  case 1:
    sprintf(string,"Emiss (+) Pmiss");
    break;
  case 2:
    sprintf(string,"Mmiss (lower cut)");
    break;
  case 3:
    sprintf(string,"Mmiss (upper cut)");
    break;
  case 4:
    sprintf(string,"Pstar");
    break;
  case 5:
    sprintf(string,"#chi^{2} kinematic fit ");
    break;
  case 6:
    sprintf(string,"Pull kinematic fit");
    break;
  case 7:
    sprintf(string,"cos[#theta_{#pi#pi}]");
    break;
  case 8:
    sprintf(string,"4#pi bkg");
    break;
  case 9:
    sprintf(string,"regeneration");
    break;
  default:
    cout << "FillCutNames ERROR! -> NCUT "<< ncut <<" not allowed." << endl;
    break;
  }
}

void FillCutSteps(int ncut,char *string)
{
  switch(ncut){
  case 0: // Mpipi
    sprintf(string,"step: #sigma = %.1f MeV",TRKMASSRES);
    break;
  case 1: // Emiss (+) Pmiss
    sprintf(string,"step: #sigma = %.1f MeV",EPMISSRES);
    break;
  case 2:// Mmiss (lower cut)
    sprintf(string,"step: #sigma = %.1f MeV^{2}",MISSMASSRES_LOW);
    break;
  case 3:// Mmiss (upper cut)
    sprintf(string,"step: #sigma = %.1f MeV^{2}",MISSMASSRES_HIGH);
    break;
  case 4:// P*
    sprintf(string,"step: #sigma = %.1f MeV",PSTARRES);
    break;
  case 5:// chi^{2} kinematic fit
    sprintf(string,"step: #sigma = %.1f",GLOBALFITCHISQUARERES);
    break;
  case 6:// Pull kinematic fit
    sprintf(string,"step: #sigma = %.1f",GLOBALFITPULLRES);
    break;
  case 7:// cos[#theta_{#pi#pi}]
    sprintf(string,"step: #sigma = %.3f",COSPIONRES);
    break;
  case 8:// 4#pi bkg
    sprintf(string,"step: #sigma = %.1f",RES4PI);
    break;
  case 9:// regeneration
    sprintf(string,"step: #alpha = %d %%",12);
    break;
  default:
    cout << "FillCutNames ERROR! -> NCUT "<< ncut <<" not allowed." << endl;
    break;
  }
}

void event::fill_data()
{
  int j;
  double tmp[4];
  // Kaon's 4-momenta from pion tracks -- p4KS p4KL
  for(j=0;j<4;j++)tmp[j]=0;
  for(j=0;j<3;j++)
    {
      p4KS[j]=pKS[0][j]+pKS[1][j];
      p4KL[j]=pKL[0][j]+pKL[1][j];
      tmp[0]+=pKS[0][j]*pKS[0][j];
      tmp[1]+=pKS[1][j]*pKS[1][j];
      tmp[2]+=pKL[0][j]*pKL[0][j];
      tmp[3]+=pKL[1][j]*pKL[1][j];
    }
  p4KS[3]=sqrt(MPION*MPION+tmp[0])+sqrt(MPION*MPION+tmp[1]);
  p4KL[3]=sqrt(MPION*MPION+tmp[2])+sqrt(MPION*MPION+tmp[3]);  
  // Pion's momenta from tracks -- p3KS p3KL
  for(j=0;j<2;j++)
    {
      p3KS[j].SetXYZ(pKS[j][0],pKS[j][1],pKS[j][2]);
      p3KL[j].SetXYZ(pKL[j][0],pKL[j][1],pKL[j][2]);
    }
  // Kaon's betagamma from p4K*
  BetaGammaS=p4KS.Beta()*p4KS.Gamma();
  BetaGammaL=p4KL.Beta()*p4KL.Gamma();
  // Flight distance from fit
  DistKS = FOut[0][0];
  DistKL = FOut[1][0];
  // DeltaT
  DeltaT = abs( DistKL/BetaGammaL - DistKS/BetaGammaS )/(C*TAUS);
  DeltaT12 =( DistKS/BetaGammaS - DistKL/BetaGammaL)/(C*TAUS);
  // need to define DistKS(L)reco variables
  // DeltaT12reco = ( DistKSreco/BetaGammaL - DistKLreco/BetaGammaS )/(C*TAUS);
  if(p4KL[2] > p4KS[2]){// if it is the other way around...
    DeltaT12*=-1.;
    //DeltaT12reco*=-1;
  }
  // Phi 4-momentum evaluation from K -- p4Phi (liblorentz.f)
  tmp[0]=0.;
  tmp[1]=0.;
  for(j=0;j<3;j++){
    p4Phi[j]=pKStot[j]+pKLtag[j];
    tmp[0]+=pKStot[j]*pKStot[j];
    tmp[1]+=pKLtag[j]*pKLtag[j];
  }
  p4Phi[3]=sqrt(tmp[0]+MK0*MK0)+sqrt(tmp[1]+MK0*MK0);
  // Kaon's missing momenta
  p4KSmiss = (p4Phi - p4KL) - p4KS; // 2 var uguali!!
  p4KLmiss = (p4Phi - p4KS) - p4KL;
  // track invariant mass
  if(p4KS.Mag2()<0. || p4KL.Mag2()<0.)cout<<"WARNING: invariant kaon mass < 0"<< endl;
  myKSLMass[0]=p4KS.Mag()-MK0;
  myKSLMass[1]=p4KL.Mag()-MK0;
}

void event::fill_mc()
{
  int j;
  double tmp[4];
  // Kaon's 4-momenta from pion tracks -- p4KSMC p4KLMC
  for(j=0;j<4;j++)tmp[j]=0;
  for(j=0;j<3;j++)
    {
      p4KSMC[j]=pKSMC[0][j]+pKSMC[1][j];
      p4KLMC[j]=pKLMC[0][j]+pKLMC[1][j];
      tmp[0]+=pKSMC[0][j]*pKSMC[0][j];
      tmp[1]+=pKSMC[1][j]*pKSMC[1][j];
      tmp[2]+=pKLMC[0][j]*pKLMC[0][j];
      tmp[3]+=pKLMC[1][j]*pKLMC[1][j];
    }
  p4KSMC[3]=sqrt(MPION*MPION+tmp[0])+sqrt(MPION*MPION+tmp[1]);
  p4KLMC[3]=sqrt(MPION*MPION+tmp[2])+sqrt(MPION*MPION+tmp[3]);  
  // Pion's momenta from tracks -- p3KSMC p3KLMC
  for(j=0;j<2;j++)
    {
      p3KSMC[j].SetXYZ(pKSMC[j][0],pKSMC[j][1],pKSMC[j][2]);
      p3KLMC[j].SetXYZ(pKLMC[j][0],pKLMC[j][1],pKLMC[j][2]);
    }

  // Phi 4-momentum evaluation from K -- p4PhiMC
  tmp[0]=0.;
  tmp[1]=0.;
  for(j=0;j<3;j++){
    p4PhiMC[j]=pKStotMC[j]+pKLtotMC[j];
    tmp[0]+=pKStotMC[j]*pKStotMC[j];
    tmp[1]+=pKLtotMC[j]*pKLtotMC[j];
  }
  p4PhiMC[3]=sqrt(tmp[0]+MK0*MK0)+sqrt(tmp[1]+MK0*MK0); // GF different from liblorentz!!
  // p4PhiMC[3]=sqrt(tmp[0]+tmp[1]+2*MK0*MK0); // liblorentz

  // Kaon's missing momenta
  p4KSmissMC = (p4PhiMC - p4KLMC) - p4KSMC;
  p4KLmissMC = (p4PhiMC - p4KSMC) - p4KLMC;
  // Decay length 
  DistKSMC=0.;
  DistKLMC=0.;
  for(j=0;j<3;j++)
    {
      DistKSMC+=(xPhiMC[j]-xKSMC[j])*(xPhiMC[j]-xKSMC[j]);
      DistKLMC+=(xPhiMC[j]-xKLMC[j])*(xPhiMC[j]-xKLMC[j]);
    }
  DistKSMC=sqrt(DistKSMC);
  DistKLMC=sqrt(DistKLMC);
  // BetaGamma
  BetaGammaSMC = p4KSMC.Beta()*p4KSMC.Gamma();
  BetaGammaLMC = p4KLMC.Beta()*p4KLMC.Gamma();

  // Time distance
  DeltaTMC=abs( DistKSMC/BetaGammaSMC - DistKLMC/BetaGammaLMC )/(TAUS*C);
  DeltaT12MC = (DistKSMC/BetaGammaSMC - DistKLMC/BetaGammaLMC)/(TAUS*C);
  if(p4KLMC[2] > p4KSMC[2]){// if it is the other way around...
    DeltaT12MC*=-1.;
  }
}

// Smears pKS[i],pKL[i] by f*gauss(m=0,s=1) p -> p*(1+f*gauss)
void smear_mc(double f)
{
  Double_t gauss;
  
  int i,j;
  for(i=0;i<2;i++){
    for(j=0;j<3;j++)
      {
	// pKStot and pKLtag NOT smeared (see fortran for reference)
	
	// GF: se richiamo sta cosa non è che la ri-inizializzo ogni volta ?
	// --> e quindi mi ritorna sempre la stessa cosa?
	gauss = randgen2->Gaus();
	pKS[i][j]*=(1+f*gauss);
	gauss = randgen2->Gaus();
	pKL[i][j]*=(1+f*gauss);
      }
  }
}

int is_signal()
{
  if((KSwordMC==2312 || KSwordMC==591873) && (KLwordMC==2312 || KLwordMC==591873))return 1;
  return 0;
}

int is_signal_ADS()
{
  if( (KSwordMC==2312 || KSwordMC==591873) && (KLwordMC==2312 || KSwordMC==591873) )return 1;
  return 0;
}



int is_selected(event *e,double* cutarray)
{
  //===========================
  //   NPar cut
  //===========================
  int status=-1;
  if(NPar==0)return status; // quando serve questo taglio? NEL MC serve spesso
  //===========================
  //   tagMask cut
  //===========================
  status=-2;
  if( ((tagMask>>5)&1)==0 && ((tagMask>>6)&1)==0 )return status; // if bit5 = 0 AND bit6=0 return (vtx requirement) GF??
  //===========================
  //   m(pipi) cut
  //===========================
  status=-3;

  if(abs(e->myKSLMass[0])>cutarray[0] || abs(e->myKSLMass[1])>cutarray[0])return status;
  //===========================
  //   Pstar cut
  //===========================
  status=-4;
  double EKstar = e->p4Phi.Mag2()/4. - MK0*MK0;
  if(EKstar<0.){
    cout << "WARNING: EKstar < 0. -- sqrt(s/4) = " <<endl;
    return status;
  }
  EKstar=sqrt(EKstar);
  // boost kaons back into phi frame
  TLorentzVector p4KSstar=e->p4KS;
  TLorentzVector p4KLstar=e->p4KL;
  p4KSstar.Boost(-e->p4Phi.BoostVector());
  p4KLstar.Boost(-e->p4Phi.BoostVector());
  double tmp[2];
  tmp[0]=0.;
  tmp[1]=0.;
  int j;
  for(j=0;j<3;j++)
    {
      tmp[0]+=p4KSstar[j]*p4KSstar[j];
      tmp[1]+=p4KLstar[j]*p4KLstar[j];
    }
  tmp[0]=sqrt(tmp[0]);
  tmp[1]=sqrt(tmp[1]);
  e->myKSLEstar[0] = tmp[0]-EKstar;
  e->myKSLEstar[1] = tmp[1]-EKstar;
  if(abs(e->myKSLEstar[0]) > cutarray[4] || abs(e->myKSLEstar[1]) > cutarray[4] )return status;
  //===========================
  //   sqrt(Emiss**2+Pmiss**2) cut
  //===========================
  status=-5;
  // calculate EPmiss
  tmp[0]=0.;
  tmp[1]=0.;
  for(j=0;j<4;j++)
    {
      tmp[0]+=e->p4KSmiss[j]*e->p4KSmiss[j];
      tmp[1]+=e->p4KLmiss[j]*e->p4KLmiss[j];
    }
  e->myKSLEPMiss[0]=sqrt(tmp[0]);
  e->myKSLEPMiss[1]=sqrt(tmp[1]);
  if(e->myKSLEPMiss[0]>cutarray[1] || e->myKSLEPMiss[1]>cutarray[1])return status;
  //===========================
  //   missing mass cut
  //===========================
  status=-6;
  // calculate EPmissmass
  e->myKSLMissMass[0]=e->p4KSmiss.Mag2();
  e->myKSLMissMass[1]=e->p4KLmiss.Mag2();
  if( e->myKSLMissMass[0]<cutarray[2] || e->myKSLMissMass[0]>cutarray[3] || e->myKSLMissMass[1]<cutarray[2] || e->myKSLMissMass[1]>cutarray[3] )return status;
  //===========================
  //   kinematic fit cut (chisquare & pulls)
  //===========================  
  status=-7;
  if(FChi2>cutarray[5] || abs(FOut[0][0]-FIn[0][0])/FOut[0][1]>cutarray[6] || abs(FOut[1][0]-FIn[1][0])/FOut[1][1]>cutarray[6])return status;
  //===========================
  //   pions angle cut
  //===========================   
  status=-8;
  if(e->p3KS[0]*e->p3KS[1]/(e->p3KS[0].Mag()*e->p3KS[1].Mag()) < cutarray[7] || e->p3KL[0]*e->p3KL[1]/(e->p3KL[0].Mag()*e->p3KL[1].Mag()) < cutarray[7])return status;
	
  //===========================
  //   ordered kaon z-momentum cut
  //===========================
	/*
  status=-9; // N.B. al posto di 2. nel taglio ci andava cutarray[8] prima che diventasse il taglio dei 4 pioni...
  if( (e->p4KS[2]>e->p4KL[2] && e->p4KS[2]< 2.) || (e->p4KS[2]<e->p4KL[2] && e->p4KL[2]< 2.) )return status;
  */
  // if passed all the cuts...
  return 1;
}

int GetCutResult(event* e, double* cutarray, double* cutresult)
{
  int i;
  for(i=0;i<NCUT;i++)cutresult[i]=1; // passed =1 not passed =0
  int status=-1; // increment if any cut is NOT passed
  //===========================
  //   NPar cut
  //===========================
  if(NPar==0)status++; // quando serve questo taglio?
  //===========================
  //   tagMask cut
  //===========================
  if( ((tagMask>>5)&1)==0 && ((tagMask>>6)&1)==0 )status++; // if bit5 = 0 AND bit6=0 return (vtx requirement) GF??
  //===========================
  //   m(pipi) cut
  //===========================
  if(abs(e->myKSLMass[0])>cutarray[0] || abs(e->myKSLMass[1])>cutarray[0]){
    status++;
    cutresult[0]=0;
  }
  //===========================
  //   Pstar cut
  //===========================
  double EKstar = e->p4Phi.Mag2()/4. - MK0*MK0;
  if(EKstar<0.){
    cout << "WARNING: EKstar < 0. -- sqrt(s/4) = " <<endl;
    status++;
  }
  EKstar=sqrt(EKstar);
  // boost kaons back into phi frame
  TLorentzVector p4KSstar=e->p4KS;
  TLorentzVector p4KLstar=e->p4KL;
  p4KSstar.Boost(-e->p4Phi.BoostVector());
  p4KLstar.Boost(-e->p4Phi.BoostVector());
  double tmp[2];
  tmp[0]=0.;
  tmp[1]=0.;
  int j;
  for(j=0;j<3;j++)
    {
      tmp[0]+=p4KSstar[j]*p4KSstar[j];
      tmp[1]+=p4KLstar[j]*p4KLstar[j];
    }
  tmp[0]=sqrt(tmp[0]);
  tmp[1]=sqrt(tmp[1]);
  e->myKSLEstar[0] = tmp[0]-EKstar;
  e->myKSLEstar[1] = tmp[1]-EKstar;
  if(abs(e->myKSLEstar[0]) > cutarray[4] || abs(e->myKSLEstar[1]) > cutarray[4] ){
    status++;
    cutresult[4]=0;
  }
  //===========================
  //   sqrt(Emiss**2+Pmiss**2) cut
  //===========================
  // calculate EPmiss
  tmp[0]=0.;
  tmp[1]=0.;
  for(j=0;j<4;j++)
    {
      tmp[0]+=e->p4KSmiss[j]*e->p4KSmiss[j];
      tmp[1]+=e->p4KLmiss[j]*e->p4KLmiss[j];
    }
  e->myKSLEPMiss[0]=sqrt(tmp[0]);
  e->myKSLEPMiss[1]=sqrt(tmp[1]);
  if(e->myKSLEPMiss[0]>cutarray[1] || e->myKSLEPMiss[1]>cutarray[1]){
      status++;
      cutresult[1]=0;
    }
  //===========================
  //   missing mass cut
  //===========================
  // calculate EPmissmass
  e->myKSLMissMass[0]=e->p4KSmiss.Mag2();
  e->myKSLMissMass[1]=e->p4KLmiss.Mag2();
    if( e->myKSLMissMass[0]<cutarray[2] || e->myKSLMissMass[0]>cutarray[3] || e->myKSLMissMass[1]<cutarray[2] || e->myKSLMissMass[1]>cutarray[3] ){
      status++;
      cutresult[2]=0;
      cutresult[3]=0;
    }
  //===========================
  //   kinematic fit cut (chisquare & pulls)
  //===========================  
    if(FChi2>cutarray[5] || abs(FOut[0][0]-FIn[0][0])/FOut[0][1]>cutarray[6] || abs(FOut[1][0]-FIn[1][0])/FOut[1][1]>cutarray[6]){
      status++;
      cutresult[5]=0;
    }
  //===========================
  //   pions angle cut
  //===========================   
    if(e->p3KS[0]*e->p3KS[1]/(e->p3KS[0].Mag()*e->p3KS[1].Mag()) < cutarray[7] || e->p3KL[0]*e->p3KL[1]/(e->p3KL[0].Mag()*e->p3KL[1].Mag()) < cutarray[7]){
      status++;
      cutresult[7]=0;
    }

    cout << "GetCutResult WARNING: NCUT = " << NCUT << ". cutresult[>8] set to 1 (default). "<< endl;
    /* 
  //===========================
  //   ordered kaon z-momentum cut
  //===========================   
    if( (e->p4KS[2]>e->p4KL[2] && e->p4KS[2]<cutarray[8]) || (e->p4KS[2]<e->p4KL[2] && e->p4KL[2]<cutarray[8]) ){
      status++;
      cutresult[8]=0;
    }
    */
  
  // if passed all the cuts...
    if(status == -1){
      return 1;
    }else{
      return 0;
    }

}

int is_selected_lorentz(event e,double* cutarray)
{
  // (old) = from v 01-0
  // (new) = from lorentz selection (timedist_v03)
  //===========================
  //   NPar cut (old)
  //===========================
  int status=-1;
  if(NPar==0)return status; // quando serve questo taglio?
  //===========================
  //   tagMask cut (old)
  //===========================
  status=-2;
  if( ((tagMask>>5)&1)==0 && ((tagMask>>6)&1)==0 )return status; // if bit5 = 0 AND bit6=0 return (vtx requirement) GF??
  //===========================
  //   m(pipi) cut (new)
  //===========================
  status=-3;
  if(abs(KSLMass[0])>cutarray[0] || abs(KSLMass[1])>cutarray[0])return status;
  //===========================
  //   Pstar cut (new)
  //===========================
  //
  status=-4;
  if(abs(KSLEstar[0])>cutarray[4] || abs(KSLEstar[1])>cutarray[4])return status;
  //===========================
  //   sqrt(Emiss**2+Pmiss**2) cut (new)
  //===========================
  status=-5;
  if(KSLEPMiss[0]>cutarray[1] || KSLEPMiss[1]>cutarray[1])return status;
  //===========================
  //   missing mass cut (new)
  //===========================
  status=-6;
  if((KSLMissMass[0]<cutarray[2] || KSLMissMass[0]>cutarray[3]) || (KSLMissMass[1]<cutarray[2] || KSLMissMass[1]>cutarray[3]))return status;
  //===========================
  //   kinematic fit cut (chisquare & pulls)
  //===========================  
  status=-7;
  if(FChi2>cutarray[5] || abs(FOut[0][0]-FIn[0][0])/FOut[0][1]>cutarray[6] || abs(FOut[1][0]-FIn[1][0])/FOut[1][1]>cutarray[6])return status;
  //===========================
  //   pions angle cut
  //===========================   
  status=-8;
  if(e.p3KS[0]*e.p3KS[1]/(e.p3KS[0].Mag()*e.p3KS[1].Mag()) < cutarray[7] || e.p3KL[0]*e.p3KL[1]/(e.p3KL[0].Mag()*e.p3KL[1].Mag()) < cutarray[7])return status;
  //===========================
  //   ordered kaon z-momentum cut
  //===========================   
  status=-9;
  if( (e.p4KS[2]>e.p4KL[2] && e.p4KS[2]<cutarray[8]) || (e.p4KS[2]<e.p4KL[2] && e.p4KL[2]<cutarray[8]) )return status;
  
  // if passed all the cuts...
  return 1;
}

void SetBranchAddresses(TTree* mytree)
{ // defines link btw leaves and global variables
  mytree->SetBranchAddress("nRun",&NRun);
  mytree->SetBranchAddress("nEv",&NEv);
  mytree->SetBranchAddress("NPar",&NPar); //
  mytree->SetBranchAddress("tagMask",&tagMask);
  mytree->SetBranchAddress("pKS",pKS);
  mytree->SetBranchAddress("pKLv",pKL);
  mytree->SetBranchAddress("KSLMass",KSLMass); //
  mytree->SetBranchAddress("KSLEPMiss",KSLEPMiss);
  mytree->SetBranchAddress("KSLMIssMass",KSLMissMass);
  mytree->SetBranchAddress("KSLEstar",KSLEstar);
  mytree->SetBranchAddress("FIn",FIn); //
  mytree->SetBranchAddress("FOut",FOut);
  mytree->SetBranchAddress("FChi2",&FChi2);
  mytree->SetBranchAddress("pKStot",pKStot);
  mytree->SetBranchAddress("pKLtag",pKLtag);
  mytree->SetBranchAddress("pKSMC",pKSMC); //
  mytree->SetBranchAddress("pKLMC",pKLMC);
  mytree->SetBranchAddress("pKStotMC",pKStotMC);
  mytree->SetBranchAddress("pKLtotMC",pKLtotMC);
  mytree->SetBranchAddress("KSwordMC",&KSwordMC);
  mytree->SetBranchAddress("KLwordMC",&KLwordMC);
  mytree->SetBranchAddress("xPhiMC",xPhiMC);
  mytree->SetBranchAddress("xKSMC",xKSMC);
  mytree->SetBranchAddress("xKLMC",xKLMC);
}

void SetBranchAddresses_Patch(TTree* mytree)
{ // defines link btw leaves and global variables
  mytree->SetBranchAddress("nrun",&NRun);
  mytree->SetBranchAddress("nev",&NEv);
  mytree->SetBranchAddress("Npar",&NPar); 
  mytree->SetBranchAddress("tagmask",&tagMask);
  mytree->SetBranchAddress("pks",pKS);
  mytree->SetBranchAddress("pklv",pKL);
  mytree->SetBranchAddress("Kslmass",KSLMass); 
  mytree->SetBranchAddress("Kslepmiss",KSLEPMiss);
  mytree->SetBranchAddress("Kslmissmass",KSLMissMass);
  mytree->SetBranchAddress("Kslestar",KSLEstar);
  mytree->SetBranchAddress("Fin",FIn); 
  mytree->SetBranchAddress("Fout",FOut);
  mytree->SetBranchAddress("Fchi2",&FChi2);
  mytree->SetBranchAddress("pkstot",pKStot);
  mytree->SetBranchAddress("pkltag",pKLtag);
  mytree->SetBranchAddress("pksmc",pKSMC);
  mytree->SetBranchAddress("pklmc",pKLMC);
  mytree->SetBranchAddress("pkstotmc",pKStotMC);
  mytree->SetBranchAddress("pkltotmc",pKLtotMC);
  mytree->SetBranchAddress("Kswordmc",&KSwordMC);
  mytree->SetBranchAddress("Klwordmc",&KLwordMC);
  mytree->SetBranchAddress("xphimc",xPhiMC);
  mytree->SetBranchAddress("xksmc",xKSMC);
  mytree->SetBranchAddress("xklmc",xKLMC);
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
int vector2file(const char* filename, int nbin,int* input)
{
  ofstream myfile;
  int i;
  myfile.open(filename);
  for(i=0;i<nbin;i++)myfile << input[i] << endl;
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
//=============================================================
//     GFloop methods
//-------------------------------------------------------------
GFloop::GFloop(char mode){
  nfile=0;
  nev=0;
  if(mode == 'D'){
    DATmode=1;
    MCmode=0;
  }else if(mode == 'M'){
    DATmode=0;
    MCmode=1;
  }else{
    cout << "GFloop: invalid call. mode = 'M'/'D' for MC / DATA" << endl;
  }
}

int GFloop::LoadNext(int verbose)
{
  if(nfile == 419 && nev==eventTOT)
    {
      myfile->Close();
      delete myfile;
      return 0; // END of files
    }
  if(nev==eventTOT || nfile == 0){ // end of file reached
    nfile++;
    if(nfile==129){nfile=133;} //BUG skip corrupted files
    nev=0;
    // open rootfile, get into directory, get TTree
    if(MCmode==1){// /home/guido or /media/guido/Elements/
      /*
      if(nfile<=100){
      sprintf(filename,"/media/guido/A8147BA8147B77E0/montecarlo/mk0_lorentz_l4rkmc_%d.root",nfile);
      }else{
      sprintf(filename,"/media/guido/Elements/ntu/mc/mk0_lorentz_l4rkmc_%d.root",nfile);
      }
      */
      sprintf(filename,"/Volumes/Elements//ntu/mc/mk0_lorentz_l4rkmc_%d.root",nfile);
    }else if(DATmode==1){
      sprintf(filename,"/Volumes/Elements//ntu/data/dk0_lorentz_l4dsys_%d.root",nfile);  
    }
    if(nfile!=1){
      myfile->Close();
      delete myfile;
    }
    myfile = new TFile(filename,"READ");
    mytree = (TTree*)myfile->Get("KCP_TAG/h1");   
    eventTOT= mytree->GetEntries();
    if(eventTOT == 0){
      cout << "GFloop: FATAL ERROR!" << endl;
      cout << filename << endl;
      cout << "eventTOT = 0" << endl;
      cout << "loop terminates here." << endl;
      return 0;
    }
    if(verbose>=2)cout << "eventTOT set to " << eventTOT << endl;
    if(verbose>=1)std::cout << "Processing " << filename << endl;
    // define link btw leaves and variables (gflib)
      if(nfile==101 && DATmode==1){ // anomalous tree structure patch (upper/lower case)
	SetBranchAddresses_Patch(mytree);
      }else{
	SetBranchAddresses(mytree);
      }   
  }
  if(verbose>=2)cout << "GFloop: nev / eventTOT = " << nev << " / " << eventTOT << endl;
  mytree->GetEntry(nev);// fetch variables 
  nev++;
  return 1;
}

//---------------------------------------------------------------------------
class fitresult
{
public:
  int id;
  int npar; // new (version > 08)
    Double_t value; // compatibility for < v08 val[0]
    Double_t error; // compatibility for < v08 err[0]
  Double_t val[10]; // up to 10 parameters are supported
  Double_t err[10];
  Double_t chisquare;
  Double_t dat[12];
  Double_t fit[12];
  Double_t errDAT[12]; // loadERR
  Double_t errMC[12];  // loadERR
  void LoadID(const char* inputfile,int setID,int verbosity); // old method compatibility (up to v07)
  void LoadID(const char* inputfile,int setID,int setNPAR,int verbosity);
  void LoadERR(const char* path);
  void Dump();

};
void fitresult::LoadID(const char* inputfile,int setID,int verbosity){
    if(verbosity == 2)cout << "[fitresult] Running compatibility mode gflib version < v08" << endl;
    this->LoadID(inputfile,setID,1,verbosity);
    value = val[0]; // compatibility
    error = err[0]; // compatibility
}

void fitresult::LoadID(const char* inputfile,int setID,int setNPAR,int verbosity){
  fstream input;
  input.open(inputfile,ios::in);
  if(!input.is_open())cerr << "[fitresult::LoadID] ERROR: could not open " << inputfile << endl;
  int i,j,k;
  int exitFLAG=0;
  if(setNPAR>0)npar=setNPAR;
  i=0;
  while( i < (NCUT*4+1) && exitFLAG != 1){ // loop over IDs = "folders + default folder" (i.e. lines in the inputfiles)
	  if(input.eof())
	  {
		  cerr << "[fitresult::LoadID] ERROR: End of file reached! " << endl;
		  cerr << "[fitresult::LoadID] val[k] and err[k] set to 0." << endl;
		  for(k=0;k<npar;k++)
		  {
			  val[k] = 0.;
			  err[k] = 0.;
		  }
		  break;
	  }
	  input >> id;
	  for(k=0;k<npar;k++)input >> val[k] >> err[k]; // input as many pairs of <val> <err> as set by setNPAR
	  input >> chisquare;
	  if(verbosity == 2)
	  {
		  cout << "id/val[0]/err[0]/.../chisquare " << id ;
		  for(k=0;k<npar;k++)cout << " / " << val[k] << " / " << err[k];
		  cout << " / " << chisquare << endl;
	  }
	  for(j=0;j<12;j++)input >> dat[j] >> fit[j];
	  if(id == setID)exitFLAG=1; // exit loop
  }
}
void fitresult::LoadERR(const char* path){
  char filename[100];
  sprintf(filename,"%serrDAT.dat",path);
  file2vector(filename,12,errDAT);
  sprintf(filename,"%serrMC.dat",path);
  file2vector(filename,12,errMC);
}
void fitresult::Dump(){
  cout << "fitresult::Dump()" << endl;
  cout << "Welcome to the DEBUG of class fitresult. Yes, this is verbose mode!" << endl;
  cout << "id = " << id << endl;
  cout << "npar = " << npar << endl;
    int k;
    for(k=0;k<npar;k++)
    {
        cout << "val[" << k << "] = " << val[k] << endl;
	cout << "err[" << k << "] = " << err[k] << endl;
    }
  cout << "chisquare = " << chisquare << endl;
  cout << "----------------------------------" << endl;
  int i;
  for(i=0;i<12;i++){
    cout << "dat / fit = " << dat[i] << " / " << fit[i] << endl;
  }
  cout << "----------------------------------" << endl;
  for(i=0;i<12;i++){
    cout << "errDAT / errMC = " << errDAT[i] << " / " << errMC[i] << endl;
  }
  cout << "==================================" << endl;
}





#endif
