/*
v04 (sviluppo in corso...)
    30-05-2016
    13-06-2016 update
    [ok] Function to smear MC reco p(trk) by f*gauss(m=0,sigma=1)
         Inserted random generator as global variable randgen2
    [ok] Array of cut names for plots available with FillCutNames
    [ok] Inserted file2vector and vector2file functions from root2dat_v03
    [ok] Inserted constant FSMEARMC = 2 x 10**-3
    [ok] Included DeltaT12 (data and MC)
    09-07-2016 update
    [ok] Inserted non-zero resolution for pull cut

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

#define NCUT 9// number of cuts
#define TRKMASSCUT 5 // 5 MeV cut on |m(pipi)-mK0|
#define EPMISSCUT 10 // 10 MeV cut on EPmiss
#define MISSMASSCUT_LOW -50 // MeV**2
#define MISSMASSCUT_HIGH 10 // MeV**2
#define PSTARCUT 10 // MeV
#define GLOBALFITCHISQUARE 30 //-logL > 30 discard
#define GLOBALFITPULL 10 // reject |in-out|/sigma(out) > 10
#define COSPIONCUT -0.975 // reject too open pions
#define PZCUT 2 // regect K (dal pz-con-segno più alto) if pz<2 MeV

#define TRKMASSRES 0.9
#define EPMISSRES 1.3
#define MISSMASSRES_LOW 4.0
#define MISSMASSRES_HIGH 4.0
#define PSTARRES 1.3
#define GLOBALFITCHISQUARERES 2.0
#define GLOBALFITPULLRES 1.0 // <------------------- this cut was not varied
#define COSPIONRES 0.001
#define PZRES 2.0

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
  Double_t DeltaT12; // t1 - t2 [ns] pz(K1) > pz(K2)
  
  // true variables end by MC
  TVector3 p3KSMC[2],p3KLMC[2];         // pions 3-momenta (MC)
  TLorentzVector p4KSMC,p4KLMC;         // kaons 4-momenta (MC)
  TLorentzVector p4KSmissMC,p4KLmissMC; //
  TLorentzVector p4PhiMC;               // phi 4-momentum (MC)
  
  Double_t BetaGammaSMC,BetaGammaLMC;   
  Float_t DistKSMC,DistKLMC;
  Float_t DeltaTMC;
  Double_t DeltaT12MC;
  
  double myKSLMass[2];
  double myKSLEPMiss[2];
  double myKSLEstar[2];
  double myKSLMissMass[2];
  
 public:
  void fill_data(); // fill reco variables
  void fill_mc();   // fill true variables
};

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
  cutarray[8] = PZCUT + nsigma[8]*PZRES;
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
  cutarray[8] = PZCUT;
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
  sprintf(string[8],"MeV");  
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
  resarray[8] = PZRES;
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
    sprintf(string,"cos[#theta(#pi+#pi-)]");
    break;
  case 8:
    sprintf(string,"pz");
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
  if(p4KL[2] > p4KS[2]){// if it is the other way around...
    DeltaT12*=-1.;
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
  if(NPar==0)return status; // quando serve questo taglio?
  //===========================
  //   tagMask cut
  //===========================
  status=-2;
  if( ((tagMask>>5)&1)==0 && ((tagMask>>6)&1)==0 )return status; // if bit5 = 0 AND bit6=0 return (vtx requirement) GF??
  //===========================
  //   m(pipi) cut
  //===========================
  status=-3;
  if(e->p4KS.Mag2()<0. || e->p4KL.Mag2()<0.)cout<<"WARNING: invariant kaon mass < 0"<< endl;
  e->myKSLMass[0]=e->p4KS.Mag()-MK0;
  e->myKSLMass[1]=e->p4KL.Mag()-MK0;
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
  status=-9;
  if( (e->p4KS[2]>e->p4KL[2] && e->p4KS[2]<cutarray[8]) || (e->p4KS[2]<e->p4KL[2] && e->p4KL[2]<cutarray[8]) )return status;
  
  // if passed all the cuts...
  return 1;
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
int file2vector(const char* filename,int nbin,double* output)
{
  ifstream myfile;
  int i;
  myfile.open(filename,ios::in);
  for(i=0;i<nbin;i++)myfile >> output[i];
  myfile.close();
  return 1;
}

#endif
