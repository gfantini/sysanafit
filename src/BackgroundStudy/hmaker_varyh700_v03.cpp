/*
  This programme makes histograms
  h1000 data histo
  h100 h200 efficiency histo
  h700 smearing matrix histo
  "cfg" TVectorD(8)
 "4pi" TVectorD(1)
 "kregen" TVectorD(1)
  "kres" TVectorD(1)
 
 
  _varyh700_v03 Guido Fantini 08-08-2017
  -- use updated gflib_v08. This does NOT include the pz cut (included in the gflib_v04!)
  -- remove configfile requirement (useless -> store to root file TVectorD(8) as "cfg" the cut configuration )
 
  _varyh700_v02 Guido Fantini 08-08-2017
  -- This modifies the resolution enlarging(K>0) / shrinking(K<0) by a factor K 

  v02 Guido Fantini 01-06-2016
  -- Including gflib_v04 I include smearing of MC momenta

  v01 Guido Fantini 19-05-2016
  -- configfile needs to have exactly NCUT rows filled with "nsigma[i]" (default 0.)

*/

#include "../inc/phys_const.inc"
#include "../inc/gflib_v08.h"
const int NselectionCuts = 8;

void hmaker_varyh700_v03(const char* rootfile_output,const int cut_to_vary=0,const int nsigma=0,const int nSigma4pi = 0,const double Kregeneration=0,const double Kresolution = 0.)
{
	// setting cuts
	double cutarray[NselectionCuts];
	double ArrayNsigma[NselectionCuts];
	if(cut_to_vary < 0 || cut_to_vary > NCUT){
		cerr << "FATAL: cut_to_vary = " << cut_to_vary << " where NCUT (max) = " << NCUT << endl;
		exit(0);
	}
	for(int i=0;i<NselectionCuts;i++)ArrayNsigma[i]=0.;
	ArrayNsigma[cut_to_vary] = (double)nsigma;
	
	// setup cut array [NselectionCuts]
	cutarray[0] = TRKMASSCUT + ArrayNsigma[0]*TRKMASSRES;
	cutarray[1] = EPMISSCUT + ArrayNsigma[1]*EPMISSRES;
	cutarray[2] = MISSMASSCUT_LOW + ArrayNsigma[2]*MISSMASSRES_LOW;
	cutarray[3] = MISSMASSCUT_HIGH + ArrayNsigma[3]*MISSMASSRES_HIGH;
	cutarray[4] = PSTARCUT + ArrayNsigma[4]*PSTARRES;
	cutarray[5] = GLOBALFITCHISQUARE + ArrayNsigma[5]*GLOBALFITCHISQUARERES;
	cutarray[6] = GLOBALFITPULL + ArrayNsigma[6]*GLOBALFITPULLRES;
	cutarray[7] = COSPIONCUT + ArrayNsigma[7]*COSPIONRES;
	//	cutarray[8] = PZCUT + nsigma[8]*PZRES; // DEPRECATED
	TVectorD Tcutarray(NselectionCuts);
	for(int i=0;i<NselectionCuts;i++)Tcutarray[i]=cutarray[i];
	// 4pi bkg variation
	TVectorD T4pi(1);
	T4pi[0] = nSigma4pi;
	// regen variation
	TVectorD Tregen(1);
	Tregen[0] = Kregeneration;
	// resolution variation
	TVectorD Tres(1);
	Tres[0] = Kresolution;
	
  // definition of root variables
  TFile* myfile;
  TTree* mytree;
  char dat_filename[100];
  char mc_filename[100];  
  char* filename;
  // histogram <---- note that they are FLOAT
  TH1F* DeltaT = new TH1F("h1000","Time distribution (data);Delta t / tau_s;Events;",70,0,70);
  TH1F* InvMass = new TH1F("hInvMass","Invariant Mass;m_{#pi#pi} - m_{K^{0}} [MeV]; Events;",60,-6,6);
  //TH1F* InvMassL = new TH1F("hInvMassKL","InvariantMass KL ;MeV; Events;",100,-15,15);
  TH1F* DeltaTall = new TH1F("h100","Time distribution (MC-signal);Delta tMC/tau_s;Events;",280,0,70);
  TH1F* DeltaTsel = new TH1F("h200","Time distribution (MC-signal-sel);Delta tMC/tau_s;Events;",280,0,70);
  TH2F* SmearingMatrix = new TH2F("h700","Smearing Matrix;DeltaT_MC;DeltaT_reco;",280,0,70,280,0,70);


  event ev;
  Int_t EventTOT;        //number of events in rootfile
  // <-- possibly include into gflib (next version)

  
 
  // definition of auxiliary variables
  int i,j,k;
  double tmp[6];
  k=1; // <<------------------------------ BEGINNING FROM FILE #...
  tmp[5]=0.;

  // loop on files
  while(k<=419)
    {
      // open rootfile, get into directory, get TTree
      filename = dat_filename; // <<------ MODIFY HERE TO SWITCH FROM DATA TO MC, NOT ANYWHERE ELSE
      //sprintf(dat_filename,"/home/guido/ntu/data/dk0_lorentz_l4dsys_%d.root",k);
      //sprintf(mc_filename,"/home/guido/ntu/mc/mk0_lorentz_l4rkmc_%d.root",k);
      sprintf(dat_filename,"../ntp/data/dk0_lorentz_l4dsys_%d.root",k);
      sprintf(mc_filename,"../ntp/mc/mk0_lorentz_l4rkmc_%d.root",k);   
      
      myfile = new TFile(filename,"READ");
      mytree = (TTree*)myfile->Get("KCP_TAG/h1"); 
      EventTOT= mytree->GetEntries();

      cout << "Processing " << filename << endl;
      cout << "Number of events written to file: " << tmp[5] << endl;
      
      // define link btw leaves and variables (gflib)
      if(k==101 && filename == dat_filename){ // anomalous tree structure patch (upper/lower case)
	  SetBranchAddresses_Patch(mytree);
	}else{
	  SetBranchAddresses(mytree);
	}   
      // ----------------------------------------------------------------------------
      
      // loop on events
      for(i=0;i<EventTOT;i++)
	{ 
	  tmp[5]+=1.; //update event counter
	  // fetch variables
	  mytree->GetEntry(i);
	  ev.fill_data();
	  tmp[0]=is_selected(&ev,cutarray); // apply cuts
	  if(tmp[0]==1){ // if selected
	    DeltaT->Fill(ev.DeltaT);
	    if(ev.DeltaT <= 1){
	      InvMass -> Fill(ev.myKSLMass[0]);
	      InvMass -> Fill(ev.myKSLMass[1]);
	    }
	  }
	  
	}     
      k++;
      myfile->Close();
      delete myfile;
    }
	for(auto ijk=0; ijk<DeltaT->GetNbinsX(); ijk++)cout << "DEBUG " << ijk << " <-> " << DeltaT->GetBinContent(ijk) << endl; // DEBUG
  for(i=0;i<6;i++)tmp[i]=0.;
  k=1;/*
  while(k<=419) // loop on files again
    {
      // open rootfile, get into directory, get TTree
      filename = mc_filename; // <<------ MODIFY HERE TO SWITCH FROM DATA TO MC, NOT ANYWHERE ELSE
      //      sprintf(dat_filename,"/home/guido/ntu/data/dk0_lorentz_l4dsys_%d.root",k);
      //      sprintf(mc_filename,"/home/guido/ntu/mc/mk0_lorentz_l4rkmc_%d.root",k);
      sprintf(dat_filename,"../ntp/data/dk0_lorentz_l4dsys_%d.root",k);
      sprintf(mc_filename,"../ntp/mc/mk0_lorentz_l4rkmc_%d.root",k);

      myfile = new TFile(filename,"READ");
      mytree = (TTree*)myfile->Get("KCP_TAG/h1"); 
      EventTOT= mytree->GetEntries();

      cout << "Processing " << filename << endl;
      cout << "Number of events written to file: " << tmp[5] << endl;

      // define link btw leaves and variables (gflib)      
      if(k==101 && filename == dat_filename){ // anomalous tree structure patch (upper/lower case)
	  SetBranchAddresses_Patch(mytree);
	}else{
	  SetBranchAddresses(mytree);
      }
      //---------------------------------------------------------------------------------

      // loop on events
      for(i=0;i<EventTOT;i++)
	{ 
	  tmp[5]+=1.; //update event counter
	  // fetch variables
	  mytree->GetEntry(i);
	  smear_mc(FSMEARMC);// smear momenta before using them for calculations
	  ev.fill_data();
	  ev.fill_mc();
	  if(is_signal())
	    { // fill mc-all
	      DeltaTall->Fill(ev.DeltaTMC);
	      tmp[0]=is_selected(&ev,cutarray); // apply selection
	      if(tmp[0]==1)//if selected...
		{
		  DeltaTsel->Fill(ev.DeltaTMC);
		  SmearingMatrix->Fill(ev.DeltaTMC,ev.DeltaT+Kresolution*(ev.DeltaTMC-ev.DeltaT));
		  // DeltaTMC == dt_trueMC   DeltaT == dt_recoMC
		}
	    }  
	}
      
      k++;
      myfile->Close();
      delete myfile;
    }
      */
  /*
  // draw distributions
  TCanvas* c0=new TCanvas();
  gStyle->SetHistMinimumZero();
  DeltaT->Draw();
  TCanvas* c1=new TCanvas();
  gStyle->SetHistMinimumZero();
  DeltaTall->Draw();
  TCanvas* c2=new TCanvas();
  DeltaTsel->Draw();
  TCanvas* c3=new TCanvas();
  SmearingMatrix->Draw();
  */

  // write histograms to file
  TFile* output = new TFile(rootfile_output,"RECREATE");
  cout << "Writing " << rootfile_output << endl;
  DeltaT->Write();
  DeltaTall->Write();
  DeltaTsel->Write();
  InvMass -> Write();
  // InvMassL -> Write();
  SmearingMatrix->Write();
  Tcutarray.Write("cfg");
	T4pi.Write("4pi");
	Tregen.Write("kregen");
  Tres.Write("kres");
  output->Write();
  output->Close();
  delete output;
}
