/*
  This programme makes histograms
  h1000 data histo
  h100 h200 efficiency histo
  h700 smearing matrix histo

  v02 Guido Fantini 01-06-2016
  -- Including gflib_v04 I include smearing of MC momenta
  -- Included new function that makes histogram of data without kinfit

  v01 Guido Fantini 19-05-2016
  -- configfile needs to have exactly NCUT rows filled with "nsigma[i]" (default 0.)

*/

#include "../inc/phys_const.inc"
#include "../inc/gflib_v04.h"


void hmaker_v02(const char* rootfile_output,const char* configfile)
{
  // setting cuts
  double cutarray[NCUT];
  SetCutArray(configfile,cutarray);

  // definition of root variables
  TFile* myfile;
  TTree* mytree;
  char dat_filename[100];
  char mc_filename[100];  
  char* filename;
  // histogram <---- note that they are FLOAT
  TH1F* DeltaT = new TH1F("h1000","Time distribution (data);Delta t / tau_s;Events;",70,0,70);
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
      /* FIXME -- this was removed to test the sw
      sprintf(dat_filename,"/home/guido/ntu/data/dk0_lorentz_l4dsys_%d.root",k);
      sprintf(mc_filename,"/home/guido/ntu/mc/mk0_lorentz_l4rkmc_%d.root",k);
      */
      // test with data in the repo
      sprintf(dat_filename,"./ntp/data/dk0_lorentz_l4dsys_%d.root",k);
      if( k > 3 )
	break;

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
	  }	  
	}     
      k++;
      myfile->Close();
      delete myfile;
    }

  for(i=0;i<6;i++)tmp[i]=0.;
  k=1;
  while(k<=419) // loop on files again
    {
      // open rootfile, get into directory, get TTree
      filename = mc_filename; // <<------ MODIFY HERE TO SWITCH FROM DATA TO MC, NOT ANYWHERE ELSE
      sprintf(dat_filename,"/home/guido/ntu/data/dk0_lorentz_l4dsys_%d.root",k);
      sprintf(mc_filename,"/home/guido/ntu/mc/mk0_lorentz_l4rkmc_%d.root",k);
      // FIXME -- this is just to test it works
      sprintf(mc_filename,"./ntp/mc/mk0_lorentz_l4rkmc_%d.root",k);
      if( k > 3 )
	break;

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
		  SmearingMatrix->Fill(ev.DeltaTMC,ev.DeltaT);
		}
	    }  
	}
      
      k++;
      myfile->Close();
      delete myfile;
    }
  
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
  SmearingMatrix->Write();
  output->Write();
  output->Close();
  delete output;
}


// This version just disregards the result of the kin fit while making the h100 dat file
// WARNING: events are filled on the histo based also on kinfit chi2 and pull cuts even if the Dt is computed without kinfit
void hmaker_v02_NoKinFit(const char* rootfile_output,const char* configfile)
{
	// setting cuts
	double cutarray[NCUT];
	SetCutArray(configfile,cutarray);
	
	// definition of root variables
	TFile* myfile;
	TTree* mytree;
	char dat_filename[100];
	char mc_filename[100];
	char* filename;
	// histogram <---- note that they are FLOAT
	TH1F* DeltaT = new TH1F("h1000","Time distribution (data);Delta t / tau_s;Events;",70,0,70);
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
	double DtNoKinFit;
	
	// loop on files
	while(k<=419)
	{
		// open rootfile, get into directory, get TTree
		filename = dat_filename; // <<------ MODIFY HERE TO SWITCH FROM DATA TO MC, NOT ANYWHERE ELSE
		/*
		 sprintf(dat_filename,"/home/guido/ntu/data/dk0_lorentz_l4dsys_%d.root",k);
		 sprintf(mc_filename,"/home/guido/ntu/mc/mk0_lorentz_l4rkmc_%d.root",k);
		 */
	 sprintf(dat_filename,"/Volumes/Elements/ntu/data/dk0_lorentz_l4dsys_%d.root",k);
	 sprintf(mc_filename,"/Volumes/Elements/ntu/mc/mk0_lorentz_l4rkmc_%d.root",k);
		
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
		  DtNoKinFit = abs( FIn[1][0]/ev.BetaGammaL - FIn[0][0]/ev.BetaGammaS )/(C*TAUS); // <-----
		  DeltaT->Fill(DtNoKinFit); // <-----
	  }
		}
		k++;
		myfile->Close();
		delete myfile;
	}
	
	for(i=0;i<6;i++)tmp[i]=0.;
	k=1;
	while(k<=419) // loop on files again
	{
		// open rootfile, get into directory, get TTree
		filename = mc_filename; // <<------ MODIFY HERE TO SWITCH FROM DATA TO MC, NOT ANYWHERE ELSE
		/*
		 sprintf(dat_filename,"/home/guido/ntu/data/dk0_lorentz_l4dsys_%d.root",k);
		 sprintf(mc_filename,"/home/guido/ntu/mc/mk0_lorentz_l4rkmc_%d.root",k);
		 */
	 sprintf(dat_filename,"/Volumes/Elements/ntu/data/dk0_lorentz_l4dsys_%d.root",k);
	 sprintf(mc_filename,"/Volumes/Elements/ntu/mc/mk0_lorentz_l4rkmc_%d.root",k);
		
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
			  SmearingMatrix->Fill(ev.DeltaTMC,ev.DeltaT);
		  }
	  }
		}
		
		k++;
		myfile->Close();
		delete myfile;
	}
	
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
	SmearingMatrix->Write();
	output->Write();
	output->Close();
	delete output;
}
