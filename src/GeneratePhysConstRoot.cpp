// 13.10.2017 FINALIZED
//            Appende un TTree che si chiama PhysicalConstante all' outputRootfile

#include "../inc/phys_const_v00.cc"

// will generate files in output_path/GeneratePhysConst_##.inc where ## = 1 .. Ngen
// Will generate a TTree of Ngen entries and append it to the rootfile outputRootfile
void GeneratePhysConstRoot(const int Ngen,const char* outputRootfile)
{
	// define the variables
	PhysConst Tau_S("Tau_S",0.89564*1E-10,0.00033*1E-10,"s","w/o CPT KS review PDG (2016) and 2017 update");
	PhysConst Tau_L("Tau_L",5.116*1E-8,0.021*1E-8,"s","KL review PDG (2016) and 2017 update");
	PhysConst DeltaM("DeltaM",0.5289*1E10,0.0010*1E10,"hbar s-1","not assuming CPT KL review PDG (2016) and 2017 update");
	// |eta+-| = A(KL -> pi+pi-) / A(KS -> pi+pi-) rev. KL p. 35/50 [PDG2016]
	PhysConst modeta("modeta",2.232e-3,0.011e-3," PDG (2016) and 2017 update");
	PhysConst phieta("phieta",0.7574728954,0.0087,"rad","(43.4° w/o CPT), Err 0.5° PDG (2016) and 2017 update");
	
	// open the rootfile of output, define the branches and the variables to store in there
	TFile* pOutput = new TFile(outputRootfile,"UPDATE");
	if(pOutput->IsZombie()){
		cerr << "FATAL! Cannot open " << outputRootfile << endl;
		exit(0);
	}
	TTree* pTree = new TTree("PhysicalConstants","Tree of random generated physical constants");
	Double_t d_Tau_S,d_Tau_L,d_DeltaM,d_modeta,d_phieta;
	pTree->Branch("TauS",&d_Tau_S,"TauS/D");
	pTree->Branch("TauL",&d_Tau_L,"TauL/D");
	pTree->Branch("DeltaM",&d_DeltaM,"DeltaM/D");
	pTree->Branch("modeta",&d_modeta,"modeta/D");
	pTree->Branch("phieta",&d_phieta,"phieta/D");
	
	// loop Ngen times to fill the tree with randomly generated values of the physical constants
	for(int i=0;i<Ngen; i++)
	{
		d_Tau_S = Tau_S.GetRandom();
		d_Tau_L = Tau_L.GetRandom();
		d_DeltaM = DeltaM.GetRandom();
		d_modeta = modeta.GetRandom();
		d_phieta = phieta.GetRandom();
		pTree->Fill();
	}
	
	// write tree, close file, goodbye
	pTree->Write("PhysicalConstants");
	pOutput->Close();
}
