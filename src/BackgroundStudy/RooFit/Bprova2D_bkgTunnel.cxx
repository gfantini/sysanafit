void Bprova2D_bkgTunnel() {
  TFile *rootfile= new TFile("Prova.root" , "RECREATE");
  TTree *tree = new TTree("Prova" , "Prova");
  double mks;
  double mkl;
  tree -> Branch("InvMassKS" , &mks);
  tree -> Branch("InvMassKL" , &mkl);
  const int nbkg = 10000;
  const int nsig = 20000;
  // generate bkg events with hit-or-miss approach
  // gauss distributed in mks + mkl variable
  // flat distributed in mks - mkl variable
  int nBkgAcc = 0;
  double mkSum = 0.;
  double mkDiff = 0.;
  const double mksSymmetricCut = 10.; // MeV
  const double mklSymmetricCut = 10.; // MeV
  while( nBkgAcc < nbkg )
    {
      mkSum  = gRandom->Rndm();
      mkDiff = (-0.5+gRandom->Rndm())*40;
      mks    = (mkSum + mkDiff)/2.;
      mkl    = (mkSum - mkDiff)/2.;
      if( abs(mks) > mksSymmetricCut ||
	  abs(mkl) > mklSymmetricCut )
	continue;
      // otherwise ...
      tree->Fill();
      nBkgAcc++;
    }

  for(int i = 0; i < nsig ; i++){
    mks = gRandom->Gaus(0,2);  
    mkl = gRandom->Gaus(0,2);
    tree -> Fill();
  }

  tree ->Print();
  tree -> Write();
  rootfile -> Close();

}
