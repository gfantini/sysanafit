void KolmogorovProva() {

  TRandom3 r;
  r.SetSeed(0);
  double meanA  =  0.0;
  double widthA =  1.0;
  double meanB  =  0.0;
  double widthB =  1.0;
  const int Nexp = 1;
  const int Nevents = 1000;
  TH1F* h  = new TH1F("Hist_ProbKS" , "Hist_ProbKS",   30,  0.0, 1.0);
  TH1F* h1  = new TH1F("gau1" , "Gauss1",   30,  -10, 10);
   TH1F* h2  = new TH1F("gau2" , "Gauss2",   30,  -10, 10);
   vector<double> x_gaussA;
   for (int i=0; i < Nevents; i++) {
    double xA = r.Gaus(meanA,widthA);
     x_gaussA.push_back(xA);
  }
  for (int j = 0; j < 5000 ; j++){
    vector<double> x_gaussB;
  for (int i=0; i < Nevents; i++) {
     double xB = r.Gaus(meanB,widthB);
      x_gaussB.push_back(xB);
      if(j == 0){
      h1 -> Fill(x_gaussA[i]);
      h2 -> Fill(x_gaussB[i]);
      }
  }
  double x_gaussA_array[Nevents];
  double x_gaussB_array[Nevents];
  sort(x_gaussA.begin(), x_gaussA.end());
  sort(x_gaussB.begin(), x_gaussB.end());
  copy(x_gaussA.begin(), x_gaussA.end(), x_gaussA_array);
  copy(x_gaussB.begin(), x_gaussB.end(), x_gaussB_array);
  double pKS = TMath::KolmogorovTest(Nevents, x_gaussA_array, Nevents, x_gaussB_array, "");
  h -> Fill(pKS);
  }
  h -> Draw();
  //  h1 -> Draw();
  //h2 -> Draw("ESame");
}
