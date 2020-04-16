

Double_t GausPol(Double_t *x , Double_t *par){

  double y = par[2]*TMath::Gaus(*x , par[0] , par[1]) + par[3];
  return y;

}

Double_t Pol(Double_t *x , Double_t *par){
  double y;
  if (*x < -5 || *x > 5){
    y = par[0];
    return y;
  } else {
    return 0;
  }

}

void BkgFit() {
  gStyle -> SetOptFit(11111111);

  TFile *rootfile = new TFile("Histograms_cutvar00_p000.root");
  TH1F *h = (TH1F*) rootfile-> Get("hInvMass");
  TF1 *f = new TF1("f" , GausPol , -5 , 5 , 4);
  TF1 *f1 = new TF1("f1" ,Pol , -5 , 5 , 1);
  f1 -> SetParName(0 , "Bkg");
  f1 -> SetParameter(0 , 1);
  f -> SetParName(0 , "Mean");
  f -> SetParName(1 , "Std Dev");
  f -> SetParName(2 , "Norm");
  f -> SetParName(3 , "Bkg");  
  //TF1 *f = new TF1("f" , "[0]*gaus(1) + pol0(4)" , -10 , 10);
   f -> SetParameters (-1 , 1 , 40 , 3);
   f -> SetNpx(300);
   h -> Fit(f, "R");
  h-> Draw();







}
