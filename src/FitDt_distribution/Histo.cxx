void Histo() {
  auto c1 = new TCanvas("c1","c1",600,500);
  gStyle->SetOptStat(0);
  TFile *rootfile = new TFile("../../ntp/Histograms_cutvar00_p000.root");
  TH1F *hDataBkg = (TH1F*) rootfile -> Get("h1000");
  TH1F *hBkg = (TH1F*) rootfile -> Get("h5000");
  TH1F *hrege = (TH1F*) rootfile -> Get("hrege");
  hDataBkg -> GetYaxis() -> SetRangeUser(0 , 1100);
  hDataBkg -> GetXaxis() -> SetRangeUser(0 , 12);
  TH1F *prova = new TH1F();

  prova -> SetFillColor(12);
  hDataBkg -> SetMarkerStyle(8);
  hDataBkg -> Draw("E");
  //prova -> Draw("Same");
  TFile *rfile = new TFile("Aggiornamento18062020/fitSL.root");
  TDirectory* dir = (TDirectory*)rfile->Get("SL");
  dir -> cd();
  TH1F *hDataNoBkg = (TH1F*) dir -> Get("hData");
  TH1F *hFit = (TH1F*) dir -> Get("hFit");
  prova -> SetBinContent(1 , hFit -> GetBinContent(1) + hBkg -> GetBinContent(1));
  prova -> Draw("same");
  hDataNoBkg -> SetMarkerStyle(8);
  hDataNoBkg -> GetYaxis() -> SetRangeUser(0 , 1300);
  hDataNoBkg -> GetXaxis() -> SetRangeUser(0 , 12);
  hDataNoBkg -> GetYaxis() -> SetTitle("Events");
  hDataNoBkg -> GetXaxis() -> SetTitle("#Deltat/#tau_{s}");


  hFit -> SetFillColor(14);  
  hFit -> SetLineColor(kRed);
  hFit -> Draw("Same");
   double j = 1;
  
  TH1F *hrege1 = new TH1F("boh" , "mah" , 12 , 0 , 12);
  for (int i = 0 ; i < 12 ; i++){
    hrege1 -> SetBinContent(j , hrege -> GetBinContent(j));
    j = j+1;
  }
  TH1F *hFitNoRege =(TH1F*) hFit -> Clone();
  hFitNoRege -> Divide(hrege1);
  hFitNoRege -> SetLineColor(kGreen);
  hFitNoRege -> SetFillColor(kGray);
  hFitNoRege -> Draw("Same");
  hDataBkg -> Draw("ESAME");
  auto legend = new TLegend(0.1,0.75,0.295,0.9);
   legend->SetHeader("Time Distribution","C"); // option "C" allows to center the header
   legend->AddEntry(hDataBkg,"Data","p");
   legend->AddEntry(hFit ,"Regeneration","f");
   legend->AddEntry(prova,"e^{+}e^{-} #rightarrow 4#pi","f");
   legend->AddEntry(hFitNoRege,"Theory Histogram","f");
   legend->Draw();
   c1 -> Modified();
   c1 -> Update();
   c1 -> SaveAs("HistSL.png");
   return 0;
}
