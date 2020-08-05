using namespace std;
void SL(){
   TFile *rfile = new TFile("Fit/Aggiornamento18062020/fitSL.root");
  TDirectory* dir = (TDirectory*)rfile->Get("SL");
  dir -> cd();
  TH1D *gr1 = (TH1D*) dir -> Get("hFit");
  gr1 -> SetLineColor(kRed);
  gr1 -> SetMarkerColor(kRed);
  //gr1 -> Draw("E");
  //rfile -> Close();
  // delete rfile;
  TFile *file = new TFile("fitSL.root");
  TDirectory* dir2 = (TDirectory*)file->Get("SL");
  dir2 -> cd();
  TH1D *gr2 = (TH1D*) dir2 -> Get("hFit");
  gr2 -> SetMarkerStyle(4);
  TH1D *gr3 = new TH1D("Difference between histograms","h1" ,  12 , 0 , 12);
  for(int i = 0; i < 12 ; i++){
    cout <<"zeta < 0 " <<  gr1 -> GetBinContent(i+1) << endl;
    cout << "zeta = 0 "  <<  gr2 -> GetBinContent(i+1) <<endl;
    gr3 -> SetBinContent(i+1 , abs(gr1 -> GetBinContent(i+1) - gr2 -> GetBinContent(i+1)));

    }
  // gr2 -> Draw("ESame");
  gr2 -> Add(gr1,-1);
  gr3 -> SetTitle("Difference between histograms");
  gr3 -> Draw();

}
