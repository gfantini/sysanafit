using namespace RooFit ;
using namespace RooStats;
void PvalueTest () {
  double pKs , pKl;
  TFile *KS = new TFile("pValue.root" , "RECREATE");
  TTree *ptree = new TTree("t" , "pvalues");
  ptree -> Branch("pks" , &pKs);
  ptree -> Branch("pkl" , &pKl);
  for(int j = 0 ; j < 200 ; j++){
    gROOT->ProcessLine(".x Bprova2D_bkgTunnel.cxx");
    //TFile *rootfile = new TFile("../src/Histograms_cutvar00_p000.root");
  TFile *rootfile = new TFile("Prova.root");
  //TTree *tree = (TTree*) rootfile -> Get("TreeIMass2D");
  TTree *tree = (TTree*) rootfile -> Get("Prova");
  RooRealVar x("InvMassKS" , "InvMassKS" , -10 , 10);
  RooRealVar y("InvMassKL" , "InvMassKL" , -10 , 10);
  RooCategory c("c","c") ;
  c.defineType("accept",1) ;
  c.defineType("reject",0) ;
  RooDataSet data("data","data", tree ,RooArgSet(x,y,c));
  RooRealVar sigmax("sigmax","sigmax",1,0.,10.);
  RooRealVar meanx1("meanx1","meanx1",0 , -10. , 10.) ;
  RooRealVar sigmax1("sigmax1","sigmax1",1 , 0 , 10.);
  RooGaussian gaussx1("gaussx1","gaussx1",x,meanx1,sigmax1);
  RooRealVar meany1("meany1","meany1",0 , -10. , 10.) ;
  RooRealVar sigmay1("sigmay1","sigmay1",1 , 0. , 10.) ;
  RooGaussian gaussy1("gaussy1","gaussy1",y,meany1,sigmay1);
  RooProdPdf gaussxy1("gaussxy1","gaussxy1",RooArgSet(gaussx1,gaussy1));
  RooGenericPdf gauss("g" , "g" , "1/(sigmax*sqrt(2*pi))*exp(-(InvMassKS*cos(-pi/4) - InvMassKL*sin(-pi/4))^2/(2*sigmax^2))" , RooArgList(x , y , sigmax));
  RooRealVar nbkg("nbkg","background fraction", 1000 , 0 , 7000) ;
  RooExtendPdf ebkg("ebkg","ebkg",gauss,nbkg);
  RooRealVar nsig("nsig","signal fraction",2000 , 0 , 7000) ;
  RooExtendPdf esig("esig","esig",gaussxy1,nsig);
  RooAddPdf model("model" , "model" , RooArgList(ebkg , esig));
  RooDataSet *Dat = model.generate(RooArgSet(x,y) , 3000);
  RooMsgService::instance().setSilentMode(true);
  model.fitTo(data , Extended(kTRUE));
  double ks , kl;
  double ks1 , kl1;
  tree->SetBranchAddress("InvMassKS", &ks);
  tree->SetBranchAddress("InvMassKL", &kl);
  double KS[tree->GetEntries()], KL[tree -> GetEntries()];
  vector<double>Ks;
  vector<double>Ks2;
  vector<double>Kl;
  vector<double>Kl2;
  TFile *rfile = new TFile("P.root" , "RECREATE");
  TTree *tree2 = (TTree*)GetAsTTree("boh" , "boh" , *Dat);
  tree2 -> Write();
  tree2->SetBranchAddress("InvMassKS", &ks1);
  tree2->SetBranchAddress("InvMassKL", &kl1);
  for (int i = 0,  N = tree->GetEntries(); i < N; i++) {
    tree->GetEntry(i);
    tree2 -> GetEntry(i);
    Ks.push_back(ks);
    Ks2.push_back(ks1);
    Kl.push_back(kl);
    Kl2.push_back(kl1);
  }
  std::sort(Ks.begin(), Ks.end());
  std::sort(Ks2.begin(), Ks2.end());
  std::sort(Kl.begin(), Kl.end());
  std::sort(Kl2.begin(), Kl2.end());
 
  pKs =  TMath::KolmogorovTest(tree->GetEntries(), &Ks[0], tree2->GetEntries(), &Ks2[0] , "");
  pKl = TMath::KolmogorovTest(tree->GetEntries(), &Kl[0], tree2->GetEntries(), &Kl2[0] , "");
  ptree-> Fill();
  rfile->Close();
  rootfile -> Close();
  delete rfile;
  delete rootfile;
  cout << "Iteration: " << j << endl;
  }
  KS-> Write();
  KS -> Close();
}
