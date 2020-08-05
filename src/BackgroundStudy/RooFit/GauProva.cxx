using namespace RooFit ;
using namespace RooStats;
void GauProva () {
  TFile *rootfile = new TFile("../src/Histograms_cutvar00_p000.root");
  //TFile *rootfile = new TFile("Prova.root");
  TTree *tree = (TTree*) rootfile -> Get("TreeIMass2D");
  //TTree *tree = (TTree*) rootfile -> Get("Prova");
  RooRealVar x("InvMassKS" , "InvMassKS" , -10 , 10);
  RooRealVar y("InvMassKL" , "InvMassKL" , -10 , 10);
  RooCategory c("c","c") ;
  c.defineType("accept",1) ;
  c.defineType("reject",0) ;
  RooDataSet data("data","data", tree ,RooArgSet(x,y,c));
  RooRealVar meanx("meanx","meanx",0 , -10 , 10) ;
  RooRealVar sigmax("sigmax","sigmax",1,0.,10.);
  RooGaussian gaussx("gaussx","gaussx",x,meanx,sigmax);
  RooRealVar meany("meany","meany",0 , -10 , 10) ;
  RooRealVar sigmay("sigmay","sigmay" , 20 , 0 , 100) ;
  RooGaussian gaussy("gaussy","gaussy",y,meany,sigmay);
  RooProdPdf gaussxy("gaussxy","gaussxy",RooArgSet(gaussx,gaussy));
  RooRealVar meanx1("meanx1","meanx1",0 , -1 , 1) ;
  RooRealVar sigmax1("sigmax1","sigmax1",1 , 0 , 10.);
  RooGaussian gaussx1("gaussx1","gaussx1",x,meanx1,sigmax1);
  RooRealVar meany1("meany1","meany1",0 , -1 , 1) ;
  RooRealVar sigmay1("sigmay1","sigmay1",1 , 0. , 10.) ;
  RooGaussian gaussy1("gaussy1","gaussy1",y,meany1,sigmay1);
  RooProdPdf gaussxy1("gaussxy1","gaussxy1",RooArgSet(gaussx1,gaussy1));
  RooGenericPdf gauss("g" , "g" , "1/(sigmax*sqrt(2*pi))*exp(-(InvMassKS*cos(-pi/4) - InvMassKL*sin(-pi/4))^2/(2*sigmax^2))" , RooArgList(x , y , sigmax));
  RooPolynomial pol("pol0","pol0",y ,RooArgList());
  RooProdPdf gaussxy2("gaussxy2","gaussxy2",RooArgSet(gauss));
  RooRealVar nbkg("nbkg","background fraction", 100 , 0 , 50000) ;
  RooExtendPdf ebkg("ebkg","ebkg",gauss,nbkg);
  RooRealVar nsig("nsig","signal fraction",200 , 0 , 50000) ;
  RooExtendPdf esig("esig","esig",gaussxy1,nsig);
  RooAddPdf model("model" , "model" , RooArgList(ebkg , esig));
  RooDataSet *Dat = model.generate(RooArgSet(x,y) , 328);
RooFitResult *fr =  model.fitTo(data , Extended(kTRUE) ,Save());
  RooArgList pars(* ebkg.getParameters(RooArgSet(x,y)));
  RooArgList obs(* ebkg.getObservables(RooArgSet(x,y)));
  //  TCanvas * canvas = new TCanvas ("canvas","canvas") ;
   cout << "Size of the observables: " << obs.getSize() << endl;
  TF2 *f =(TF2*) ebkg.asTF(obs , pars , RooArgSet(x,y)) ;
  TF12 *f12 = new TF12("f12",f,0,"x");;
  TH2D* dh2 = (TH2D*)data.createHistogram("Ks vs Kl invariant mass" , x , Binning(15) , YVar(y, Binning(15)));
  dh2 -> SetMarkerStyle(4);
  dh2 -> SetLineColor(kRed);
  
  dh2 -> Draw("Lego");
  
  TH2D* mh2 =(TH2D*) model.createHistogram("Ks vs Kl invariant mass model" , x , Binning(15) , YVar(y, Binning(15)));
  mh2 -> Draw("surfSame");
  //canvas -> SaveAs("InvMass10MeV.png");
  double ks , kl;
  tree->SetBranchAddress("InvMassKS", &ks);
  tree->SetBranchAddress("InvMassKL", &kl);
  double KS[tree->GetEntries()], KL[tree -> GetEntries()];
  TH1F *h = new  TH1F("c" , "c" , 100 , -10 , 10);

  cout << "BKG pdf integral in the range x[-4,4] , y [-4,4]: " << f -> Integral(-4 , 4 , -4 , 4) << endl;
  cout << "Kolmogorov-Smirnov test: " <<  dh2 -> KolmogorovTest(mh2) << endl;
  /* TFile *rfile = new TFile("P.root" , "RECREATE");
  TTree *tree2 = (TTree*)GetAsTTree("boh" , "boh" , *Dat);
  tree2 -> Write();
  for (int i = 0,  N = tree->GetEntries(); i < N; i++) {
    tree->GetEntry(i);
    h -> Fill(ks);
  }
  rfile->Close();*/
  rootfile -> Close();
  //  delete rfile;
  delete rootfile;
  RooMCStudy mgr(model , RooArgSet(x,y) , Silence());
  mgr.addFitResult(*fr);
  mgr.generateAndFit(1000 , 328);
  TCanvas *c1 = new TCanvas("c" , "ParameterPlots" , 800, 800 );
  mgr.plotParam(nbkg) -> Draw();
  c1 -> SaveAs("AggiornamentoNew/nbkg.png");
  mgr.plotParam(nsig) -> Draw();
  c1 -> Modified();
  c1 -> Update();
  c1 -> SaveAs("AggiornamentoNew/nsig.png");
  mgr.plotPull(nbkg , FitGauss(true)) -> Draw();
  c1 -> Modified();
  c1 -> Update();
  c1 -> SaveAs("AggiornamentoNew/pull.png");
  mgr.plotNLL() -> Draw();
  c1 -> Modified();
  c1 -> Update();
  c1 -> SaveAs("AggiornamentoNew/LogLikelihood.png");
}
/*  RooDataHist hist ("boh", "mah" , RooArgSet(x,y) , dh2);
  RooChi2Var chi2 ("chi2" , "chi2" , model, hist);
  double chi = chi2.getVal();
  cout << "Value of the Chi2 for the fit: " << chi/(dh2->GetNbinsX()*dh2->GetNbinsY()-2.0) << endl;
  cout << model.createIntegral(RooArgSet(x,y) , NormSet(RooArgSet(x,y))) -> getVal() <<endl;*/
