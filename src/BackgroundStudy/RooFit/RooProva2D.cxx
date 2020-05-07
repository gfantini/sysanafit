  using namespace RooFit ;
void RooProva2D() {
  TFile *rootfile = new TFile("../src/Histograms_cutvar00_p000.root");
  //TFile *rootfile = new TFile("Prova.root");
  TTree *tree = (TTree*) rootfile -> Get("TreeIMass2D");
  RooRealVar x("InvMassKS" , "InvMassKS" , -5 , 5);
  RooRealVar y("InvMassKL" , "InvMassKL" , -5 , 5);
  RooCategory c("c","c") ;
  c.defineType("accept",1) ;
  c.defineType("reject",0) ;
  RooDataSet data("data","data", tree ,RooArgSet(x,y,c));
  RooRealVar meanx("meanx","meanx",2,-10,10) ;
  RooRealVar sigmax("sigmax","sigmax",1,0.,10.);
  RooGaussian gaussx("gaussx","gaussx",x,meanx,sigmax);
  RooRealVar meany("meany","meany",-2,-10,10) ;
  RooRealVar sigmay("sigmay","sigmay",5,0.,10.) ;
  RooGaussian gaussy("gaussy","gaussy",y,meany,sigmay);
  RooPolynomial bkg("pol0","pol0",x,RooArgList()) ;
  RooProdPdf gaussxy("gaussxy","gaussxy",RooArgSet(gaussx,gaussy)) ;
  RooRealVar nsig("nsig","signal fraction",500,0.,50000.) ;
  RooRealVar nbkg("nbkg","background fraction",40,0.,50000) ;
  RooExtendPdf ebkg("ebkg","ebkg",bkg,nbkg) ;
  RooExtendPdf esig("esig","esig",gaussxy,nsig);
  RooAddPdf model("model" , "model" , RooArgList(esig, ebkg));
  RooPlot* framex = x.frame() ;
  RooPlot* framey = y.frame();
  model.fitTo(data , Extended(kTRUE));
  data.plotOn(framey , Binning(14));
  model.plotOn(framey , LineColor(kRed));
  model.plotOn(framey , Components("pol0") , LineStyle(kDashed));
  framey -> Draw();
  data.plotOn(framex , Binning(14));
  model.plotOn(framex,  LineColor(kRed));
  model.plotOn(framex  , Components("pol0") , LineStyle(kDashed));
  framex -> Draw();

  

}
