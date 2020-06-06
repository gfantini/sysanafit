using namespace RooFit ;
void GauProva () {
  TFile *rootfile = new TFile("../src/Histograms_cutvar00_p000.root");
  //TFile *rootfile = new TFile("Prova.root");
 TTree *tree = (TTree*) rootfile -> Get("TreeIMass2D");
 //TTree *tree = (TTree*) rootfile -> Get("Prova");
 RooRealVar x("InvMassKS" , "InvMassKS" , -5 , 5);
 RooRealVar y("InvMassKL" , "InvMassKL" , -5 , 5);
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
 RooRealVar sigmax1("sigmax1","sigmax1",1 , 0 , 2);
 RooGaussian gaussx1("gaussx1","gaussx1",x,meanx1,sigmax1);
 RooRealVar meany1("meany1","meany1",0 , -1 , 1) ;
 RooRealVar sigmay1("sigmay1","sigmay1",1 , 0. , 2) ;
 RooGaussian gaussy1("gaussy1","gaussy1",y,meany1,sigmay1);
 RooProdPdf gaussxy1("gaussxy1","gaussxy1",RooArgSet(gaussx1,gaussy1));
 RooGenericPdf gauss("g" , "g" , "1/(sigmax*sqrt(2*pi))*exp(-(InvMassKS*cos(-pi/4) - InvMassKL*sin(-pi/4))^2/(2*sigmax^2))" , RooArgList(x , y , sigmax));
 RooPolynomial pol("pol0","pol0",y ,RooArgList()) ;
 RooProdPdf gaussxy2("gaussxy2","gaussxy2",RooArgSet(gauss));
 RooRealVar nbkg("nbkg","background fraction", 1000 , 0 , 50000) ;
 RooExtendPdf ebkg("ebkg","ebkg",gauss,nbkg);
 RooRealVar nsig("nsig","signal fraction",2000 , 0 , 50000) ;
 RooExtendPdf esig("esig","esig",gaussxy1,nsig);
 RooAddPdf model("model" , "model" , RooArgList(ebkg , esig));
 model.fitTo(data , Extended(kTRUE));
 TH2D* mh2 =(TH2D*) model.createHistogram("Ks vs Kl invariant mass model" , x , Binning(50) , YVar(y, Binning(50)));
 mh2 -> Draw("Surf");

 TH2D* dh2 = (TH2D*)data.createHistogram("Ks vs Kl invariant mass" , x , Binning(10) , YVar(y, Binning(10)));
 dh2 -> SetMarkerStyle(4);
 dh2 -> SetMarkerSize(1);
 //dh2 -> Draw("lego");
 RooPlot* framex = x.frame() ;
 RooPlot* framey = y.frame();
 model.plotOn(framey , LineColor(kRed));
 model.plotOn(framex,  LineColor(kRed));
 //gauss.plotOn(framex);
 // framey -> Draw();
 /* RooRealVar mu("meanx","meanx",2,-10,10) ;
    RooRealVar sigma("sigmax","sigmax",1,0.,10.);
    RooGaussian bkg(x  = -y , mu , sigma );
 */

}
