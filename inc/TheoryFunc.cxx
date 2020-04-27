Double_t ComputeIntegralSL2(double par,const Double_t xmin = 0.,const Double_t xmax = 12.);
Double_t ComputeIntegral002(const Double_t par,const Double_t xmin = 0.,const Double_t xmax = 12.);
Double_t ComputeIntegralGamma2(const Double_t par,const Double_t xmin = 0.,const Double_t xmax = 12.);
Double_t ComputeIntegralOmega2(const Double_t ReOmega,const Double_t ImOmega,const Double_t xmin = 0.,const Double_t xmax = 12.);

Double_t FuncSL(Double_t *dt , Double_t *par){

 complex <double> ai (0., 1.);
  Double_t dm = 5.293*0.08954;
  Double_t gs = 1.;
  Double_t gl = 0.08954/51.2;
  complex<double> lambdas (0. , -gs/2.);
  complex<double> lambdal (dm , -gl/2.);
  complex<double> ak01 , ak0b1 , ak02 , ak0b2;
  double aidt= 0;  
  double t1 = *dt;
  double t2 = *dt - par[0];
  ak01 = exp(-ai*lambdas*(t1));
  ak0b1 = exp(-ai*lambdal*(t1));
  ak02 = exp(-ai*lambdas*(t2));
  ak0b2 = exp(-ai*lambdal*(t2));
  aidt = pow(abs(ak01*ak0b2),2) + pow(abs(ak02*ak0b1), 2) - 2*(1.-par[1])*((ak01*ak0b2*conj(ak02*ak0b1))).real();
  return abs(aidt);

}

Double_t Func00(double *dt , double *par){
  complex <double> ai (0., 1.);
  Double_t dm = 5.293*0.08954;
  Double_t gs = 1.;
  Double_t gl = 0.08954/51.2;
  complex<double> lambdas (0. , -gs/2.);
  complex<double> lambdal (dm , -gl/2.);
  complex<double> ak01 , ak0b1 , ak02 , ak0b2;
  double t;
  double aidt = 0;
  Double_t tmax = 50;
  double par1 = pow(0.002232 , 2);
  t = *dt;
  double t1 = t;
  double t2 = t-par[0];
  ak01 = exp(-ai*lambdas*t1);
  ak0b1 = exp(-ai*lambdal*t1);
  ak02 = exp(-ai*lambdas*t2);
  ak0b2 = exp(-ai*lambdal*t2);
  aidt =  pow(abs(ak01*ak0b2),2) + pow(abs(ak02*ak0b1), 2) - 2*((ak01*ak0b2*conj(ak02*ak0b1))).real() + par[1]/2*(-pow(abs(ak01*ak0b2),2) - pow(abs(ak02*ak0b1), 2) + 2*((ak01*ak0b2*conj(ak02*ak0b1)).real() - (ak01*ak0b2*ak02*ak0b1).real())) + 0.5*par[1]/par1*pow(abs(ak01*ak02),2);     
  return  aidt;
}

Double_t FuncGamma(double *dt , double *par){
  complex <double> ai (0., 1.);
  Double_t dm = 5.293*0.08954;
  Double_t gs1 = 1.;
  Double_t gl1 = 0.08954/51.2;
  complex<double> lambdas (0. , -gs1/2.);
  complex<double> lambdal (dm , -gl1/2.);
  complex<double> ak01 , ak0b1 , ak02 , ak0b2;
   double t;
  double aidt;
  double gs = 1./(0.8954*1E-10)*6.58211951440*1E-25;
  double gl = 1./(5.12*1E-8)*6.58211951440*1E-25;
  double par1 = pow(0.002232 , 2);
  t = *dt;
  aidt = 0.;
  double t1 = t;
  double t2 = t-par[0];
  ak01 = exp(-ai*lambdas*t1);
  ak0b1 = exp(-ai*lambdal*t1);
  ak02 = exp(-ai*lambdas*t2);
  ak0b2 = exp(-ai*lambdal*t2);
  aidt =  (1+par[1]/((gs - gl)*par1))*(pow(abs(ak0b1*ak02) , 2) + pow(abs(ak01*ak0b2) , 2)) - 2*(ak01*ak0b2*ak02*conj(ak0b1)).real() - 2*par[1]/((gs - gl)*par1)*pow(abs(ak01*ak02),2);
  return  aidt;
}

Double_t FuncOmega(double *dt , double *par){
  complex <double> ai (0., 1.);
  Double_t dm = 5.293*0.08954;
  Double_t gs = 1.;
  Double_t gl = 0.08954/51.2;
  complex<double> lambdas (0. , -gs/2.);
  complex<double> lambdal (dm , -gl/2.);
  complex<double> ak01 , ak0b1 , ak02 , ak0b2;
  double t;
  double aidt;
  Double_t phi = 43.4*TMath::Pi()/180;
  complex<double> b = exp(-ai*phi);
  complex<double>omega (par[1], par[2]);
  double par1 = pow(0.002232 , 2);
  t = *dt;
  aidt = 0.;
  
  double t1 = t;
  double t2 = t-par[0];
  ak01 = exp(-ai*lambdas*t1);
  ak0b1 = exp(-ai*lambdal*t1);
  ak02 = exp(-ai*lambdas*t2);
  ak0b2 = exp(-ai*lambdal*t2);
  aidt = pow(abs(ak0b1*ak02) , 2) + pow(abs(ak01*ak0b2) , 2) - 2*(ak01*ak0b2*conj(ak02*ak0b1)).real() + pow(abs(omega*ak01*ak02), 2)/par1 +  2*1/TMath::Sqrt(par1)*(pow(abs(ak01) , 2)*(omega*b*ak02*conj(ak0b2)).real() - pow(abs(ak02) , 2)*(omega*b*ak01*conj(ak0b1)).real());
  return  aidt;
}


Double_t FProva1(Double_t dt , Double_t par , Double_t par1 , bool normalize){

  TF1 f("f" , FuncOmega , 0 , 12 , 3);
  f.FixParameter(0 , dt);
  f.FixParameter(1 , par);
  f.FixParameter(2 , par1);
  f.SetNpx(200);
  double aidt = 0;
  aidt = f.Integral(dt , 100);
  if(normalize) return aidt/ComputeIntegralOmega2(par , par1);  
  else return  aidt;
}
Double_t FProva(Double_t dt , Double_t par , bool normalize){

  TF1 f("f" , FuncSL , 0 , 12 , 2);
  f.FixParameter(0 , dt);
  f.FixParameter(1 , par);
  f.SetNpx(100000);
  double aidt = 0;
  aidt = f.Integral(dt , 50);
  if(normalize) return aidt/ComputeIntegralSL2(par);  
  else return  aidt;
}



void prova(){
  double par1[7] = {0.5e-4 , 1e-3 , 9e-3 , 4e-2 , 6e-2 , 8e-2 , 3e-1}; //SL
  //double par1[7] = {0.5e-6 , 1e-6, 2e-6, 3e-6 , 4e-6 , 5e-6 , 8e-6}; //00
  //double par1[9] = {1e-22 , 5e-22 , 8e-22 , 1e-21 , 2e-21 , 3e-21 , 5e-21 , 8e-21 , 1e-20}; //Gamma
  //double par1[5] = { -5e-4 , -9e-4 , -1e-3}; //ReOmega
  //double par2[5] = { -5e-4 , -9e-4, -1e-3}; //IOmega
  double par[2];
  double x3[7] , y3[7] , x4[9];
  int h = 0;
  int n = 12/0.1;
  for(int g = 0 ; g < 7 ; g++ ){
    //for(int j = 0 ; j < 3 ; j++){
      par[0] = par1[g];
      //par[1] = par2[j];
      int i = 0;
      double x[n], y[n];
      double  y1[n];
      double y2[n];
      double dt = 0;
      TCanvas *c1 = new TCanvas("c1" , "c1" , 500 , 400);
      c1 -> SetGrid();
      int npar;
      int ndim;
      while(i < n) {
	
	x[i] =dt;
	/*y[i] = FProva1(dt , par[0] , par[1], true);
	y1[i] = fOmega(dt , par[0] , par[1] , true);*/ //Omega
	y[i] = FProva(dt , par[0] , true);
	y1[i] = fZetaSL(dt , par[0] , true);
	y2[i] = abs(y[i] - y1[i])*100/((y[i] + y1[i])/2);
	if(i == 0){
	  cout << y[i] - y1[i] << endl;
	  cout << (y[i] + y1[i])/2 << endl;
	  cout << y[i] << endl;
	    cout << y1[i] << endl;
	}
	//cout << i << endl;  
	dt = dt+0.1;    
	i++;
	
	
      }
      string str = to_string(par[0]*1e4);
      string str1 = to_string(par[1]*1e4);
      TString s = "Funzione in base #zeta_{SL} con #zeta_{SL} = " + str+ "10^{-4}"; //  + " e Im#omega = " + str1 + "10^{-4}";
      TGraph *func1 = new TGraph(n , x , y1);
      TGraph *func2 = new TGraph(n , x , y);
      TGraph *func3 = new TGraph(n , x , y2);
      func1 -> SetTitle(s);
      func1 ->  GetYaxis() -> SetTitle("I(#Deltat)");
      func1 -> GetXaxis() -> SetTitle("#Deltat/#tau_{s}");
      func1 -> SetLineColor(kRed);
      func1 -> SetMarkerStyle(3);
      func1 -> SetMarkerColor(kRed);
      func1 -> SetMarkerSize(1);
      func1 -> GetHistogram() -> SetMinimum(0);
      func1 -> Draw();
      func2 -> SetLineColor(kGreen);
      func2 -> SetMarkerColor(kGreen);
      func2 -> SetMarkerStyle(4);
      func2 -> SetMarkerSize(1);
      func2->Draw("P SAME");
      TLegend  *legend = new TLegend(0.2 , 0.2 , 0.4 , 0.4);
      legend -> SetHeader("Theory curves" , "C");
      legend -> AddEntry(func1 , "GF" , "P");
      legend -> AddEntry(func2 , "RDA" , "P");
      legend -> Draw(); 
      c1 -> Modified();
      c1 -> Update();
      TString st = "PlotSL1/Gamma1_" + str + ".pdf" ;//+  "_" + str1 + ".pdf";
      c1 -> SaveAs(st);
      TString s1 = "Scarto delle funzioni in base #zeta_{SL} con #zeta_{SL} =  " + str + "10^{-4}";// + " e Im#omega = " + str1 + "10^{-4}";
      TString s2 = "PlotSL1/DGamma1_" + str + ".pdf" ;// +  "_" + str1 + ".pdf";
      func3 -> SetTitle(s1);
      
      func3 -> SetLineColor(kGreen);
      func3 -> SetMarkerColor(kRed);
      func3 ->  GetYaxis() -> SetTitle("%#Delta I(#Deltat)");
      func3 -> GetXaxis() -> SetTitle("#Deltat/#tau_{s}");
      func3 -> SetMarkerStyle(4);
      func3 -> SetMarkerSize(1);
      func3 -> Draw("AP");
      c1 -> Modified();
      c1 -> Update();
      c1 -> SaveAs(s2);
      y3[g] =   TMath::MaxElement(n , func3 -> GetY());
      x3[g] = par[0];
      // x4[h] = par[1]*1e4;
      h++;
      delete c1;
      
      
      
      // }
  }
  
  TCanvas *c2 = new TCanvas("c2" , "c2" , 500 , 400);
  TGraph *func4 = new TGraph(7 , x3 , y3);
  // func3 ->  GetXaxis() -> SetRangeUser(0, 1.5e-3);
  func4 -> SetTitle("%Difference between the functions");// ;Re#omega 10^{-4} ; Im#omega 10^{-4}; %#Delta I(#Deltat)");
  func4 ->  GetYaxis() -> SetTitle("%#Delta I(#Deltat)");
  func4 -> GetXaxis() -> SetTitle("#zeta_{SL}");
  func4 -> SetMarkerStyle(4);
  func4 -> SetMarkerSize(1);
  func4 -> SetMarkerColor(kRed);
  func4-> Draw("AP");
  c2 -> SetGrid();
  c2 -> Modified();
  c2 -> Update();
  c2 -> SaveAs("PlotSL1/Prova.pdf");

}


Double_t ComputeIntegralSL2(const Double_t par,const Double_t xmin = 0.,const Double_t xmax = 12.)
{
  auto lambda = [](double *x, double *p){ return FProva(x[0],p[0] , false); };
  TF1 fSL("f1", lambda ,0,100,1); // 1 == number of parameters
  fSL.SetParameter(0,par);
  return fSL.Integral(xmin,xmax); // relative error on integral 10^-12
}

Double_t ComputeIntegral002(const Double_t par,const Double_t xmin = 0.,const Double_t xmax = 12.)
{
  auto lambda = [](double *x, double *p){ return FProva(x[0], p[0],false); };
  TF1 f00("f1", lambda ,0,100,1); // 1 == number of parameters
  f00.SetParameter(0,par);
  return f00.Integral(xmin,xmax); // relative error on integral 10^-12
}



Double_t ComputeIntegralGamma2(const Double_t par,const Double_t xmin = 0.,const Double_t xmax = 12.)
{
  auto lambda = [](double *x, double *p){ return FProva(x[0], p[0],false); };
  TF1 fGamma("f1", lambda ,0,100,1); // 1 == number of parameters
  fGamma.SetParameter(0,par);
  return fGamma.Integral(xmin,xmax); // relative error on integral 10^-12
}

Double_t ComputeIntegralOmega2(const Double_t ReOmega,const Double_t ImOmega,const Double_t xmin = 0.,const Double_t xmax = 12.)
{
  auto lambda = [](double *x, double *p){ return FProva1(x[0],p[0],p[1],false); }; // x[0] means "|dt|"
	TF1 fOmega("f1",lambda,0,100,2); // 2 == number of parameters
	fOmega.SetParameter(0,ReOmega);
	fOmega.SetParameter(1,ImOmega);
	return fOmega.Integral(xmin,xmax); // relative error on integral 10^-12
}
