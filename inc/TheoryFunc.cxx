
Double_t ComputeIntegralSL2(double par,const Double_t xmin = 0.,const Double_t xmax = 12.);
Double_t ComputeIntegral002(const Double_t par,const Double_t xmin = 0.,const Double_t xmax = 12.);
Double_t ComputeIntegralGamma2(const Double_t par,const Double_t xmin = 0.,const Double_t xmax = 12.);
Double_t ComputeIntegralOmega2(const Double_t ReOmega,const Double_t ImOmega,const Double_t xmin = 0.,const Double_t xmax = 12.);
Double_t FuncSL(double dt , double par , const bool normalize){
  //Definizione dei parametri
  complex <double> ai (0., 1.);
  Double_t dm = 5.301*0.08935;
  Double_t gs = 1.;
  Double_t gl = 0.08935/51.7;
  complex<double> lambdas (0. , -gs/2.);
  complex<double> lambdal (dm , -gl/2.);
  complex<double> ak01 , ak0b1 , ak02 , ak0b2;
  double binwidth = 0.1;
  double t;
  double aidt;
  Double_t tmax = 50;
   t = dt;
   aidt = 0.;
   /*LOOP per l'integrazione , la variabile d'integrazione t viene incrementata di
     0.005 fino a raggiungere il valore tmax */
   while (t  < tmax){
     double t1 = t;
     double t2 = t-dt;
     ak01 = exp(-ai*lambdas*t1);
     ak0b1 = exp(-ai*lambdal*t1);
     ak02 = exp(-ai*lambdas*t2);
     ak0b2 = exp(-ai*lambdal*t2);
     aidt = aidt + pow(abs(ak01*ak0b2),2) + pow(abs(ak02*ak0b1), 2) - 2*(1.-par)*((ak01*ak0b2*conj(ak02*ak0b1))).real();     
     t = t + binwidth;
   }

   
   if(normalize) return aidt/ComputeIntegralSL2(par);
   else return  aidt;
   
}



Double_t Func00(double dt , double par , bool normalize){
  complex <double> ai (0., 1.);
  Double_t dm = 5.301*0.08935;
  Double_t gs = 1.;
  Double_t gl = 0.08935/51.7;
  complex<double> lambdas (0. , -gs/2.);
  complex<double> lambdal (dm , -gl/2.);
  complex<double> ak01 , ak0b1 , ak02 , ak0b2;
  double binwidth = 0.1;
  double t;
  double aidt;
  Double_t tmax = 50;
  double par1 = pow(0.00232 , 2);
  t = dt;
  aidt = 0.;
  /*LOOP per l'integrazione , la variabile d'integrazione t viene incrementata di
    0.005 fino a raggiungere il valore tmax */
  while (t  < tmax){
    double t1 = t;
    double t2 = t-dt;
    ak01 = exp(-ai*lambdas*t1);
    ak0b1 = exp(-ai*lambdal*t1);
     ak02 = exp(-ai*lambdas*t2);
     ak0b2 = exp(-ai*lambdal*t2);
     aidt = aidt + pow(abs(ak01*ak0b2),2) + pow(abs(ak02*ak0b1), 2) - 2*((ak01*ak0b2*conj(ak02*ak0b1))).real() + par/2*(-pow(abs(ak01*ak0b2),2) - pow(abs(ak02*ak0b1), 2) + 2*((ak01*ak0b2*conj(ak02*ak0b1)).real() - (ak01*ak0b2*ak02*ak0b1).real())) + 0.5*par/par1*pow(abs(ak01*ak02),2);     
     t = t + binwidth;
  } 
   if(normalize) return aidt/ComputeIntegral002(par);
   else return  aidt;
}


Double_t FuncGamma(double dt , double par , bool normalize){
  complex <double> ai (0., 1.);
  Double_t dm = 5.301*0.08935;
  Double_t gs1 = 1.;
  Double_t gl1 = 0.08935/51.7;
  complex<double> lambdas (0. , -gs1/2.);
  complex<double> lambdal (dm , -gl1/2.);
  complex<double> ak01 , ak0b1 , ak02 , ak0b2;
  double binwidth = 0.1;
  double t;
  double aidt;
  Double_t tmax = 50;
  double gs = 1./(0.8935*1E-10)*6.58211951440*1E-25;
  double gl = 1./(5.17*1E-8)*6.58211951440*1E-25;
  double par1 = pow(0.00232 , 2);
  t = dt;
  aidt = 0.;
  /*LOOP per l'integrazione , la variabile d'integrazione t viene incrementata di
    0.005 fino a raggiungere il valore tmax */
  while (t  < tmax){
    double t1 = t;
    double t2 = t-dt;
    ak01 = exp(-ai*lambdas*t1);
    ak0b1 = exp(-ai*lambdal*t1);
     ak02 = exp(-ai*lambdas*t2);
     ak0b2 = exp(-ai*lambdal*t2);
     aidt =  aidt + (1+par/((gs - gl)*par1))*(pow(abs(ak0b1*ak02) , 2) + pow(abs(ak01*ak0b2) , 2)) - 2*(ak01*ak0b2*ak02*conj(ak0b1)).real() - 2*par/((gs - gl)*par1)*pow(abs(ak01*ak02),2);
     t = t + binwidth;
  }
  if(normalize) return aidt/ComputeIntegralGamma2(par);
   else return  aidt;
}

Double_t FuncOmega(double dt , double par2 , double par3, bool normalize){
  complex <double> ai (0., 1.);
  Double_t dm = 5.301*0.08935;
  Double_t gs = 1.;
  Double_t gl = 0.08935/51.7;
  complex<double> lambdas (0. , -gs/2.);
  complex<double> lambdal (dm , -gl/2.);
  complex<double> ak01 , ak0b1 , ak02 , ak0b2;
  double binwidth = 0.1;
  double t;
  double aidt;
  Double_t tmax = 50;
  Double_t phi = 43.4*TMath::Pi()/180;
  complex<double> b = exp(-ai*phi);
  complex<double>omega (par2, par3);
  double par1 = pow(0.00232 , 2);
  t = dt;
  aidt = 0.;
  /*LOOP per l'integrazione , la variabile d'integrazione t viene incrementata di
    0.005 fino a raggiungere il valore tmax */
  while (t  < tmax){
    double t1 = t;
    double t2 = t-dt;
    ak01 = exp(-ai*lambdas*t1);
    ak0b1 = exp(-ai*lambdal*t1);
     ak02 = exp(-ai*lambdas*t2);
     ak0b2 = exp(-ai*lambdal*t2);
     aidt = aidt + pow(abs(ak0b1*ak02) , 2) + pow(abs(ak01*ak0b2) , 2) - 2*(ak01*ak0b2*conj(ak02*ak0b1)).real() + pow(abs(omega*ak01*ak02), 2)/par1 +  2*1/TMath::Sqrt(par1)*(pow(abs(ak01) , 2)*(omega*b*ak02*conj(ak0b2)).real() - pow(abs(ak02) , 2)*(omega*b*ak01*conj(ak0b1)).real());
     t = t + binwidth;
  } 
  if(normalize) return aidt/ComputeIntegralOmega2(par2 ,par3);
   else return  aidt;
}
void prova(){


 TCanvas *c1 = new TCanvas("c1" , "c1" , 500 , 400);
 TMultiGraph *mg = new TMultiGraph();
 //double par1[5] = {0 , 5e-3 , 4e-2 , 3e-1}; //SL
 //double par1[5] = {0 , 1e-6 , 4e-6 , 8e-6}; //00
 //double par1[5] = {0 , 1e-21 , 5e-21 , 1e-20}; //Gamma
 double par1[5] = {0 , -5e-4 , -1e-3 , -1.5e-3}; //ReOmega
 double par2[5] = {0 , -5e-4 , -1e-3 , -1.5e-3}; //IOmega
 for(int j = 0 ; j < 4 ; j++) {
  double par[4];
  par[2] = par1[j];
  par[3] = par2[j];
  double dt = 0; 
  int i = 0;
  int n = 12/0.1;
  double x[n], y[n];
   double x1[n], y1[n];
   

  while(i < n) {

    x[i] =dt;
    y[i] = fOmega(dt , par[2], par[3] , true);

    dt = dt+0.1;
    
    i++;
    
    }


  i = 0;
  dt = 0;
 

  TGraph *func1 = new TGraph(n , x , y);
  
  c1 -> SetGrid();
  func1 -> SetLineColor(kRed);
  func1 -> SetMarkerStyle(3);
  func1 -> SetMarkerColor(kRed);
  func1 -> SetMarkerSize(1);
  string str = to_string(par[2]);
  string str1 = to_string(par[3]); 
  TString s = "Funzione in base #omega con Re#omega = " + str + " e I#omega = " + str1; 
  func1 -> SetTitle(s);
  func1 ->  GetYaxis() -> SetTitle("I(#Deltat)");
  func1 -> GetXaxis() -> SetTitle("#Deltat/#tau_{s}");
  func1 -> GetHistogram() -> SetMinimum(0);
  func1->Draw("AP");
   while(i < n) {

    x1[i] =dt;
    y1[i] = FuncOmega(dt , par[2] , par[3] , true);
    dt = dt+0.1;
    i++;
    
    }

   TGraph *func2 = new TGraph(n , x1 , y1);
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
   TString st = "PlotOmega/Omega_" + str + ".pdf";
   c1 -> SaveAs(st);
  
   double x2[n] , y2[n];
   for (int i = 0 ; i < n ; i++) {
     x2[i] = x[i];
     y2[i] = abs(y[i] - y1[i]);
   }

   TGraph *func3 = new TGraph(n , x2 , y2);
   TString s1 = "Scarto delle funzioni in base #omega con Re#omega = " + str + "e I#omega = " + str1;
   TString s2 = "PlotOmega/DOmega_" + str + ".pdf";
   func3 -> SetTitle(s1);
   func3 -> SetLineColor(kGreen);
   func3 -> SetMarkerColor(kRed);
   func3 ->  GetYaxis() -> SetTitle("#Delta I(#Deltat)");
   func3 -> GetXaxis() -> SetTitle("#Deltat/#tau_{s}");
   func3 -> SetMarkerStyle(4);
   func3 -> SetMarkerSize(1);
   func3 -> Draw("AP");
   c1 -> Modified();
   c1 -> Update();
   c1 -> SaveAs(s2);
   //mg -> Add(v[j], "p");
 }

 //mg -> Draw("a");
 //c1 -> BuildLegend(); 
}



Double_t ComputeIntegralSL2(const Double_t par,const Double_t xmin = 0.,const Double_t xmax = 12.)
{
  auto lambda = [](double *x, double *p){ return FuncSL(x[0],p[0] , false); };
  TF1 fSL("f1", lambda ,0,100,1); // 1 == number of parameters
  fSL.SetParameter(0,par);
  return fSL.Integral(xmin,xmax); // relative error on integral 10^-12
}

Double_t ComputeIntegral002(const Double_t par,const Double_t xmin = 0.,const Double_t xmax = 12.)
{
  auto lambda = [](double *x, double *p){ return Func00(x[0], p[0],false); };
  TF1 f00("f1", lambda ,0,100,1); // 1 == number of parameters
  f00.SetParameter(0,par);
  return f00.Integral(xmin,xmax); // relative error on integral 10^-12
}

Double_t ComputeIntegralGamma2(const Double_t par,const Double_t xmin = 0.,const Double_t xmax = 12.)
{
  auto lambda = [](double *x, double *p){ return FuncGamma(x[0], p[0],false); };
  TF1 fGamma("f1", lambda ,0,100,1); // 1 == number of parameters
  fGamma.SetParameter(0,par);
  return fGamma.Integral(xmin,xmax); // relative error on integral 10^-12
}

Double_t ComputeIntegralOmega2(const Double_t ReOmega,const Double_t ImOmega,const Double_t xmin = 0.,const Double_t xmax = 12.)
{
  auto lambda = [](double *x, double *p){ return FuncOmega(x[0],p[0],p[1],false); }; // x[0] means "|dt|"
	TF1 fOmega("f1",lambda,0,100,2); // 2 == number of parameters
	fOmega.SetParameter(0,ReOmega);
	fOmega.SetParameter(1,ImOmega);
	return fOmega.Integral(xmin,xmax); // relative error on integral 10^-12
}
