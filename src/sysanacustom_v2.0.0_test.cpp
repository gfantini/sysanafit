/* 
 * Questa macro fitta (con impostazioni di default) tutti e quattro i parametri
 * e li stampa a schermo
 *
 * PARAMETRI IN INPUT
 * inputFileName = nome del file root che contiene tutti gli istogrammi (cfr. sysanacustom)
 */


#include "sysanacustom_v2.0.0.cpp"
TGraph gLL;

void test(const string inputFileName)
{
  // initialize
  GFitinterf fitSL(inputFileName,"SL");
  // define some graphs to store profiled likelihood
  double value,error;
  fitSL.Fit(true);// perform fit silently, save result into class member
  // compute profiled likelihood
  value = fitSL.GetParameter("ZetaSL");
  error = fitSL.GetSymmetricError("ZetaSL");
  gLL = fitSL.GetProfiledLikelihood("ZetaSL",value-2*error,value+2*error,0.05*error); // model,min,max,step
  gLL.SetMarkerStyle(kFullDotMedium);
  // FIXME -- plotting does not work 
  // FIXME gLL.DrawClone("AP");

  // get best fit histogram
  // get TFitResultPtr
  ROOT::Fit::FitResult resultSL = fitSL.GetFitResult();
  //-----------------------------------------------------------------
  /*
  TFile *fp = new TFile("fitSL.root","RECREATE");
  fp->WriteObject(&resultSL,"fitSL"); // try and write it to disk
  //  THIS IS HOW YOU CAN READ THE FIT RESULT
  //  root [8] ROOT::Fit::FitResult* pfr;
  //  root [9] fp->GetObject("fitSL;1",pfr);
  fp->Close();
  */
  
  TFile* f = new TFile( "fitSL.root","RECREATE" );
  f->mkdir("SL");
  f->cd("SL");
  gDirectory->WriteObject(&gLL,"ProfiledLikelihood","SingleKey | WriteDelete");
  f->Close();
  
  // save fit result to file (method uses UPDATE on file root)
  fitSL.Save("fitSL.root",0); // 2nd argument is a way of keeping track of different versions of the fit in the same file
  std::cout << "Check fitSL.root for fit result" << std::endl;

  exit(0); // exit here --- modified to avoid crashing

  GFitinterf fit00(inputFileName,"00");
  GFitinterf fitGamma(inputFileName,"gamma");
  GFitinterf fitOmega(inputFileName,"omega");
  // perform fit
  fit00.Fit();
  fitGamma.Fit();
  fitOmega.Fit();
  // output best fit histogram
  // output fit result pointer
  map<string,ROOT::Fit::FitResult> FitResult;
  FitResult["SL"] = fitSL.GetFitResult();
  FitResult["00"] = fit00.GetFitResult();
  FitResult["gamma"] = fitGamma.GetFitResult();
  FitResult["omega"] = fitOmega.GetFitResult();
  
  TCanvas* c3 = new TCanvas();
  // compute profiled likelihood
  value = fitGamma.GetParameter("gamma");
  error = fitGamma.GetSymmetricError("gamma");
  gLL = fitGamma.GetProfiledLikelihood("gamma",value-2*error,value+2*error,0.05*error); // model,min,max,step
  gLL.SetMarkerStyle(kFullDotMedium);
  gLL.DrawClone("AP");
  
  // dump on stdout fit results
  TH1D* hData;
  TH1D* hFit;
  for(auto model : {"SL","00","gamma","omega"}){
    FitResult[model].Print(std::cout); // print symmetric error
    // print asymmetric error
    if(FitResult[model].HasMinosError(2)){
      cout << "MINOS error [2]: LOW " << FitResult[model].LowerError(2) << endl;
      cout << "MINOS error [2]: UP  " << FitResult[model].UpperError(2) << endl;
    }
    if(FitResult[model].HasMinosError(3)){ // omega
      cout << "-- ImOmega --" << endl;
      cout << "MINOS error [3]: LOW " << FitResult[model].LowerError(3) << endl;
      cout << "MINOS error [3]: UP  " << FitResult[model].UpperError(3) << endl;
    }
  }
  cout << "*********************************" << endl;
  /*
    for (auto model : {fitSL,fit00,fitGamma,fitOmega}) {
    // print bin contents and error
    hData = model.GetFitResultHdata(true);
    hFit = model.GetFitResultHtheory(true);
    cout << "--------------------------" << endl;
    }
  */
  
  cout << "Goodbye!" << endl;
}
