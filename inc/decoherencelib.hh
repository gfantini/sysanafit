// Guido Fantini
// 13.10.2017 rese coerenti le costanti fisiche come variabili globali lette da tutte le funzioni f**** che implementano i modelli teorici
//            testati i fit wrt versione prima della modifica e ci sono piccole variazioni perché ho aggiornato le costanti da pdg vecchi (ads)
//            alla stessa versione delle costanti, riportate sul pdg 2016
// 23.08.2017 aggiunta e testata (a occhio) fZeta00
// 24.08.2017 aggiunta e testata (a occhio) fOmega
// 13.09.2017 aggiunte ComputeIntegral... 00/Gamma
//            should be implemented some tool to check the normalization is done properly

/* USAGE INSTRUCTIONS: INTERACTIVE MODE
 * root [0] .L decoherencelib.hh
 * root [1] const Double_t par[4] = {<modelID>,<norm_par>,<par1>,<par2>};
 * root [2] TH1D hTheory = GetTheoryHistogram(par,12,0.,12.);
 * root [3] hTheory.Draw();
 */

/* CODER INFO:
 * The ONLY handle users should have to get something from here is GetTheoryHistogram.
 * It assumes a precise meaning of the "par" vector, hard coded inside him and ComputeIntegral functions, explained before its prototype.
 * - - - - -
 * physical constants begin with g_ or are hard-coded inside -> should find a way to deal with this
 * All f<ParName> functions are general, not normalized
 * All ComputeIntegral<ParName> are not meant to be called from outside. They are called by GetTheoryHistogram.
 * - - - - - - -
 * KNOWN BUGS
 * - - - - - - -
 * GetTheoryHistogram produces a TH1D named "GetTheoryHistogram".
 * Multiple calls without deleting the previous one will lead to memory leak and crash.
 * This is not a problem in fitting & other stuff AS LONG AS IT IS CALLED ONLY ONCE AND IN A LAMBDA FUNCTION so that at each call the TH1D is deleted
 */
#include "phys_const_v00.cc"

#ifndef _DECOHERENCELIB_
#define _DECOHERENCELIB_

// To modify those global variables use ResetGlobalPhysicalConstants
// definition of physical constants later to be put somewhere else
//const Double_t g_PlanckReduced = 6.592119514*1E-22; // [MeV s]
const Double_t g_PlanckReduced = 6.58211951440*1E-25; // [GeV s] PDG 2016
Double_t g_TauS = 0.8954*1E-10; // [s] PDG 2015
Double_t g_TauL = 5.116*1E-8; // [s] PDG 2015
Double_t g_DeltaM = 0.5293*1E10; // [htagliato s^-1]
// |eta+-| = A(KL -> pi+pi-) / A(KS -> pi+pi-) rev. KL p. 35/50 [PDG2016]
Double_t g_modeta = 2.232*1E-3; // [adim]
Double_t g_phieta = 0.7574728954; //[rad?] pdg2013 (43.4° w/o CPT) | Err 0.5° pdg 2017


// definition of global variables
int g_GetTheoryHistogramCounter = 1;
int g_MultiplyHistogramCounter = 1;
int g_GetSmearedHistogramCounter = 1;


// ********************************************************************
/* - - - - - IMPLEMENTATION OF THEORETICAL MODELS  - - - - - - - - - -
 * For each <model> should exist, in growing order of complexity
 * f<model>                      that returns I(Dt,parameters) non-normalized, non-integrated in any bin
 * ComputeIntegral<model>        that returns the integral of f<model> from Dt(min) ... Dt(max)
 *                               in such a way that calling ONLY this function it is possible to get the histogram of data
 *
 *
 *
 *
 */
// ********************************************************************
// dt == |t1-t2|/TauS
Double_t fZetaSL(const Double_t dt,const Double_t zetaSL,const bool normalize);
// Integrate fZetaSL (non-normalized)
// to get the bin content (normalized) for fitting use
// ComputeIntegralSL(zetaSL,bin_min,bin_max) / ComputeIntegralSL(zetaSL, fit_range_min, fit_range_max)
Double_t ComputeIntegralSL(const Double_t zetaSL,const Double_t xmin = 0.,const Double_t xmax = 12.);
Double_t fGamma(const Double_t dt,const Double_t gamma , const bool normalize); // WARNING <-- global constants not implemented!
Double_t ComputeIntegralGamma(const Double_t gamma,const Double_t xmin = 0.,const Double_t xmax = 12.);
Double_t fZeta00(const Double_t dt,const Double_t zeta00 , const bool normalize); // WARNING <-- global constants not implemented!
Double_t ComputeIntegral00(const Double_t zeta00,const Double_t xmin = 0.,const Double_t xmax = 12.);
Double_t fOmega(const Double_t dt,const Double_t ReOmega,const Double_t ImOmega , const bool normalize); // WARNING <-- global constants not implemented!
Double_t ComputeIntegralOmega(const Double_t ReOmega,const Double_t ImOmega,const Double_t xmin = 0.,const Double_t xmax = 12.);

// **************************** BEGIN SECTION ABOUT HISTOGRAM BUILDING
// par[0] defines which theoretical model
// par[1] is always the normalization
// = 0    SL     par[2] = zetaSL
// = 1    00     par[2] = zeta00
// = 2    GAMMA  par[2] = gamma
// = 3    OMEGA  par[2] = Re(omega)  par[3] = Im(omega)
// **************************** NOTICE: the error evaluation is performed ONLY the end of the sequence
TH1D GetTheoryHistogram(const Double_t *par,const int nbin,const double dt_min,const double dt_max);
TH1D MultiplyHistogram(TH1D histo1,TH1D histo2,const int nbin,const double DtMin,const double DtMax, const bool debug);
TH1D GetSmearedHistogram(TH1D htrue,TH2D SmMat,const bool debug); // [debug=false]
void ResetGlobalPhysicalConstants(Double_t TauS,Double_t TauL,Double_t DeltaM,Double_t modeta,Double_t phieta);

//-----------------------------------------------------------------------------------------------
// UNDER THIS LINE IMPLEMENTATION BEGINS . . . . are you sure you really wanna dig in this sh**?
//-----------------------------------------------------------------------------------------------
Double_t fZetaSL(const Double_t dt,const Double_t zetaSL,const bool normalize)
{
	Double_t GS = 1./g_TauS;
	Double_t GL = 1./g_TauL;
	
	Double_t NonNormalized = TMath::Exp(-g_TauS*dt/g_TauS) + TMath::Exp(-g_TauS*dt/g_TauL) -2*(1.-zetaSL)*TMath::Exp(-0.5*(GS+GL)*g_TauS*dt) * TMath::Cos(g_DeltaM*g_TauS*dt);
	if(normalize) return NonNormalized/ComputeIntegralSL(zetaSL);
	else return NonNormalized;
}

Double_t ComputeIntegralSL(const Double_t zetaSL,const Double_t xmin = 0.,const Double_t xmax = 12.)
{
	auto lambda = [](double *x, double *p){ return fZetaSL(x[0],p[0],false); };
	TF1 fSL("f1",lambda,0,100,1); // 1 == number of parameters
	fSL.SetParameter(0,zetaSL);
	return fSL.Integral(xmin,xmax); // relative error on integral 10^-12
}

//------------------------------------------------------------------------------------
Double_t fGamma(Double_t dt,Double_t gamma, const bool normalize)
{// versione della funzione gamma mostrata al kaon meeting 05/2017
	// all numerical constants taken from PDG 2016
	// err(hbar) circa 6 ppb (trascurabile)
	//double hbar = 6.58211951440*1E-25; // [GeV s]
	double hbar = g_PlanckReduced;
	
	// |eta+-| = A(KL -> pi+pi-) / A(KS -> pi+pi-) rev. KL p. 35/50
	//double modeta = 2.232*1E-3; // adim. (same ad ADS)
	double modeta = g_modeta;
	
	// KS width = hbar / tauS
	// err(TauS) +- 0.00033*1E-10 [s]
	//double GammaS = 1./(0.89564*1E-10); // [hbar s^-1]
	double GammaS = 1./g_TauS;
	
	// KL width = hbar / tauL
	// err(TauL) +- 0.021*1E-8 [s]
	//double GammaL = 1./(5.099*1E-8); // [hbar s^-1]
	double GammaL = 1./g_TauL;
	
	// GammaS - GammaL
	double DG = GammaS - GammaL;
	
	// GammaS + GammaL
	double SG = GammaS + GammaL;
	
	// m(KS) - m(KL)
	// err(dm) =  +- 0.0010*1E10
	//double dm =  0.5289*1E10; // [hbar sec^-1] (same as ADS)
	double dm = g_DeltaM;
	
	// - - - - - - - - -
	// CONSTANT REDEFINITION (why?)
	/*
	GammaS = 1./(0.8953*1E-10);
	GammaL = 1./(511.6*1E-10);
	DG = GammaS - GammaL;
	SG = GammaS + GammaL;
	*/
	// - - - - - - - - -
	
	
	// convert dt from [units of taus] --> [s]
	// assuming hbar = 1.
	dt = dt/GammaS;
	// convert gamma from [GeV] --> [s^-1]
	gamma = gamma/hbar;
	Double_t result = 0.;
	result += (1.+gamma/(DG*modeta*modeta))*4.*cosh(0.5*DG*dt)*exp(-0.5*SG*dt);
	result += -4.*cos(dm*dt)*exp(-0.5*SG*dt);
	result += (-2.*SG/GammaS)*(gamma/(DG*modeta*modeta))*exp(-GammaS*dt);
	
	if(normalize ) {return result/ComputeIntegralGamma(gamma*hbar);}
	else return result;
}

Double_t ComputeIntegralGamma(const Double_t gamma,const Double_t xmin = 0.,const Double_t xmax = 12.)
{
  auto lambda = [](double *x, double *p){ return fGamma(x[0],p[0] , false); };
  TF1 fGamma("f1",lambda,0,100,1); // 1 == number of parameters
  fGamma.SetParameter(0,gamma);
  return fGamma.Integral(xmin,xmax); // relative error on integral 10^-12
}
Double_t fZeta00(const Double_t dt,const Double_t zeta00, const bool normalize)
{
	
	// all those physical constants are copied from the phys_const in FITINTERF
	//Double_t etapm_mod = 2.232e-3; // same as ADS
	Double_t etapm_mod = g_modeta;
	//Double_t etapm_phi = 0.7574728954;
	Double_t etapm_phi = g_phieta;
	//Double_t delta_mSL = 0.5289; // 10^10 [hbar sec^-1] (same as ADS)
	Double_t delta_mSL = g_DeltaM*1E-10;
	//Double_t tauS      = 0.8953; // 10^-10 [sec] (same as ADS)
	Double_t tauS = g_TauS*1E10;
	//Double_t tauL      = 511.6;  // 10^-10 [sec] (same as ADS) pdg 2013
	Double_t tauL = g_TauL*1E10;
	// some declarations and definitions in order to match ADS conventions
	TComplex aux1,aux2,v1;
	TComplex aI  = TComplex(0.,1.);
	TComplex a2I = 2.*aI;
	Double_t p_1  = etapm_phi;
	Double_t p_2  = p_1;
	Double_t m_1 = etapm_mod;
	Double_t m_2 = m_1;
	Double_t dm  = delta_mSL*tauS;
	Double_t gs  = 1.;
	Double_t gl  = tauS/tauL;
	Double_t sgk = gl + gs;
	Double_t z00 = zeta00;
	
	
	aux1 = TComplex::Exp(-aI*(p_1+p_2))*TComplex::Exp(((a2I*dm-3.*sgk)*dt)/2.);
	aux2 = -4.*TComplex::Power(aI,dm)*gl*gs*sgk*m_1*m_2*z00;
	aux2*= (exp(sgk*dt)-TComplex::Exp(a2I*(p_1 + p_2) + (-a2I*dm+sgk)*dt));
	aux2+= 4.*pow(dm,2)*(  2.*(TComplex::Exp(a2I*p_2)*exp(sgk*dt) + TComplex::Exp(a2I*p_1)*TComplex::Exp((-a2I*dm+sgk)*dt)) *gl*gs*m_1*m_2*(z00-2.) +TComplex::Exp(aI*(p_1+p_2))* TComplex::Exp(((-a2I*dm+3.*gl+gs)*dt)/2.) *gl*(gl*z00 + gs*(-2.*pow(m_2,2)*(z00-2.) + z00)) +TComplex::Exp(aI*(p_1+p_2))* TComplex::Exp(((-a2I*dm+gl+3.*gs)*dt)/2.) *gs*pow(m_1,2)*(gs*pow(m_2,2)*z00+ gl*(4.+(-2.+pow(m_2,2))*z00)));
	aux2+= pow(sgk,2)*(gl*(gl*z00 + gs*(-2.*pow(m_2,2)*(z00-2.) + z00))* TComplex::Exp(aI*(p_1+p_2))*TComplex::Exp(((-a2I*dm + 3*gl + gs)*dt)/2.)-2.*TComplex::Exp((-a2I*dm +sgk)*dt)*gl*gs*m_1*m_2*(TComplex::Exp(2.*aI*dm*dt)* (-1.*(TComplex::Exp(2.*aI*p_2)*(z00-2.))+z00) + TComplex::Exp(2.*aI*p_1)*(2. + (-1. + TComplex::Exp(2.*aI*p_2))*z00)) +TComplex::Exp(aI*(p_1 + p_2) + ((-a2I*dm+gl+3.*gs)*dt)/2.)* gs*pow(m_1,2)*(gs*pow(m_2,2)*z00 + gl*(4. + (pow(m_2,2)-2.)*z00)));
	v1 = aux1*aux2;
	if(normalize){return v1.Re()/(4.*gl*gs*sgk*(4.*pow(dm,2) + pow(sgk,2)))/(ComputeIntegral00(zeta00));}
	else {return v1.Re()/(4.*gl*gs*sgk*(4.*pow(dm,2) + pow(sgk,2)));};
}
Double_t ComputeIntegral00(const Double_t zeta00,const Double_t xmin = 0.,const Double_t xmax = 12.)
{
  auto lambda = [](double *x, double *p){ return fZeta00(x[0],p[0] , false); };
  TF1 f00("f1",lambda,0,100,1); // 1 == number of parameters
  f00.SetParameter(0,zeta00);
  return f00.Integral(xmin,xmax); // relative error on integral 10^-12
}
Double_t fOmega(const Double_t dt,const Double_t ReOmega,const Double_t ImOmega, const bool normalize)
{
	// all those physical constants are copied from the phys_const in FITINTERF
	//Double_t etapm_mod = 2.232e-3; // same as ADS
	Double_t etapm_mod = g_modeta;
	//Double_t etapm_phi = 0.7574728954;
	Double_t etapm_phi = g_phieta;
	//Double_t delta_mSL = 0.5289; // 10^10 [hbar sec^-1] (same as ADS)
	Double_t delta_mSL = g_DeltaM*1E-10;
	//Double_t tauS      = 0.8953; // 10^-10 [sec] (same as ADS)
	Double_t tauS = g_TauS*1E10;
	//Double_t tauL      = 511.6;  // 10^-10 [sec] (same as ADS) pdg 2013
	Double_t tauL = g_TauL*1E10;
	
	TComplex aux1;
	TComplex aI  = TComplex(0.,1.);
	TComplex aI2 = 2.*aI;
	Double_t p_1 = etapm_phi;
	Double_t p_2 = p_1;
	Double_t m_1 = etapm_mod;
	Double_t m_2 = m_1;
	Double_t dm  = delta_mSL*tauS;
	Double_t DG  = 1. - tauS/tauL;
	Double_t SG  = 1. + tauS/tauL;
	Double_t gs  = 1.;
	Double_t gl  = tauS/tauL;
	
	Double_t modomega = sqrt(ReOmega*ReOmega + ImOmega*ImOmega);
	Double_t phiomega = TMath::ATan2(ImOmega,ReOmega);
	
	// when you see "/**/" means breakline in old FORTRAN function
	aux1= ( modomega*( (-4.*TComplex::Exp( (aI*(2.*p_1-2.*phiomega+(-2.*dm+aI*SG)*dt))/2. )* /**/ m_1)/(aI2*dm+gl+3.*gs) + /**/ (4.*TComplex::Exp(aI*(p_2-phiomega+aI*gs*dt))*m_2)/ /**/ (aI2*dm+gl+3.*gs) + /**/ (4.*TComplex::Exp(-aI*(p_2+phiomega)-gl*dt)*pow(m_1,2)*m_2)/ /**/ (-aI2*dm+3.*gl+gs)+ /**/ (4.*TComplex::Exp(-aI*(p_1+phiomega)- /**/ ((-aI2*dm+gl+gs)*dt)/2.)* /**/ m_1*pow(m_2,2))/(aI2*dm-3.*gl-gs)+ /**/ modomega/(exp(gs*dt)*gs) /**/ -(2.*TComplex::Exp(-aI*p_1-aI*p_2-((-aI2*dm+SG)*dt)/2.)* /**/ m_1*m_2*modomega)/(-aI2*dm+SG) - /**/ (2.*TComplex::Exp((aI*(2.*p_1+2.*p_2 + /**/ (-2.*dm + aI*SG)*dt))/2.)* /**/ m_1*m_2*modomega)/(aI2*dm+SG) + /**/ (pow(m_1,2)*pow(m_2,2)*modomega)/(exp(gl*dt)*gl)))/2.;
	aux1+=( (TComplex::Exp(aI*p_1-((2.*gl+gs)*dt)/2.)*m_1)/SG - /**/ (TComplex::Exp(aI*p_2 + aI*dm*dt-((gl+2.*gs)*dt)/2.) /**/*m_2)/SG + /**/(2.*TComplex::Exp(aI*phiomega+aI*dm*dt- /**/ ((gl+2.*gs)*dt)/2.) /**/ *modomega)/ /**/ (aI2*dm-gl-3.*gs) + /**/ (2.*TComplex::Exp(aI*(p_1+p_2+phiomega) - /**/ ((2.*gl+gs)*dt)/2. )*m_1* /**/ m_2*modomega) / (aI2*dm+3.*gl+gs))* /**/ TComplex::Conjugate(TComplex::Exp(aI*p_1+(gs*dt)/2.)*m_1 - /**/ TComplex::Exp(aI*p_2+aI*dm*dt+(gl*dt)/2.)*m_2);
	if(normalize) return aux1.Re()/ComputeIntegralOmega(ReOmega , ImOmega);
	else return aux1.Re();
}
Double_t ComputeIntegralOmega(const Double_t ReOmega,const Double_t ImOmega,const Double_t xmin = 0.,const Double_t xmax = 12.)
{
  auto lambda = [](double *x, double *p){ return fOmega(x[0],p[0],p[1] ,false); }; // x[0] means "|dt|"
	TF1 fOmega("f1",lambda,0,100,2); // 2 == number of parameters
	fOmega.SetParameter(0,ReOmega);
	fOmega.SetParameter(1,ImOmega);
	return fOmega.Integral(xmin,xmax); // relative error on integral 10^-12
}

TH1D GetTheoryHistogram(const Double_t *par,const int nbin,const double DtMin,const double DtMax)
{
	if(DtMax < DtMin)cerr << "[GetTheoryHistogram] ERROR: DtMax < DtMin. Result makes no sense!" << endl;
	const double BinWidth = (DtMax - DtMin)/(double)nbin;
	TH1D h(Form("GetTheoryHistogram%d",g_GetTheoryHistogramCounter),"Theoretical function",nbin,DtMin,DtMax);
	g_GetTheoryHistogramCounter++;
	Double_t normalization = 0.;
	if(par[0] == 0.) // ----------------------------------------- zetaSL
	{
		normalization = ComputeIntegralSL(par[2],0.,12.);
		for(int ibin = 0; ibin < nbin; ibin++)
		{ // fill histogram with integral of theoretical model
			h.SetBinContent(ibin+1,ComputeIntegralSL(par[2],DtMin + ibin*BinWidth,DtMin + (ibin+1)*BinWidth)/normalization );
		}
	}
	else if(par[0] == 1.) // ----------------------------------------- zeta00
	{
		normalization = ComputeIntegral00(par[2],0.,12.);
		for(int ibin = 0.; ibin < nbin; ibin++)
		{// fill histogram h with integral of theoretical model over the bin (normalized)
			h.SetBinContent(ibin+1,ComputeIntegral00(par[2],DtMin + ibin*BinWidth,DtMin+ (ibin+1)*BinWidth)/normalization);
		}
	}
	else if(par[0] == 2.) // ----------------------------------------- gamma
	{
		normalization = ComputeIntegralGamma(par[2],0.,12.);
		for(int ibin = 0; ibin < nbin; ibin++)
		{ // fill histogram with integral of theoretical model
			h.SetBinContent(ibin+1,ComputeIntegralGamma(par[2],DtMin + ibin*BinWidth,DtMin + (ibin+1)*BinWidth)/normalization );
		}
	}else if(par[0] == 3.)// ----------------------------------------- omega
	{
		normalization = ComputeIntegralOmega(par[2],par[3],0.,12.);
		for(int ibin = 0; ibin < nbin; ibin++)
		{ // fill histogram with integral of theoretical model
			h.SetBinContent(ibin+1,ComputeIntegralOmega(par[2],par[3],DtMin + ibin*BinWidth,DtMin + (ibin+1)*BinWidth)/normalization );
		}
	}else{
		cerr << "[GetTheoryHistogram] ERROR: par[0] = " << par[0] << " -> no such model." << endl;
	}
	return h;
}
TH1D MultiplyHistogram(TH1D histo1,TH1D histo2,const int nbin,const double DtMin,const double DtMax, const bool debug=false)
{
	if(debug)
	{
		cerr << "[MultiplyHistogram] histo1 nbin / xmin / xmax = ";
		cerr << histo1.GetNbinsX() << " / " << histo1.GetBinLowEdge(1) << " / " << histo1.GetBinLowEdge(histo1.GetNbinsX()+1) << endl;
		cerr << "[MultiplyHistogram] histo2 nbin / xmin / xmax = ";
		cerr << histo2.GetNbinsX() << " / " << histo2.GetBinLowEdge(1) << " / " << histo2.GetBinLowEdge(histo2.GetNbinsX()+1) << endl;
		if(nbin != histo1.GetNbinsX() || nbin != histo2.GetNbinsX())cerr << "[MultiplyHistogram] ERROR: nbin mismatch! " << endl;
		if(DtMin != histo1.GetBinLowEdge(1) || DtMin != histo2.GetBinLowEdge(1))cerr << "[[MultiplyHistogram] ERROR: Low edge of bin 1 mismatch! " << endl;
		if( (DtMax-DtMin)/(double)nbin != histo1.GetBinLowEdge(2) - histo1.GetBinLowEdge(1) ) cerr <<"[MultiplyHistogram] ERROR: BinWidth mismatch histo 1! " << endl;
		if( (DtMax-DtMin)/(double)nbin != histo2.GetBinLowEdge(2) - histo1.GetBinLowEdge(1) ) cerr <<"[MultiplyHistogram] ERROR: BinWidth mismatch histo 2! " << endl;
	}
	TH1D hresult(Form("hresult%d",g_MultiplyHistogramCounter),"Histogram Product",nbin,DtMin,DtMax);
	g_MultiplyHistogramCounter++; // this is needed in order not to overwrite histograms if multiple calls happen
	// blindly performs multiplication bin-by-bin
	for(int ibin = 1; ibin <= nbin; ibin++)hresult.SetBinContent(ibin,histo1.GetBinContent(ibin)*histo2.GetBinContent(ibin));
	return hresult;
}
TH1D GetSmearedHistogram(TH1D htrue,TH2D SmMat,const bool debug = false)
{
	const int nbins = htrue.GetNbinsX();
	const double DtMin = htrue.GetBinLowEdge(1);
	const double DtMax = htrue.GetBinLowEdge(nbins + 1);
	if(debug){
		cerr << "[GetSmearedHistogram] DEBUG: nbins = " << nbins << endl;
		cerr << "[GetSmearedHistogram] DEBUG: DtMin = " << DtMin << endl;
		cerr << "[GetSmearedHistogram] DEBUG: DtMax = " << DtMax << endl;
		if(htrue.GetNbinsX() != SmMat.GetNbinsX() || htrue.GetNbinsX() != SmMat.GetNbinsY())cerr << "[GetSmearedHistogram] ERROR: bin mismatch!" << endl;
	}
	TH2D *SmMatNorm = (TH2D*) SmMat.Clone();
	Double_t tmp_norm_true;
	for(int ibintrue = 1; ibintrue <= nbins; ibintrue++)
	{ // normalize smearing matrix
		tmp_norm_true = 0.;
		for(int ibinreco = 1; ibinreco <= nbins; ibinreco++)tmp_norm_true += SmMatNorm->GetBinContent(ibintrue,ibinreco);
		for(int ibinreco = 1; ibinreco <= nbins; ibinreco++)SmMatNorm->SetBinContent(ibintrue,ibinreco,SmMatNorm->GetBinContent(ibintrue,ibinreco)/tmp_norm_true);
	}
	// print to stderr the whole normalized smearing matrix
	if(debug)for(int ibintrue = 1; ibintrue <= nbins; ibintrue++)for(int ibinreco = 1; ibinreco <= nbins; ibinreco++)cerr << "[GetSmearedHistogram::DEBUG] SmMatNorm.GetBinContent(" << ibintrue << ","<<ibinreco<<")="<< SmMatNorm->GetBinContent(ibintrue,ibinreco) << endl;
	
	// compute smeared histogram
	TH1D hSmerd(Form("hSmeared%d",g_GetSmearedHistogramCounter),"Smeared Histogram",nbins,DtMin,DtMax); // define histogram for smeared data
	g_GetSmearedHistogramCounter++; // this is needed in order not to overwrite histograms if multiple calls happen
	Double_t tmp_add;
	for(int ibinreco = 1; ibinreco <= nbins; ibinreco++)
	{
		for(int ibintrue = 1; ibintrue <= nbins; ibintrue++)
		{
			tmp_add = SmMatNorm->GetBinContent(ibintrue,ibinreco) * htrue.GetBinContent(ibintrue);
			hSmerd.SetBinContent(ibinreco, hSmerd.GetBinContent(ibinreco) + tmp_add );
		}
	}
	SmMatNorm->Delete();
	return hSmerd;
}

void ResetGlobalPhysicalConstants(Double_t TauS,Double_t TauL,Double_t DeltaM,Double_t modeta,Double_t phieta){
	g_TauS = TauS;
	g_TauL = TauL;
	g_DeltaM = DeltaM;
	g_modeta = modeta;
	g_phieta = phieta;
	cerr << "[ResetGlobalPhysicalConstants] RESET Physical constants global" << endl;
}

#endif
