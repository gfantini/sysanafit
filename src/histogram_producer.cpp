/*
 *  Questa macro unifica i tool di produzione istogrammi per la nuova era in cui diciamo addio al FORTRAN
 *  ATTENZIONE: gflib_v06 (caricata qui) elimina il taglio il pz che c'era in lorentz e in gflib_v04
 *  CUTVARNUM = 0...NCUT variabile di cui voglio variare il taglio
 *  CUTNSIGMA = <int>    numero di sigma di cui lo voglio variare
 *  NSIGMA4PI = <int>    numero di sigma (RES4PI, hard coded in inc/gflib_v08) di cui lo voglio variare (NON chiamare insieme al taglio in mpp)
 *  KREGEN    = <double> variazione alla rigenerazione: regen = default*(1.00+Kregen). Normal values are 0.00 | +/- 0.12 | +/- 0.24
 *  KRES      = <double> costante di variazione percentuale della risoluzione (0 = risoluzione nominale)
 */

#include "hmaker_varyh700_v03.cpp"
#include "root2dat_v06.cpp"

void MakeHistograms(const string folder = "./", const string name = "Histograms_cutvar00_p000.root",const int cutVarNum = 0,const int cutNsigma = 0,const int nsigma4pi=0,const double kRegen=0,const double kRes=0.)
{
  string filename=folder+name;
  // info printout
  cerr << "Welcome to [MakeHistograms]" << endl;
  cerr << "folder \t= " << folder << endl;
  cerr << "name \t= " << name << endl;
  cerr << "cutVarNum \t=" << cutVarNum << endl;
  cerr << "cutNsigma \t=" << cutNsigma << endl;
  cerr << "kRes \t=" << kRes << endl;
  /* production of basic histograms: adds to root file
   * TH1F h1000 data histo
   * TH1F h100 h200 efficiency histo
   * TH1F h700 smearing matrix histo
   * "cfg" TVectorD(8)
   * "4pi" TVectorD(1)
   * "kregen" TVectorD(1)
   * "kres" TVectorD(1)
   */
  cerr << "Calling hmaker to produce histograms." << endl;
  hmaker_varyh700_v03(filename.c_str(),cutVarNum,cutNsigma,nsigma4pi,kRegen,kRes);

  /* prouction of all the other needed histograms
   * TH1D
   *
   *
   * NOTICE: 1st argument of root2dat is not needed
   * WARNING: root2dat takes input information from some .dat files stored in /fit/dat/
   */
  cerr << "Calling root2dat_v06" << endl;
  root2dat_v06(nullptr,folder.c_str(),name.c_str(),0,"./fit/dat/regen_default.dat","./fit/dat/3000v.dat","./fit/dat/3000e.dat",cutVarNum,cutNsigma,nsigma4pi,kRegen,kRes);
  cerr << "Program terminated. Produced"<< folder+name << endl;
  TFile* pFile = new TFile(filename.c_str(),"READ");
  pFile->ls();
  pFile->Close();
}

void MakeCutvar(const string folder)
{
  int part,sigma;
  for(int cut=0; cut<8; cut++){// loop over variables (cuts)
    cerr << "[MakeCutvar] ********************" << endl;
    cerr << "[MakeCutvar] cut = " << cut << endl;
    
    part = 0;
    sigma = -3;
    while(sigma <= 3){// loop over variations of the single cut
      cerr << "[MakeCutVar] sigma = " << sigma << endl;
      MakeHistograms(folder,Form("Histograms_cutvar%02d_p%03d.root",cut,part),cut,sigma); // 4pi/regen/res == default
      part++;
      sigma++;
    }
  }
}

void Make4pivar(const string folder){
	
}

void MakeRegenvar(const string folder){
}

// vary resolution
void MakeResvar(const string folder)
{
  double kres[4] = {-0.015,-0.0075,0.0075,0.015};
  for(int part=0;part<4;part++){
    MakeHistograms(folder,Form("Histograms_resvar_p%03d.root",part),0,0,kres[part]);
  }
  
}
