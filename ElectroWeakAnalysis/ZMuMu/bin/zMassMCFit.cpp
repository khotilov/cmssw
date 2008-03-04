#include "PhysicsTools/Utilities/interface/Parameter.h"
#include "PhysicsTools/Utilities/interface/BreitWigner.h"
#include "PhysicsTools/Utilities/interface/Constant.h"
#include "PhysicsTools/Utilities/interface/Product.h"
#include "PhysicsTools/Utilities/interface/Sum.h"
#include "PhysicsTools/Utilities/interface/Gaussian.h"
#include "PhysicsTools/Utilities/interface/ZLineShape.h"
#include "PhysicsTools/Utilities/interface/HistoChiSquare.h"
#include "PhysicsTools/Utilities/interface/RootMinuit.h"
#include "PhysicsTools/Utilities/interface/RootFunctionAdapter.h"
#include "TFile.h"
#include "TH1.h"
#include "TF1.h"
#include "TMath.h"
#include "TCanvas.h"
#include "TROOT.h"
#include <boost/shared_ptr.hpp>
#include <boost/program_options.hpp>
using namespace boost;
namespace po = boost::program_options;

#include <iostream>
#include <algorithm> 
#include <iterator>
#include <string>
#include <vector>
using namespace std;
using namespace ::function;

// A helper function to simplify the main part.
template<class T>
ostream& operator<<(ostream& os, const vector<T>& v)
{
    copy(v.begin(), v.end(), ostream_iterator<T>(cout, " ")); 
    return os;
}

int main(int ac, char *av[]) { 
  try {
      double fMin, fMax;
      string ext;
      po::options_description desc("Allowed options");
      desc.add_options()
	("help", "produce help message")
	("include-path,I", po::value< vector<string> >(), 
	 "include path")
	("input-file", po::value< vector<string> >(), "input file")
	("min,m", po::value<double>(&fMin)->default_value(80), "minimum value for fit range")
	("max,M", po::value<double>(&fMax)->default_value(120), "maximum value for fit range")
	("breitwigner", "fit to a breit-wigner")
	("gauss", "fit to a gaussian")
	("bwint", "fit to a breit-wigner plus Z/photon interference term")
	("bwintgam", 
	 "fit to a breit-wigner plus Z/photon interference term and photon propagator")
	("output-file,O", po::value<string>(&ext)->default_value(".ps"), 
	 "output file format")
	;
      
      po::positional_options_description p;
      p.add("input-file", -1);
      
      po::variables_map vm;
      po::store(po::command_line_parser(ac, av).
	    options(desc).positional(p).run(), vm);
      po::notify(vm);
           
      if (vm.count("help")) {
	cout << "Usage: options_description [options]\n";
	cout << desc;
	return 0;
      }
      
      if (vm.count("include-path"))
        {
	  cout << "Include paths are: " 
	       << vm["include-path"].as< vector<string> >() << "\n";
        }
      
      vector<string> v_file;
      vector<TH1D*> v_ZMCMassHistos;
      vector<string> v_eps;
      
      if (vm.count("input-file"))
	{
	  cout << "Input files are: " 
	       << vm["input-file"].as< vector<string> >() << "\n";
	  v_file = vm["input-file"].as< vector<string> >();
	  for(vector<string>::const_iterator it = v_file.begin(); 
	      it != v_file.end(); ++it) { 
	     TFile * root_file = new TFile(it->c_str(),"read");
	     TDirectory *Histos = (TDirectory*) root_file->GetDirectory("ZHisto");
	     TDirectory *MCHistos = (TDirectory*) Histos->GetDirectory("ZMCHisto");
	     TH1D * zMass = (TH1D*) MCHistos->Get("ZMCMass");
	     zMass->GetXaxis()->SetTitle("Mass (GeV/c^{2})");
	     v_ZMCMassHistos.push_back(zMass);
	     gROOT->SetStyle("Plain");
	     string f_string = *it;
	     replace(f_string.begin(), f_string.end(), '.', '_');
	     string eps_string = f_string + ext;
	     v_eps.push_back(eps_string);
	     cout << ">>> histogram loaded\n";
	  }
	  cout << v_file.size() << ", " << v_ZMCMassHistos.size() << ", " << v_eps.size() << endl;
	  cout <<">>> Input files loaded\n";
	} 
      //PDG values for Z mass and width
      Parameter mass("Mass", 91.1876);
      Parameter gamma("Gamma", 2.4952);
      Parameter dmass("Mass Error", 0.0021);
      Parameter dgamma("Gamma Error", 0.0023);
      //Parameters for Z Line Shape
      Parameter f_gamma("Photon factor", 0);
      Parameter f_int("Interference factor", 0.001);
      Parameter df_gamma("Photon factor Error", 0);
      Parameter df_int("Interference factor Error", 0);
      //Parameters for fits with one gaussian
      Parameter yield("Yield", 1000);
      Parameter mean("Mean", 90);
      Parameter sigma("Sigma", 1.);
      Parameter dyield("Yield Error", 0);
      Parameter dmean("Mean Error", 0);
      Parameter dsigma("Sigma Error", 0);
      
      if (vm.count("breitwigner"))
	{
	  cout << "Fitting histograms in input file to a Breit-Wigner\n";
	  cout << ">>> set pars: " << endl;
	  cout << yield << " ; " << dyield << endl; 
	  cout << mass << " ; " << dmass << endl; 
	  cout << gamma << " ; " << dgamma << endl; 
	  for(size_t i = 0; i < v_ZMCMassHistos.size(); ++i) { 
	    cout << ">>> load histogram\n";
	    TH1D * zMass = v_ZMCMassHistos[i];
	    cout << ">>> histogram loaded\n";
	    BreitWigner bw(mass, gamma);
	    Constant c(yield);
	    typedef Product<Constant, BreitWigner> FitFunction;
	    FitFunction f = c * bw;
	    typedef fit::HistoChiSquare<FitFunction> ChiSquared;
	    ChiSquared chi2(f, zMass, fMin, fMax);
	    int fullBins = chi2.degreesOfFreedom();
	    cout << "N. deg. of freedom: " << fullBins << endl;
	    fit::RootMinuit<ChiSquared> minuit(3, chi2, true);
	    minuit.setParameter(0, yield, 10, 100, 100000);
	    minuit.setParameter(1, mass, .1, 70., 110);
	    minuit.setParameter(2, gamma, 1, 1, 10);
	    double amin = minuit.minimize();
	    cout << "fullBins = " << fullBins 
		 << "; free pars = " << minuit.getNumberOfFreeParameters() 
		 << endl;
	    unsigned int ndof = fullBins - minuit.getNumberOfFreeParameters();
	    cout << "Chi^2 = " << amin << "/" << ndof << " = " << amin/ndof 
		 << "; prob: " << TMath::Prob( amin, ndof )
		 << endl;
	    dyield = minuit.getParameterError(0);
	    cout << yield << " ; " << dyield << endl;
	    dmass = minuit.getParameterError(1);
	    cout << mass << " ; " << dmass << endl;
	    dgamma = minuit.getParameterError(2);
	    cout << gamma << " ; " << dgamma << endl;
	    TF1 fun = root::tf1("fun", f, fMin, fMax, yield, mass, gamma);
	    fun.SetParNames(yield.name().c_str(), mass.name().c_str(), gamma.name().c_str());
	    fun.SetLineColor(kRed);
	    TCanvas *canvas = new TCanvas("canvas");
	    zMass->Draw("e");
	    fun.Draw("same");	
	    string epsFilename = "ZMassMCFitBW_" + v_eps[i];
	    canvas->SaveAs(epsFilename.c_str());
	    canvas->SetLogy();
	    string epsLogFilename = "ZMassMCFitBW_Log_" + v_eps[i];
	    canvas->SaveAs(epsLogFilename.c_str());
	  }
	}
      
      if (vm.count("gauss"))
	{
	  cout << "Fitting histograms in input files to a Gaussian\n"; 
	  cout << ">>> set pars: " << endl;
	  cout << yield << " ; " << dyield << endl; 
	  cout << mean << " ; " << dmean << endl; 
	  cout << sigma << " ; " << dsigma << endl;
	  for(size_t i = 0; i < v_ZMCMassHistos.size(); ++i) { 
	    TH1D * zMass = v_ZMCMassHistos[i]; 
	    Gaussian gaus(mean, sigma);
	    Constant c(yield);
	    typedef Product<Constant, Gaussian> FitFunction;
	    FitFunction f = c * gaus;
	    typedef fit::HistoChiSquare<FitFunction> ChiSquared;
	    ChiSquared chi2(f, zMass, fMin, fMax);
	    int fullBins = chi2.degreesOfFreedom();
	    cout << "N. deg. of freedom: " << fullBins << endl;
	    fit::RootMinuit<ChiSquared> minuit(3, chi2, true);
	    minuit.setParameter(0, yield, 10, 100, 100000);
	    minuit.setParameter(1, mean, 0.001, 80, 100);
	    minuit.setParameter(2, sigma, 0.1, -5., 5.);
	    double amin = minuit.minimize();
	    cout << "fullBins = " << fullBins 
		 << "; free pars = " << minuit.getNumberOfFreeParameters() 
		 << endl;
	    unsigned int ndof = fullBins - minuit.getNumberOfFreeParameters();
	    cout << "Chi^2 = " << amin << "/" << ndof << " = " << amin/ndof 
		 << "; prob: " << TMath::Prob( amin, ndof )
		 << endl;
	    dyield = minuit.getParameterError(0);
	    cout << yield << " ; " << dyield << endl;
	    dmean = minuit.getParameterError(1);
	    cout << mean << " ; " << dmean << endl;
	    dsigma = minuit.getParameterError(2);
	    cout << sigma << " ; " << dsigma << endl;
	    TF1 fun = root::tf1("fun", f, fMin, fMax, yield, mean, sigma);
	    fun.SetParNames(yield.name().c_str(), mean.name().c_str(), sigma.name().c_str());
	    fun.SetLineColor(kRed);
	    TCanvas *canvas = new TCanvas("canvas");
	    zMass->Draw("e");
	    fun.Draw("same");	
	    string epsFilename = "ZMassMCFitG_" + v_eps[i];
	    canvas->SaveAs(epsFilename.c_str());
	    canvas->SetLogy();
	    string epsLogFilename = "ZMassMCFitG_Log_" + v_eps[i];
	    canvas->SaveAs(epsLogFilename.c_str());
	  }
	}
      
      if (vm.count("bwint"))
	{
	  cout << "Fitting histograms in input files to the Breit-Wigner plus Z/photon interference term\n";
	  cout << ">>> set pars: " << endl;
	  cout << yield << " ; " << dyield << endl; 
	  cout << mass << " ; " << dmass << endl; 
	  cout << gamma << " ; " << dgamma << endl; 
	  cout << f_gamma << " ; " << df_gamma << endl; 
	  cout << f_int << " ; " << df_int << endl; 
	  for(size_t i = 0; i < v_ZMCMassHistos.size(); ++i) { 
	    TH1D * zMass = v_ZMCMassHistos[i]; 
	    ZLineShape zls(mass, gamma, f_gamma, f_int);
	    Constant c(yield);
	    typedef Product<Constant, ZLineShape> FitFunction;
	    FitFunction f = c * zls;
	    cout << "set functions" << endl;
	    vector<shared_ptr<double> > pars;
	    pars.push_back(yield.ptr());
	    pars.push_back(mass.ptr());
	    pars.push_back(gamma.ptr());
	    pars.push_back(f_gamma.ptr());
	    pars.push_back(f_int.ptr());
	    typedef fit::HistoChiSquare<FitFunction> ChiSquared;
	    ChiSquared chi2(f, zMass, fMin, fMax);
	    int fullBins = chi2.degreesOfFreedom();
	    cout << "N. deg. of freedom: " << fullBins << endl;
	    fit::RootMinuit<ChiSquared> minuit(5, chi2, true);
	    minuit.setParameter(0, yield, 10, 100, 100000);
	    minuit.setParameter(1, mass, .1, 70., 110);
	    minuit.setParameter(2, gamma, 1, 1, 10);
	    minuit.setParameter(3, f_gamma, 0.1, -100, 1000);
	    minuit.fixParameter(3);
	    minuit.setParameter(4, f_int, .0001, -1000000, 1000000);
	    double amin = minuit.minimize();
	    cout << "fullBins = " << fullBins 
		 << "; free pars = " << minuit.getNumberOfFreeParameters() 
		 << endl;
	    unsigned int ndof = fullBins - minuit.getNumberOfFreeParameters();
	    cout << "Chi^2 = " << amin << "/" << ndof << " = " << amin/ndof 
		 << "; prob: " << TMath::Prob( amin, ndof )
		 << endl;
	    dyield = minuit.getParameterError(0);
	    cout << yield << " ; " << dyield << endl;
	    dmass = minuit.getParameterError(1);
	    cout << mass << " ; " << dmass << endl;
	    dgamma = minuit.getParameterError(2);
	    cout << gamma << " ; " << dgamma << endl;
	    df_gamma = minuit.getParameterError(3);
	    cout << f_gamma << " ; " << df_gamma << endl;
	    df_int = minuit.getParameterError(4);
	    cout << f_int << " ; " << df_int << endl;
	    TF1 fun = root::tf1("fun", f, fMin, fMax, pars); 
	    fun.SetParNames(yield.name().c_str(), mass.name().c_str(), gamma.name().c_str(), 
	                    f_gamma.name().c_str(), f_int.name().c_str());
	    fun.SetLineColor(kRed); 
	    TCanvas *canvas = new TCanvas("canvas");
	    zMass->Draw("e");
	    fun.Draw("same");
	    string epsFilename = "ZMassMCFitBwIn_" + v_eps[i];
	    canvas->SaveAs(epsFilename.c_str());
	    canvas->SetLogy();
	    string epsLogFilename = "ZMassMCFitBwIn_Log_" + v_eps[i];
	    canvas->SaveAs(epsLogFilename.c_str());
	  }
	}
      
      if (vm.count("bwintgam"))
	{
	  cout << "Fitting histograms in input files to the Breit-Wigner plus Z/photon interference term and photon propagator\n";
	  cout << ">>> set pars: " << endl;
	  cout << yield << " ; " << dyield << endl; 
	  cout << mass << " ; " << dmass << endl; 
	  cout << gamma << " ; " << dgamma << endl; 
	  cout << f_gamma << " ; " << df_gamma << endl; 
	  cout << f_int << " ; " << df_int << endl; 
	  for(size_t i = 0; i < v_ZMCMassHistos.size(); ++i) { 
	    TH1D * zMass = v_ZMCMassHistos[i]; 
	    ZLineShape zls(mass, gamma, f_gamma, f_int);
	    Constant c(yield);
	    typedef Product<Constant, ZLineShape> FitFunction;
	    FitFunction f = c * zls;
	    cout << "set functions" << endl;
	    vector<shared_ptr<double> > pars;
	    pars.push_back(yield.ptr());
	    pars.push_back(mass.ptr());
	    pars.push_back(gamma.ptr());
	    pars.push_back(f_gamma.ptr());
	    pars.push_back(f_int.ptr());
	    typedef fit::HistoChiSquare<FitFunction> ChiSquared;
	    ChiSquared chi2(f, zMass, fMin, fMax);
	    int fullBins = chi2.degreesOfFreedom();
	    cout << "N. deg. of freedom: " << fullBins << endl;
	    fit::RootMinuit<ChiSquared> minuit(5, chi2, true);
	    minuit.setParameter(0, yield, 10, 100, 100000);
	    minuit.setParameter(1, mass, .1, 70., 110);
	    minuit.setParameter(2, gamma, 1, 1, 10);
	    minuit.setParameter(3, f_gamma, 0.1, -100, 1000);
	    minuit.setParameter(4, f_int, .0001, -1000000, 1000000);
	    double amin = minuit.minimize();
	    cout << "fullBins = " << fullBins 
		 << "; free pars = " << minuit.getNumberOfFreeParameters() 
		 << endl;
	    unsigned int ndof = fullBins - minuit.getNumberOfFreeParameters();
	    cout << "Chi^2 = " << amin << "/" << ndof << " = " << amin/ndof 
		 << "; prob: " << TMath::Prob( amin, ndof )
		 << endl;
	    dyield = minuit.getParameterError(0);
	    cout << yield << " ; " << dyield << endl;
	    dmass = minuit.getParameterError(1);
	    cout << mass << " ; " << dmass << endl;
	    dgamma = minuit.getParameterError(2);
	    cout << gamma << " ; " << dgamma << endl;
	    df_gamma = minuit.getParameterError(3);
	    cout << f_gamma << " ; " << df_gamma << endl;
	    df_int = minuit.getParameterError(4);
	    cout << f_int << " ; " << df_int << endl;
	    TF1 fun = root::tf1("fun", f, fMin, fMax, pars); 
	    fun.SetParNames(yield.name().c_str(), mass.name().c_str(), gamma.name().c_str(), 
	                    f_gamma.name().c_str(), f_int.name().c_str());
	    fun.SetLineColor(kRed); 
	    TCanvas *canvas = new TCanvas("canvas");
	    zMass->Draw("e");
	    fun.Draw("same");
	    string epsFilename = "ZMassMCFitBwInGa_" + v_eps[i];
	    canvas->SaveAs(epsFilename.c_str());
	    canvas->SetLogy();
	    string epsLogFilename = "ZMassMCFitBwInGa_Log_" + v_eps[i];
	    canvas->SaveAs(epsLogFilename.c_str());
	  }
	}
      
      cout << "It works!\n";
  }
  catch(exception& e) {
    cerr << "error: " << e.what() << "\n";
    return 1;
  }
  catch(...) {
    cerr << "Exception of unknown type!\n";
  }
  return 0;
}
