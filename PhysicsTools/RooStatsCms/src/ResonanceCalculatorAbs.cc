#include "PhysicsTools/RooStatsCms/interface/ResonanceCalculatorAbs.hh"
#include "PhysicsTools/RooStatsCms/interface/FeldmanCousinsBinomialInterval.h"

#include <limits>

#include "TH1.h"
#include "TMath.h"
#include "TLatex.h"
#include "TFile.h"
#include "TTree.h"

#include "RooWorkspace.h"
#include "RooAbsPdf.h"
#include "RooRealVar.h"
#include "RooAbsData.h"
#include "RooDataHist.h"
#include "RooDataSet.h"
#include "RooPlot.h"
#include "RooFitResult.h"
#include "RooRandom.h"
#include "RooStats/SimpleLikelihoodRatioTestStat.h"
#include "RooStats/ModelConfig.h"

////////////////////////////////////////////////////////////////////////////////
// implementation
////////////////////////////////////////////////////////////////////////////////

int ResonanceCalculatorAbs::printLevel_=1;

ResonanceCalculatorAbs::ResonanceCalculatorAbs()
  : numPEs_(100), nBinsToDraw_(100), whichTestStatistic_(0), searchStepSize_(1.0), fixedSearchStepSize_(-999), randomSeed_(1), fitStrategy_(2), controlMin_(1.), controlMax_(0.), fixNbkg_(false), dataIntegral_(-1), ws_(0)
{
  // default constructor
}

ResonanceCalculatorAbs::ResonanceCalculatorAbs(const ResonanceCalculatorAbs& other)
  : numPEs_(other.numPEs_), nBinsToDraw_(other.nBinsToDraw_), whichTestStatistic_(other.whichTestStatistic_),
    searchStepSize_(other.searchStepSize_), fixedSearchStepSize_(other.fixedSearchStepSize_), randomSeed_(other.randomSeed_), fitStrategy_(other.fitStrategy_), controlMin_(other.controlMin_), controlMax_(other.controlMax_), fixNbkg_(other.fixNbkg_),
    dataBinning_(other.dataBinning_), dataIntegral_(other.dataIntegral_)
{
  // copy constructor

  ws_ = new RooWorkspace(*other.ws_);
  data_ = ws_->data(other.data_->GetName());
  obs_ = ws_->var(other.obs_->GetName());
  model_ = ws_->pdf(other.model_->GetName());
  signal_ = ws_->pdf(other.signal_->GetName());
  background_ = ws_->pdf(other.background_->GetName());
  sigmass_ = ws_->var(other.sigmass_->GetName());
  sigwidth_ = ws_->function(other.sigwidth_->GetName());
  nsig_ = ws_->var(other.nsig_->GetName());  
  nbkg_ = ws_->var(other.nbkg_->GetName());

  findSignalAndBackgroundParams();

}

ResonanceCalculatorAbs::~ResonanceCalculatorAbs()
{
  // destructor

  if(ws_) delete ws_;
  if(sigparams_) delete sigparams_;
  if(bkgparams_) delete bkgparams_;
  
}

double ResonanceCalculatorAbs::calculate(const char* rootfilename, bool fastAlgo)
{
  double bumpMass, bumpTestStat;
  std::vector<std::pair<double, double> > teststatsfloat, teststatsfixed;
  return calculate(rootfilename, teststatsfloat, teststatsfixed, bumpMass, bumpTestStat, fastAlgo);
}

double ResonanceCalculatorAbs::calculate(const char* rootfilename, std::vector<std::pair<double, double> >& teststatsfloat_, std::vector<std::pair<double, double> >& teststatsfixed_, double& bumpMass_, double& bumpTestStat_, bool fastAlgo)
{
  // find the bump, and put the results and diagnostic plots into rootfilename
  
  // clear out the vector, just in case
  teststatsfloat_.clear();
  teststatsfixed_.clear();
  
  // set the random seed
  RooRandom::randomGenerator()->SetSeed(randomSeed_);
  
  // open the file and create a directory based on the random seed
  TFile* rootfile=new TFile(rootfilename, "RECREATE");
  rootfile->cd();
  TString dirname="seed";
  dirname+=randomSeed_;
  
  // setup a tree to store the test statistic values
  TTree* tree=new TTree("rescalctree","Resonance Calculator Tree");
  int randomSeed;
  float bumpMass, bumpTestStat;
  int numPEs;
  float *testStatFixed, *testStatFloat, *testStatMass;
  testStatFixed = new float[numPEs_];
  testStatFloat = new float[numPEs_];
  testStatMass  = new float[numPEs_];
  tree->Branch("randomSeed",&randomSeed,"randomSeed/I");                 // random seed
  tree->Branch("bumpMass",&bumpMass,"bumpMass/F");                       // mass of the most significant bump in the data
  tree->Branch("bumpTestStat",&bumpTestStat,"bumpTestStat/F");           // test statistic value of the most sigificant bump in the data
  tree->Branch("numPEs",&numPEs,"numPEs/I");                             // number of pseudoexperiments run
  tree->Branch("testStatFixed",testStatFixed,"testStatFixed[numPEs]/F"); // test statistic value for a given PE and with a mass fixed at bumpMass
  tree->Branch("testStatFloat",testStatFloat,"testStatFloat[numPEs]/F"); // test statistic value for a given PE with a floating mass
  tree->Branch("testStatMass",testStatMass,"testStatMass[numPEs]/F");    // mass of the most significant bump for a given PE
  randomSeed=randomSeed_;
  numPEs=numPEs_;

  // find a bump
  rootfile->mkdir(dirname+"_data");
  rootfile->cd(dirname+"_data");
  if(controlMin_>controlMax_) {
    if(fastAlgo) scanFastestForBump("data");
    else         scanForBump("data");
  } else {
    if(fastAlgo) scanFastestForBumpWithControl("data");
    else         scanForBumpWithControl("data");
  }

  if(nsig_->getVal()==0) {
    if(printLevel_>=0) std::cout << "There was no bump found." << std::endl;
    numPEs_=0; // don't run any PEs, that would be silly
  }
  
  bumpMass=sigmass_->getVal();
  bumpTestStat=evaluateTestStatistic();
  bumpTestStat_=bumpTestStat;
  bumpMass_=bumpMass;
  
  // keep the original parameters and for later reference
  RooArgList originalBackgroundParams;
  originalBackgroundParams.addClone(*bkgparams_);
  RooAbsData* originalData = data_;

  //  make PEs based on the background
  for(int i=0; i<numPEs_; i++) {
    // copy the original background parameters from the data fit
    copyValuesToBkgParams(&originalBackgroundParams);
    
    // make pseudodata
    if(dynamic_cast<RooDataHist*>(originalData))
      data_ = generateBinned(background_, dynamic_cast<RooDataHist*>(originalData), dataIntegral_);
    else if(dynamic_cast<RooDataSet*>(originalData))
      data_ = generateUnbinned(background_, dataIntegral_);
    else
      assert(1);
    
    // set the label and descend into the directory
    TString label("PE");
    label+=i;
    rootfile->cd("/");
    rootfile->mkdir(dirname+"_"+label);
    rootfile->cd(dirname+"_"+label);

    // look for a bump in the pseudodata with the same mass as the bump in the data
    sigmass_->setVal(bumpMass);
    doBkgOnlyExcludeWindowFit(label+"fixedMassToData");
    doSigOnlyFixMassFit(label+"fixedMassToData");
    testStatFixed[i]=evaluateTestStatistic();

    // scan for a bump in the pseudodata
    if(controlMin_>controlMax_) {
      if(fastAlgo) scanFastestForBump(label);
      else         scanForBump(label);
    } else {
      if(fastAlgo) scanFastestForBumpWithControl(label);
      else         scanForBumpWithControl(label);
    }
    testStatFloat[i]=evaluateTestStatistic();
    testStatMass[i]=sigmass_->getVal();

    teststatsfloat_.push_back(std::pair<double, double>(testStatFloat[i], 1.0));
    teststatsfixed_.push_back(std::pair<double, double>(testStatFixed[i], 1.0));
    
    // delete the pseudodata when we're through with it
    delete data_;
  }
  
  // fill the tree
  tree->Fill();

  // reset the data back to the original
  data_=originalData;
  
  // write the workspace and close the file
  rootfile->cd(dirname+"_data");
  ws_->Write();
  rootfile->cd("/");
  tree->Write();
  rootfile->Close();
  
  delete[] testStatFloat;
  delete[] testStatFixed;
  delete[] testStatMass;

  // if we aren't running any PEs, just return the bump test statistic now
  if(numPEs_<=0) return bumpTestStat_;

  return getZScore(teststatsfloat_, bumpTestStat_);
}


void ResonanceCalculatorAbs::setBinnedData(TH1* dataHist, int minBin, int maxBin)
{
  // Creates a dataset based on the histogram and imports it into the workspace.
  // The minimum and maximum values of the observable are set here, being determined by
  // the minBin and maxBin
  
  if(!dataHist) {
    std::cerr << "ResonanceCalculatorAbs::setBinnedData(): Bad Histogram input" << std::endl;
    return;
  }

  // use the name of the histogram to give the dataset a name
  TString dataname = dataHist->GetName();

  // if there is already a dataset around, get a new name
  if(data_) dataname = TString(data_->GetName())+"_";

  // set the range of the observable to the minimum and maximum values found in the histogram
  obs_->setRange(dataHist->GetXaxis()->GetBinLowEdge(minBin), dataHist->GetXaxis()->GetBinUpEdge(maxBin));
  obs_->setVal(dataHist->GetXaxis()->GetBinLowEdge(minBin));

  // store the data binning (amazingly, this is needed to throw PEs)
  if(dataHist->GetXaxis()->GetXbins()->GetArray()) {
    std::cerr << "ResonanceCalculatorAbs::setBinnedData(): You cannot use variable binning!!" << std::endl;
    assert(0);
    //    RooBinning binning(dataHist->GetNbinsX(), dataHist->GetXaxis()->GetXbins()->GetArray());
  } else {
    RooBinning binning(dataHist->GetNbinsX(), dataHist->GetXaxis()->GetBinLowEdge(minBin), dataHist->GetXaxis()->GetBinUpEdge(maxBin));
    dataBinning_ = binning;
  }

  // create the datahist
  data_ = new RooDataHist(dataname, "binned data", RooArgSet(*obs_), RooFit::Import(*dataHist, kFALSE));
  ws_->import(*data_);

  // delete the dataset now that we've imported it into the workspace
  delete data_;

  // get back a pointer to the imported dataset
  data_ = dynamic_cast<RooDataHist*>(ws_->data(dataname));  

  // store the integral of the data
  dataIntegral_ = data_->sumEntries();
  nbkg_->setVal(dataIntegral_);

  // determine the default range for the minimum and maximum mass to search for
  findMinMaxMass();

  return;
}

void ResonanceCalculatorAbs::setBinnedData(TH1* dataHist)
{
  // Creates a dataset based on the histogram and imports it into the workspace.
  // The minimum and maximum values of the observable are set here, being determined by
  // the histogram boundaries.

  if(!dataHist) {
    std::cerr << "ResonanceCalculatorAbs::setBinnedData(): Bad Histogram input" << std::endl;
    return;
  }

  // use the name of the histogram to give the dataset a name
  TString dataname = dataHist->GetName();

  // if there is already a dataset around, get a new name
  if(data_) dataname = TString(data_->GetName())+"_";

  // set the range of the observable to the minimum and maximum values found in the histogram
  obs_->setRange(dataHist->GetXaxis()->GetXmin(), dataHist->GetXaxis()->GetXmax());
  obs_->setVal(dataHist->GetXaxis()->GetXmin());

  // store the data binning (amazingly, this is needed to throw PEs)
  if(dataHist->GetXaxis()->GetXbins()->GetArray()) {
    std::cerr << "ResonanceCalculatorAbs::setBinnedData(): You cannot use variable binning!!" << std::endl;
    assert(0);
    //    RooBinning binning(dataHist->GetNbinsX(), dataHist->GetXaxis()->GetXbins()->GetArray());
  } else {
    RooBinning binning(dataHist->GetNbinsX(), dataHist->GetXaxis()->GetXmin(), dataHist->GetXaxis()->GetXmax());
    dataBinning_ = binning;
  }

  // create the datahist
  data_ = new RooDataHist(dataname, "binned data", RooArgSet(*obs_), RooFit::Import(*dataHist, kFALSE));
  ws_->import(*data_);

  // delete the dataset now that we've imported it into the workspace
  delete data_;

  // get back a pointer to the imported dataset
  data_ = dynamic_cast<RooDataHist*>(ws_->data(dataname));  

  // store the integral of the data
  dataIntegral_ = data_->sumEntries();
  nbkg_->setVal(dataIntegral_);

  // determine the default range for the minimum and maximum mass to search for
  findMinMaxMass();

  return;
}

void ResonanceCalculatorAbs::setBinnedData(const char* filename, int nbins, double minx, double maxx)
{
  // input is a text file with a list of data points, but the data gets binned
  // the number of bins, as well as the minimum and maximum bin values must be specified

  // open the file
  ifstream is;
  is.open(filename, ifstream::in);
  if(!is) {
    std::cerr << "ResonanceCalculatorAbs::setBinnedData(): Could not open file " << filename << std::endl;
    return;
  }

  // create a temporary histogram
  TH1D* dataHist=new TH1D("dataHist","data histogram",nbins, minx, maxx);

  // fill the histogram
  while(is.good()) {
    double temp;
    is >> temp;
    if(!is.good()) break;
    dataHist->Fill(temp);
  }
  is.close();

  // use the name of the histogram to give the dataset a name
  TString dataname = dataHist->GetName();

  // if there is already a dataset around, get a new name
  if(data_) dataname = TString(data_->GetName())+"_";

  // set the range of the observable to the minimum and maximum values found in the histogram
  obs_->setRange(minx, maxx);
  obs_->setVal(minx);

  // store the data binning (amazingly, this is needed to throw PEs)
  if(dataHist->GetXaxis()->GetXbins()->GetArray()) {
    RooBinning binning(dataHist->GetNbinsX(), dataHist->GetXaxis()->GetXbins()->GetArray());
    dataBinning_ = binning;
  } else {
    RooBinning binning(dataHist->GetNbinsX(), dataHist->GetXaxis()->GetXmin(), dataHist->GetXaxis()->GetXmax());
    dataBinning_ = binning;
  }

  // create the datahist
  data_ = new RooDataHist(dataname, "binned data", RooArgSet(*obs_), RooFit::Import(*dataHist, kFALSE));
  ws_->import(*data_);

  // delete the dataset and histogram now that we've imported it into the workspace
  delete data_;
  delete dataHist;

  // get back a pointer to the imported dataset
  data_ = dynamic_cast<RooDataHist*>(ws_->data(dataname));  

  // store the integral of the data
  dataIntegral_ = data_->sumEntries();
  nbkg_->setVal(dataIntegral_);

  // determine the default range for the minimum and maximum mass to search for
  findMinMaxMass();
  
  return;
}

void ResonanceCalculatorAbs::setUnbinnedData(const char* filename, double minx, double maxx)
{
  // input is a text file with a list of data points.  The data remains unbinned, but
  // the minimum and maximum bin values must be specified

  // open the file
  ifstream is;
  is.open(filename, ifstream::in);
  if(!is) {
    std::cerr << "ResonanceCalculatorAbs::setUnbinnedData(): Could not open file " << filename << std::endl;
    return;
  }

  // fill the temporary dataset (recording the min/max data value)
  std::vector<double> tempDataset;
  while(is.good()) {
    double temp;
    is >> temp;
    if(!is.good()) break;
    if(temp>=minx && temp<=maxx) tempDataset.push_back(temp);
  }
  is.close();

  // set the range of the observable
  obs_->setRange(minx, maxx);
  obs_->setVal(minx);

  // create the dataset
  TString dataname = "dataSet";
  if(data_) dataname = TString(data_->GetName())+"_";
  data_ = new RooDataSet(dataname, "unbinned data", RooArgSet(*obs_));
  for(unsigned int i=0; i<tempDataset.size(); i++) {
    obs_->setVal(tempDataset[i]);
    data_->add(RooArgSet(*obs_));
  }

  // import the dataset
  ws_->import(*data_);

  // store the integral of the data
  dataIntegral_ = data_->sumEntries();
  nbkg_->setVal(dataIntegral_);

  // delete the dataset and histogram now that we've imported it into the workspace
  delete data_;

  // get back a pointer to the imported dataset
  data_ = ws_->data(dataname);

  // determine the default range for the minimum and maximum mass to search for
  findMinMaxMass();

  return;
}

void ResonanceCalculatorAbs::setupWorkspaceViaFactory(const char* sigexpr, const char* bkgexpr, const char* widthexpr)
{
  // setup internal workspace via a factory mechanism

  // create the workspace
  ws_ = new RooWorkspace("ws");

  // setup internal observable
  ws_->factory("obs[0, 0, 100]");
  obs_ = ws_->var("obs");
  obs_->setConstant(false);

  // setup signal mass parameter
  ws_->factory("signalmass[1000.]");
  sigmass_ = ws_->var("signalmass");

  // setup background pdf
  ws_->factory(bkgexpr);
  background_ = ws_->pdf("background");
  if(!background_) {
    std::cerr << "ResonanceCalculatorAbs::setupWorkspaceViaFactory(): Failed to import the background pdf" << std::endl;
    assert(0);
  }

  // setup signal width parameter
  ws_->factory(widthexpr);
  sigwidth_ = ws_->var("signalwidth");
  if(!sigwidth_) sigwidth_ = ws_->function("signalwidth");
  if(!sigwidth_) {
    std::cerr << "ResonanceCalculatorAbs::setupWorkspaceViaFactory(): Failed to import the signal width expression" << std::endl;
    assert(0);
  }

  // setup background pdf
  ws_->factory(sigexpr);
  signal_ = ws_->pdf("signal");
  if(!signal_) {
    std::cerr << "ResonanceCalculatorAbs::setupWorkspaceViaFactory(): Failed to import the signal pdf" << std::endl;
    assert(0);
  }

  // setup model
  ws_->factory(TString("SUM::model(nsig[1]*")+signal_->GetName()+", nbkg[1]*"+background_->GetName()+")");
  model_=ws_->pdf("model");
  nsig_=ws_->var("nsig");
  nbkg_=ws_->var("nbkg");
  nsig_->setRange(0.0, numeric_limits<double>::max());
  nbkg_->setRange(0.0, numeric_limits<double>::max());

  // find the signal and background params and set them to the pointers, accordingly
  findSignalAndBackgroundParams();

  // set other defaults
  data_ = 0;

  if(printLevel_>=1)
    ws_->Print();

  return;
}

void ResonanceCalculatorAbs::setupWorkspace(void)
{
  // setup the internal workspace (called from the constructors)

  // create the workspace
  ws_ = new RooWorkspace("ws");

  // setup internal observable
  ws_->factory("obs[0, 0, 100]");
  obs_ = ws_->var("obs");
  obs_->setConstant(false);

  // setup signal mass parameter
  ws_->factory("signalmass[1000.]");
  sigmass_ = ws_->var("signalmass");

  // setup background pdf
  RooAbsPdf* backgroundPdf=setupBackgroundPdf();
  if(ws_->import(*backgroundPdf)) {
    std::cerr << "ResonanceCalculatorAbs::setupWorkspace(): Failed to import the background pdf named " << backgroundPdf->GetName() << std::endl;
    assert(0);
  }
  background_ = ws_->pdf(backgroundPdf->GetName());
  delete backgroundPdf;

  // setup signal width parameter
  RooAbsReal* signalWidth=setupSignalWidth();
  if(ws_->import(*signalWidth)) {
    std::cerr << "ResonanceCalculatorAbs::setupWorkspace(): Failed to import the signal width named " << signalWidth->GetName() << std::endl;
    assert(0);
  }
  sigwidth_ = ws_->var(signalWidth->GetName());
  if(!sigwidth_) sigwidth_ = ws_->function(signalWidth->GetName());
  delete signalWidth;

  // setup signal pdf
  RooAbsPdf* signalPdf=setupSignalPdf();
  if(ws_->import(*signalPdf)) {
    std::cerr << "ResonanceCalculatorAbs::setupWorkspace(): Failed to import the signal pdf named " << signalPdf->GetName() << std::endl;
    assert(0);
  }
  signal_ = ws_->pdf(signalPdf->GetName());
  delete signalPdf;

  // setup model
  ws_->factory(TString("SUM::model(nsig[1]*")+signal_->GetName()+", nbkg[1]*"+background_->GetName()+")");
  model_=ws_->pdf("model");
  nsig_=ws_->var("nsig");
  nbkg_=ws_->var("nbkg");
  nsig_->setRange(0.0, numeric_limits<double>::max());
  nbkg_->setRange(0.0, numeric_limits<double>::max());

  // find the signal and background params and set them to the pointers, accordingly
  findSignalAndBackgroundParams();

  // set other defaults
  data_ = 0;

  if(printLevel_>=1)
    ws_->Print();

  return;
}

void ResonanceCalculatorAbs::scanFastestForBump(const char* label)
{
  // do B-only fit, first
  doBkgOnlyFit(label);
  
  // setup a SCAN for the best bump
  nbkg_->setConstant(true);
  setBkgParamsConst(true);
  nsig_->setConstant(false);
  nsig_->setRange(0.01,nbkg_->getVal());
  nsig_->setVal(nbkg_->getVal()/2.0);
  sigmass_->setConstant(false);
  setSigParamsConst(false);

  // SCAN
  obs_->setRange("full",obs_->getMin(), obs_->getMax());
    RooFitResult* fit=model_->fitTo(*data_, RooFit::Save(kTRUE), RooFit::Extended(kTRUE), RooFit::Strategy(0), RooFit::PrintLevel(printLevel_), RooFit::Minimizer("GSLSimAn"));
  if(printLevel_>=2) fit->Print();
  drawFitResult(TString(label)+"_siman","full");

  // do things more carefully
  doBkgOnlyExcludeWindowFit(label);
  double minmassrange=std::max(sigmass_->getVal()-1.1*getStepSize(), sigmass_->getMin());
  double maxmassrange=std::min(sigmass_->getVal()+1.1*getStepSize(), sigmass_->getMax());
  doSigOnlyFloatMassFit(label, minmassrange, maxmassrange);
  
  return;
}

void ResonanceCalculatorAbs::scanForBump(const char* label)
{
  // do B-only fit, first
  doBkgOnlyFit(label);

  // scan the mass range starting at the minimum
  sigmass_->setVal(sigmass_->getMin());
  double currmass=sigmass_->getVal();
  double nextmass=currmass+getStepSize();
  double bestSig=0.0, bestSigMass=0.0;
  RooArgList bestSignalParams;
  int counter=0;
  while(currmass<sigmass_->getMax()) {
    ++counter;

    // do a fit
    double minmassrange=std::max(currmass-1.1*getStepSize(), sigmass_->getMin());
    double maxmassrange=std::min(currmass+1.1*getStepSize(), sigmass_->getMax());
    char newlabel[1000];
    sprintf(newlabel, "%sRegion%d",label,counter);
    if(doBkgPlusSigFloatMassFit(newlabel, minmassrange, maxmassrange)->status()!=0) continue;

    // if there is a bump, fit again, more carefully
    if(nsig_->getVal()>0.0) {

      // fit for bkg-only outside the mass region
      if(doBkgOnlyExcludeWindowFit(newlabel)->status()!=0) continue;

      // fit for the signal, with fixed background
      if(doSigOnlyFloatMassFit(newlabel, minmassrange, maxmassrange)->status()!=0) continue;
      //      if(doBkgPlusSigFloatMassFit(newlabel, minmassrange, maxmassrange)->status()!=0) continue;
    }

    // see if the bump is significant
    double significance=evaluateTestStatistic();
    if(significance>bestSig) {
      bestSignalParams.removeAll();
      bestSignalParams.addClone(*sigparams_);
      bestSig=significance;
      bestSigMass=sigmass_->getVal();
    }

    if(printLevel_>=1) {
      std::cout << "ResonanceCalculatorAbs::scanForBump(): "
		<< "In the range [" << minmassrange << ", " << maxmassrange << "]"
		<< " we found a bump at m=" << sigmass_->getVal() << " with significance="
		<< significance << ".  The best significance so far is "
		<< bestSig << " at m=" << bestSigMass << std::endl;
    }

    // step
    sigmass_->setVal(nextmass);
    currmass=sigmass_->getVal();
    nextmass=currmass+getStepSize();
  }

  if(bestSig<=0.0) {
    if(printLevel_>=0) std::cout << "ResonanceCalculatorAbs::scanForBump(): "
				 << "No significant excess found!" << std::endl;
    return;
  }

  // redo fit
  copyValuesToSigParams(&bestSignalParams);
  sigmass_->setVal(bestSigMass);
  doBkgOnlyExcludeWindowFit(label);
  doSigOnlyFixMassFit(label);


  return;
}

void ResonanceCalculatorAbs::scanFastestForBumpWithControl(const char* label)
{
  // check for an overlap between the signal and control regions
  if((sigmass_->getMin()>controlMin_ && sigmass_->getMin()<controlMax_) ||
     (sigmass_->getMax()>controlMin_ && sigmass_->getMax()<controlMax_) ||
     (sigmass_->getMin()<controlMin_ && sigmass_->getMax()>controlMax_)) {
    if(printLevel_>=0) std::cout << "The signal and control regions overlap.  Are you sure you want to do that?" << std::endl;
  }
  
  // do B-only fit, first
  double oldMin=obs_->getMin();
  double oldMax=obs_->getMax();
  obs_->setRange(controlMin_, controlMax_);
  doBkgOnlyFit(label);
  obs_->setRange(oldMin, oldMax);
  
  // setup a SCAN for the best bump
  nbkg_->setConstant(true);
  setBkgParamsConst(true);
  nsig_->setConstant(false);
  nsig_->setRange(0.01,nbkg_->getVal());
  nsig_->setVal(nbkg_->getVal()/2.0);
  sigmass_->setConstant(false);
  setSigParamsConst(false);

  // SCAN
  obs_->setRange("full",obs_->getMin(), obs_->getMax());
    RooFitResult* fit=model_->fitTo(*data_, RooFit::Save(kTRUE), RooFit::Extended(kTRUE), RooFit::Strategy(0), RooFit::PrintLevel(printLevel_), RooFit::Minimizer("GSLSimAn"));
  if(printLevel_>=2) fit->Print();
  drawFitResult(TString(label)+"_siman","full");

  // do things more carefully
  double minmassrange=std::max(sigmass_->getVal()-1.1*getStepSize(), sigmass_->getMin());
  double maxmassrange=std::min(sigmass_->getVal()+1.1*getStepSize(), sigmass_->getMax());
  doSigOnlyFloatMassFit(label, minmassrange, maxmassrange);
  
  return;
}

void ResonanceCalculatorAbs::scanForBumpWithControl(const char* label)
{
  // check for an overlap between the signal and control regions
  if((sigmass_->getMin()>controlMin_ && sigmass_->getMin()<controlMax_) ||
     (sigmass_->getMax()>controlMin_ && sigmass_->getMax()<controlMax_) ||
     (sigmass_->getMin()<controlMin_ && sigmass_->getMax()>controlMax_)) {
    if(printLevel_>=0) std::cout << "The signal and control regions overlap.  Are you sure you want to do that?" << std::endl;
  }
  
  // do B-only fit, first
  double oldMin=obs_->getMin();
  double oldMax=obs_->getMax();
  obs_->setRange(controlMin_, controlMax_);
  doBkgOnlyFit(label);
  obs_->setRange(oldMin, oldMax);

  // scan the mass range starting at the minimum
  sigmass_->setVal(sigmass_->getMin());
  double currmass=sigmass_->getVal();
  double nextmass=currmass+getStepSize();
  double bestSig=0.0, bestSigMass=0.0;
  int counter=0;
  while(currmass<sigmass_->getMax()) {
    ++counter;

    // do a fit
    double minmassrange=std::max(currmass-1.1*getStepSize(), sigmass_->getMin());
    double maxmassrange=std::min(currmass+1.1*getStepSize(), sigmass_->getMax());
    char newlabel[1000];
    sprintf(newlabel, "%sRegion%d",label,counter);
    if(doSigOnlyFloatMassFit(newlabel, minmassrange, maxmassrange)->status()!=0) continue;

    // see if the bump is significant
    double significance=evaluateTestStatistic();
    if(significance>bestSig) {
      bestSig=significance;
      bestSigMass=sigmass_->getVal();
    }

    if(printLevel_>=1) {
      std::cout << "ResonanceCalculatorAbs::scanForBumpWithControl(): "
		<< "In the range [" << minmassrange << ", " << maxmassrange << "]"
		<< " we found a bump at m=" << sigmass_->getVal() << " with significance="
		<< significance << ".  The best significance so far is "
		<< bestSig << " at m=" << bestSigMass << std::endl;
    }

    // step
    sigmass_->setVal(nextmass);
    currmass=sigmass_->getVal();
    nextmass=currmass+getStepSize();
  }

  if(bestSig<=0.0) {
    if(printLevel_>=0) std::cout << "ResonanceCalculatorAbs::scanForBumpWithControl(): "
				 << "No significant excess found!" << std::endl;
    return;
  }

  // redo fit
  sigmass_->setVal(bestSigMass);
  doSigOnlyFixMassFit(label);

  return;
}


RooFitResult* ResonanceCalculatorAbs::doBkgOnlyFit(const char* label)
{
  // do a fit to the background only
  if(!fixNbkg_) nbkg_->setConstant(false);
  setBkgParamsConst(false);
  nsig_->setMin(0.0);
  nsig_->setVal(0.0);
  nsig_->setConstant(true);
  sigmass_->setConstant(true);
  setSigParamsConst(true);

  // set range
  obs_->setRange("full",obs_->getMin(), obs_->getMax());

  // do the fit
  return doFit(TString("BkgOnly_")+label, "full");
}

RooFitResult* ResonanceCalculatorAbs::doBkgOnlyExcludeWindowFit(const char* label)
{
  // do a fit to the background only excluding a window +/- 3 units in width about the resonance mass
  if(!fixNbkg_) nbkg_->setConstant(false);
  setBkgParamsConst(false);
  nsig_->setMin(0.0);
  nsig_->setVal(0);
  nsig_->setConstant(true);
  sigmass_->setConstant(true);
  setSigParamsConst(true);

  // exclude window +/- 3 units in width
  double minobs=obs_->getMin();
  double maxobs=obs_->getMax();
  double minexclude = std::max(minobs, sigmass_->getVal()-3.0*sigwidth_->getVal());
  double maxexclude = std::min(maxobs, sigmass_->getVal()+3.0*sigwidth_->getVal());

  // do the fit
  RooFitResult* fit;
  if(maxexclude>=maxobs) {
    obs_->setRange("SB1", minobs, minexclude);
    if(printLevel_>=2) std::cout << "ResonanceCalculatorAbs::doBkgOnlyExcludeWindowFit(): "
				 << "Fitting the background in the range: [" << minobs << ", " << minexclude << "]" << std::endl;
    fit=doFit(TString("BkgOnlyExcludeWindow_")+label, "SB1");
    obs_->removeRange("SB1");

  } else if(minexclude<=minobs) {
    obs_->setRange("SB2", maxexclude, maxobs);
    if(printLevel_>=2) std::cout << "ResonanceCalculatorAbs::doBkgOnlyExcludeWindowFit(): "
				 << "Fitting the background in the range: [" << maxexclude << ", " << maxobs << "]" << std::endl;
    fit=doFit(TString("BkgOnlyExcludeWindow_")+label, "SB2");
    obs_->removeRange("SB2");

  } else {
    obs_->setRange("SB1", minobs, minexclude);
    obs_->setRange("SB2", maxexclude, maxobs);
    if(printLevel_>=2) std::cout << "ResonanceCalculatorAbs::doBkgOnlyExcludeWindowFit(): "
				 << "Fitting the background in the range: [" << minobs << ", " << minexclude << "], [" << maxexclude << ", " << maxobs << "]" << std::endl;
    fit=doFit(TString("BkgOnlyExcludeWindow_")+label, "SB1,SB2");
    obs_->removeRange("SB1");
    obs_->removeRange("SB2");
  }
  
  // set back to the original full range
  obs_->setRange("full",minobs, maxobs);

  return fit;
}

RooFitResult* ResonanceCalculatorAbs::doBkgPlusSigFixMassFit(const char* label)
{
  // do a fit to signal+background with a fixed mass
  if(!fixNbkg_) nbkg_->setConstant(false);
  setBkgParamsConst(false);
  nsig_->setRange(0.01,nbkg_->getVal());
  nsig_->setVal(nbkg_->getVal()/2.0); // set seed away from 0
  nsig_->setConstant(false);
  sigmass_->setConstant(true);
  setSigParamsConst(false);

  // set range
  obs_->setRange("full",obs_->getMin(), obs_->getMax());

  // do the fit
  return doFit(TString("BkgPlusSigFixMass_")+label, "full");
}

RooFitResult* ResonanceCalculatorAbs::doSigOnlyFixMassFit(const char* label)
{
  // do a fit to signal-only with a fixed mass
  nbkg_->setConstant(true);
  setBkgParamsConst(true);
  nsig_->setRange(0.01,nbkg_->getVal());
  nsig_->setVal(nbkg_->getVal()/2.0); // set seed away from 0
  nsig_->setConstant(false);
  sigmass_->setConstant(true);
  setSigParamsConst(false);

  // set range
  obs_->setRange("full",obs_->getMin(), obs_->getMax());

  // do the fit
  return doFit(TString("SigOnlyFixMass_")+label, "full");
}

RooFitResult* ResonanceCalculatorAbs::doSigOnlyFloatMassFit(const char* label, double minMass, double maxMass)
{
  // do a fit to signal only with a floating mass
  nbkg_->setConstant(true);
  setBkgParamsConst(true);
  nsig_->setRange(0.01,nbkg_->getVal());
  nsig_->setVal(nbkg_->getVal()/2.0); // set seed away from 0
  nsig_->setConstant(false);
  sigmass_->setConstant(false);
  setSigParamsConst(false);

  // save the resmass range from the old setting
  double oldMassMin = sigmass_->getMin();
  double oldMassMax = sigmass_->getMax();
  sigmass_->setRange(minMass, maxMass);
  sigmass_->setVal((minMass+maxMass)/2.0);
  
  // set observable range
  obs_->setRange("full",obs_->getMin(), obs_->getMax());

  RooFitResult* fit=doFit(TString("SigOnlyFloatMass_")+label, "full");

  // set the range back to the original
  sigmass_->setRange(oldMassMin, oldMassMax);

  return fit;
}


RooFitResult* ResonanceCalculatorAbs::doBkgPlusSigFloatMassFit(const char* label, double minMass, double maxMass)
{
  // do a fit to signal+background with a floating mass
  if(!fixNbkg_) nbkg_->setConstant(false);
  setBkgParamsConst(false);
  nsig_->setRange(0.01,nbkg_->getVal());
  nsig_->setVal(nbkg_->getVal()/2.0); // set seed away from 0
  nsig_->setConstant(false);
  sigmass_->setConstant(false);
  setSigParamsConst(false);

  // save the resmass range from the old setting
  double oldMassMin = sigmass_->getMin();
  double oldMassMax = sigmass_->getMax();
  sigmass_->setRange(minMass, maxMass);
  sigmass_->setVal((minMass+maxMass)/2.0);
  
  // set observable range
  obs_->setRange("full",obs_->getMin(), obs_->getMax());

  RooFitResult* fit=doFit(TString("BkgPlusSigFloatMass_")+label, "full");

  // set the range back to the original
  sigmass_->setRange(oldMassMin, oldMassMax);

  return fit;
}

RooFitResult* ResonanceCalculatorAbs::doFit(const char* label, const char* range)
{
  // do a fit of the model to the data

  RooFitResult* fit;
  if(range)
    fit=model_->fitTo(*data_, RooFit::Save(kTRUE), RooFit::Extended(kTRUE), RooFit::Strategy(fitStrategy_), RooFit::PrintLevel(printLevel_), RooFit::Range(range));
  else
    fit=model_->fitTo(*data_, RooFit::Save(kTRUE), RooFit::Extended(kTRUE), RooFit::Strategy(fitStrategy_), RooFit::PrintLevel(printLevel_));

  if(printLevel_>=2) fit->Print();

  drawFitResult(label, range);

  if(fit->status()!=0)
    std::cerr << "ResonanceCalculatorAbs::doFit(" << label << ") failed to converge properly." << std::endl;

  if(printLevel_>=0) std::cout << "ResonanceCalculatorAbs::doFit(" << label << ") storing RooFitResult as " << label << "_result" << std::endl;
  fit->Write(TString(label)+"_result");

  return fit;
}

void ResonanceCalculatorAbs::drawFitResult(const char* label, const char* range) {
  // draw fit
  RooPlot* plot = obs_->frame(RooFit::Name(label));
  data_->plotOn(plot, RooFit::Binning(RooBinning(nBinsToDraw_, obs_->getMin(), obs_->getMax())));
  if(range) {
    model_->plotOn(plot, RooFit::LineColor(kBlue+2), RooFit::NormRange(range), RooFit::Range(range));
    if(nsig_->getVal()>0) {
      model_->plotOn(plot, RooFit::Components("signal"), RooFit::LineColor(kRed+2), RooFit::NormRange(range), RooFit::Range(range));
      model_->plotOn(plot, RooFit::Components("background"), RooFit::LineColor(kBlue+2), RooFit::LineStyle(kDotted), RooFit::NormRange(range), RooFit::Range(range));
    }
  } else {
    model_->plotOn(plot, RooFit::LineColor(kBlue+2));
    if(nsig_->getVal()>0) {
      model_->plotOn(plot, RooFit::Components("signal"), RooFit::LineColor(kRed+2));
      model_->plotOn(plot, RooFit::Components("background"), RooFit::LineColor(kBlue+2), RooFit::LineStyle(kDotted));
    }
  }
  
  model_->paramOn(plot, RooFit::Layout(0.43, 0.88, 0.92), RooFit::Format("NEU",RooFit::AutoPrecision(1)) );
  plot->GetYaxis()->SetRangeUser(0.1,plot->GetMaximum()*2.0);
  plot->GetYaxis()->SetTitleOffset(1.2);
  plot->Write();
  delete plot;

  // calculate pull and diff
  RooArgSet observables(*obs_);
  double integral = nsig_->getVal()+nbkg_->getVal();
  TH1* dataHist=data_->createHistogram(TString("dataHist")+label,*obs_,RooFit::Binning(RooBinning(nBinsToDraw_, obs_->getMin(), obs_->getMax())));
  TH1D* pullHist=(TH1D*)dataHist->Clone(TString(label)+"_pull");
  TH1D* diffHist=(TH1D*)dataHist->Clone(TString(label)+"_diff");
  //  double chi2=0;
  for(int i=1; i<=dataHist->GetNbinsX(); i++) {
    double content = dataHist->GetBinContent(i);
    double error = content>25 ? sqrt(content) :
      (TMath::ChisquareQuantile(0.95,2*(content+1)-TMath::ChisquareQuantile(0.05,2*content)))/4.0;
    
    // evaluate the pdf
    obs_->setVal(dataHist->GetBinCenter(i));
    double fitval = model_->getVal(&observables)*integral*dataHist->GetBinWidth(i);
    
    // calculate the pull and diff
    double pullval=error==0 ? 0 : (content-fitval)/error;
    double diffval=(content-fitval)/fitval;
    double differr=error/fitval;
    //    chi2+=fabs(pullval);
    pullHist->SetBinContent(i, pullval);
    pullHist->SetBinError(i, 1.0);
    diffHist->SetBinContent(i, diffval);
    diffHist->SetBinError(i, differr);
  }
  
  // draw pull
  pullHist->SetMarkerStyle(20);
  //  pullHist->Draw("E");
  pullHist->GetYaxis()->SetTitle("Pull");
  pullHist->Write();
  
  // draw diff
  diffHist->SetMarkerStyle(20);

  diffHist->SetMaximum(0.3);
  diffHist->SetMinimum(-0.3);
  diffHist->GetYaxis()->SetTitle("(data-fit)/fit");
  diffHist->Write();

  delete dataHist;
  delete pullHist;
  delete diffHist;

  return;
}

void ResonanceCalculatorAbs::findSignalAndBackgroundParams(void)
{
  // get the signal and background parameters
  // RooArgSet objects must be deallocated later

  sigparams_ = signal_->getParameters(RooArgSet(*obs_));
  bkgparams_ = background_->getParameters(RooArgSet(*obs_));

  // remove the signal mass from the list of signal parameters
  sigparams_->remove(*sigmass_);

  if(printLevel_>=2) {
    std::cout << "ResonanceCalculatorAbs::findSignalAndBackgroundParams(): Signal parameters: ";
    sigparams_->Print();
    std::cout << "ResonanceCalculatorAbs::findSignalAndBackgroundParams(): Background parameters: ";
    bkgparams_->Print();
    std::cout << std::endl;
  }

  return;
}

void ResonanceCalculatorAbs::setSigParamsConst(bool isConst)
{
  // set all of the signal parameters (except the mass) to a constant

  TIterator *iter = sigparams_->createIterator();
  TObject* obj ;
  while((obj=iter->Next())) {
    RooRealVar* var=dynamic_cast<RooRealVar*>(obj);
    if(var) var->setConstant(isConst);
  }
  delete iter;

  return;
}

void ResonanceCalculatorAbs::setBkgParamsConst(bool isConst)
{
  // set all of the background parameters to a constant

  TIterator *iter = bkgparams_->createIterator();
  TObject* obj ;
  while((obj=iter->Next())) {
    RooRealVar* var=dynamic_cast<RooRealVar*>(obj);
    if(var) var->setConstant(isConst);
  }
  delete iter;

  return;
}

void ResonanceCalculatorAbs::copyValuesToBkgParams(RooArgList* params)
{
  // set all of the background parameters to the values of the parameters in the argument

  TIterator *iter = bkgparams_->createIterator();
  TObject* obj;
  while((obj=iter->Next())) {
    RooRealVar* var=dynamic_cast<RooRealVar*>(obj);
    if(var) {
      RooRealVar* varcopy=dynamic_cast<RooRealVar*>(params->find(var->GetName()));
      if(varcopy) var->setVal(varcopy->getVal());
    }
  }

  delete iter;

  return;
}

void ResonanceCalculatorAbs::copyValuesToSigParams(RooArgList* params)
{
  // set all of the background parameters to the values of the parameters in the argument

  TIterator *iter = sigparams_->createIterator();
  TObject* obj;
  while((obj=iter->Next())) {
    RooRealVar* var=dynamic_cast<RooRealVar*>(obj);
    if(var) {
      RooRealVar* varcopy=dynamic_cast<RooRealVar*>(params->find(var->GetName()));
      if(varcopy) var->setVal(varcopy->getVal());
    }
  }

  delete iter;

  return;
}


void ResonanceCalculatorAbs::findMinMaxMass(void)
{
  // find the mass range based on +/-3 units in width from the data boundaries

  // set the range so that we can move the value of sigmass
  sigmass_->setRange(obs_->getMin()-1.0, obs_->getMax()+1.0);

  sigmass_->setVal(obs_->getMin());
  sigmass_->setMin(sigmass_->getVal()+1.0*sigwidth_->getVal());

  sigmass_->setVal(obs_->getMax());
  sigmass_->setMax(sigmass_->getVal()-1.0*sigwidth_->getVal());

  if(printLevel_>=2) {
    std::cout << "ResonanceCalculatorAbs::findMinMaxMass(): Chose the range ["
	      << sigmass_->getMin() << ", " << sigmass_->getMax() << "] for sigmass" << std::endl;
  }
  return;
}

double ResonanceCalculatorAbs::getStepSize(void) const
{
  return fixedSearchStepSize_>0 ? fixedSearchStepSize_ : searchStepSize_*sigwidth_->getVal();
}

double ResonanceCalculatorAbs::evaluateTestStatistic(void)
{
  if(whichTestStatistic_==0) {
    nsig_->setMin(0.0);
    double originalVal=nsig_->getVal();
    nsig_->setConstant(true);
    RooArgSet poi(*nsig_);
    RooArgSet *sbParams=dynamic_cast<RooArgSet*>(poi.snapshot());
    nsig_->setVal(0.0);
    RooArgSet *bParams=dynamic_cast<RooArgSet*>(poi.snapshot());
    RooStats::SimpleLikelihoodRatioTestStat slrts(*model_, *model_, *bParams, *sbParams);
    double sig=slrts.Evaluate(*data_, poi);
    delete bParams;
    delete sbParams;
    nsig_->setVal(originalVal);
    return sig>0 ? sqrt(2*sig) : 0.0;
  }
  else if(whichTestStatistic_==1) {
    return nsig_->getVal()/nsig_->getError();
  } else {
    assert(1);
  }
  return 0.0;
}

RooDataHist* ResonanceCalculatorAbs::generateBinned(RooAbsPdf* pdf, RooDataHist* templateDataHist, int numEntries)
{
  // generate a pseudo dataset
  RooDataSet* dataset = pdf->generate(RooArgSet(*obs_), RooFit::NumEvents(numEntries), RooFit::Name("tempdata"));

  // create a histogram out of the template data histogram and empty it
  obs_->setBinning(dataBinning_);
  TH1* hist=templateDataHist->createHistogram("temphist", *obs_);
  for(int bin=0; bin<=hist->GetNbinsX()+1; bin++) {
    hist->SetBinContent(bin, 0.0);
    hist->SetBinError(bin, 0.0);
  }

  // fill the histogram with the unbinned data
  hist = dataset->fillHistogram(hist, RooArgList(*obs_));

  // create the new histogrammed data
  RooDataHist* datahist=new RooDataHist("BinnedPseudoData", "binned data", RooArgSet(*obs_), RooFit::Import(*hist, kFALSE));

  // delete temporaries
  delete dataset;
  delete hist;

  return datahist;
}

RooDataSet* ResonanceCalculatorAbs::generateUnbinned(RooAbsPdf* pdf, int numEntries)
{
  return pdf->generate(RooArgSet(*obs_), RooFit::NumEvents(numEntries), RooFit::Name("pseudodata"));
}



double ResonanceCalculatorAbs::pValueToZScore(double pvalue)
{
  if(pvalue<=0.0 && numeric_limits<double>::has_infinity) {
    return numeric_limits<double>::infinity();
  } else if(pvalue<=0.0) {
    return numeric_limits<double>::max();
  } else if(pvalue>=1.0 && numeric_limits<double>::has_infinity) {
    return -numeric_limits<double>::infinity();
  } else if(pvalue>=1.0) {
    return -numeric_limits<double>::max();
  } else {
    return TMath::Sqrt(2.0)*TMath::ErfcInverse(2*pvalue);
  }
}

double ResonanceCalculatorAbs::zScoreToPValue(double zscore)
{
  if(numeric_limits<double>::has_infinity && zscore==numeric_limits<double>::infinity()) {
    return 0.0;
  } else if(!numeric_limits<double>::has_infinity && zscore==numeric_limits<double>::max()) {
    return 0.0;
  } else if(numeric_limits<double>::has_infinity && zscore==-numeric_limits<double>::infinity()) {
    return 1.0;
  } else if(!numeric_limits<double>::has_infinity && zscore==-numeric_limits<double>::max()) {
    return 1.0;
  } else {
    return TMath::Erfc(zscore/TMath::Sqrt(2.0))/2.0;
  }
}

double ResonanceCalculatorAbs::getPValue(std::vector<std::pair<double, double> >& teststats, double bumpTestStat)
{
  double integral=0.0;
  double tailintegral=0.0;
  for(unsigned int i=0; i<teststats.size(); i++) {
    integral += teststats[i].second;
    if(teststats[i].first>=bumpTestStat) tailintegral+= teststats[i].second;
  }
  if(integral<=0.0) return 0.0;
  return tailintegral/integral;
}
 
double ResonanceCalculatorAbs::getZScore(std::vector<std::pair<double, double> >& teststats, double bumpTestStat)
{
  return pValueToZScore(getPValue(teststats, bumpTestStat));
}

std::pair<double, double> ResonanceCalculatorAbs::getPValueRange(std::vector<std::pair<double, double> >& teststats, double bumpTestStat, double alpha)
{
  if(teststats.size()==0) {
    std::cout << "test stats does not have an entries.  The range is invalid." << std::endl;
    return std::pair<double, double>(0.0, 0.0);
  }

  double integral=0.0;
  double tailintegral=0.0;
  int numtail=0;
  for(unsigned int i=0; i<teststats.size(); i++) {
    integral += teststats[i].second;
    if(teststats[i].first>=bumpTestStat) {
      tailintegral+=teststats[i].second;
      ++numtail;
    }
  }
  int a=numtail;
  int b=static_cast<int>(teststats.size())-numtail;
  double x= a>0 ? tailintegral/a : 1.0;
  double y= b>0 ? (integral-tailintegral)/b : 1.0;

  // p-value, p = x*a/(x*a+y*b)
  // f=a/(a+b)
  // p=x*f/(x*f+y-y*f)
  // fc computes the interval on f, not p
  // we then translate the interval on f to an interval on p
  
  FeldmanCousinsBinomialInterval fc;
  fc.init(alpha);
  fc.calculate(a,a+b);
  double upper=b>0 ? fc.upper() : 1.0;
  double lower=a>0 ? fc.lower() : 0.0;
  double r1=x*upper/(x*upper+y-y*upper);
  double r2=x*lower/(x*lower+y-y*lower);
  std::pair<double, double> range(std::min(r1, r2), std::max(r1, r2));
  return range;
}

std::pair<double, double> ResonanceCalculatorAbs::getZScoreRange(std::vector<std::pair<double, double> >& teststats, double bumpTestStat, double alpha)
{
  std::pair<double, double> range=getPValueRange(teststats, bumpTestStat, alpha);
  double z1 = pValueToZScore(range.first);
  double z2 = pValueToZScore(range.second);
  range.first=std::min(z1, z2);
  range.second=std::max(z1, z2);
  return range;
}
