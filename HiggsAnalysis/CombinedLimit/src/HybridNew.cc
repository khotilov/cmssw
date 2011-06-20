#include <stdexcept>
#include <cstdio>
#include <unistd.h>
#include <sys/types.h>
#include <sys/wait.h>
#include <errno.h>

#include "../interface/HybridNew.h"
#include <TFile.h>
#include <TF1.h>
#include <TKey.h>
#include <TLine.h>
#include <TCanvas.h>
#include <TGraphErrors.h>
#include <TStopwatch.h>
#include "RooRealVar.h"
#include "RooArgSet.h"
#include "RooAbsPdf.h"
#include "RooFitResult.h"
#include "RooRandom.h"
#include "RooAddPdf.h"
#include "RooConstVar.h"
#include <RooStats/ModelConfig.h>
#include <RooStats/HybridCalculator.h>
#include <RooStats/SimpleLikelihoodRatioTestStat.h>
#include <RooStats/RatioOfProfiledLikelihoodsTestStat.h>
#include <RooStats/ProfileLikelihoodTestStat.h>
#include <RooStats/ToyMCSampler.h>
#include <RooStats/HypoTestPlot.h>
#include "../interface/Combine.h"
#include "../interface/CloseCoutSentry.h"
#include "../interface/RooFitGlobalKillSentry.h"
#include "../interface/SimplerLikelihoodRatioTestStat.h"
#include "../interface/ProfiledLikelihoodRatioTestStat.h"
#include "../interface/SimplerLikelihoodRatioTestStatExt.h"
#include "../interface/ProfiledLikelihoodRatioTestStatExt.h"
#include "../interface/utils.h"
#include "../interface/ProfileLikelihood.h"

using namespace RooStats;

HybridNew::WorkingMode HybridNew::workingMode_ = MakeLimit;
unsigned int HybridNew::nToys_ = 500;
double HybridNew::clsAccuracy_ = 0.005;
double HybridNew::rAbsAccuracy_ = 0.1;
double HybridNew::rRelAccuracy_ = 0.05;
std::string HybridNew::rule_ = "CLs";
std::string HybridNew::testStat_ = "LEP";
bool HybridNew::genNuisances_ = true;
bool HybridNew::genGlobalObs_ = false;
bool HybridNew::fitNuisances_ = false;
unsigned int HybridNew::iterations_ = 1;
unsigned int HybridNew::nCpu_ = 0; // proof-lite mode
unsigned int HybridNew::fork_ = 1; // fork mode
double HybridNew::rValue_  = 1.0;
bool HybridNew::CLs_ = false;
bool HybridNew::saveHybridResult_  = false;
bool HybridNew::readHybridResults_ = false; 
bool  HybridNew::expectedFromGrid_ = false; 
float HybridNew::quantileForExpectedFromGrid_ = 0.5; 
std::string HybridNew::gridFile_ = "";
bool HybridNew::importanceSamplingNull_ = false;
bool HybridNew::importanceSamplingAlt_  = false;
std::string HybridNew::algo_ = "logSecant";
bool HybridNew::optimizeProductPdf_     = true;
bool HybridNew::optimizeTestStatistics_ = true;
std::string HybridNew::plot_;
std::string HybridNew::minimizerAlgo_ = "Minuit2";
float       HybridNew::minimizerTolerance_ = 1e-2;

HybridNew::HybridNew() : 
LimitAlgo("HybridNew specific options") {
    options_.add_options()
        ("rule",    boost::program_options::value<std::string>(&rule_)->default_value(rule_),            "Rule to use: CLs, CLsplusb")
        ("testStat",boost::program_options::value<std::string>(&testStat_)->default_value(testStat_),    "Test statistics: LEP, TEV, LHC (previously known as Atlas), Profile.")
        ("singlePoint",  boost::program_options::value<float>(),  "Just compute CLs for the given value of r")
        ("onlyTestStat", "Just compute test statistics, not actual p-values (works only with --singlePoint)")
        ("generateNuisances",            boost::program_options::value<bool>(&genNuisances_)->default_value(genNuisances_), "Generate nuisance parameters for each toy")
        ("generateExternalMeasurements", boost::program_options::value<bool>(&genGlobalObs_)->default_value(genGlobalObs_), "Generate external measurements for each toy, taken from the GlobalObservables of the ModelConfig")
        ("fitNuisances", boost::program_options::value<bool>(&fitNuisances_)->default_value(fitNuisances_), "Fit the nuisances, when not generating them.")
        ("searchAlgo", boost::program_options::value<std::string>(&algo_)->default_value(algo_),         "Algorithm to use to search for the limit (bisection, logSecant)")
        ("toysH,T", boost::program_options::value<unsigned int>(&nToys_)->default_value(nToys_),         "Number of Toy MC extractions to compute CLs+b, CLb and CLs")
        ("clsAcc",  boost::program_options::value<double>(&clsAccuracy_ )->default_value(clsAccuracy_),  "Absolute accuracy on CLs to reach to terminate the scan")
        ("rAbsAcc", boost::program_options::value<double>(&rAbsAccuracy_)->default_value(rAbsAccuracy_), "Absolute accuracy on r to reach to terminate the scan")
        ("rRelAcc", boost::program_options::value<double>(&rRelAccuracy_)->default_value(rRelAccuracy_), "Relative accuracy on r to reach to terminate the scan")
        ("iterations,i", boost::program_options::value<unsigned int>(&iterations_)->default_value(iterations_), "Number of times to throw 'toysH' toys to compute the p-values (for --singlePoint if clsAcc is set to zero disabling adaptive generation)")
        ("fork",    boost::program_options::value<unsigned int>(&fork_)->default_value(fork_),           "Fork to N processes before running the toys (set to 0 for debugging)")
        ("nCPU",    boost::program_options::value<unsigned int>(&nCpu_)->default_value(nCpu_),           "Use N CPUs with PROOF Lite (experimental!)")
        ("saveHybridResult",  "Save result in the output file")
        ("readHybridResults", "Read and merge results from file (requires 'toysFile' or 'grid')")
        ("grid",    boost::program_options::value<std::string>(&gridFile_),            "Use the specified file containing a grid of SamplingDistributions for the limit (implies readHybridResults).\n For --singlePoint or --signif use --toysFile=x.root --readHybridResult instead of this.")
        ("expectedFromGrid", boost::program_options::value<float>(&quantileForExpectedFromGrid_)->default_value(0.5), "Use the grid to compute the expected limit for this quantile")
        ("importanceSamplingNull", boost::program_options::value<bool>(&importanceSamplingNull_)->default_value(importanceSamplingNull_),  
                                   "Enable importance sampling for null hypothesis (background only)") 
        ("importanceSamplingAlt",  boost::program_options::value<bool>(&importanceSamplingAlt_)->default_value(importanceSamplingAlt_),    
                                   "Enable importance sampling for alternative hypothesis (signal plus background)") 
        ("optimizeTestStatistics", boost::program_options::value<bool>(&optimizeTestStatistics_)->default_value(optimizeTestStatistics_), 
                                   "Use optimized test statistics if the likelihood is not extended (works for LEP and TEV test statistics).")
        ("optimizeProductPdf",     boost::program_options::value<bool>(&optimizeProductPdf_)->default_value(optimizeProductPdf_),      
                                   "Optimize the code factorizing pdfs")
        ("minimizerAlgo",      boost::program_options::value<std::string>(&minimizerAlgo_)->default_value(minimizerAlgo_), "Choice of minimizer used for profiling (Minuit vs Minuit2)")
        ("minimizerTolerance", boost::program_options::value<float>(&minimizerTolerance_)->default_value(minimizerTolerance_),  "Tolerance for minimizer used for profiling")
        ("plot",   boost::program_options::value<std::string>(&plot_), "Save a plot of the result (test statistics distributions or limit scan)")
        ("frequentist", "Shortcut to switch to Frequentist mode (--generateNuisances=0 --generateExternalMeasurements=1 --fitNuisances=1)")
    ;
}

void HybridNew::applyOptions(const boost::program_options::variables_map &vm) {
    if (vm.count("expectedFromGrid") && !vm["expectedFromGrid"].defaulted()) {
        if (!vm.count("grid")) throw std::invalid_argument("HybridNew: Can't use --expectedFromGrid without --grid!");
        if (quantileForExpectedFromGrid_ <= 0 || quantileForExpectedFromGrid_ >= 1.0) throw std::invalid_argument("HybridNew: the quantile for the expected limit must be between 0 and 1");
        expectedFromGrid_ = true;
    } 
    if (vm.count("frequentist")) {
        genNuisances_ = 0; genGlobalObs_ = withSystematics; fitNuisances_ = withSystematics;
        if (vm["testStat"].defaulted()) testStat_ = "LHC";
    }
    if (vm.count("singlePoint")) {
        if (doSignificance_) throw std::invalid_argument("HybridNew: Can't use --significance and --singlePoint at the same time");
        rValue_ = vm["singlePoint"].as<float>();
        workingMode_ = ( vm.count("onlyTestStat") ? MakeTestStatistics : MakePValues );
        rValue_ = vm["singlePoint"].as<float>();
    } else if (vm.count("onlyTestStat")) {
        throw std::invalid_argument("HybridNew: --onlyTestStat works only with --singlePoint");
    } else if (doSignificance_) {
        workingMode_ = MakeSignificance;
    } else {
        workingMode_ = MakeLimit;
    }
    saveHybridResult_ = vm.count("saveHybridResult");
    readHybridResults_ = vm.count("readHybridResults") || vm.count("grid");
    if (readHybridResults_ && !(vm.count("toysFile") || vm.count("grid")))     throw std::invalid_argument("HybridNew: must have 'toysFile' or 'grid' option to have 'readHybridResults'\n");
    validateOptions(); 
}

void HybridNew::applyDefaultOptions() { 
    workingMode_ = MakeLimit;
    validateOptions(); 
}

void HybridNew::validateOptions() {
    if (fork_ > 1) nToys_ /= fork_; // makes more sense
    if (rule_ == "CLs") {
        CLs_ = true;
    } else if (rule_ == "CLsplusb") {
        CLs_ = false;
    } else {
        throw std::invalid_argument("HybridNew: Rule must be either 'CLs' or 'CLsplusb'");
    }
    if (testStat_ == "SimpleLikelihoodRatio"      || testStat_ == "SLR" ) { testStat_ = "LEP";     }
    if (testStat_ == "RatioOfProfiledLikelihoods" || testStat_ == "ROPL") { testStat_ = "TEV";     }
    if (testStat_ == "ProfileLikelihood"          || testStat_ == "PL")   { testStat_ = "Profile"; }
    if (testStat_ == "ModifiedProfileLikelihood"  || testStat_ == "MPL")  { testStat_ = "LHC";     }
    if (testStat_ == "Atlas") { testStat_ = "LHC"; std::cout << "Note: the Atlas test statistics is now known as LHC test statistics.\n" << std::endl; }
    if (testStat_ != "LEP" && testStat_ != "TEV" && testStat_ != "LHC" && testStat_ != "Profile") {
        throw std::invalid_argument("HybridNew: Test statistics should be one of 'LEP' or 'TEV' or 'LHC' (previously known as 'Atlas') or 'Profile'");
    }
    if (verbose) {
        if (testStat_ == "LEP")     std::cout << ">>> using the Simple Likelihood Ratio test statistics (Q_LEP)" << std::endl;
        if (testStat_ == "TEV")     std::cout << ">>> using the Ratio of Profiled Likelihoods test statistics (Q_TEV)" << std::endl;
        if (testStat_ == "LHC")     std::cout << ">>> using the Profile Likelihood test statistics modified for upper limits (Q_LHC)" << std::endl;
        if (testStat_ == "Profile") std::cout << ">>> using the Profile Likelihood test statistics not modified for upper limits (Q_Profile)" << std::endl;
    }
}

bool HybridNew::run(RooWorkspace *w, RooStats::ModelConfig *mc_s, RooStats::ModelConfig *mc_b, RooAbsData &data, double &limit, double &limitErr, const double *hint) {
    RooFitGlobalKillSentry silence(verbose <= 1 ? RooFit::WARNING : RooFit::DEBUG);
    ProfileLikelihood::MinimizerSentry minimizerConfig(minimizerAlgo_, minimizerTolerance_);
    perf_totalToysRun_ = 0; // reset performance counter
    switch (workingMode_) {
        case MakeLimit:            return runLimit(w, mc_s, mc_b, data, limit, limitErr, hint);
        case MakeSignificance:     return runSignificance(w, mc_s, mc_b, data, limit, limitErr, hint);
        case MakePValues:          return runSinglePoint(w, mc_s, mc_b, data, limit, limitErr, hint);
        case MakeTestStatistics:   return runTestStatistics(w, mc_s, mc_b, data, limit, limitErr, hint);
    }
    assert("Shouldn't get here" == 0);
}

bool HybridNew::runSignificance(RooWorkspace *w, RooStats::ModelConfig *mc_s, RooStats::ModelConfig *mc_b, RooAbsData &data, double &limit, double &limitErr, const double *hint) {
    HybridNew::Setup setup;
    std::auto_ptr<RooStats::HybridCalculator> hc(create(w, mc_s, mc_b, data, 1.0, setup));
    std::auto_ptr<HypoTestResult> hcResult;
    if (readHybridResults_) {
        hcResult.reset(readToysFromFile());
    } else {
        hcResult.reset(evalGeneric(*hc));
        for (unsigned int i = 1; i < iterations_; ++i) {
            std::auto_ptr<HypoTestResult> more(evalGeneric(*hc));
            hcResult->Append(more.get());
            if (verbose) std::cout << "\tCLb = " << hcResult->CLb() << " +/- " << hcResult->CLbError() << std::endl;
        }
    }
    if (hcResult.get() == 0) {
        std::cerr << "Hypotest failed" << std::endl;
        return false;
    }
    if (saveHybridResult_) {
        TString name = TString::Format("HypoTestResult_r%g_%u", 0., RooRandom::integer(std::numeric_limits<UInt_t>::max() - 1));
        writeToysHere->WriteTObject(new HypoTestResult(*hcResult), name);
        if (verbose) std::cout << "Hybrid result saved as " << name << " in " << writeToysHere->GetFile()->GetName() << " : " << writeToysHere->GetPath() << std::endl;
    }
    if (testStat_ == "LHC" || testStat_ == "Profile") {
        // I don't need to flip the P-values for significances, only for limits
        // hcResult->SetPValueIsRightTail(!hcResult->GetPValueIsRightTail());
        hcResult->SetTestStatisticData(hcResult->GetTestStatisticData()-1e-9); // issue with < vs <= in discrete models
    } else {
        hcResult->SetTestStatisticData(hcResult->GetTestStatisticData()+1e-9); // issue with < vs <= in discrete models
    }
    limit = hcResult->Significance();
    double sigHi = RooStats::PValueToSignificance( 1 - (hcResult->CLb() + hcResult->CLbError()) ) - limit;
    double sigLo = RooStats::PValueToSignificance( 1 - (hcResult->CLb() - hcResult->CLbError()) ) - limit;
    limitErr = std::max(sigHi,-sigLo);
    std::cout << "\n -- Hybrid New -- \n";
    std::cout << "Significance: " << limit << "  " << sigLo << "/+" << sigHi << " (CLb " << hcResult->CLb() << " +/- " << hcResult->CLbError() << ")\n";
    return isfinite(limit);
}

bool HybridNew::runLimit(RooWorkspace *w, RooStats::ModelConfig *mc_s, RooStats::ModelConfig *mc_b, RooAbsData &data, double &limit, double &limitErr, const double *hint) {
  RooRealVar *r = dynamic_cast<RooRealVar *>(mc_s->GetParametersOfInterest()->first()); r->setConstant(true);
  w->loadSnapshot("clean");

  if ((hint != 0) && (*hint > r->getMin())) {
    r->setMax(std::min<double>(3.0 * (*hint), r->getMax()));
    r->setMin(std::max<double>(0.3 * (*hint), r->getMin()));
  }
  
  typedef std::pair<double,double> CLs_t;

  double clsTarget = 1 - cl; 
  CLs_t clsMin(1,0), clsMax(0,0), clsMid(0,0);
  double rMin = r->getMin(), rMax = r->getMax(); 
  limit    = 0.5*(rMax + rMin);
  limitErr = 0.5*(rMax - rMin);
  bool done = false;

  TF1 expoFit("expoFit","[0]*exp([1]*(x-[2]))", rMin, rMax);

  if (readHybridResults_) { 
      if (verbose > 0) std::cout << "Search for upper limit using pre-computed grid of p-values" << std::endl;

      if (!gridFile_.empty()) {
        if (grid_.empty()) {
            std::auto_ptr<TFile> gridFile(TFile::Open(gridFile_.c_str()));
            if (gridFile.get() == 0) throw std::runtime_error(("Can't open grid file "+gridFile_).c_str());
            TDirectory *toyDir = gridFile->GetDirectory("toys");
            if (!toyDir) throw std::logic_error("Cannot use readHypoTestResult: empty toy dir in input file empty");
            readGrid(toyDir);
        }
        updateGridData(w, mc_s, mc_b, data, true, clsTarget);
      } else readAllToysFromFile(); 

      useGrid();

      double minDist=1e3;
      rMin = limitPlot_->GetX()[0]; 
      rMax = limitPlot_->GetX()[limitPlot_->GetN()-1];
      for (int i = 0, n = limitPlot_->GetN(); i < n; ++i) {
          double x = limitPlot_->GetX()[i], y = limitPlot_->GetY()[i], ey = limitPlot_->GetErrorY(i);
          if (verbose > 0) printf("  r %.2f: %s = %6.4f +/- %6.4f\n", x, CLs_ ? ", CLs = " : ", CLsplusb = ", y, ey);
          if (y-3*max(ey,0.01) >= clsTarget) { rMin = x; clsMin = CLs_t(y,ey); }
          if (fabs(y-clsTarget) < minDist) { limit = x; minDist = fabs(y-clsTarget); }
          rMax = x; clsMax = CLs_t(y,ey);    
          if (y+3*max(ey,0.01) <= clsTarget) break; 
      }
      if (verbose > 0) std::cout << " after scan x ~ " << limit << ", bounds [ " << rMin << ", " << rMax << "]" << std::endl;
      limitErr = std::max(limit-rMin, rMax-limit);
      expoFit.SetRange(rMin,rMax);

      if (limitErr < std::max(rAbsAccuracy_, rRelAccuracy_ * limit)) {
          if (verbose > 1) std::cout << "  reached accuracy " << limitErr << " below " << std::max(rAbsAccuracy_, rRelAccuracy_ * limit) << std::endl;
          done = true; 
      }
  } else {
      limitPlot_.reset(new TGraphErrors());

      if (verbose > 0) std::cout << "Search for upper limit to the limit" << std::endl;
      for (int tries = 0; tries < 6; ++tries) {
          clsMax = eval(w, mc_s, mc_b, data, rMax);
          if (clsMax.first == 0 || clsMax.first + 3 * fabs(clsMax.second) < clsTarget ) break;
          rMax += rMax;
          if (tries == 5) { 
              std::cerr << "Cannot set higher limit: at " << r->GetName() << " = " << rMax << " still get " << (CLs_ ? "CLs" : "CLsplusb") << " = " << clsMax.first << std::endl;
              return false;
          }
      }
      if (verbose > 0) std::cout << "Search for lower limit to the limit" << std::endl;
      clsMin = (CLs_ && rMin == 0 ? CLs_t(1,0) : eval(w, mc_s, mc_b, data, rMin));
      if (clsMin.first != 1 && clsMin.first - 3 * fabs(clsMin.second) < clsTarget) {
          if (CLs_) { 
              rMin = 0;
              clsMin = CLs_t(1,0); // this is always true for CLs
          } else {
              rMin = -rMax / 4;
              for (int tries = 0; tries < 6; ++tries) {
                  clsMin = eval(w, mc_s, mc_b, data, rMin);
                  if (clsMin.first == 1 || clsMin.first - 3 * fabs(clsMin.second) > clsTarget) break;
                  rMin += rMin;
                  if (tries == 5) { 
                      std::cerr << "Cannot set lower limit: at " << r->GetName() << " = " << rMin << " still get " << (CLs_ ? "CLs" : "CLsplusb") << " = " << clsMin.first << std::endl;
                      return false;
                  }
              }
          }
      }

      if (verbose > 0) std::cout << "Now doing proper bracketing & bisection" << std::endl;
      do {
          // determine point by bisection or interpolation
          limit = 0.5*(rMin+rMax); limitErr = 0.5*(rMax-rMin);
          if (algo_ == "logSecant" && clsMax.first != 0) {
              double logMin = log(clsMin.first), logMax = log(clsMax.first), logTarget = log(clsTarget);
              limit = rMin + (rMax-rMin) * (logTarget - logMin)/(logMax - logMin);
              if (clsMax.second != 0 && clsMin.second != 0) {
                  limitErr = hypot((logTarget-logMax) * (clsMin.second/clsMin.first), (logTarget-logMin) * (clsMax.second/clsMax.first));
                  limitErr *= (rMax-rMin)/((logMax-logMin)*(logMax-logMin));
              }
          }
          r->setError(limitErr);

          // exit if reached accuracy on r 
          if (limitErr < std::max(rAbsAccuracy_, rRelAccuracy_ * limit)) {
              if (verbose > 1) std::cout << "  reached accuracy " << limitErr << " below " << std::max(rAbsAccuracy_, rRelAccuracy_ * limit) << std::endl;
              done = true; break;
          }

          // evaluate point 
          clsMid = eval(w, mc_s, mc_b, data, limit, true, clsTarget);
          if (clsMid.second == -1) {
              std::cerr << "Hypotest failed" << std::endl;
              return false;
          }

          // if sufficiently far away, drop one of the points
          if (fabs(clsMid.first-clsTarget) >= 2*clsMid.second) {
              if ((clsMid.first>clsTarget) == (clsMax.first>clsTarget)) {
                  rMax = limit; clsMax = clsMid;
              } else {
                  rMin = limit; clsMin = clsMid;
              }
          } else {
              if (verbose > 0) std::cout << "Trying to move the interval edges closer" << std::endl;
              double rMinBound = rMin, rMaxBound = rMax;
              // try to reduce the size of the interval 
              while (clsMin.second == 0 || fabs(rMin-limit) > std::max(rAbsAccuracy_, rRelAccuracy_ * limit)) {
                  rMin = 0.5*(rMin+limit); 
                  clsMin = eval(w, mc_s, mc_b, data, rMin, true, clsTarget); 
                  if (fabs(clsMin.first-clsTarget) <= 2*clsMin.second) break;
                  rMinBound = rMin;
              } 
              while (clsMax.second == 0 || fabs(rMax-limit) > std::max(rAbsAccuracy_, rRelAccuracy_ * limit)) {
                  rMax = 0.5*(rMax+limit); 
                  clsMax = eval(w, mc_s, mc_b, data, rMax, true, clsTarget); 
                  if (fabs(clsMax.first-clsTarget) <= 2*clsMax.second) break;
                  rMaxBound = rMax;
              } 
              expoFit.SetRange(rMinBound,rMaxBound);
              break;
          }
      } while (true);

  }

  if (!done) { // didn't reach accuracy with scan, now do fit
      if (verbose) {
          std::cout << "\n -- HybridNew, before fit -- \n";
          std::cout << "Limit: " << r->GetName() << " < " << limit << " +/- " << limitErr << " [" << rMin << ", " << rMax << "]\n";
      }

      expoFit.FixParameter(0,clsTarget);
      expoFit.SetParameter(1,log(clsMax.first/clsMin.first)/(rMax-rMin));
      expoFit.SetParameter(2,limit);
      double rMinBound, rMaxBound; expoFit.GetRange(rMinBound, rMaxBound);
      limitErr = std::max(fabs(rMinBound-limit), fabs(rMaxBound-limit));
      int npoints = 0; 
      for (int j = 0; j < limitPlot_->GetN(); ++j) { 
          if (limitPlot_->GetX()[j] >= rMinBound && limitPlot_->GetX()[j] <= rMaxBound) npoints++; 
      }
      for (int i = 0, imax = (readHybridResults_ ? 0 : 8); i <= imax; ++i, ++npoints) {
          limitPlot_->Sort();
          limitPlot_->Fit(&expoFit,(verbose <= 1 ? "QNR EX0" : "NR EXO"));
          if (verbose) {
              std::cout << "Fit to " << npoints << " points: " << expoFit.GetParameter(2) << " +/- " << expoFit.GetParError(2) << std::endl;
          }
          if ((rMin < expoFit.GetParameter(2))  && (expoFit.GetParameter(2) < rMax) && (expoFit.GetParError(2) < 0.5*(rMaxBound-rMinBound))) { 
              // sanity check fit result
              limit = expoFit.GetParameter(2);
              limitErr = expoFit.GetParError(2);
              if (limitErr < std::max(rAbsAccuracy_, rRelAccuracy_ * limit)) break;
          }
          // add one point in the interval. 
          double rTry = RooRandom::uniform()*(rMaxBound-rMinBound)+rMinBound; 
          if (i != imax) eval(w, mc_s, mc_b, data, rTry, true, clsTarget);
      } 
  }

  if (!plot_.empty() && limitPlot_.get()) {
      TCanvas *c1 = new TCanvas("c1","c1");
      limitPlot_->Sort();
      limitPlot_->SetLineWidth(2);
      double xmin = r->getMin(), xmax = r->getMax();
      for (int j = 0; j < limitPlot_->GetN(); ++j) {
        if (limitPlot_->GetY()[j] > 1.4*clsTarget || limitPlot_->GetY()[j] < 0.6*clsTarget) continue;
        xmin = std::min(limitPlot_->GetX()[j], xmin);
        xmax = std::max(limitPlot_->GetX()[j], xmax);
      }
      limitPlot_->GetXaxis()->SetRangeUser(xmin,xmax);
      limitPlot_->GetYaxis()->SetRangeUser(0.5*clsTarget, 1.5*clsTarget);
      limitPlot_->Draw("AP");
      expoFit.Draw("SAME");
      TLine line(limitPlot_->GetX()[0], clsTarget, limitPlot_->GetX()[limitPlot_->GetN()-1], clsTarget);
      line.SetLineColor(kRed); line.SetLineWidth(2); line.Draw();
      line.DrawLine(limit, 0, limit, limitPlot_->GetY()[0]);
      line.SetLineWidth(1); line.SetLineStyle(2);
      line.DrawLine(limit-limitErr, 0, limit-limitErr, limitPlot_->GetY()[0]);
      line.DrawLine(limit+limitErr, 0, limit+limitErr, limitPlot_->GetY()[0]);
      c1->Print(plot_.c_str());
  }

  std::cout << "\n -- Hybrid New -- \n";
  std::cout << "Limit: " << r->GetName() << " < " << limit << " +/- " << limitErr << " @ " << cl * 100 << "% CL\n";
  if (verbose > 1) std::cout << "Total toys: " << perf_totalToysRun_ << std::endl;
  return true;
}

bool HybridNew::runSinglePoint(RooWorkspace *w, RooStats::ModelConfig *mc_s, RooStats::ModelConfig *mc_b, RooAbsData &data, double &limit, double &limitErr, const double *hint) {
    RooRealVar *r = dynamic_cast<RooRealVar *>(mc_s->GetParametersOfInterest()->first()); r->setConstant(true);
    std::pair<double, double> result = eval(w, mc_s, mc_b, data, rValue_, clsAccuracy_ != 0);
    std::cout << "\n -- Hybrid New -- \n";
    std::cout << (CLs_ ? "CLs = " : "CLsplusb = ") << result.first << " +/- " << result.second << std::endl;
    if (verbose > 1) std::cout << "Total toys: " << perf_totalToysRun_ << std::endl;
    limit = result.first;
    limitErr = result.second;
    return true;
}

bool HybridNew::runTestStatistics(RooWorkspace *w, RooStats::ModelConfig *mc_s, RooStats::ModelConfig *mc_b, RooAbsData &data, double &limit, double &limitErr, const double *hint) {
    RooRealVar *r = dynamic_cast<RooRealVar *>(mc_s->GetParametersOfInterest()->first());
    HybridNew::Setup setup;
    std::auto_ptr<RooStats::HybridCalculator> hc(create(w, mc_s, mc_b, data, rValue_, setup));
    RooArgSet nullPOI(*setup.modelConfig_bonly.GetSnapshot());
    if (testStat_ == "LHC" || testStat_ == "Profile") nullPOI.setRealValue(r->GetName(), rValue_);
    limit = -2 * setup.qvar->Evaluate(data, nullPOI);
    if (testStat_ == "LHC" || testStat_ == "Profile") limit = -limit; // there's a sign flip for these two
    std::cout << "\n -- Hybrid New -- \n";
    std::cout << "-2 ln Q_{"<< testStat_<<"} = " << limit << std::endl;
    return true;
}

std::pair<double, double> HybridNew::eval(RooWorkspace *w, RooStats::ModelConfig *mc_s, RooStats::ModelConfig *mc_b, RooAbsData &data, double rVal, bool adaptive, double clsTarget) {
    if (readHybridResults_) {
        std::auto_ptr<RooStats::HypoTestResult> result(readToysFromFile(rVal));
        std::pair<double, double> ret(-1,-1);
        if (result.get() == 0) { 
            std::cerr << "HypoTestResults for r = " << rVal << " not found in file" << std::endl;
        } else {
            ret.first  = CLs_ ? result->CLs()      : result->CLsplusb();
            ret.second = CLs_ ? result->CLsError() : result->CLsplusbError();
        }
        return ret;
    }

    HybridNew::Setup setup;
    RooRealVar *r = dynamic_cast<RooRealVar *>(mc_s->GetParametersOfInterest()->first());
    r->setVal(rVal);
    if (verbose) std::cout << "  " << r->GetName() << " = " << rVal << " +/- " << r->getError() << std::endl;
    std::auto_ptr<RooStats::HybridCalculator> hc(create(w, mc_s, mc_b, data, rVal, setup));
    std::pair<double, double> ret = eval(*hc, rVal, adaptive, clsTarget);

    // add to plot 
    if (limitPlot_.get()) { 
        limitPlot_->Set(limitPlot_->GetN()+1);
        limitPlot_->SetPoint(limitPlot_->GetN()-1, rVal, ret.first); 
        limitPlot_->SetPointError(limitPlot_->GetN()-1, 0, ret.second);
    }

    return ret;
}



std::auto_ptr<RooStats::HybridCalculator> HybridNew::create(RooWorkspace *w, RooStats::ModelConfig *mc_s, RooStats::ModelConfig *mc_b, RooAbsData &data, double rVal, HybridNew::Setup &setup) {
  using namespace RooStats;
  
  w->loadSnapshot("clean");
  // realData_ = &data;  

  RooArgSet  poi(*mc_s->GetParametersOfInterest());
  RooRealVar *r = dynamic_cast<RooRealVar *>(poi.first());
 
  r->setMax(rVal); r->setVal(rVal); 
  if (testStat_ == "LHC" || testStat_ == "Profile") {
    r->setConstant(false); r->setMin(0); 
    if (workingMode_ == MakeSignificance) {
        rVal = 0;
        r->setVal(0);
        r->removeMax(); 
    }
  } else {
    r->setConstant(true);
  }


  std::auto_ptr<RooFitResult> fitMu, fitZero;
  if (fitNuisances_ && mc_s->GetNuisanceParameters() && withSystematics) {
    TStopwatch timer;
    r->setConstant(true); r->setVal(0);
    timer.Start();
    {
        CloseCoutSentry sentry(verbose < 3);
        fitZero.reset(mc_s->GetPdf()->fitTo(data, RooFit::Save(1), RooFit::Minimizer("Minuit2","minimize"), RooFit::Strategy(0), RooFit::Hesse(0), RooFit::Constrain(*mc_s->GetNuisanceParameters())));
    }
    if (verbose > 1) { std::cout << "Zero signal fit" << std::endl; fitZero->Print("V"); }
    if (verbose > 1) { std::cout << "Fitting of the background hypothesis done in " << timer.RealTime() << " s" << std::endl; }
    r->setConstant(true); r->setVal(rVal);
    timer.Start();
    {
       CloseCoutSentry sentry(verbose < 3);
       fitMu.reset(mc_s->GetPdf()->fitTo(data, RooFit::Save(1), RooFit::Minimizer("Minuit2","minimize"), RooFit::Strategy(0), RooFit::Hesse(0), RooFit::Constrain(*mc_s->GetNuisanceParameters())));
    }
    if (verbose > 1) { std::cout << "Reference signal fit (r = " << rVal << ")" << std::endl; fitMu->Print("V"); }
    if (verbose > 1) { std::cout << "Fitting of the signal-plus-background hypothesis done in " << timer.RealTime() << " s" << std::endl; }
  } else { fitNuisances_ = false; }

  // since ModelConfig cannot allow re-setting sets, we have to re-make everything 
  setup.modelConfig = ModelConfig("HybridNew_mc_s", w);
  setup.modelConfig.SetPdf(*mc_s->GetPdf());
  setup.modelConfig.SetObservables(*mc_s->GetObservables());
  setup.modelConfig.SetParametersOfInterest(*mc_s->GetParametersOfInterest());
  if (withSystematics) {
    if (genNuisances_ && mc_s->GetNuisanceParameters()) setup.modelConfig.SetNuisanceParameters(*mc_s->GetNuisanceParameters());
    if (genGlobalObs_ && mc_s->GetGlobalObservables())  setup.modelConfig.SetGlobalObservables(*mc_s->GetGlobalObservables());
    // if (genGlobalObs_ && mc_s->GetGlobalObservables())  snapGlobalObs_.reset(mc_s->GetGlobalObservables()->snapshot()); 
  }

  setup.modelConfig_bonly = ModelConfig("HybridNew_mc_b", w);
  setup.modelConfig_bonly.SetPdf(*mc_b->GetPdf());
  setup.modelConfig_bonly.SetObservables(*mc_b->GetObservables());
  setup.modelConfig_bonly.SetParametersOfInterest(*mc_b->GetParametersOfInterest());
  if (withSystematics) {
    if (genNuisances_ && mc_b->GetNuisanceParameters()) setup.modelConfig_bonly.SetNuisanceParameters(*mc_b->GetNuisanceParameters());
    if (genGlobalObs_ && mc_b->GetGlobalObservables())  setup.modelConfig_bonly.SetGlobalObservables(*mc_b->GetGlobalObservables());
  }

  if (withSystematics && !genNuisances_) {
      // The pdf will contain non-const parameters which are not observables
      // and the HybridCalculator will assume they're nuisances and try to generate them
      // to avoid this, we force him to generate a fake nuisance instead
      if (w->var("__HybridNew_fake_nuis__") == 0) { 
        w->factory("__HybridNew_fake_nuis__[0.5,0,1]");
        w->factory("Uniform::__HybridNew_fake_nuisPdf__(__HybridNew_fake_nuis__)");
      }
      setup.modelConfig.SetNuisanceParameters(RooArgSet(*w->var("__HybridNew_fake_nuis__")));
      setup.modelConfig_bonly.SetNuisanceParameters(RooArgSet(*w->var("__HybridNew_fake_nuis__")));
  }
 
  // create snapshots
  RooArgSet poiZero; 
  poiZero.addClone(*r); 
  poiZero.setRealValue(r->GetName(), 0);
  if (fitNuisances_) poi.add(fitMu->floatParsFinal());
  if (fitNuisances_) poiZero.addClone(fitZero->floatParsFinal());
  setup.modelConfig.SetSnapshot(poi);
  setup.modelConfig_bonly.SetSnapshot(poiZero);
  (const_cast<RooArgSet&>(*setup.modelConfig.GetSnapshot())) = poi;           // make sure we reset the values
  (const_cast<RooArgSet&>(*setup.modelConfig_bonly.GetSnapshot())) = poiZero; // make sure we reset the values

  // Create pdfs without nusiances terms, can be used for LEP tests statistics and
  // for generating toys when not generating global observables
  RooAbsPdf *factorizedPdf_s = setup.modelConfig.GetPdf(), *factorizedPdf_b = setup.modelConfig_bonly.GetPdf();
  if (withSystematics && optimizeProductPdf_ && !genGlobalObs_) {
        RooArgList constraints;
        RooAbsPdf *factorizedPdf_s = utils::factorizePdf(*mc_s->GetObservables(), *mc_s->GetPdf(), constraints);
        RooAbsPdf *factorizedPdf_b = utils::factorizePdf(*mc_b->GetObservables(), *mc_b->GetPdf(), constraints);
        setup.cleanupList.addOwned(*factorizedPdf_s);
        setup.cleanupList.addOwned(*factorizedPdf_b);
        setup.modelConfig.SetPdf(*factorizedPdf_s);
        setup.modelConfig_bonly.SetPdf(*factorizedPdf_b);
  }

  if (testStat_ == "LEP") {
      //SLR is evaluated using the central value of the nuisance parameters, so I believe we have to put them in the snapshots
      RooArgSet snapS; snapS.addClone(poi); 
      if (withSystematics) snapS.addClone(*mc_s->GetNuisanceParameters());
      RooArgSet snapB; snapB.addClone(snapS);
      snapS.setRealValue(r->GetName(), rVal);
      snapB.setRealValue(r->GetName(),    0);
      if (optimizeTestStatistics_) {
          if (!mc_s->GetPdf()->canBeExtended()) {
              setup.qvar.reset(new SimplerLikelihoodRatioTestStat(*factorizedPdf_s,*factorizedPdf_s, snapB, snapS));
          } else {
              setup.qvar.reset(new SimplerLikelihoodRatioTestStatOpt(*mc_s->GetObservables(), *factorizedPdf_s, *factorizedPdf_s, snapB, snapS));
          }
      } else {
          setup.qvar.reset(new SimpleLikelihoodRatioTestStat(*factorizedPdf_b,*factorizedPdf_s));
          ((SimpleLikelihoodRatioTestStat&)*setup.qvar).SetNullParameters(snapB); // Null is B
          ((SimpleLikelihoodRatioTestStat&)*setup.qvar).SetAltParameters(snapS);
      }
  } else if (testStat_ == "TEV") {
      if (optimizeTestStatistics_) {
          setup.qvar.reset(new ProfiledLikelihoodRatioTestStatOpt(*mc_s->GetObservables(), *mc_s->GetPdf(), *mc_s->GetPdf(), mc_s->GetNuisanceParameters(), poiZero, poi));
          ((ProfiledLikelihoodRatioTestStatOpt&)*setup.qvar).setPrintLevel(verbose-1);
      } else {   
          setup.qvar.reset(new RatioOfProfiledLikelihoodsTestStat(*mc_s->GetPdf(), *mc_s->GetPdf(), setup.modelConfig.GetSnapshot()));
          ((RatioOfProfiledLikelihoodsTestStat&)*setup.qvar).SetSubtractMLE(false);
      }
  } else if (testStat_ == "LHC" || testStat_ == "Profile") {
    setup.qvar.reset(new ProfileLikelihoodTestStat(*mc_s->GetPdf()));
    if (testStat_ == "LHC") {
       ((ProfileLikelihoodTestStat&)*setup.qvar).SetOneSided(true);
        if (optimizeTestStatistics_) {
            setup.qvar.reset(new ProfiledLikelihoodTestStatOpt(*mc_s->GetObservables(), *mc_s->GetPdf(), mc_s->GetNuisanceParameters(),  *setup.modelConfig.GetSnapshot(), verbose-1));
        }
    }
  }
  
  setup.toymcsampler.reset(new ToyMCSampler(*setup.qvar, nToys_));

  if (!mc_b->GetPdf()->canBeExtended()) setup.toymcsampler->SetNEventsPerToy(1);
  
  if (nCpu_ > 0) {
    if (verbose > 1) std::cout << "  Will use " << nCpu_ << " CPUs." << std::endl;
    setup.pc.reset(new ProofConfig(*w, nCpu_, "", kFALSE)); 
    setup.toymcsampler->SetProofConfig(setup.pc.get());
  }   

  std::auto_ptr<HybridCalculator> hc(new HybridCalculator(data,setup.modelConfig, setup.modelConfig_bonly, setup.toymcsampler.get()));
  if (genNuisances_ && withSystematics) {
    RooAbsPdf *nuisancePdf = utils::makeNuisancePdf(*mc_s);
    hc->ForcePriorNuisanceNull(*nuisancePdf);
    hc->ForcePriorNuisanceAlt(*nuisancePdf);
    setup.cleanupList.addOwned(*nuisancePdf);
  } else if (genGlobalObs_ && !genNuisances_) {
      hc->ForcePriorNuisanceNull(*w->pdf("__HybridNew_fake_nuisPdf__"));
      hc->ForcePriorNuisanceAlt(*w->pdf("__HybridNew_fake_nuisPdf__"));
  }

  // we need less B toys than S toys
  if (workingMode_ == MakeSignificance) {
      // need only B toys. just keep a few S+B ones to avoid possible divide-by-zero errors somewhere
      hc->SetToys(nToys_, int(0.01*nToys_)+1);
  } else if (!CLs_) {
      // we need only S+B toys to compute CLs+b
      hc->SetToys(int(0.01*nToys_)+1, nToys_);
  } else {
      // need both, but more S+B than B 
      hc->SetToys(int(0.25*nToys_)+1, nToys_);
  }

  static const char * istr = "__HybridNew__importanceSamplingDensity";
  if(importanceSamplingNull_) {
    if(verbose > 1) std::cout << ">>> Enabling importance sampling for null hyp." << std::endl;
    if(!withSystematics) {
      throw std::invalid_argument("Importance sampling is not available without systematics");
    }
    RooArgSet importanceSnapshot;
    importanceSnapshot.addClone(poi);
    importanceSnapshot.addClone(*mc_s->GetNuisanceParameters());
    if (verbose > 2) importanceSnapshot.Print("V");
    hc->SetNullImportanceDensity(mc_b->GetPdf(), &importanceSnapshot);
  }
  if(importanceSamplingAlt_) {
    if(verbose > 1) std::cout << ">>> Enabling importance sampling for alt. hyp." << std::endl;
    if(!withSystematics) {
      throw std::invalid_argument("Importance sampling is not available without systematics");
    }
    if (w->pdf(istr) == 0) {
      w->factory("__oneHalf__[0.5]");
      RooAddPdf *sum = new RooAddPdf(istr, "fifty-fifty", *mc_s->GetPdf(), *mc_b->GetPdf(), *w->var("__oneHalf__"));
      w->import(*sum); 
    }
    RooArgSet importanceSnapshot;
    importanceSnapshot.addClone(poi);
    importanceSnapshot.addClone(*mc_s->GetNuisanceParameters());
    if (verbose > 2) importanceSnapshot.Print("V");
    hc->SetAltImportanceDensity(w->pdf(istr), &importanceSnapshot);
  }

  return hc;
}

std::pair<double,double> 
HybridNew::eval(RooStats::HybridCalculator &hc, double rVal, bool adaptive, double clsTarget) {
    std::auto_ptr<HypoTestResult> hcResult(evalGeneric(hc));
    if (hcResult.get() == 0) {
        std::cerr << "Hypotest failed" << std::endl;
        return std::pair<double, double>(-1,-1);
    }
    if (testStat_ == "LHC" || testStat_ == "Profile") {
        // I need to flip the P-values
        hcResult->SetPValueIsRightTail(!hcResult->GetPValueIsRightTail());
        hcResult->SetTestStatisticData(hcResult->GetTestStatisticData()-1e-9); // issue with < vs <= in discrete models
    } else {
        hcResult->SetTestStatisticData(hcResult->GetTestStatisticData()+1e-9); // issue with < vs <= in discrete models
    }
    double clsMid    = (CLs_ ? hcResult->CLs()      : hcResult->CLsplusb());
    double clsMidErr = (CLs_ ? hcResult->CLsError() : hcResult->CLsplusbError());
    if (verbose) std::cout << (CLs_ ? "\tCLs = " : "\tCLsplusb = ") << clsMid << " +/- " << clsMidErr << std::endl;
    if (adaptive) {
        hc.SetToys(CLs_ ? nToys_ : 1, 4*nToys_);
        while (clsMidErr >= clsAccuracy_ && (clsTarget == -1 || fabs(clsMid-clsTarget) < 3*clsMidErr) ) {
            std::auto_ptr<HypoTestResult> more(evalGeneric(hc));
            if (testStat_ == "LHC" || testStat_ == "Profile") more->SetPValueIsRightTail(!more->GetPValueIsRightTail());
            hcResult->Append(more.get());
            clsMid    = (CLs_ ? hcResult->CLs()      : hcResult->CLsplusb());
            clsMidErr = (CLs_ ? hcResult->CLsError() : hcResult->CLsplusbError());
            if (verbose) std::cout << (CLs_ ? "\tCLs = " : "\tCLsplusb = ") << clsMid << " +/- " << clsMidErr << std::endl;
        }
    } else if (iterations_ > 1) {
        for (unsigned int i = 1; i < iterations_; ++i) {
            std::auto_ptr<HypoTestResult> more(evalGeneric(hc));
            if (testStat_ == "LHC" || testStat_ == "Profile") more->SetPValueIsRightTail(!more->GetPValueIsRightTail());
            hcResult->Append(more.get());
            clsMid    = (CLs_ ? hcResult->CLs()      : hcResult->CLsplusb());
            clsMidErr = (CLs_ ? hcResult->CLsError() : hcResult->CLsplusbError());
            if (verbose) std::cout << (CLs_ ? "\tCLs = " : "\tCLsplusb = ") << clsMid << " +/- " << clsMidErr << std::endl;
        }
    }
    if (verbose > 0) {
        std::cout <<
            "\tCLs      = " << hcResult->CLs()      << " +/- " << hcResult->CLsError()      << "\n" <<
            "\tCLb      = " << hcResult->CLb()      << " +/- " << hcResult->CLbError()      << "\n" <<
            "\tCLsplusb = " << hcResult->CLsplusb() << " +/- " << hcResult->CLsplusbError() << "\n" <<
            std::endl;
    }
    perf_totalToysRun_ += (hcResult->GetAltDistribution()->GetSize() + hcResult->GetNullDistribution()->GetSize());

    if (!plot_.empty() && workingMode_ != MakeLimit) {
        HypoTestPlot plot(*hcResult, 30);
        TCanvas *c1 = new TCanvas("c1","c1");
        plot.Draw();
        c1->Print(plot_.c_str());
        delete c1;
    }
    if (saveHybridResult_) {
        TString name = TString::Format("HypoTestResult_r%g_%u", rVal, RooRandom::integer(std::numeric_limits<UInt_t>::max() - 1));
        writeToysHere->WriteTObject(new HypoTestResult(*hcResult), name);
        if (verbose) std::cout << "Hybrid result saved as " << name << " in " << writeToysHere->GetFile()->GetName() << " : " << writeToysHere->GetPath() << std::endl;
    }

    return std::pair<double, double>(clsMid, clsMidErr);
} 

RooStats::HypoTestResult * HybridNew::evalGeneric(RooStats::HybridCalculator &hc, bool noFork) {
    if (fork_ && !noFork) return evalWithFork(hc);
    return hc.GetHypoTest();
}

RooStats::HypoTestResult * HybridNew::evalWithFork(RooStats::HybridCalculator &hc) {
    TStopwatch timer;
    std::auto_ptr<RooStats::HypoTestResult> result(0);
    char *tmpfile = tempnam(NULL,"rstat");
    unsigned int ich = 0;
    std::vector<UInt_t> newSeeds(fork_);
    for (ich = 0; ich < fork_; ++ich) {
        newSeeds[ich] = RooRandom::integer(std::numeric_limits<UInt_t>::max()-1);
        if (!fork()) break; // spawn children (but only in the parent thread)
    }
    if (ich == fork_) { // if i'm the parent
        int cstatus, ret;
        do {
            do { ret = waitpid(-1, &cstatus, 0); } while (ret == -1 && errno == EINTR);
        } while (ret != -1);
        if (ret == -1 && errno != ECHILD) throw std::runtime_error("Didn't wait for child");
        for (ich = 0; ich < fork_; ++ich) {
            TFile *f = TFile::Open(TString::Format("%s.%d.root", tmpfile, ich));
            if (f == 0) throw std::runtime_error(TString::Format("Child didn't leave output file %s.%d.root", tmpfile, ich).Data());
            RooStats::HypoTestResult *res = dynamic_cast<RooStats::HypoTestResult *>(f->Get("result"));
            if (res == 0)  throw std::runtime_error(TString::Format("Child output file %s.%d.root is corrupted", tmpfile, ich).Data());
            if (result.get()) result->Append(res); else result.reset(dynamic_cast<RooStats::HypoTestResult *>(res->Clone()));
            f->Close();
            unlink(TString::Format("%s.%d.root",    tmpfile, ich).Data());
            unlink(TString::Format("%s.%d.out.txt", tmpfile, ich).Data());
            unlink(TString::Format("%s.%d.err.txt", tmpfile, ich).Data());
        }
    } else {
        RooRandom::randomGenerator()->SetSeed(newSeeds[ich]); 
        freopen(TString::Format("%s.%d.out.txt", tmpfile, ich).Data(), "w", stdout);
        freopen(TString::Format("%s.%d.err.txt", tmpfile, ich).Data(), "w", stderr);
        std::cout << " I'm child " << ich << ", seed " << newSeeds[ich] << std::endl;
        RooStats::HypoTestResult *hcResult = evalGeneric(hc, /*noFork=*/true);
        TFile *f = TFile::Open(TString::Format("%s.%d.root", tmpfile, ich), "RECREATE");
        f->WriteTObject(hcResult, "result");
        f->ls();
        f->Close();
        fflush(stdout); fflush(stderr);
        std::cout << "And I'm done" << std::endl;
        throw std::runtime_error("done"); // I have to throw instead of exiting, otherwise there's no proper stack unwinding
                                          // and deleting of intermediate objects, and when the statics get deleted it crashes
                                          // in 5.27.06 (but not in 5.28)
    }
    free(tmpfile);
    if (verbose > 1) { std::cout << "      Evaluation of p-values done in  " << timer.RealTime() << " s" << std::endl; }
    return result.release();
}

#if 0
/// Another implementation of frequentist toy tossing without RooStats.
/// Can use as a cross-check if needed

RooStats::HypoTestResult * HybridNew::evalFrequentist(RooStats::HybridCalculator &hc) {
    int toysSB = (workingMode_ == MakeSignificance ? 1 : nToys_);
    int toysB  = (workingMode_ == MakeSignificance ? nToys_ : (CLs_ ? nToys_/4+1 : 1));
    RooArgSet obs(*hc.GetAlternateModel()->GetObservables());
    RooArgSet nuis(*hc.GetAlternateModel()->GetNuisanceParameters());
    RooArgSet gobs(*hc.GetAlternateModel()->GetGlobalObservables());
    RooArgSet nullPoi(*hc.GetNullModel()->GetSnapshot());
    std::auto_ptr<RooAbsCollection> parS(hc.GetAlternateModel()->GetPdf()->getParameters(obs));
    std::auto_ptr<RooAbsCollection> parB(hc.GetNullModel()->GetPdf()->getParameters(obs));
    RooArgList constraintsS, constraintsB;
    RooAbsPdf *factorS = hc.GetAlternateModel()->GetPdf();
    RooAbsPdf *factorB = hc.GetNullModel()->GetPdf();
    //std::auto_ptr<RooAbsPdf> factorS(utils::factorizePdf(obs, *hc.GetAlternateModel()->GetPdf(),  constraintsS));
    //std::auto_ptr<RooAbsPdf> factorB(utils::factorizePdf(obs, *hc.GetNullModel()->GetPdf(), constraintsB));
    std::auto_ptr<RooAbsPdf> nuisPdf(utils::makeNuisancePdf(const_cast<RooStats::ModelConfig&>(*hc.GetAlternateModel())));
    std::vector<Double_t> distSB, distB;
    *parS = *snapGlobalObs_;
    Double_t tsData = hc.GetTestStatSampler()->GetTestStatistic()->Evaluate(*realData_, nullPoi);
    if (verbose > 2) std::cout << "Test statistics on data: " << tsData << std::endl;
    for (int i = 0; i < toysSB; ++i) {
       // Initialize parameters to snapshot
       *parS = *hc.GetAlternateModel()->GetSnapshot(); 
       // Throw global observables, and set them
       if (verbose > 2) { std::cout << "Generating global obs starting from " << std::endl; parS->Print("V"); }
       std::auto_ptr<RooDataSet> gdata(nuisPdf->generate(gobs, 1));
       *parS = *gdata->get(0);
       if (verbose > 2) { std::cout << "Generated global obs" << std::endl; utils::printRAD(&*gdata); }
       // Throw observables
       if (verbose > 2) { std::cout << "Generating obs starting from " << std::endl; parS->Print("V"); }
       std::auto_ptr<RooDataSet> data(factorS->generate(obs, RooFit::Extended()));
       if (verbose > 2) { std::cout << "Generated obs" << std::endl; utils::printRAD(&*data); }
       // Evaluate T.S.
       distSB.push_back(hc.GetTestStatSampler()->GetTestStatistic()->Evaluate(*data, nullPoi));
       //if std::cout << "Test statistics on S+B : " << distSB.back() << std::endl;
    }
    for (int i = 0; i < toysB; ++i) {
       // Initialize parameters to snapshot
       *parB = *hc.GetNullModel()->GetSnapshot();
       //*parB = *hc.GetAlternateModel()->GetSnapshot(); 
       // Throw global observables, and set them
       if (verbose > 2) { std::cout << "Generating global obs starting from " << std::endl; parB->Print("V"); }
       std::auto_ptr<RooDataSet> gdata(nuisPdf->generate(gobs, 1));
       *parB = *gdata->get(0);
       if (verbose > 2) { std::cout << "Generated global obs" << std::endl; utils::printRAD(&*gdata); }
       // Throw observables
       if (verbose > 2) { std::cout << "Generating obs starting from " << std::endl; parB->Print("V"); }
       std::auto_ptr<RooDataSet> data(factorB->generate(obs, RooFit::Extended()));
       if (verbose > 2) { std::cout << "Generated obs" << std::endl; utils::printRAD(&*data); }
       // Evaluate T.S.
       distB.push_back(hc.GetTestStatSampler()->GetTestStatistic()->Evaluate(*data, nullPoi));
       //std::cout << "Test statistics on B   : " << distB.back() << std::endl;
    }
    // Load reference global observables
    RooStats::HypoTestResult *ret = new RooStats::HypoTestResult();
    ret->SetTestStatisticData(tsData);
    ret->SetAltDistribution(new RooStats::SamplingDistribution("sb","sb",distSB));
    ret->SetNullDistribution(new RooStats::SamplingDistribution("b","b",distB));
    return ret;
}
#endif

RooStats::HypoTestResult * HybridNew::readToysFromFile(double rValue) {
    if (!readToysFromHere) throw std::logic_error("Cannot use readHypoTestResult: option toysFile not specified, or input file empty");
    TDirectory *toyDir = readToysFromHere->GetDirectory("toys");
    if (!toyDir) throw std::logic_error("Cannot use readHypoTestResult: empty toy dir in input file empty");
    if (verbose) std::cout << "Reading toys for r = " << rValue << std::endl;
    TString prefix = TString::Format("HypoTestResult_r%g_",rValue);
    std::auto_ptr<RooStats::HypoTestResult> ret;
    TIter next(toyDir->GetListOfKeys()); TKey *k;
    while ((k = (TKey *) next()) != 0) {
        if (TString(k->GetName()).Index(prefix) != 0) continue;
        RooStats::HypoTestResult *toy = dynamic_cast<RooStats::HypoTestResult *>(toyDir->Get(k->GetName()));
        if (toy == 0) continue;
        if (verbose) std::cout << " - " << k->GetName() << std::endl;
        if (ret.get() == 0) {
            ret.reset(new RooStats::HypoTestResult(*toy));
        } else {
            ret->Append(toy);
        }
    }

    if (verbose > 0) {
        std::cout <<
            "\tCLs      = " << ret->CLs()      << " +/- " << ret->CLsError()      << "\n" <<
            "\tCLb      = " << ret->CLb()      << " +/- " << ret->CLbError()      << "\n" <<
            "\tCLsplusb = " << ret->CLsplusb() << " +/- " << ret->CLsplusbError() << "\n" <<
            std::endl;
    }
    return ret.release();
}

void HybridNew::readAllToysFromFile() {
    if (!readToysFromHere) throw std::logic_error("Cannot use readHypoTestResult without grid and without option toysFile");
    TDirectory *toyDir = readToysFromHere->GetDirectory("toys");
    if (!toyDir) throw std::logic_error("Cannot use readHypoTestResult: empty toy dir in input file empty");
    readGrid(toyDir);
}

void HybridNew::readGrid(TDirectory *toyDir) {
    clearGrid();

    TIter next(toyDir->GetListOfKeys()); TKey *k;
    while ((k = (TKey *) next()) != 0) {
        TString name(k->GetName());
        if (name.Index("HypoTestResult_r") != 0 || name.Index("_", name.Index("_")+1) == -1) continue;
        name.ReplaceAll("HypoTestResult_r","");
        name.Remove(name.Index("_"),name.Length());
        double rVal = atof(name.Data());
        if (verbose > 2) std::cout << "  Do " << k->GetName() << " -> " << name << " --> " << rVal << std::endl;
        RooStats::HypoTestResult *toy = dynamic_cast<RooStats::HypoTestResult *>(toyDir->Get(k->GetName()));
        RooStats::HypoTestResult *&merge = grid_[rVal];
        if (merge == 0) merge = new RooStats::HypoTestResult(*toy);
        else merge->Append(toy);
        merge->ResetBit(1);
    }
}
void HybridNew::updateGridData(RooWorkspace *w, RooStats::ModelConfig *mc_s, RooStats::ModelConfig *mc_b, RooAbsData &data, bool smart, double clsTarget_) {
    typedef std::map<double, RooStats::HypoTestResult *>::iterator point;
    if (!smart) {
        for (point it = grid_.begin(), ed = grid_.end(); it != ed; ++it) {
            it->second->ResetBit(1);
            updateGridPoint(w, mc_s, mc_b, data, it);
        }
    } else {
        typedef std::pair<double,double> CLs_t;
        std::vector<point> points; points.reserve(grid_.size()); 
        std::vector<CLs_t> values; values.reserve(grid_.size());
        for (point it = grid_.begin(), ed = grid_.end(); it != ed; ++it) { points.push_back(it); values.push_back(CLs_t(-99, -99)); }
        int iMin = 0, iMax = points.size()-1;
        while (iMax-iMin > 3) {
            if (verbose > 1) std::cout << "Bisecting range [" << iMin << ", " << iMax << "]" << std::endl; 
            int iMid = (iMin+iMax)/2;
            CLs_t clsMid = values[iMid] = updateGridPoint(w, mc_s, mc_b, data, points[iMid]);
            if (verbose > 1) std::cout << "    Midpoint " << iMid << " value " << clsMid.first << " +/- " << clsMid.second << std::endl; 
            if (clsMid.first - 3*max(clsMid.second,0.01) > clsTarget_) { 
                if (verbose > 1) std::cout << "    Replacing Min" << std::endl; 
                iMin = iMid; continue;
            } else if (clsMid.first + 3*max(clsMid.second,0.01) < clsTarget_) {
                if (verbose > 1) std::cout << "    Replacing Max" << std::endl; 
                iMax = iMid; continue;
            } else {
                if (verbose > 1) std::cout << "    Tightening Range" << std::endl; 
                while (iMin < iMid-1) {
                    int iLo = (iMin+iMid)/2;
                    CLs_t clsLo = values[iLo] = updateGridPoint(w, mc_s, mc_b, data, points[iLo]);
                    if (verbose > 1) std::cout << "        Lowpoint " << iLo << " value " << clsLo.first << " +/- " << clsLo.second << std::endl; 
                    if (clsLo.first - 3*max(clsLo.second,0.01) > clsTarget_) iMin = iLo; 
                    else break;
                }
                while (iMax > iMid+1) {
                    int iHi = (iMax+iMid)/2;
                    CLs_t clsHi = values[iHi] = updateGridPoint(w, mc_s, mc_b, data, points[iHi]);
                    if (verbose > 1) std::cout << "        Highpoint " << iHi << " value " << clsHi.first << " +/- " << clsHi.second << std::endl; 
                    if (clsHi.first + 3*max(clsHi.second,0.01) < clsTarget_) iMax = iHi; 
                    else break;
                }
                break;
            }
        }
        if (verbose > 1) std::cout << "Final range [" << iMin << ", " << iMax << "]" << std::endl; 
        for (int i = 0; i < iMin; ++i) {
            points[i]->second->SetBit(1);
            if (verbose > 1) std::cout << "  Will not use point " << i << " (r " << points[i]->first << ")" << std::endl;
        }
        for (int i = iMin; i <= iMax; ++i) {
            points[i]->second->ResetBit(1);
            if (values[i].first < -2) {
                if (verbose > 1) std::cout << "   Updaing point " << i << " (r " << points[i]->first << ")" << std::endl; 
                updateGridPoint(w, mc_s, mc_b, data, points[i]);
            }
            else if (verbose > 1) std::cout << "   Point " << i << " (r " << points[i]->first << ") was already updated during search." << std::endl; 
        }
        for (int i = iMax+1, n = points.size(); i < n; ++i) {
            points[i]->second->SetBit(1);
            if (verbose > 1) std::cout << "  Will not use point " << i << " (r " << points[i]->first << ")" << std::endl;
        }
    }
}
std::pair<double,double> HybridNew::updateGridPoint(RooWorkspace *w, RooStats::ModelConfig *mc_s, RooStats::ModelConfig *mc_b, RooAbsData &data, std::map<double, RooStats::HypoTestResult *>::iterator point) {
    typedef std::pair<double,double> CLs_t;
    bool isProfile = (testStat_ == "LHC" || testStat_ == "Profile");
    if (point->first == 0 && CLs_) return std::pair<double,double>(1,0);
    if (expectedFromGrid_) {
        std::vector<Double_t> btoys = point->second->GetNullDistribution()->GetSamplingDistribution();
        std::sort(btoys.begin(), btoys.end());
        Double_t testStat = btoys[std::min<int>(floor((1.-quantileForExpectedFromGrid_) * btoys.size()+0.5), btoys.size())];
        point->second->SetTestStatisticData(testStat + (isProfile ? -1e-9 : 1e-9));
    } else {
        RooArgSet  poi(*mc_s->GetParametersOfInterest());
        RooRealVar *r = dynamic_cast<RooRealVar *>(poi.first());
        Setup setup;
        std::auto_ptr<RooStats::HybridCalculator> hc = create(w, mc_s, mc_b, data, point->first, setup);
        RooArgSet nullPOI(*setup.modelConfig_bonly.GetSnapshot());
        if (isProfile) nullPOI.setRealValue(r->GetName(), rValue_);
        double testStat = setup.qvar->Evaluate(data, nullPOI);
        point->second->SetTestStatisticData(testStat + (isProfile ? -1e-9 : 1e-9));
    }
    return CLs_ ? CLs_t(point->second->CLs(), point->second->CLsError()) : CLs_t(point->second->CLsplusb(), point->second->CLsplusbError());
}
void HybridNew::useGrid() {
    typedef std::pair<double,double> CLs_t;
    int i = 0, n = 0;
    limitPlot_.reset(new TGraphErrors(1));
    for (std::map<double, RooStats::HypoTestResult *>::iterator itg = grid_.begin(), edg = grid_.end(); itg != edg; ++itg, ++i) {
        if (itg->second->TestBit(1)) continue;
        CLs_t val(1,0);
        if (CLs_) {
            if (itg->first > 0) val =CLs_t(itg->second->CLs(), itg->second->CLsError());
        } else {
            val = CLs_t(itg->second->CLsplusb(), itg->second->CLsplusbError());
        }
        limitPlot_->Set(n+1);
        limitPlot_->SetPoint(     n, itg->first, val.first); 
        limitPlot_->SetPointError(n, 0,          val.second);
        n++;
    }
}
void HybridNew::clearGrid() {
    for (std::map<double, RooStats::HypoTestResult *>::iterator it = grid_.begin(), ed = grid_.end(); it != ed; ++it) {
        delete it->second;
    }
    grid_.clear();
}
