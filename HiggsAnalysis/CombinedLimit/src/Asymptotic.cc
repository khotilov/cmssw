#include <stdexcept>

#include "../interface/Asymptotic.h"
#include <RooRealVar.h>
#include <RooArgSet.h>
#include <RooAbsPdf.h>
#include <RooFitResult.h>
#include <RooRandom.h>
#include <RooMinimizer.h>
#include <RooStats/ModelConfig.h>
#include <Math/DistFuncMathCore.h>
#include "../interface/Combine.h"
#include "../interface/CloseCoutSentry.h"
#include "../interface/RooFitGlobalKillSentry.h"
#include "../interface/ProfiledLikelihoodRatioTestStatExt.h"
#include "../interface/ToyMCSamplerOpt.h"
#include "../interface/ProfileLikelihood.h"
#include "../interface/utils.h"

using namespace RooStats;

double Asymptotic::rAbsAccuracy_ = 0.0005;
double Asymptotic::rRelAccuracy_ = 0.005;
std::string Asymptotic::what_ = "both"; 
bool  Asymptotic::qtilde_ = true; 
float Asymptotic::quantileForExpected_ = 0.5; 
std::string Asymptotic::minimizerAlgo_ = "Minuit2";
float       Asymptotic::minimizerTolerance_ = 1e-3;
int         Asymptotic::minimizerStrategy_  = 1;

Asymptotic::Asymptotic() : 
LimitAlgo("Asymptotic specific options") {
    options_.add_options()
        ("rAbsAcc", boost::program_options::value<double>(&rAbsAccuracy_)->default_value(rAbsAccuracy_), "Absolute accuracy on r to reach to terminate the scan")
        ("rRelAcc", boost::program_options::value<double>(&rRelAccuracy_)->default_value(rRelAccuracy_), "Relative accuracy on r to reach to terminate the scan")
        ("run", boost::program_options::value<std::string>(&what_)->default_value(what_), "What to run: both (default), observed, expected.")
        ("minimizerAlgo",      boost::program_options::value<std::string>(&minimizerAlgo_)->default_value(minimizerAlgo_), "Choice of minimizer used for profiling (Minuit vs Minuit2)")
        ("minimizerTolerance", boost::program_options::value<float>(&minimizerTolerance_)->default_value(minimizerTolerance_),  "Tolerance for minimizer used for profiling")
        ("minimizerStrategy",  boost::program_options::value<int>(&minimizerStrategy_)->default_value(minimizerStrategy_),      "Stragegy for minimizer")
        ("qtilde", boost::program_options::value<bool>(&qtilde_)->default_value(qtilde_),  "Allow only non-negative signal strengths (default is true).")
    ;
}

void Asymptotic::applyOptions(const boost::program_options::variables_map &vm) {
    if (what_ != "observed" && what_ != "expected" && what_ != "both") 
        throw std::invalid_argument("Asymptotic: option 'run' can only be 'observed', 'expected' or 'both' (the default)");
}

void Asymptotic::applyDefaultOptions() { 
    what_ = "observed";
}

bool Asymptotic::run(RooWorkspace *w, RooStats::ModelConfig *mc_s, RooStats::ModelConfig *mc_b, RooAbsData &data, double &limit, double &limitErr, const double *hint) {
    RooFitGlobalKillSentry silence(verbose <= 1 ? RooFit::WARNING : RooFit::DEBUG);
    ProfileLikelihood::MinimizerSentry minimizerConfig(minimizerAlgo_, minimizerTolerance_);
    if (verbose > 0) std::cout << "Will use minimizer " << minimizerAlgo_ << " with strategy " << minimizerStrategy_ << " and tolerance " << minimizerTolerance_ << std::endl;

    bool ret = false;
    std::vector<std::pair<float,float> > expected;
    if (what_ != "observed") expected = runLimitExpected(w, mc_s, mc_b, data, limit, limitErr, hint);
    if (what_ != "expected") ret = runLimit(w, mc_s, mc_b, data, limit, limitErr, hint);

    if (verbose >= 0) {
        const char *rname = mc_s->GetParametersOfInterest()->first()->GetName();
        std::cout << "\n -- Asymptotic -- " << "\n";
        if (ret && what_ != "expected") printf("Observed Limit: %s < %6.4f\n", rname, limit);
        for (std::vector<std::pair<float,float> >::const_iterator it = expected.begin(), ed = expected.end(); it != ed; ++it) {
            printf("Expected %4.1f%%: %s < %6.4f\n", it->first*100, rname, it->second);
        }
        std::cout << std::endl;
    }
    return ret;
}

bool Asymptotic::runLimit(RooWorkspace *w, RooStats::ModelConfig *mc_s, RooStats::ModelConfig *mc_b, RooAbsData &data, double &limit, double &limitErr, const double *hint) {
  RooRealVar *r = dynamic_cast<RooRealVar *>(mc_s->GetParametersOfInterest()->first()); 
  w->loadSnapshot("clean");
  RooAbsData &asimov = *asimovDataset(w, mc_s, mc_b, data);

  r->setConstant(false);
  r->setMin(qtilde_ ? 0 : -r->getMax());
 
  if (params_.get() == 0) params_.reset(mc_s->GetPdf()->getParameters(data));
  RooArgSet constraints; if (withSystematics) constraints.add(*mc_s->GetNuisanceParameters());
  nllD_.reset(mc_s->GetPdf()->createNLL(data,   RooFit::Constrain(constraints)));
  nllA_.reset(mc_s->GetPdf()->createNLL(asimov, RooFit::Constrain(constraints)));

  if (verbose > 0) std::cout << (qtilde_ ? "Restricting" : "Not restricting") << " " << r->GetName() << " to positive values." << std::endl;
  if (verbose > 1) params_->Print("V");
 
  if (verbose > 0) std::cout << "\nMake global fit of real data" << std::endl;
  {
    CloseCoutSentry sentry(verbose < 3);
    *params_ = snapGlobalObsData;
    RooMinimizer minim(*nllD_);
    minim.setStrategy(minimizerStrategy_);
    minim.setPrintLevel(-1);
    nllutils::robustMinimize(*nllD_, minim, verbose-1);
    fitFreeD_.reset(minim.save());
  }
  if (verbose > 0) std::cout << "NLL at global minimum of data: " << fitFreeD_->minNll() << " (" << r->GetName() << " = " << r->getVal() << ")" << std::endl;

  r->setMin(0);

  if (verbose > 1) fitFreeD_->Print("V");
  if (verbose > 0) std::cout << "\nMake global fit of asimov data" << std::endl;
  {
    CloseCoutSentry sentry(verbose < 3);
    *params_ = snapGlobalObsAsimov;
    RooMinimizer minim(*nllA_);
    minim.setStrategy(minimizerStrategy_);
    minim.setPrintLevel(-1);
    nllutils::robustMinimize(*nllA_, minim, verbose-1);
    fitFreeA_.reset(minim.save());
  }
  if (verbose > 0) std::cout << "NLL at global minimum of asimov: " << fitFreeA_->minNll() << " (" << r->GetName() << " = " << r->getVal() << ")" << std::endl;
  if (verbose > 1) fitFreeA_->Print("V");

  *params_ = fitFreeD_->floatParsFinal();
  r->setConstant(true);

  double clsTarget = 1-cl;
  double rMin = std::max<double>(0, r->getVal()), rMax = rMin + 3 * ((RooRealVar*)fitFreeD_->floatParsFinal().find(r->GetName()))->getError();
  for (int tries = 0; tries < 5; ++tries) {
    double cls = getCLs(*r, rMax);
    if (cls == -999) { std::cerr << "Minimization failed in an unrecoverable way" << std::endl; return false; }
    if (cls < clsTarget) break;
    rMax *= 2;
  }
  do {
    limit = 0.5*(rMin + rMax); limitErr = 0.5*(rMax - rMin);
    double cls = getCLs(*r, limit);
    if (cls == -999) { std::cerr << "Minimization failed in an unrecoverable way" << std::endl; break; }
    if (cls > clsTarget) {
        rMin = limit;
    } else {
        rMax = limit;
    }
  } while (limitErr > std::max(rRelAccuracy_ * limit, rAbsAccuracy_));

  return true;
}

double Asymptotic::getCLs(RooRealVar &r, double rVal) {
  r.setMax(1.1 * rVal);
  r.setConstant(true);

  CloseCoutSentry sentry(verbose < 2);
  RooMinimizer minimD(*nllD_), minimA(*nllA_);
  minimD.setStrategy(minimizerStrategy_); minimD.setPrintLevel(-1); 
  minimA.setStrategy(minimizerStrategy_); minimA.setPrintLevel(-1);

  *params_ = fitFixD_.get() ? fitFixD_->floatParsFinal() : fitFreeD_->floatParsFinal();
  *params_ = snapGlobalObsData;
  r.setVal(rVal);
  r.setConstant(true);
  if (!nllutils::robustMinimize(*nllD_, minimD, verbose-1)) return -999;
  fitFixD_.reset(minimD.save());
  if (verbose >= 2) fitFixD_->Print("V");
  double qmu = 2*(nllD_->getVal() - fitFreeD_->minNll()); if (qmu < 0) qmu = 0;

  *params_ = fitFixA_.get() ? fitFixA_->floatParsFinal() : fitFreeA_->floatParsFinal();
  *params_ = snapGlobalObsAsimov;
  r.setVal(rVal);
  r.setConstant(true);
  if (!nllutils::robustMinimize(*nllA_, minimA, verbose-1)) return -999;
  fitFixA_.reset(minimA.save());
  if (verbose >= 2) fitFixA_->Print("V");
  double qA  = 2*(nllA_->getVal() - fitFreeA_->minNll()); if (qA < 0) qA = 0;

  double CLsb = ROOT::Math::normal_cdf_c(sqrt(qmu));
  double CLb  = ROOT::Math::normal_cdf(sqrt(qA)-sqrt(qmu));
  double CLs  = (CLb == 0 ? 0 : CLsb/CLb);
  if (qtilde_ && qmu > qA) {
    // In this region, things are tricky
    double mos = sqrt(2*qA); // mu/sigma
    CLsb = ROOT::Math::normal_cdf_c( (qmu + qA)/(2*mos) );
    CLb  = ROOT::Math::normal_cdf  ( (qmu - qA)/(2*mos) );
  }
  sentry.clear();
  if (verbose > 0) printf("At %s = %f:\tq_mu = %.5f\tq_A  = %.5f\tCLsb = %7.5f\tCLb  = %7.5f\tCLs  = %7.5f\n", r.GetName(), rVal, qmu, qA, CLsb, CLb, CLs);
  return CLs; 
}   

std::vector<std::pair<float,float> > Asymptotic::runLimitExpected(RooWorkspace *w, RooStats::ModelConfig *mc_s, RooStats::ModelConfig *mc_b, RooAbsData &data, double &limit, double &limitErr, const double *hint) {
    // See equation 35-38 of AN 2011/298 and references cited therein
    //
    //   (35)  sigma^2   = mu^2 / q_mu(Asimov)
    //   (38)  mu_median = sigma * normal_quantile(1-0.5*(1-cl))
    //
    // -->  q_mu(Asimov) = pow(normal_quantile(1-0.5*(1-cl)), 2)
    //      can be solved to find mu_median
    //
    // --> then (38) gives sigma, and the quantiles are given by (37)
    //      mu_N = sigma * (normal_quantile(1 - quantile*(1-cl), 1.0) + normal_quantile(quantile));
    //
    // 1) get parameter of interest
    RooArgSet  poi(*mc_s->GetParametersOfInterest());
    RooRealVar *r = dynamic_cast<RooRealVar *>(poi.first());

    // 2) get asimov dataset
    RooAbsData *asimov = asimovDataset(w, mc_s, mc_b, data);

    // 2b) load asimov global observables
    if (params_.get() == 0) params_.reset(mc_s->GetPdf()->getParameters(data));
    *params_ = snapGlobalObsAsimov;

    // 3) solve for q_mu
    CloseCoutSentry sentry(verbose < 3);
    r->setConstant(false);
    
    std::auto_ptr<RooAbsReal> nll(mc_s->GetPdf()->createNLL(*asimov, RooFit::Constrain(*mc_s->GetNuisanceParameters())));
    RooMinimizer minim(*nll);
    minim.setStrategy(0);
    minim.setPrintLevel(-1);
    minim.setErrorLevel(0.5*pow(ROOT::Math::normal_quantile(1-0.5*(1-cl),1.0), 2)); // the 0.5 is because qmu is -2*NLL
                        // eventually if cl = 0.95 this is the usual 1.92!
    nllutils::robustMinimize(*nll, minim, verbose-1);
    int minosStat = -1;
    for (int tries = 0; tries < 3; ++tries) {
        minosStat = minim.minos(RooArgSet(*r));
        if (minosStat != -1) break;
        minim.setStrategy(2);
        if (tries == 1) { 
            if (minimizerAlgo_.find("Minuit2") != std::string::npos) {
                minim.minimize("Minuit","minimize");
            } else {
                minim.minimize("Minuit2","minmize");
            }
        }
    }
    sentry.clear();
    if (minosStat == -1) {
        std::cerr << "Minos did not converge. No expected limit available" << std::endl;
        return std::vector<std::pair<float,float> >(); 
    }
    
    // 3) get ingredients for equation 37
    double median = r->getAsymErrorHi();
    double sigma  = median / ROOT::Math::normal_quantile(1-0.5*(1-cl),1.0);
    double alpha = 1-cl;

    std::vector<std::pair<float,float> > expected;
    const double quantiles[5] = { 0.025, 0.16, 0.50, 0.84, 0.975 };
    if (verbose >= 0) std::cout << "\n -- Expected Asymptotic -- " << "\n";
    for (int iq = 0; iq < 5; ++iq) {
        double N     = ROOT::Math::normal_quantile(quantiles[iq], 1.0);
        limit = sigma*(ROOT::Math::normal_quantile(1 - alpha * quantiles[iq], 1.0) + N);
        limitErr = 0;
        Combine::commitPoint(true, quantiles[iq]);
        expected.push_back(std::pair<float,float>(quantiles[iq], limit));
    }
    return expected;
}

RooAbsData * Asymptotic::asimovDataset(RooWorkspace *w, RooStats::ModelConfig *mc_s, RooStats::ModelConfig *mc_b, RooAbsData &data) {
    if (w->data("_Asymptotic_asimovDataset_") == 0) {
        // CREATE THE DATASET

        std::cout << "Generating Asimov dataset" << std::endl;
        RooArgSet  poi(*mc_s->GetParametersOfInterest());
        RooRealVar *r = dynamic_cast<RooRealVar *>(poi.first());
        r->setConstant(true); r->setVal(0);
        CloseCoutSentry sentry(verbose < 3);
        mc_s->GetPdf()->fitTo(data, RooFit::Minimizer("Minuit2","minimize"), RooFit::Strategy(1), RooFit::Constrain(*mc_s->GetNuisanceParameters()));
        toymcoptutils::SimPdfGenInfo newToyMC(*mc_b->GetPdf(), *mc_s->GetObservables(), false); 
        RooRealVar *weightVar = 0;
        RooAbsData *asimov = newToyMC.generateAsimov(weightVar); // as simple as that
        asimov->SetName("_Asymptotic_asimovDataset_");
        w->import(*asimov);
        delete weightVar;

        // NOW SNAPSHOT THE GLOBAL OBSERVABLES
        if (withSystematics && mc_s->GetGlobalObservables()) {
            RooArgSet gobs(*mc_s->GetGlobalObservables());
            // snapshot data global observables
            snapGlobalObsData.removeAll();
            utils::setAllConstant(gobs, true);
            gobs.snapshot(snapGlobalObsData);
            // now get the ones for the asimov dataset.
            // we do this by fitting the nuisance pdf with floating globa observables but fixed nuisances (which we call observables)
            // we clone the pdf, to avoid messing up the status of the nuisance parameters
            if (params_.get() == 0) params_.reset(mc_s->GetPdf()->getParameters(data));
            RooArgSet paramsSetToConstants;
            std::auto_ptr<TIterator> iter(params_->createIterator());
            for (RooAbsArg *a = (RooAbsArg *) iter->Next(); a != 0; a = (RooAbsArg *) iter->Next()) {
                RooRealVar *rrv = dynamic_cast<RooRealVar *>(a);
                if (rrv) {
                    if (gobs.find(rrv->GetName())) {
                        rrv->setConstant(false);
                    } else {
                        if (!rrv->isConstant()) paramsSetToConstants.add(*rrv);
                        rrv->setConstant(true);
                    }
                }
            }
            mc_s->GetPdf()->fitTo(*w->data("_Asymptotic_asimovDataset_"), RooFit::Minimizer("Minuit2","minimize"), RooFit::Strategy(1), RooFit::Constrain(gobs));
            // snapshot
            snapGlobalObsAsimov.removeAll();
            utils::setAllConstant(gobs, true);
            gobs.snapshot(snapGlobalObsAsimov);
            // revert things to normal
            gobs = snapGlobalObsData;
            utils::setAllConstant(paramsSetToConstants, false);
            utils::setAllConstant(gobs, true);
            if (verbose > 2) {
                sentry.clear();
                std::cout <<   "Global observables in data\n";   snapGlobalObsData.Print("V");
                std::cout << "\nGlobal observables in asimov\n"; snapGlobalObsAsimov.Print("V");
            }
        }
    }
    return w->data("_Asymptotic_asimovDataset_");
}
