#include "HiggsAnalysis/HiggsToTauTau/interface/PlotLimits.h"

PlotLimits::PlotLimits(const char* output, const edm::ParameterSet& cfg) : 
  output_(output),
  dataset_(cfg.getParameter<std::string>("dataset")),
  xaxis_(cfg.getParameter<std::string>("xaxis")),
  yaxis_(cfg.getParameter<std::string>("yaxis")),
  min_  (cfg.existsAs<double>("min" ) ? cfg.getParameter<double>("min" ) :    0.),
  max_  (cfg.existsAs<double>("max" ) ? cfg.getParameter<double>("max" ) :   -1.),
  log_  (cfg.existsAs<int   >("log" ) ? cfg.getParameter<int   >("log" ) :     0),
  png_  (cfg.existsAs<bool  >("png" ) ? cfg.getParameter<bool  >("png" ) : false),
  pdf_  (cfg.existsAs<bool  >("pdf" ) ? cfg.getParameter<bool  >("pdf" ) : false),
  txt_  (cfg.existsAs<bool  >("txt" ) ? cfg.getParameter<bool  >("txt" ) : false),
  root_ (cfg.existsAs<bool  >("root") ? cfg.getParameter<bool  >("root") : false),
  mssm_ (cfg.existsAs<bool  >("mssm") ? cfg.getParameter<bool  >("mssm") : false)
{
  // mass-points
  bins_=cfg.getParameter<std::vector<double> >("masspoints");
  for(unsigned int i=0; i<bins_.size(); ++i){ valid_.push_back(true); }

  // minima as function of mass (if available as such)
  if(cfg.existsAs<std::vector<edm::ParameterSet> >("xmins")){
    std::vector<edm::ParameterSet> mins = cfg.getParameter<std::vector<edm::ParameterSet> >("xmins");
    for(std::vector<edm::ParameterSet>::const_iterator m=mins.begin(); m!=mins.end(); ++m){
      double key = m->getParameter<double>("mass");
      double val = m->getParameter<double>("min");
      xmins_[key] = val;
    }
  }
  if(cfg.existsAs<std::vector<edm::ParameterSet> >("ymins")){
    std::vector<edm::ParameterSet> mins = cfg.getParameter<std::vector<edm::ParameterSet> >("ymins");
    for(std::vector<edm::ParameterSet>::const_iterator m=mins.begin(); m!=mins.end(); ++m){
      double key = m->getParameter<double>("mass");
      double val = m->getParameter<double>("min");
      ymins_[key] = val;
    }
  }
  // maxima as function of mass (if available as such)
  if(cfg.existsAs<std::vector<edm::ParameterSet> >("xmaxs")){
    std::vector<edm::ParameterSet> maxs = cfg.getParameter<std::vector<edm::ParameterSet> >("xmaxs");
    for(std::vector<edm::ParameterSet>::const_iterator m=maxs.begin(); m!=maxs.end(); ++m){
      double key = m->getParameter<double>("mass");
      double val = m->getParameter<double>("max");
      xmaxs_[key] = val;
    }
  }
  if(cfg.existsAs<std::vector<edm::ParameterSet> >("ymaxs")){
    std::vector<edm::ParameterSet> maxs = cfg.getParameter<std::vector<edm::ParameterSet> >("ymaxs");
    for(std::vector<edm::ParameterSet>::const_iterator m=maxs.begin(); m!=maxs.end(); ++m){
      double key = m->getParameter<double>("mass");
      double val = m->getParameter<double>("max");
      ymaxs_[key] = val;
    }
  }

  // common configs
  label_ = cfg.existsAs<std::string>("outputLabel") ? cfg.getParameter<std::string>("outputLabel") : std::string();
  verbosity_ = cfg.existsAs<unsigned int>("verbosity") ? cfg.getParameter<unsigned int>("verbosity") : 0,
  // specifics to plot signal strength
  signal_ = cfg.existsAs<std::string>("signal") ? cfg.getParameter<std::string>("signal") : std::string();
  // specifics to plot scan-2d
  model_ = cfg.existsAs<std::string>("model") ? cfg.getParameter<std::string>("model") : std::string();
  temp_ = cfg.existsAs<bool>("temp") ? cfg.getParameter<bool>("temp") : false;
  // specifics to plot signal strength
  bestfit_ = cfg.existsAs<bool>("bestfit") ? cfg.getParameter<bool>("bestfit") : false;
  // specifics to plot xsec limits
  injected_ = cfg.existsAs<bool>("injected") ? cfg.getParameter<bool>("injected") : false;
  BG_Higgs_ = cfg.existsAs<bool>("BG_Higgs") ? cfg.getParameter<bool>("BG_Higgs") : false;
  // specifics to plot MSSM mA-tanb limits
  higgs125_ =cfg.existsAs<bool>("higgs125" ) ? cfg.getParameter<bool>("higgs125" ) : false;
  outerband_=cfg.existsAs<bool>("outerband") ? cfg.getParameter<bool>("outerband") : false;
}

TGraph*
PlotLimits::fillCentral(const char* directory, TGraph* plot, const char* filename)
{
  std::vector<double> central;
  // fill pre-defined values from previous results
  if (std::string(filename).find("HIG")!=std::string::npos){
    bool initial = masses_.empty();
    for(unsigned int imass=0; imass<bins_.size(); ++imass){
      if(std::string(filename)==std::string("HIG-11-020-obs")){
	prepareHIG_11_020(central, "observed", false, bins_[imass], initial);
      }
      else if(std::string(filename)==std::string("HIG-11-020-exp")){
	prepareHIG_11_020(central, "expected", false, bins_[imass], initial);
      }
      else if(std::string(filename)==std::string("HIG-11-029-obs")){
	prepareHIG_11_029(central, "observed", bins_[imass], initial);
      }
      else if(std::string(filename)==std::string("HIG-11-029-exp")){
	prepareHIG_11_029(central, "expected", bins_[imass], initial);
      }
      else if(std::string(filename)==std::string("HIG-12-018-obs")){
	prepareHIG_12_018(central, "observed", bins_[imass], initial);
      }
      else if(std::string(filename)==std::string("HIG-12-018-exp")){
	prepareHIG_12_018(central, "expected", bins_[imass], initial);
      }
      else if(std::string(filename)==std::string("HIG-12-032-obs")){
	prepareHIG_12_032(central, "observed", bins_[imass], initial);
      }
      else if(std::string(filename)==std::string("HIG-12-032-exp")){
	prepareHIG_12_032(central, "expected", bins_[imass], initial);
      }
      else if(std::string(filename)==std::string("HIG-12-043-obs")){
	prepareHIG_12_043(central, "observed", bins_[imass], initial);
      }
      else if(std::string(filename)==std::string("HIG-12-043-exp")){
	prepareHIG_12_043(central, "expected", bins_[imass], initial);
      }
      else if(std::string(filename)==std::string("HIG-12-050-obs")){
	prepareHIG_12_050(central, "observed", bins_[imass], initial);
      }
      else if(std::string(filename)==std::string("HIG-12-050-exp")){
	prepareHIG_12_050(central, "expected", bins_[imass], initial);
      }
    }
  }
  else{
    if(std::string(filename)==std::string("MEDIAN") || std::string(filename)==std::string("MEAN")){
      prepareByToy(directory, central, filename);
    }
    else if(std::string(filename).find("obs.Asymptotic") != std::string::npos){
      prepareByValue(directory, central, "higgsCombine-obs.Asymptotic.mH$MASS", -1.);
    }
    else if(std::string(filename).find("exp.Asymptotic") != std::string::npos){
      prepareByValue(directory, central, "higgsCombine-exp.Asymptotic.mH$MASS", 0.5);
    }
    else if(std::string(filename).find("MaxLikelihood") != std::string::npos){
      prepareByValue(directory, central, filename, -1.);
    }
    else{
      prepareByFile(directory, central, filename);
    }
  }
  for(unsigned int imass=0, ipoint=0; imass<bins_.size(); ++imass){
    if(valid_[imass] && std::string(filename).find("HIG")==std::string::npos){
      plot->SetPoint(ipoint++, bins_[imass], central[imass]);
      if(verbosity_>1){ std::cout << "INFO: central [" << bins_[imass] << "] = " << central[imass] << "[" << (valid_[imass] ? "OK]" : "FAILED]") << std::endl; }
    }
    else{
      if(valid_[imass]){
	for(unsigned int jmass=0; jmass<masses_.size(); ++jmass){
	  if(masses_[jmass]==bins_[imass]){
	    plot->SetPoint(ipoint++, masses_[jmass], central[jmass]);
	    if(verbosity_>1){ std::cout << "INFO: central [" << masses_[jmass] << "] = " << central[jmass] << "[" << (valid_[imass] ? "OK]" : "FAILED]") << std::endl; }
	    break;
	  }
	}
      }
    }
  }
  return plot;
}

TGraphAsymmErrors*
PlotLimits::fillBand(const char* directory, TGraphAsymmErrors* plot, const char* method, bool innerBand)
{
  std::vector<double> upper, lower, expected;

  if (std::string(method).find("HIG")!=std::string::npos){
    bool initial = masses_.empty();
    for(unsigned int imass=0; imass<bins_.size(); ++imass){
      if(std::string(method) == std::string("HIG-11-020")){
	bool xsec = false;
	prepareHIG_11_020(expected, "expected"                       , xsec, bins_[imass], initial);
	prepareHIG_11_020(upper   , innerBand ? "+1sigma" : "+2sigma", xsec, bins_[imass], initial);
	prepareHIG_11_020(lower   , innerBand ? "-1sigma" : "-2sigma", xsec, bins_[imass], initial);
      }
      else if(std::string(method) == std::string("HIG-11-029")){
	prepareHIG_11_029(expected, "expected"                       , bins_[imass], initial);
	prepareHIG_11_029(upper   , innerBand ? "+1sigma" : "+2sigma", bins_[imass], initial);
	prepareHIG_11_029(lower   , innerBand ? "-1sigma" : "-2sigma", bins_[imass], initial);
      }
      else if(std::string(method) == std::string("HIG-12-018")){
	prepareHIG_12_018(expected, "expected"                       , bins_[imass], initial);
	prepareHIG_12_018(upper   , innerBand ? "+1sigma" : "+2sigma", bins_[imass], initial);
	prepareHIG_12_018(lower   , innerBand ? "-1sigma" : "-2sigma", bins_[imass], initial);
      }
      else if(std::string(method) == std::string("HIG-12-032")){
	prepareHIG_12_032(expected, "expected"                       , bins_[imass], initial);
	prepareHIG_12_032(upper   , innerBand ? "+1sigma" : "+2sigma", bins_[imass], initial);
	prepareHIG_12_032(lower   , innerBand ? "-1sigma" : "-2sigma", bins_[imass], initial);
      }
      else if(std::string(method) == std::string("HIG-12-043")){
	prepareHIG_12_043(expected, "expected"                       , bins_[imass], initial);
	prepareHIG_12_043(upper   , innerBand ? "+1sigma" : "+2sigma", bins_[imass], initial);
	prepareHIG_12_043(lower   , innerBand ? "-1sigma" : "-2sigma", bins_[imass], initial);
      }
      else if(std::string(method) == std::string("HIG-12-050")){
	prepareHIG_12_050(expected, "expected"                       , bins_[imass], initial);
	prepareHIG_12_050(upper   , innerBand ? "+1sigma" : "+2sigma", bins_[imass], initial);
	prepareHIG_12_050(lower   , innerBand ? "-1sigma" : "-2sigma", bins_[imass], initial);
      }
    }
  }
  else{
    if(std::string(method) == std::string("TOYBASED")){
      prepareByToy(directory, expected, "MEDIAN");
      prepareByToy(directory, upper   , innerBand ? "+1SIGMA" : "+2SIGMA");
      prepareByToy(directory, lower   , innerBand ? "-1SIGMA" : "-2SIGMA");
    }
    else if(std::string(method) == std::string("CLs")){
      prepareCLs(directory, expected, ".quant0.500");
      prepareCLs(directory, upper   , innerBand ? ".quant0.840" : ".quant0.975");
      prepareCLs(directory, lower   , innerBand ? ".quant0.160" : ".quant0.027");
    }
    else if(std::string(method).find("Asymptotic") != std::string::npos){
      prepareByValue(directory, expected, method, 0.50);
      prepareByValue(directory, upper   , method, innerBand ? 0.84 : 0.975);
      prepareByValue(directory, lower   , method, innerBand ? 0.16 : 0.025);
    }
    else if(std::string(method).find("MaxLikelihood") != std::string::npos){
      prepareByValue(directory, expected, method, -1);
      prepareByValue(directory, upper, method, 0.84);
      prepareByValue(directory, lower, method, 0.16);
    }
    else{
      std::cout << "ERROR: chose wrong method to fill uncertainty band: " << method << std::endl;
      exit(1);
    }
  }
  for(unsigned int imass=0, ipoint=0; imass<bins_.size(); ++imass){
    if(valid_[imass] && std::string(method).find("HIG")==std::string::npos){
      plot->SetPoint(ipoint, bins_[imass], expected[imass]);
      plot->SetPointEYhigh(ipoint, upper[imass] - expected[imass]);
      plot->SetPointEYlow (ipoint, expected[imass] - lower[imass]);
      ++ipoint; // only add valid mass points to the TGraph
      if(verbosity_>1){
	std::cout << "INFO: Calculating " << (innerBand ? "inner" : "outer") << " uncertainty band" << std::endl;
	std::cout << "INFO: upper    [" << bins_[imass] << "] = " << " (" << upper[imass] << ") " << "[" << (valid_[imass] ? "OK]" : "FAILED]") << std::endl;
	std::cout << "INFO: expected [" << bins_[imass] << "] = " << expected[imass]          << "[" << (valid_[imass] ? "OK]" : "FAILED]") << std::endl;
	std::cout << "INFO: lower    [" << bins_[imass] << "] = " << " (" << lower[imass] << ") " << "[" << (valid_[imass] ? "OK]" : "FAILED]") << std::endl;
      }
    }
    else{
      if(valid_[imass]){
	for(unsigned int jmass=0; jmass<masses_.size(); ++jmass){
	  if(masses_[jmass]==bins_[imass]){
	    plot->SetPoint(ipoint, masses_[jmass], expected[jmass]);
	    plot->SetPointEYhigh(ipoint, upper[jmass] - expected[jmass]);
	    plot->SetPointEYlow (ipoint, expected[jmass] - lower[jmass]);
	    ++ipoint; // only add valid mass points to the TGraph
	    if(verbosity_>1){
	      std::cout << "INFO: Calculating " << (innerBand ? "inner" : "outer") << " uncertainty band" << std::endl;
	      std::cout << "INFO: upper    [" << masses_[jmass] << "] = " << " (" << upper[jmass] << ") " << "[" << (valid_[imass] ? "OK]" : "FAILED]") << std::endl;
	      std::cout << "INFO: expected [" << masses_[jmass] << "] = " << expected[jmass]          << "[" << (valid_[imass] ? "OK]" : "FAILED]") << std::endl;
	      std::cout << "INFO: lower    [" << masses_[jmass] << "] = " << " (" << lower[jmass] << ") " << "[" << (valid_[imass] ? "OK]" : "FAILED]") << std::endl;
	    }
	    break;
	  }
	}
      }
    }
  }
  return plot;
}

float
PlotLimits::maximum(TGraph* graph)
{
  float maximum = -1;
  if(max_<0){
    for(int ibin=0; ibin<graph->GetN(); ++ibin){
      if(maximum<0 || graph->GetY()[ibin]>maximum){
	maximum=graph->GetY()[ibin];
      }
    }
    maximum*=2.5;
  }
  else{
    maximum=max_;
  }
  std::cout << "M A X I M U M   I S : " << maximum << std::endl;
  return maximum;
}
