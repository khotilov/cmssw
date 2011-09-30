#ifndef HiggsAnalysis_CombinedLimit_ChannelCompatibilityCheck_h
#define HiggsAnalysis_CombinedLimit_ChannelCompatibilityCheck_h
/** \class ChannelCompatibilityCheck
 *
 * Do a ML fit of the data with background and signal+background hypothesis and print out diagnostics plots 
 *
 * \author Giovanni Petrucciani (UCSD)
 *
 *
 */
#include "../interface/LimitAlgo.h"
#include "../interface/ProfileLikelihood.h"

class ChannelCompatibilityCheck : public LimitAlgo {
public:
  ChannelCompatibilityCheck() ;
  virtual bool run(RooWorkspace *w, RooStats::ModelConfig *mc_s, RooStats::ModelConfig *mc_b, RooAbsData &data, double &limit, double &limitErr, const double *hint);
  virtual const std::string & name() const {
    static const std::string name("ChannelCompatibilityCheck");
    return name;
  }
  virtual void applyOptions(const boost::program_options::variables_map &vm) ;

protected:
  std::string nameForLabel(const char *label) ;

  static std::string minimizerAlgo_;
  static float       minimizerTolerance_;
  static int         minimizerStrategy_;

  static float mu_;
  static bool  fixedMu_;

  static bool saveFitResult_;

  static std::vector<std::string> groups_;
};


#endif
