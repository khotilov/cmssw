#ifndef BackgroundHandler_h
#define BackgroundHandler_h

#include "MuonAnalysis/MomentumScaleCalibration/interface/Functions.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include <vector>
#include "MuonAnalysis/MomentumScaleCalibration/interface/MuScleFitUtils.h"

/**
 * This class is used to handle the different background functions for the
 * different regions. <br>
 * It uses the backgroundFunctions defined in Functions.h and the
 * backgroundFunctionService defined in Functions.cc. <br>
 * More details are in the description of backgroundFunctionBase in Functions.h. <br>
 * <br>
 * A bool selects if to use the regions functions or the resonances functions.
 */
class BackgroundHandler
{
public:
  BackgroundHandler( const std::vector<int> & identifiers,
                     const std::vector<double> & leftWindowFactors,
                     const std::vector<double> & rightWindowFactors,
                     const double * ResMass );
  ~BackgroundHandler();

  /// Returns the total number of parameters used for the regions
  inline int regionsParNum()
  {
    return parNumsResonances_[0];
  }

  void countEventsInBackgroundWindows(const std::vector<std::pair<reco::Particle::LorentzVector,reco::Particle::LorentzVector> > & muonPairs, const double & weight);

  /// Sets initial parameters for all the functions
  void setParameters(double* Start, double* Step, double* Mini, double* Maxi, int* ind, TString* parname, const std::vector<double> & parBgr, const std::vector<int> & parBgrOrder, const int muonType);

  /// returns true if the parameter is to be unlocked
  bool unlockParameter(const std::vector<int> & resfind, const unsigned int ipar);

  /// Returns the appropriate window factors depending on whether the background is being fitted and on the resonance
  std::pair<double, double> windowFactors( const bool doBackgroundFit, const int ires );

  /**
   * Returns the appropriate resMass value depending on whether the background is being fitted and on the resonance. <br>
   * The resMass used for the region is the mean of the mass of the corresponding resonances, so for the Z is the same Z mass,
   * for the Upsilons is the arithmetic mean of the Upsilon masses and the same for the J/Psi and Psi2S region.
   */
  double resMass( const bool doBackgroundFit, const int ires );

  /**
   * Computes the rescaled parameters from the regions functions to the
   * resonances functions. It takes into account the difference in intervals
   * and rescales the parameters so that the fraction of events is correctly accounter for. <br>
   * It uses the list of all muon pairs to compute the number of events in each resonance window.
   */
  void rescale( std::vector<double> & parBgr, const double * ResMass, const double massWindowHalfWidth[][3], const int muonType,
                const std::vector<std::pair<reco::Particle::LorentzVector,reco::Particle::LorentzVector> > & muonPairs, const double & weight = 1. );

  /**
   * Returns the background fraction parameter (parBgr[0], but shifted to the correct function) and
   * the value returned by the background function. <br>
   * Depending on the value of doBackgroundFit it returns the values for the regions or the resonances.
   */
  std::pair<double, double> backgroundFunction( const bool doBackgroundFit,
                                                const double * parval, const int resTotNum, const int ires,
                                                const bool * resConsidered, const double * ResMass, const double ResHalfWidth[],
                                                const int MuonType, const double & mass, const int nbins );
private:
  /// Performs the rescaling of parameters for a single resonance
  double applyRescale( TF1* backgroundFunctionForIntegral, const double backgroundWindowEvents, const double resonanceWindowEvents,
                       const double & leftRegionWidth, const double & rightRegionWidth,
                       const double & leftResonanceWidth, const double & rightResonanceWidth ) const;

  /// Used to check the consistency of passed parameters
  void consistencyCheck( const std::vector<int> & identifiers,
                         const std::vector<double> & leftWindowFactors,
                         const std::vector<double> & rightWindowFactors ) const throw(cms::Exception);

  // Correspondence between regions and halfWidths used:
  // - for the Upsilons region we use the Upsilon
  // - for the J/Psi and Psi2S region we use the J/Psi
  int regToResHW_[3];
  // Correspondence between resonances and regions:
  // - Z -> region 0
  // - Uspilons -> region 1
  // - J/Psi and Psi2S -> region 2
  int resToReg_[6];

  // Used in the shifts of the parval
  // Contains 0 as the first number and Sum_0^(ires-1)(parNum(i)) for the rest,
  // so that calling parNums[ires] returns the sub of the number of parameters
  // of the previous functions (0 if none) and allows to shift the parval to the
  // parameters of the actual function.
  int parNumsRegions_[3];
  // These start from the correct parameters (take into account that the parRegions are
  // before the parResonances).
  int parNumsResonances_[6];

  // Holds the mass values used as the center of each region.
  double resMassForRegion_[3];
  // Holds the mass values for the resonances
  double resMassForResonance_[6];

  std::vector<backgroundFunctionBase*> backgroundFunctionsForRegions_;
  std::vector<backgroundFunctionBase*> backgroundFunctionsForResonances_;

  // Using double because weights are taken into account and the sum of event_i*weight is a double
  std::vector<double> regionWindowEvents_;
  std::vector<double> resonanceWindowEvents_;
  std::vector<double> leftWindowFactors_;
  std::vector<double> rightWindowFactors_;
};

#endif // BackgroundHandler_h
