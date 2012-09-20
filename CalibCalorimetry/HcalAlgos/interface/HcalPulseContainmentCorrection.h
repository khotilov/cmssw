#ifndef CALIBCALORIMETRY_HCALALGOS_HCALPULSECONTAINMENTCORRECTION_H
#define CALIBCALORIMETRY_HCALALGOS_HCALPULSECONTAINMENTCORRECTION_H 1

#include <map>

/** \class HcalPulseContainmentCorrection
  *
  * Amplitude correction for pulse containment in time.
  * Currently only for HPD pulse shape.
  *  
  * $Date: $
  * $Revision: $
  * \author P. Dudero - Minnesota
  */
class HcalPulseContainmentCorrection {
public:
  HcalPulseContainmentCorrection(int num_samples,
                                 float fixedphase_ns,
                                 float max_fracerror);

  double getCorrection(double fc_ampl) const;
  double fractionContained(double fc_ampl) const { return 1.0/this->getCorrection(fc_ampl); }

private:
  std::map<double,double> mCorFactors_;
};

#endif
