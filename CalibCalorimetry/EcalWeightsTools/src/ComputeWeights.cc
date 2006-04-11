/* \file ComputeWeights.cc
 *
 *  $Date: 2005/09/05 10:13:49 $
 *  $Revision: 1.1 $
 *  \author R. Bruneliere - CERN
 *
 * $Date: 2006/10/04 15:07:31 $
 * $Revision: 1.2 $
 * Updated by Alex Zabi.
 */

#include "CalibCalorimetry/EcalWeightsTools/interface/ComputeWeights.h"

#include <iostream>
#include <iomanip>

// Constructor
ComputeWeights::ComputeWeights(int verbosity, 
			       bool doFitBaseline, bool doFitTime, 
			       int nPulseSamples, int nPrePulseSamples) :
  verbosity_(verbosity), doFitBaseline_(doFitBaseline),
  doFitTime_(doFitTime), nPulseSamples_(nPulseSamples),
  nPrePulseSamples_(nPrePulseSamples)
{
  if (verbosity_) {
    std::cout << "ComputeWeights::ComputeWeights: Constructing with setup:"
	      << std::endl;
    if (doFitBaseline_)
      std::cout << "  - baseline weights are computed" << std::endl;
    if (doFitTime_)
      std::cout << "  - time jitter weights are computed" << std::endl;
    std::cout << "  - the number of samples used to extract amplitude in the"
	      << " pulse is " << nPulseSamples_ << std::endl; 
    if (doFitBaseline_)
      std::cout << "  - the number of samples used to extract baseline in the"
		<< " pre-pulse is " << nPrePulseSamples_ << std::endl; 
    std::cout << std::endl;
  }
}//CONSTRUCTOR

// Destructor
ComputeWeights::~ComputeWeights()
{
  if (verbosity_)
    std::cout << "ComputeWeights::~ComputeWeights: Destructing ComputeWeights"
	      << std::endl;
}

// Compute weigths from an input pulse shape
bool ComputeWeights::compute(const std::vector<double>& pulseShape,
			     const std::vector<double>& pulseShapeDerivative)
{
  int nSamples = pulseShape.size();
  int nParams = 1 + int(doFitBaseline_) + int(doFitTime_);

  // Check if nSamples is large enough
  if (nSamples < nPulseSamples_ || doFitBaseline_ && 
      nSamples < (nPulseSamples_ + nPrePulseSamples_)) {
    std::cout << "ComputeWeights::compute: Error: nSamples = "
	      << nSamples << " is too small" << std::endl;
    return false;
  }//check samples

  // INITIALIZE WEIGHTS MATRICES
  if (weights_.num_row() != nSamples) {
    weights_ = HepMatrix(nSamples, nSamples, 0); // Fill matrices with zeros
    chi2_ = HepSymMatrix(nSamples, 0);
  } else {
    for (int iColumn = 0; iColumn < nSamples; iColumn++) {
      for (int iRow = 0; iRow < nParams; iRow++)
	weights_[iRow][iColumn] = 0.;
      for (int iRow = 0; iRow < nSamples; iRow++)
	chi2_[iRow][iColumn] = 0.;
    }
  }

  // LOOK FOR MAXIMUM
  double pulseMax = 0.;
  int sampleMax = 0;
  for (int iSample = 0; iSample < nSamples; iSample++)
    if (pulseShape[iSample] > pulseMax) {
      pulseMax = pulseShape[iSample];
      sampleMax = iSample;
    }//check max
  if (pulseMax == 0.) {
    std::cout << "ComputeWeights::compute: Warning: could not compute"
	      << " weights because max = 0." << std::endl;
    return false;
  }//check max

  //LOOKING FOR SAMPLE FIRST
  int firstSample = sampleMax - 2;
  double sumffMax = 0.;
  for (int sampleFirst = sampleMax - 2; sampleFirst < nPulseSamples_ - nSamples; sampleFirst++) {
    double sumff = 0.;
    for (int iSample = sampleFirst; iSample < sampleFirst + nSamples; iSample++)
      sumff += pulseShape[iSample]*pulseShape[iSample];
    if (sumff > sumffMax) {
      sumffMax = sumff;
      firstSample = sampleFirst;
    }
  }// loop sample

  if (firstSample + nPulseSamples_ > nSamples) {
    if (verbosity_)
      std::cout << "ComputeWeights::compute: Warning: firstSample cannot be "
		<< firstSample << " because they are too few samples beyond."
		<< std::endl << "firstSample is set to "
		<< nSamples - nPulseSamples_ << std::endl;
    firstSample = nSamples - nPulseSamples_;
  }//check max samples considered

  // Fill coef matrix
  int size = nPulseSamples_;
  if (doFitBaseline_) size += nPrePulseSamples_;
  HepMatrix coef(size, nParams);
  for (int iRow = 0; iRow < nPulseSamples_; iRow++)
    for (int iColumn = 0; iColumn < nParams; iColumn++) {
      if (iColumn == 0)
	coef[iRow][iColumn] = pulseShape[firstSample + iRow];
      else if (iColumn == 1 && doFitBaseline_)
	coef[iRow][iColumn] = 1.;
      else // doFitTime_ || nParams == 3
	coef[iRow][iColumn] = pulseShapeDerivative[firstSample + iRow];
    }
  for (int iRow = nPulseSamples_; iRow < size; iRow++)
    for (int iColumn = 0; iColumn < nParams; iColumn++) {
      if (iColumn == 1)
	coef[iRow][iColumn] = 1.;
      else
	coef[iRow][iColumn] = 0.;
    }
  HepMatrix tCoef = coef.T(); // transpose coef

  // Covariance matrix
  HepSymMatrix invCov(size, 1); // By default, set it to identity (1)

  // Variance matrix = [tCoef * invCov * coef]^-1
  HepMatrix tCoeffInvCov = tCoef*invCov;
  HepMatrix variance = tCoeffInvCov*coef;
  int ierr;
  variance.invert(ierr);
  if (ierr) {
    std::cout << "ComputeWeights::compute: Error: impossible to invert "
	      << "variance matrix." << std::endl;
    std::cout << variance;
    return false;
  }//check inversion

  // Weights matrix = variance * tCoef * invCov
  HepMatrix variancetCoef = variance*tCoef;
  HepMatrix weights = variancetCoef*invCov;

  // Chi2 matrix = (1 - coef * weights)^T * invCov * (1 - coef * weights)
  HepMatrix delta = coef*weights;
  delta *= -1.;
  HepMatrix unit(size, size, 1);
  delta += unit;
  HepMatrix tDelta = delta.T();
  HepMatrix tDeltaInvCov = tDelta*invCov;
  HepMatrix chi2 = tDeltaInvCov*delta;

  // Copy matrices into class members
  for (int iColumn = 0; iColumn < nPulseSamples_; iColumn++) {
    for (int iRow = 0; iRow < nParams; iRow++)
      weights_[iRow][firstSample + iColumn] = weights[iRow][iColumn];
    for (int iRow = 0; iRow < nPulseSamples_; iRow++)
      chi2_[firstSample + iRow][firstSample + iColumn] = chi2[iRow][iColumn];
  }
  if (doFitBaseline_) {
    for (int iColumn = 0; iColumn < nPrePulseSamples_; iColumn++) {
      for (int iRow = 0; iRow < nParams; iRow++)
	weights_[iRow][iColumn] = weights[iRow][iColumn + nPulseSamples_];
      for (int iRow = 0; iRow < nPrePulseSamples_; iRow++)
	chi2_[iRow][iColumn] = chi2[iRow + nPulseSamples_] 
				   [iColumn + nPulseSamples_];
    }
  }

  return true;
}

// Get weight used to compute amplitude
double ComputeWeights::getAmpWeight(int iSample) const
{
  if (!weights_.num_row()) {
    if (verbosity_)
      std::cout << "ComputeWeights::getAmpWeight: Warning: you should call"
		<< " the method ComputeWeights::compute() before asking a" 
		<< " weight." << std::endl;
    return 0.;
  }
  if (iSample < 0 || iSample >= weights_.num_col()) {
    if (verbosity_)
      std::cout << "ComputeWeights::getAmpWeight: Warning: iSample = "
		<< iSample << " should inside the range [0, " 
		<< weights_.num_col() - 1 << "]" << std::endl;
    return 0.;
  }
  return weights_[0][iSample];
}//Get Amplitude Weights

// Get weight used to compute dynamic pedestal
double ComputeWeights::getPedWeight(int iSample) const
{
  if (weights_.num_row() < 2 || !doFitBaseline_) {
    if (verbosity_)
      std::cout << "ComputeWeights::getPedWeight: Warning: pedestal weights"
		<< " are computed only if doFitBaseline = true" << std::endl;
    return 0.;
  }
  if (iSample < 0 || iSample >= weights_.num_col()) {
    if (verbosity_)
      std::cout << "ComputeWeights::getPedWeight: Warning: iSample = "
		<< iSample << " should inside the range [0, " 
		<< weights_.num_col() - 1 << "]" << std::endl;
    return 0.;
  }
  return weights_[1][iSample];
}//Get Pedestal Weights

// Get weight used to compute time jitter
double ComputeWeights::getTimeWeight(int iSample) const
{
  if (weights_.num_row() < 3 && doFitBaseline_ || !doFitTime_) {
    if (verbosity_)
      std::cout << "ComputeWeights::getTimeWeight: Warning: time weights"
		<< " are computed only if doFitTime = true" << std::endl;
    return 0.;
  }
  if (iSample < 0 || iSample >= weights_.num_col()) {
    if (verbosity_)
      std::cout << "ComputeWeights::getTimeWeight: Warning: iSample = "
		<< iSample << " should inside the range [0, " 
		<< weights_.num_col() - 1 << "]" << std::endl;
    return 0.;
  }
  if (doFitBaseline_)
    return weights_[2][iSample];
  return weights_[1][iSample];
}//Get Time Weights

// Get chi2 matrix
double ComputeWeights::getChi2Matrix(int iSample1, int iSample2) const
{
  if (iSample1 < 0 || iSample1 >= chi2_.num_row() || 
      iSample2 < 0 || iSample2 >= chi2_.num_col()) {
    if (verbosity_)
      std::cout << "ComputeWeights::getChi2Matrix: Warning: iSample1 and "
		<< "iSample2 should inside the range [0, " 
		<< weights_.num_col() - 1 << "]" << std::endl;
    return 0.;
  }
  return chi2_[iSample1][iSample2];
}//Get chi2
