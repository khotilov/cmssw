#include "RecoEgamma/EgammaElectronAlgos/interface/ElectronMomentumCorrector.h"

#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/GsfTrackReco/interface/GsfTrack.h"
#include "DataFormats/Math/interface/LorentzVector.h"

#include "TrackingTools/GsfTools/interface/MultiGaussianState1D.h"
#include "TrackingTools/GsfTools/interface/GaussianSumUtilities1D.h"
#include "TrackingTools/GsfTools/interface/MultiGaussianStateTransform.h"
#include "TrackingTools/TrajectoryState/interface/TrajectoryStateOnSurface.h"


#include "FWCore/MessageLogger/interface/MessageLogger.h"

/****************************************************************************
 *
 * Class based E-p combination for the final electron momentum. It relies on
 * the electron classification and on the dtermination of track momentum and ecal
 * supercluster energy errors. The track momentum error is taken from the gsf fit. 
 * The ecal supercluster energy error is taken from a class dependant parametrisation 
 * of the energy resolution. 
 *
 *
 * \author Federico Ferri - INFN Milano, Bicocca university
 * \author Ivica Puljak - FESB, Split
 * \author Stephanie Baffioni - Laboratoire Leprince-Ringuet - �cole polytechnique, CNRS/IN2P3
 *
 * \version $Id: ElectronMomentumCorrector.cc,v 1.11 2009/03/25 02:15:43 charlot Exp $
 *
 ****************************************************************************/


/** The electron classification.
   barrel  :   0: golden,  10: bigbrem,  20: narrow, 30-34: showering,
             (30: showering nbrem=0, 31: showering nbrem=1, 32: showering nbrem=2 ,33: showering nbrem=3, 34: showering nbrem>=4)
              40: crack, 41: eta gaps, 42: phi gaps
   endcaps : 100: golden, 110: bigbrem, 120: narrow, 130-134: showering
            (130: showering nbrem=0, 131: showering nbrem=1, 132: showering nbrem=2 ,133: showering nbrem=3, 134: showering nbrem>=4)
             140: crack
*/
void ElectronMomentumCorrector::correct(reco::GsfElectron &electron, TrajectoryStateOnSurface & vtxTsos) {

  if (electron.isMomentumCorrected())
   {
    edm::LogWarning("ElectronMomentumCorrector::correct")<<"already done" ;
	return ;
   }

  newMomentum_ = electron.p4() ; // default
  int elClass = electron.classification() ;

  // irrelevant classification
  if ( (elClass <= reco::GsfElectron::UNKNOWN) ||
	   (elClass>reco::GsfElectron::GAP) )
   {
	edm::LogWarning("ElectronMomentumCorrector::correct")<<"unexpected classification" ;
	return ;
   }

  float scEnergy = electron.ecalEnergy() ;
  errorEnergy_ = electron.ecalEnergyError() ;

  float trackMomentum  = electron.trackMomentumAtVtx().R() ;
  errorTrackMomentum_ = 999. ;

  // retreive momentum error 
  MultiGaussianState1D qpState(MultiGaussianStateTransform::multiState1D(vtxTsos,0));
  GaussianSumUtilities1D qpUtils(qpState);
  errorTrackMomentum_ = trackMomentum*trackMomentum*sqrt(qpUtils.mode().variance());

  float finalMomentum = electron.p4().t(); // initial
  float finalMomentumError = 999.;

  // calculate E/p and corresponding error
  float eOverP = scEnergy / trackMomentum;
  float errorEOverP = sqrt(
		   (errorEnergy_/trackMomentum)*(errorEnergy_/trackMomentum) +
		   (scEnergy*errorTrackMomentum_/trackMomentum/trackMomentum)*
		   (scEnergy*errorTrackMomentum_/trackMomentum/trackMomentum));

  if ( eOverP  > 1 + 2.5*errorEOverP )
   {
    finalMomentum = scEnergy; finalMomentumError = errorEnergy_;
    if ((elClass==reco::GsfElectron::GOLDEN) && electron.isEB() && (eOverP<1.15))
     {
	  if (scEnergy<15) {finalMomentum = trackMomentum ; finalMomentumError = errorTrackMomentum_;}
     }
   }
  else if ( eOverP < 1 - 2.5*errorEOverP )
   {
    finalMomentum = scEnergy; finalMomentumError = errorEnergy_;
    if (elClass==reco::GsfElectron::SHOWERING)
     {
      if (electron.isEB())
       {
	    if(scEnergy<18) {finalMomentum = trackMomentum; finalMomentumError = errorTrackMomentum_;}
       }
      else if (electron.isEE())
       {
	    if(scEnergy<13) {finalMomentum = trackMomentum; finalMomentumError = errorTrackMomentum_;}
       }
      else
       { edm::LogWarning("ElectronMomentumCorrector::correct")<<"nor barrel neither endcap electron ?!" ; }
     }
    else if (electron.isGap())
     {
	  if(scEnergy<60) {finalMomentum = trackMomentum; finalMomentumError = errorTrackMomentum_;}
     }
   }
  else 
   {
    // combination
    finalMomentum = (scEnergy/errorEnergy_/errorEnergy_ + trackMomentum/errorTrackMomentum_/errorTrackMomentum_) /
                       (1/errorEnergy_/errorEnergy_ + 1/errorTrackMomentum_/errorTrackMomentum_);
    float finalMomentumVariance = 1 / (1/errorEnergy_/errorEnergy_ + 1/errorTrackMomentum_/errorTrackMomentum_);
    finalMomentumError = sqrt(finalMomentumVariance);
   } 

  math::XYZTLorentzVector oldMomentum = electron.p4() ;
  newMomentum_ = math::XYZTLorentzVector
   ( oldMomentum.x()*finalMomentum/oldMomentum.t(),
     oldMomentum.y()*finalMomentum/oldMomentum.t(),
     oldMomentum.z()*finalMomentum/oldMomentum.t(),
     finalMomentum ) ;

  // final set
  electron.correctMomentum(newMomentum_,errorTrackMomentum_,finalMomentumError);
 }

