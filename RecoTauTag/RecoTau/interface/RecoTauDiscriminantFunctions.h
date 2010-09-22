#ifndef RecoTauTag_RecoTau_RecoTauDiscriminantFunctions_h
#define RecoTauTag_RecoTau_RecoTauDiscriminantFunctions_h

/*
 * RecoTauDiscriminantFunctions
 *
 * Collection of unary functions used to compute tau discriminant values.
 * Each function here (may be) used in an MVA discriminator.
 *
 * The functions all have the form 
 *      ReturnType Function(const PFTau& tau)
 * where ReturnType is either vector<double> or double.
 *
 * Author: Evan K. Friis, UC Davis
 *
 * $Id $
 */

#include "DataFormats/TauReco/interface/PFTau.h"
#include <vector>

namespace reco { namespace tau { namespace disc {

// Save typing
typedef const PFTau& Tau;
typedef std::vector<double> VDouble;

/// For three prong events, take the track that has charge opposite to the
/// composite charge.
PFCandidateRef mainTrack(const PFTau& tau);

double Pt(Tau tau) { return tau.pt(); }
double Eta(Tau tau) { return tau.eta(); }
double Mass(Tau tau) { return tau.mass(); }
double DecayMode(Tau tau) { return tau.decayMode(); } 

// Number of objects in isolation cone
double OutlierN(Tau);  

// Number of charged objects in isolation cone
double OutlierNCharged(Tau);

double OutlierSumPt(Tau);
double ChargedOutlierSumPt(Tau);
double NeutralOutlierSumPt(Tau);

// Pt of the main track
double MainTrackPt(Tau);
// Eta of the main track
double MainTrackEta(Tau);
// Angle of main track to tau axis
double MainTrackAngle(Tau);

// Exactly the same as "Mass", needed for backwards compatability
double InvariantMassOfSignal(Tau tau) { return tau.mass(); }

// Quanitites of tracks 
VDouble TrackPt(Tau);
VDouble TrackAngle(Tau);
VDouble TrackEta(Tau);

// Quanitites of PiZeros 
VDouble PiZeroPt(Tau);
VDouble PiZeroAngle(Tau);
VDouble PiZeroEta(Tau);

// Isolation quantities
VDouble OutlierPt(Tau);
VDouble OutlierAngle(Tau);
VDouble ChargedOutlierPt(Tau);
VDouble ChargedOutlierAngle(Tau);
VDouble NeutralOutlierPt(Tau);
VDouble NeutralOutlierAngle(Tau);

// Dalitz for three prongs
VDouble Dalitz(Tau);

// Deprecated functions needed for backwards compatability
VDouble FilteredObjectPt(Tau);
VDouble GammaOccupancy(Tau);
VDouble GammaPt(Tau);
VDouble InvariantMassOfSignalWithFiltered(Tau);
VDouble InvariantMass(Tau);
VDouble OutlierMass(Tau);

}}} // end namespace reco::tau::disc
#endif 
