#include "FWCore/Framework/interface/MakerMacros.h"

#include "PhysicsTools/RecoAlgos/interface/TrackFullCloneSelectorBase.h"
#include "DiffractiveForwardAnalysis/SingleDiffractiveWAnalysis/interface/TrackAssociatedWithPVSelector.h"

namespace reco { 
  namespace modules {

    typedef TrackFullCloneSelectorBase< ::TrackAssociatedWithPVSelector > TrackAssociatedWithPVSelector;

    DEFINE_FWK_MODULE(TrackAssociatedWithPVSelector);
  }
}
