#ifndef BTauReco_JetCrystalsAssociation_h
#define BTauReco_JetCrystalsAssociation_h
// \class JetCrystalsAssociation
// 
// \short association of Ecal Crystals to jet 
// 
//

#include "DataFormats/JetReco/interface/CaloJetCollection.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/Common/interface/Ref.h"
#include "DataFormats/Common/interface/RefProd.h"
#include "DataFormats/Common/interface/RefVector.h"
#include "DataFormats/Math/interface/LorentzVector.h"
#include "DataFormats/Math/interface/LorentzVectorFwd.h"
#include <vector>

namespace reco {
  typedef math::PtEtaPhiELorentzVector           EMLorentzVector;
  typedef math::PtEtaPhiELorentzVectorCollection EMLorentzVectorCollection;
  typedef math::PtEtaPhiELorentzVectorRef        EMLorentzVectorRef;
  typedef math::PtEtaPhiELorentzVectorRefProd    EMLorentzVectorRefProd;
  typedef math::PtEtaPhiELorentzVectorRefVector  EMLorentzVectorRefVector;
        
  typedef 
  std::vector<std::pair<edm::RefToBase<Jet>, EMLorentzVectorRefVector>  >
    JetCrystalsAssociationCollection;
    
  typedef
  JetCrystalsAssociationCollection::value_type JetCrystalsAssociation;
  
  typedef
  edm::Ref<JetCrystalsAssociationCollection> JetCrystalsAssociationRef;
  
  typedef
  edm::RefProd<JetCrystalsAssociationCollection> JetCrystalsAssociationRefProd;
  
  typedef
  edm::RefVector<JetCrystalsAssociationCollection> JetCrystalsAssociationRefVector; 
}
#endif
