#include "DataFormats/MuonReco/interface/CaloMuon.h"
using namespace reco;

CaloMuon::CaloMuon() {
   energyValid_  = false;
   caloCompatibility_ = -9999.;
}

float CaloMuon::getCaloCompatibility() const 
{ 
   return caloCompatibility_; 
}
