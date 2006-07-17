#ifndef L1Trigger_L1JetParticleFwd_h
#define L1Trigger_L1JetParticleFwd_h
// -*- C++ -*-
//
// Package:     L1Trigger
// Class  :     L1JetParticleFwd
// 
/**\class L1JetParticleFwd L1JetParticleFwd.h DataFormats/L1Trigger/interface/L1JetParticleFwd.h

 Description: <one line class summary>

 Usage:
    <usage>

*/
//
// Original Author:  Werner Sun
//         Created:  Sat Jul 15 14:28:43 EDT 2006
// $Id$
//

// system include files

// user include files

// forward declarations
#include <vector>
#include "DataFormats/Common/interface/Ref.h"
#include "DataFormats/Common/interface/RefVector.h"


namespace level1 {

   class L1JetParticle ;

   typedef std::vector< L1JetParticle > L1JetParticleCollection ;

   typedef edm::Ref< L1JetParticleCollection > L1JetParticleRef ;
   typedef edm::RefVector< L1JetParticleCollection > L1JetParticleRefVector ;
}

#endif
