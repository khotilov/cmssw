#ifndef Fireworks_Electrons_makeSuperCluster_h
#define Fireworks_Electrons_makeSuperCluster_h
// -*- C++ -*-
//
// Package:     Electrons
// Class  :     makeSuperCluster
//
/**\class makeSuperCluster makeSuperCluster.h Fireworks/Electrons/interface/makeSuperCluster.h

   Description: <one line class summary>

   Usage:
    <usage>

 */
//
// Original Author:  Chris Jones
//         Created:  Fri Dec  5 15:32:20 EST 2008
// $Id: makeSuperCluster.h,v 1.1 2008/12/05 20:57:15 chrjones Exp $
//

// system include files
#include "DataFormats/EgammaReco/interface/SuperCluster.h"

// user include files

// forward declarations
class FWEventItem;
class TEveElement;

namespace fireworks {
   bool makeRhoPhiSuperCluster(const FWEventItem&,
                               const reco::SuperClusterRef& iCluster,
                               float iPhi,
                               TEveElement& oItemHolder);
   bool makeRhoZSuperCluster(const FWEventItem&,
                             const reco::SuperClusterRef& iCluster,
                             float iPhi,
                             TEveElement& oItemHolder);
}
#endif
