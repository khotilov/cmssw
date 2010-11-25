// -*- C++ -*-
//
// Package:     Candidates
// Class  :     FWCandidate3DProxyBuilder
// 
// Implementation:
//     <Notes on implementation>
//
// Original Author:  Colin Bernet
//         Created:  Fri May 28 15:58:19 CEST 2010
// $Id: FWPFCandidate3DProxyBuilder.cc,v 1.1 2010/05/28 16:22:30 amraktad Exp $
//

// system include files
#include "TEveTrack.h"
#include "TEveTrackPropagator.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"

// user include files
#include "Fireworks/Core/interface/FWSimpleProxyBuilderTemplate.h"
#include "Fireworks/Core/interface/FWEvePtr.h"
#include "Fireworks/Core/src/CmsShowMain.h"
#include "Fireworks/Core/interface/FWEventItem.h"
#include "Fireworks/ParticleFlow/interface/setTrackTypePF.h"

class FWPFCandidate3DProxyBuilder : public FWSimpleProxyBuilderTemplate<reco::PFCandidate>  {
      
public:
   FWPFCandidate3DProxyBuilder() {}
   virtual ~FWPFCandidate3DProxyBuilder();
   
   // ---------- const member functions ---------------------
   
   // ---------- static member functions --------------------
   
   // ---------- member functions ---------------------------
   REGISTER_PROXYBUILDER_METHODS();

private:
   FWPFCandidate3DProxyBuilder(const FWPFCandidate3DProxyBuilder&); // stop default
   
   const FWPFCandidate3DProxyBuilder& operator=(const FWPFCandidate3DProxyBuilder&); // stop default

   void build(const reco::PFCandidate& iData, unsigned int iIndex, TEveElement& oItemHolder, const FWViewContext*);

};


FWPFCandidate3DProxyBuilder::~FWPFCandidate3DProxyBuilder()
{
}

void 
FWPFCandidate3DProxyBuilder::build(const reco::PFCandidate& iData, unsigned int iIndex, TEveElement& oItemHolder, const FWViewContext*) 
{
  TEveRecTrack t;
  t.fBeta = 1.;
  t.fP = TEveVector( iData.px(), iData.py(), iData.pz() );
  t.fV = TEveVector( iData.vertex().x(), iData.vertex().y(), iData.vertex().z() );
  t.fSign = iData.charge();
  TEveTrack* trk = new TEveTrack(&t, context().getTrackPropagator());

  trk->MakeTrack();

  fireworks::setTrackTypePF( iData, trk ); 
  setupAddElement(trk, &oItemHolder);
}

//
// static member functions
//
REGISTER_FWPROXYBUILDER(FWPFCandidate3DProxyBuilder, reco::PFCandidate,"PF Candidates", FWViewType::kAll3DBits | FWViewType::kAllRPZBits );
