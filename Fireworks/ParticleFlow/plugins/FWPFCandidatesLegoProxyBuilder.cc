// -*- C++ -*-
//
// Package: Fireworks
// Class  : FWPFCandidatesLegoProxyBuilder

/*

 Description: [one line class summary]

 Usage:
    <usage>

*/
//
// Original Author: Colin Bernet
//         Created: Fri May 28 14:54:08 2010 
//
//


#include "Fireworks/Core/interface/FWSimpleProxyBuilderTemplate.h"
#include "Fireworks/Core/interface/Context.h"
#include "Fireworks/Core/interface/FWEventItem.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "Fireworks/ParticleFlow/interface/FWLegoEvePFCandidate.h"
#include "Fireworks/ParticleFlow/src/FWPFScale.h"
#include "Fireworks/ParticleFlow/interface/setTrackTypePF.h"

// forward declarations

class FWPFCandidatesLegoProxyBuilder : public FWSimpleProxyBuilderTemplate<reco::PFCandidate> {
public:
   FWPFCandidatesLegoProxyBuilder();
   virtual ~FWPFCandidatesLegoProxyBuilder();

   // ---------- const member functions ---------------------

   // ---------- static member functions --------------------

   // ---------- member functions ---------------------------
   virtual bool havePerViewProduct(FWViewType::EType) const { return true; }

   virtual void scaleProduct(TEveElementList* parent, FWViewType::EType, const FWViewContext* vc);

   virtual void localModelChanges(const FWModelId& iId, TEveElement* iCompound,
                                  FWViewType::EType viewType, const FWViewContext* vc);
   REGISTER_PROXYBUILDER_METHODS();
private:
   FWPFCandidatesLegoProxyBuilder(const FWPFCandidatesLegoProxyBuilder&); // stop default
   const FWPFCandidatesLegoProxyBuilder& operator=(const FWPFCandidatesLegoProxyBuilder&); // stop default
   
   void build(const reco::PFCandidate&, unsigned int, TEveElement&, const FWViewContext*);

   // ---------- member data --------------------------------
};

//
// constructors and destructor
//
FWPFCandidatesLegoProxyBuilder::FWPFCandidatesLegoProxyBuilder()
{
}

FWPFCandidatesLegoProxyBuilder::~FWPFCandidatesLegoProxyBuilder()
{
}

//
// member functions
//
void 
FWPFCandidatesLegoProxyBuilder::build(const reco::PFCandidate& iData, unsigned int iIndex, TEveElement& oItemHolder, const FWViewContext* vc)
{
   FWLegoEvePFCandidate* evePFCandidate = new FWLegoEvePFCandidate( iData , vc, context());
   evePFCandidate->SetMarkerColor(item()->defaultDisplayProperties().color());
   fireworks::setTrackTypePF( iData,  evePFCandidate);
   setupAddElement( evePFCandidate, &oItemHolder );
}

void
FWPFCandidatesLegoProxyBuilder::scaleProduct(TEveElementList* parent, FWViewType::EType type, const FWViewContext* vc)
{
   for (TEveElement::List_i i = parent->BeginChildren(); i!= parent->EndChildren(); ++i)
   {
      if ((*i)->HasChildren())
      {
         TEveElement* el = (*i)->FirstChild();  // there is only one child added in this proxy builder
         FWLegoEvePFCandidate* cand = dynamic_cast<FWLegoEvePFCandidate*> (el);  
         cand->updateScale(vc, context());
      }
   }
}

void
FWPFCandidatesLegoProxyBuilder::localModelChanges(const FWModelId& iId, TEveElement* parent,
                                                  FWViewType::EType viewType, const FWViewContext* vc)
{
   // line set marker is not same color as line, have to fix it here
   if ((parent)->HasChildren())
   {
      TEveElement* el = (parent)->FirstChild();  // we know there is only one child added in this proxy builder
      FWLegoEvePFCandidate* cand = dynamic_cast<FWLegoEvePFCandidate*> (el); 
      const FWDisplayProperties& dp = item()->modelInfo(iId.index()).displayProperties();
      cand->SetMarkerColor( dp.color());
      cand->ElementChanged();
   }  
}

REGISTER_FWPROXYBUILDER(FWPFCandidatesLegoProxyBuilder, reco::PFCandidate, "PFCandidates", FWViewType::kLegoBit);
