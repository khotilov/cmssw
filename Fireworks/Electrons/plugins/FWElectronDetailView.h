// -*- C++ -*-
#ifndef Fireworks_Electrons_FWElectronDetailView_h
#define Fireworks_Electrons_FWElectronDetailView_h

//
// Package:     Electrons
// Class  :     FWElectronDetailView
//
// Implementation:
//     <Notes on implementation>
//
// Original Author:
//         Created:  Sun Jan  6 23:57:00 EST 2008
// $Id: FWElectronDetailView.h,v 1.7 2009/03/31 23:27:20 jmuelmen Exp $
//


// user include files
#include "FWECALDetailView.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectronFwd.h"

class FWElectronDetailView : public FWECALDetailView<reco::GsfElectron> {

public:
   FWElectronDetailView();
   virtual ~FWElectronDetailView();

   virtual TEveElement* build (const FWModelId &id, const reco::GsfElectron*);
   void showInterestingHits(TGLViewerBase*);

protected:
   virtual bool	drawTrack () { return true; }
   virtual math::XYZPoint trackPositionAtCalo (const reco::GsfElectron &);
   virtual double deltaEtaSuperClusterTrackAtVtx (const reco::GsfElectron &);
   virtual double deltaPhiSuperClusterTrackAtVtx (const reco::GsfElectron &);
   TEveElement* build_projected (const FWModelId &id, const reco::GsfElectron*);
   class TEveElementList *makeLabels (const reco::GsfElectron &);

   void fillReducedData (const std::vector<DetId> &detids,
			 TEveCaloDataVec *data);

private:
   FWElectronDetailView(const FWElectronDetailView&); // stop default
   const FWElectronDetailView& operator=(const FWElectronDetailView&); // stop default
};

#endif
