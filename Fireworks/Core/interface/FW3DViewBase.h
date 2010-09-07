#ifndef Fireworks_Core_FW3DViewBase_h
#define Fireworks_Core_FW3DViewBase_h
// -*- C++ -*-
//
// Package:     Core
// Class  :     FW3DViewBase
//
/**\class FW3DViewBase FW3DViewBase.h Fireworks/Core/interface/FW3DViewBase.h

   Description: <one line class summary>

   Usage:
    <usage>

 */
//
// Original Author:  Chris Jones
//         Created:  Thu Feb 21 11:22:37 EST 2008
// $Id: FW3DViewBase.h,v 1.6 2010/08/30 15:42:32 amraktad Exp $
//

// system include files

// user include files
#include "Rtypes.h"
#include "Fireworks/Core/interface/FWEveView.h"
#include "Fireworks/Core/interface/FWLongParameter.h"
#include "Fireworks/Core/interface/FWBoolParameter.h"

// forward declarations
class TEveElementList;
class TEveGeoShape;
class TEveWindowSlot;

class FW3DViewGeometry;
class FWColorManager;

class FW3DViewBase : public FWEveView
{
public:
   FW3DViewBase(TEveWindowSlot*, FWViewType::EType);
   virtual ~FW3DViewBase();

   // ---------- const member functions ---------------------

   virtual void addTo(FWConfiguration&) const;
   virtual void setFrom(const FWConfiguration&);

   // ---------- static member functions --------------------

   // ---------- member functions ---------------------------
   virtual void setContext(const fireworks::Context&);

   // To be fixed.
   void updateGlobalSceneScaleParameters();

private:
   FW3DViewBase(const FW3DViewBase&);    // stop default

   const FW3DViewBase& operator=(const FW3DViewBase&);    // stop default

   // ---------- member data --------------------------------
   FW3DViewGeometry*  m_geometry;

   // parameters
   FWBoolParameter m_showMuonBarrel;
   FWBoolParameter m_showMuonEndcap;
   FWBoolParameter m_showPixelBarrel;
   FWBoolParameter m_showPixelEndcap;
   FWBoolParameter m_showTrackerBarrel;
   FWBoolParameter m_showTrackerEndcap;
   FWBoolParameter m_showWireFrame;
   FWLongParameter m_geomTransparency;
   
   void showWireFrame( bool );
};


#endif
