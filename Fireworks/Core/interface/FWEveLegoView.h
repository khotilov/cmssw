#ifndef Fireworks_Core_FWEveLegoView_h
#define Fireworks_Core_FWEveLegoView_h
// -*- C++ -*-
//
// Package:     Core
// Class  :     FWEveLegoView
//
/**\class FWEveLegoView FWEveLegoView.h Fireworks/Core/interface/FWEveLegoView.h

   Description: <one line class summary>

   Usage:
    <usage>

 */
//
// Original Author:  Chris Jones
//         Created:  Thu Feb 21 11:22:37 EST 2008
// $Id: FWEveLegoView.h,v 1.17 2009/09/18 09:49:43 amraktad Exp $
//

// system include files
#include "Rtypes.h"

// user include files
#include "Fireworks/Core/interface/FWViewBase.h"
#include "Fireworks/Core/interface/FWBoolParameter.h"
#include "Fireworks/Core/interface/FWEvePtr.h"

// forward declarations
class TGFrame;
class TGLEmbeddedViewer;
class TEvePad;
class TEveViewer;
class TEveScene;
class TEveElementList;
class TEveCaloLegoOverlay;
class TGLMatrix;
class TEvwWindowSlot;
class FWViewContextMenuHandlerGL;

class FWEveLegoView : public FWViewBase
{

public:
   FWEveLegoView(TEveWindowSlot*, TEveElementList*);
   virtual ~FWEveLegoView();

   // ---------- const member functions ---------------------
   TGFrame* frame() const;
   const std::string& typeName() const;
   virtual void addTo(FWConfiguration&) const;

   virtual void setFrom(const FWConfiguration&);
   virtual void saveImageTo(const std::string& iName) const;
   virtual FWViewContextMenuHandlerBase* contextMenuHandler() const;

   // ---------- static member functions --------------------
   static const std::string& staticTypeName();

   // ---------- member functions ---------------------------
   void finishSetup();
   // set energy thresholds from the parameters
   void setMinEnergy();
   void setBackgroundColor(Color_t);


private:
   FWEveLegoView(const FWEveLegoView&);    // stop default

   const FWEveLegoView& operator=(const FWEveLegoView&);    // stop default

   void setMinEcalEnergy(double);
   void setMinHcalEnergy(double);
   void setCameras();
   void setAutoRebin();
   void showScales();
   
   // ---------- member data --------------------------------
   FWEvePtr<TEveViewer> m_viewer;
   TGLEmbeddedViewer* m_embeddedViewer;
   FWEvePtr<TEveScene> m_scene;
   boost::shared_ptr<FWViewContextMenuHandlerGL>   m_viewContextMenu;

   TEveCaloLego* m_lego;
   TEveCaloLegoOverlay* m_overlay;
   // FWLongParameter m_range;
   //FWDoubleParameter m_minEcalEnergy;
   //FWDoubleParameter m_minHcalEnergy;
   //double m_minEcalEnergyInit;
   //double m_minHcalEnergyInit;
   FWBoolParameter m_autoRebin;
   FWBoolParameter m_showScales;

   TGLMatrix*  m_cameraMatrix;
   TGLMatrix*  m_cameraMatrixBase;
   TGLMatrix*  m_cameraMatrixRef;
   TGLMatrix*  m_cameraMatrixBaseRef;
   double m_orthoCameraZoom;
   TGLMatrix*  m_orthoCameraMatrix;
   double*     m_orthoCameraZoomRef;
   TGLMatrix*  m_orthoCameraMatrixRef;
   bool m_topView;
   bool m_cameraSet;
};


#endif
