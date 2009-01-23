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
// $Id: FWEveLegoView.h,v 1.11 2009/01/22 16:05:12 amraktad Exp $
//

// system include files
#include "Rtypes.h"

// user include files
#include "Fireworks/Core/interface/FWViewBase.h"
#include "Fireworks/Core/interface/FWLongParameter.h"
#include "Fireworks/Core/interface/FWDoubleParameter.h"
#include "Fireworks/Core/interface/FWEvePtr.h"

// forward declarations
class TGFrame;
class TGLEmbeddedViewer;
class TEvePad;
class TEveViewer;
class TEveScene;
class TEveElementList;
class TGLMatrix;

class FWEveLegoView : public FWViewBase
{

public:
   FWEveLegoView(TGFrame*, TEveElementList*);
   virtual ~FWEveLegoView();

   // ---------- const member functions ---------------------
   TGFrame* frame() const;
   const std::string& typeName() const;
   virtual void addTo(FWConfiguration&) const;

   virtual void saveImageTo(const std::string& iName) const;

   // ---------- static member functions --------------------
   static const std::string& staticTypeName();

   // ---------- member functions ---------------------------
   void finishSetup();
   virtual void setFrom(const FWConfiguration&);
   // set energy thresholds from the parameters
   void setMinEnergy();

private:
   FWEveLegoView(const FWEveLegoView&);    // stop default

   const FWEveLegoView& operator=(const FWEveLegoView&);    // stop default

   void setMinEcalEnergy(double);
   void setMinHcalEnergy(double);
   void setCameras();

   // ---------- member data --------------------------------
   TEvePad* m_pad;
   FWEvePtr<TEveViewer> m_viewer;
   TGLEmbeddedViewer* m_embeddedViewer;
   FWEvePtr<TEveScene> m_scene;

   // FWLongParameter m_range;
   //FWDoubleParameter m_minEcalEnergy;
   //FWDoubleParameter m_minHcalEnergy;
   //double m_minEcalEnergyInit;
   //double m_minHcalEnergyInit;


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
