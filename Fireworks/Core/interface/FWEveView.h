// -*- C++ -*-
//
// Package:     Core
// Class  :     FWEveView
// 
// Implementation:
//     [Notes on implementation]
//
// Original Author:  Alja Mrak-Tadel
//         Created:  Thu Mar 16 14:11:32 CET 2010
// $Id: FWEveView.h,v 1.20 2010/09/24 16:22:26 amraktad Exp $
//


#ifndef Fireworks_Core_FWEveView_h
#define Fireworks_Core_FWEveView_h

// user include files
#include "Fireworks/Core/interface/FWViewBase.h"
#include "Fireworks/Core/interface/FWViewType.h"
#include "Fireworks/Core/interface/FWDoubleParameter.h"
#include "Fireworks/Core/interface/FWBoolParameter.h"
#include "Fireworks/Core/interface/FWLongParameter.h"
#include "Fireworks/Core/interface/FWEnumParameter.h"
#include "Fireworks/Core/interface/FWEvePtr.h"

// forward declarations
class TGLViewer;
class TGLOrthoCamera;
class TGLPerspectiveCamera;
class TGLCameraGuide;
class TGLAnnotation;
class TEveViewer;
class TEveElementList;
class TEveScene;
class TEveWindowSlot;
class TEveCaloViz;

class FWEventAnnotation;
class CmsAnnotation;
class FWViewContextMenuHandlerGL;
class FWColorManager;
class FWViewContext;
class ViewerParameterGUI;

namespace fireworks
{
   class Context;
}


class FWEveView : public FWViewBase
{
public:
   enum EScaleMode { kFixedScale, kAutoScale, kCombinedScale, kNone };

   FWEveView(TEveWindowSlot*, FWViewType::EType, unsigned int version = 5);
   virtual ~FWEveView();

   virtual void setFrom(const FWConfiguration&);
   virtual void setBackgroundColor(Color_t);
   virtual void eventEnd();
   virtual void eventBegin();

   virtual void setContext(const fireworks::Context& x) { m_context = &x ;}
   const fireworks::Context& context()  { return *m_context; } 

   // ---------- const member functions --------------------- 

   virtual void addTo(FWConfiguration&) const;
   virtual FWViewContextMenuHandlerBase* contextMenuHandler() const;
   virtual void saveImageTo(const std::string& iName) const;
   virtual void populateController(ViewerParameterGUI&) const;
   virtual const std::string& typeName() const;

   bool  useGlobalScales() const { return m_useGlobalScales.value(); }

   TGLViewer*  viewerGL() const;
   TEveViewer* viewer()      { return m_viewer; }
   TEveScene*  eventScene()  { return m_eventScene;}
   TEveScene*  geoScene()    { return m_geoScene; }

   TEveElement*   ownedProducts()  { return m_ownedProducts; }
   FWViewContext* viewContext() { return m_viewContext.get(); }
   

   // ---------- static member functions --------------------
   
   FWViewType::EType typeId() const { return m_type.id(); }
   virtual void updateEnergyScales();
   virtual void setUseGlobalEnergyScales(bool);
protected:
   virtual void resetCamera();
   virtual void pointLineScalesChanged();

   // scales
   virtual TEveCaloViz* getEveCalo() const { return 0;}

   long   energyScaleMode();
   double energyMaxAbsVal();
   double energyMaxTowerHeight();

   virtual void energyScalesChanged();
   virtual void setMaxTowerHeight();

   // config
   void addToOrthoCamera(TGLOrthoCamera*, FWConfiguration&) const;
   void setFromOrthoCamera(TGLOrthoCamera*, const FWConfiguration&);
   void addToPerspectiveCamera(TGLPerspectiveCamera*, const std::string&, FWConfiguration&) const;
   void setFromPerspectiveCamera(TGLPerspectiveCamera*,  const std::string&, const FWConfiguration&);


private:
   FWEveView(const FWEveView&);    // stop default
   const FWEveView& operator=(const FWEveView&);    // stop default
  

   // ---------- member data --------------------------------

   FWViewType           m_type;
   TEveViewer*          m_viewer;
   TEveScene*           m_eventScene;
   TEveElement*         m_ownedProducts;
   TEveScene*           m_geoScene;

   FWEventAnnotation*   m_overlayEventInfo;  
   CmsAnnotation*       m_overlayLogo;
   TGLAnnotation*       m_energyMaxValAnnotation;
   TGLCameraGuide*      m_cameraGuide;

   const fireworks::Context*  m_context;

   // style parameters
#if ROOT_VERSION_CODE >= ROOT_VERSION(5,26,0)
   FWDoubleParameter   m_imageScale;
#endif
   FWEnumParameter   m_eventInfoLevel;
   FWBoolParameter   m_drawCMSLogo;

   FWBoolParameter   m_pointSmooth;
   FWDoubleParameter m_pointSize;
   FWBoolParameter   m_lineSmooth; 
   FWDoubleParameter m_lineWidth;
   FWDoubleParameter m_lineOutlineScale;
   FWDoubleParameter m_lineWireframeScale;

   FWBoolParameter   m_showCameraGuide;

   //protected:
   // scale parameters are  protected for use in configuration backward copatibility
   // later when is more material,  members will be moved new FWEveCaloView class or oder helper strucutre
   FWBoolParameter    m_useGlobalScales;
   FWEnumParameter    m_energyScaleMode;
   FWDoubleParameter  m_energyMaxAbsVal;
   FWDoubleParameter  m_energyMaxTowerHeight;

private:
   boost::shared_ptr<FWViewContextMenuHandlerGL>   m_viewContextMenu;
   std::auto_ptr<FWViewContext> m_viewContext;
};


#endif
