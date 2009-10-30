#ifndef Fireworks_Core_FWDetailViewBase_h
#define Fireworks_Core_FWDetailViewBase_h
// -*- C++ -*-
//
// Package:     Core
// Class  :     FWDetailViewBase
//
/**\class FWDetailViewBase FWDetailViewBase.h Fireworks/Core/interface/FWDetailViewBase.h

   Description: Base class for detailed views

   Usage:
    <usage>

 */
//
// Original Author:  Chris Jones
//         Created:  Fri Jan  9 13:35:52 EST 2009
// $Id: FWDetailViewBase.h,v 1.5 2009/06/18 16:03:49 amraktad Exp $
//

// system include files
#include "Rtypes.h"

// user include files
#include "Fireworks/Core/interface/FWSimpleProxyHelper.h"

// forward declarations
class TEveElement;
class TGLViewer;
class FWModelId;
class TCanvas;

class FWDetailViewBase
{
public:
   virtual ~FWDetailViewBase ();

   ///the calling code takes ownership of the returned object
   TEveElement* build (const FWModelId &);
   virtual void  clearOverlayElements() {}

   void         setViewer (TGLViewer *v) {
      m_viewer = v;
   }

   void         setTextCanvas (TCanvas *v) {
      m_textCanvas = v;
   }

   void         setViewCanvas (TCanvas *c) {
      m_viewCanvas = c;
   }


   void  setUseGL (Bool_t x) {
      m_useGL = x;
   }

   Bool_t useGL() const {
      return m_useGL;
   }

protected:
   FWDetailViewBase(const std::type_info&);

   TCanvas*      textCanvas() const {
      return m_textCanvas;
   }

   TGLViewer*       viewer() const {
      return m_viewer;
   }

   TCanvas*       viewCanvas () const {
      return m_viewCanvas;
   }

   const Double_t*  rotationCenter() const {
      return m_rotationCenter;
   }
   Double_t*        rotationCenter() {
      return m_rotationCenter;
   }

   void getCenter( Double_t* vars )
   {
      vars[0] = rotationCenter()[0];
      vars[1] = rotationCenter()[1];
      vars[2] = rotationCenter()[2];
   }

   void resetCenter() {
      rotationCenter()[0] = 0;
      rotationCenter()[1] = 0;
      rotationCenter()[2] = 0;
   }

private:
   FWDetailViewBase(const FWDetailViewBase&); // stop default

   const FWDetailViewBase& operator=(const FWDetailViewBase&); // stop default

   virtual TEveElement* build(const FWModelId&, const void*) = 0;

   Bool_t           m_useGL;

   TGLViewer       *m_viewer;
   TCanvas         *m_viewCanvas;
   TCanvas         *m_textCanvas;
   Double_t         m_rotationCenter[3];


   FWSimpleProxyHelper m_helper;

   // ---------- member data --------------------------------

};


#endif
