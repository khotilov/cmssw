#ifndef Fireworks_Core_FWEveViewManager_h
#define Fireworks_Core_FWEveViewManager_h
// -*- C++ -*-
//
// Package:     Core
// Class  :     FWEveViewManager
// 
/**\class FWEveViewManager FWEveViewManager.h Fireworks/Core/interface/FWEveViewManager.h

 Description: [one line class summary]

 Usage:
    <usage>

*/
//
// Original Author:  Chris Jones, Alja Mrak-Tadel
//         Created:  Thu Mar 18 14:12:45 CET 2010
// $Id: FWEveViewManager.h,v 1.1 2010/04/06 20:00:35 amraktad Exp $
//

// system include files
#include <vector>
#include <map>
#include <boost/shared_ptr.hpp>

// user include files
#include "Fireworks/Core/interface/FWViewManagerBase.h"
#include "Fireworks/Core/interface/FWViewType.h"
#include "Fireworks/Core/interface/FWProxyBuilderBase.h"

// forward declarations
class TEveCompund;
class TEveScene;
class TEveWindowSlot;
class FWViewBase;
class FWEveView;
class FWProxyBuilderBase;
class FWGUIManager;
class FWSelectionManager;

class FWEveViewManager : public FWViewManagerBase
{
private:
   struct BuilderInfo
   {
      std::string m_name;
      int         m_viewBit;

      BuilderInfo(std::string name, int viewBit) :
         m_name(name),
         m_viewBit(viewBit)
      {}
   };

public:
   FWEveViewManager(FWGUIManager*);
   virtual ~FWEveViewManager();

   // ---------- const member functions ---------------------

   // ---------- static member functions --------------------

   // ---------- member functions ---------------------------
   virtual void newItem(const FWEventItem*);
   virtual void eventBegin();
   virtual void eventEnd();

   void selectionAdded(TEveElement*);
   void selectionRemoved(TEveElement*);
   void selectionCleared();

   FWTypeToRepresentations supportedTypesAndRepresentations() const;

protected:
   virtual void modelChangesComing();
   virtual void modelChangesDone();
   virtual void colorsChanged();

private:
   FWEveViewManager(const FWEveViewManager&); // stop default
   const FWEveViewManager& operator=(const FWEveViewManager&); // stop default

   FWViewBase* create3DRecHitView  (TEveWindowSlot* iParent);
   FWViewBase* create3DEView       (TEveWindowSlot* iParent);
   FWViewBase* createLegoView      (TEveWindowSlot* iParent);
   FWViewBase* createGlimpseView   (TEveWindowSlot* iParent);
   FWViewBase* createRhoPhiView    (TEveWindowSlot* iParent);
   FWViewBase* createRhoZView      (TEveWindowSlot* iParent);
   FWEveView*  finishViewCreate    (boost::shared_ptr<FWEveView>);

   void makeProxyBuilderFor(const FWEventItem* iItem);
   void beingDestroyed(const FWViewBase*);

   bool haveViewForBit (int) const;

   // ---------- member data --------------------------------
   
   FWSelectionManager*     m_selectionManager; // set via FWEventItem

   typedef  std::map<std::string,  std::vector<BuilderInfo> >  TypeToBuilder;
   
   typedef  std::vector<boost::shared_ptr<FWProxyBuilderBase> >  BuilderVec;   
   typedef std::vector<boost::shared_ptr<FWProxyBuilderBase> >::iterator BuilderVec_it;
   
   typedef std::vector<boost::shared_ptr<FWEveView > >::iterator EveViewVec_it;
   
   TypeToBuilder            m_typeToBuilder;

   std::map<int, BuilderVec> m_builders;

   std::vector< std::vector<boost::shared_ptr<FWEveView> > >  m_views;
   std::vector<TEveScene*>  m_scenes;

   // TODO ...
   // std::map<const FWEventItem*, std::vector<TEveCompund*> >  m_interactions;

};


#endif
