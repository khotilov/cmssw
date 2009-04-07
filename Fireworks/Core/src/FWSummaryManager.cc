// -*- C++ -*-
//
// Package:     Core
// Class  :     FWSummaryManager
//
// Implementation:
//     <Notes on implementation>
//
// Original Author:  Chris Jones
//         Created:  Tue Mar  4 09:35:32 EST 2008
// $Id: FWSummaryManager.cc,v 1.9 2009/03/18 15:42:00 chrjones Exp $
//

// system include files
#include <boost/bind.hpp>
//#include "TGPack.h"
#include "TGFrame.h"

// user include files
#include "Fireworks/Core/interface/FWSummaryManager.h"
#include "Fireworks/Core/interface/FWSelectionManager.h"
#include "Fireworks/Core/interface/FWEventItemsManager.h"
#include "Fireworks/Core/interface/FWModelChangeManager.h"
#include "Fireworks/Core/interface/FWGUIManager.h"
#include "Fireworks/Core/interface/FWEventItem.h"

#include "Fireworks/Core/src/FWCollectionSummaryWidget.h"
#include "Fireworks/Core/interface/FWDataCategories.h"

#include "Fireworks/Core/src/FWCompactVerticalLayout.h"

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
FWSummaryManager::FWSummaryManager(TGFrame* iParent,
                                   FWSelectionManager* sm,
                                   FWEventItemsManager* eim,
                                   FWGUIManager* gm,
                                   FWModelChangeManager* cm) :
   m_guiManager(gm),
m_itemChanged(false)
{
   sm->selectionChanged_.connect(boost::bind(&FWSummaryManager::selectionChanged,this,_1));
   eim->newItem_.connect(boost::bind(&FWSummaryManager::newItem,
                                     this, _1) );
   eim->goingToClearItems_.connect(boost::bind(&FWSummaryManager::removeAllItems,this));


   m_pack = new TGVerticalFrame(iParent);
   m_pack->SetLayoutManager( new FWCompactVerticalLayout(m_pack));
   const unsigned int backgroundColor=0x2f2f2f;
   m_pack->SetBackgroundColor(backgroundColor);
   /*m_eventObjects =  new TEveElementList("Physics Objects");
      m_listTree->OpenItem(m_eventObjects->AddIntoListTree(m_listTree,
                                                        reinterpret_cast<TGListTreeItem*>(0))
                        );*/
   cm->changeSignalsAreDone_.connect(boost::bind(&FWSummaryManager::changesDone,this));
}

// FWSummaryManager::FWSummaryManager(const FWSummaryManager& rhs)
// {
//    // do actual copying here;
// }

FWSummaryManager::~FWSummaryManager()
{
}

//
// assignment operators
//
// const FWSummaryManager& FWSummaryManager::operator=(const FWSummaryManager& rhs)
// {
//   //An exception safe implementation is
//   FWSummaryManager temp(rhs);
//   swap(rhs);
//
//   return *this;
// }

//
// member functions
//
void
FWSummaryManager::newItem(FWEventItem* iItem)
{
   //TEveElementList* lst = new TEveElementList(iItem->name().c_str(),"",kTRUE);
   //lst->SetMainColor(iItem->defaultDisplayProperties().color());
   TGLayoutHints* hints = new TGLayoutHints(kLHintsExpandX);
   FWCollectionSummaryWidget* lst = new FWCollectionSummaryWidget(m_pack,*iItem,hints);
   m_pack->AddFrame(lst, hints);
   m_collectionWidgets.push_back(lst);
   iItem->goingToBeDestroyed_.connect(boost::bind(&FWSummaryManager::itemDestroyed,this,_1));
   iItem->itemChanged_.connect(boost::bind(&FWSummaryManager::itemChanged,this,_1));
   lst->Connect("requestForInfo(FWEventItem*)","FWSummaryManager",this,"requestForInfo(FWEventItem*)");
   lst->Connect("requestForFilter(FWEventItem*)","FWSummaryManager",this,"requestForFilter(FWEventItem*)");
   lst->Connect("requestForErrorInfo(FWEventItem*)","FWSummaryManager",this,"requestForError(FWEventItem*)");
   //lst->AddIntoListTree(m_listTree,m_eventObjects);
   //NOTE: Why didn't I call AddElement of m_eventObjects???  Because it will go in the wrong list tree?
   //Need to hand this lst into some container so we can get rid of it, since it doesn't seem to go into
   // m_eventObjects
}

void 
FWSummaryManager::itemDestroyed(const FWEventItem* iItem)
{
   m_pack->HideFrame(m_collectionWidgets[iItem->id()]);
   m_pack->RemoveFrame(m_collectionWidgets[iItem->id()]);
   delete m_collectionWidgets[iItem->id()];
   m_collectionWidgets[iItem->id()]=0;
   m_pack->Layout();
   gClient->NeedRedraw(m_pack);
}

void
FWSummaryManager::itemChanged(const FWEventItem*)
{
   m_itemChanged = true;
}

void
FWSummaryManager::removeAllItems()
{
   for(std::vector<FWCollectionSummaryWidget*>::iterator it = m_collectionWidgets.begin(), 
       itEnd = m_collectionWidgets.end();
       it != itEnd;
       ++it) {
      if(0!=*it) {
         m_pack->RemoveFrame(*it);
         delete *it;
         *it=0;
      }
   }
   m_collectionWidgets.clear();
}

void
FWSummaryManager::selectionChanged(const FWSelectionManager& iSM)
{
}

void
FWSummaryManager::changesDone()
{
   if(m_itemChanged) {
      m_pack->Layout();
      m_itemChanged=false;
   }
}

void 
FWSummaryManager::requestForInfo(FWEventItem* iItem)
{
   m_guiManager->updateEDI(iItem);
   m_guiManager->showEDIFrame(kData);
}
void 
FWSummaryManager::requestForFilter(FWEventItem* iItem)
{
   m_guiManager->updateEDI(iItem);
   m_guiManager->showEDIFrame(kFilter);
}
void 
FWSummaryManager::requestForError(FWEventItem* iItem)
{
   m_guiManager->updateEDI(iItem);
   m_guiManager->showEDIFrame();
}

//
// const member functions
//
TGFrame*
FWSummaryManager::widget() const
{
   return m_pack;
}

//
// static member functions
//
