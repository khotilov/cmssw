#ifndef Fireworks_Core_FWSummaryManager_h
#define Fireworks_Core_FWSummaryManager_h
// -*- C++ -*-
//
// Package:     Core
// Class  :     FWSummaryManager
//
/**\class FWSummaryManager FWSummaryManager.h Fireworks/Core/interface/FWSummaryManager.h

   Description: <one line class summary>

   Usage:
    <usage>

 */
//
// Original Author:  Chris Jones
//         Created:  Tue Mar  4 09:35:58 EST 2008
// $Id: FWSummaryManager.h,v 1.6 2009/03/04 17:07:27 chrjones Exp $
//

// system include files
#include <vector>
#include "Rtypes.h"

// user include files

// forward declarations
class TGPack;
class TGFrame;
class TGCompositeFrame;

class FWEventItem;

class FWSelectionManager;
class FWEventItemsManager;
class FWGUIManager;
class FWModelChangeManager;
class FWCollectionSummaryWidget;

class FWSummaryManager
{

public:
   FWSummaryManager(TGFrame* iParent,
                    FWSelectionManager*,
                    FWEventItemsManager*,
                    FWGUIManager*,
                    FWModelChangeManager*);
   virtual ~FWSummaryManager();

   // ---------- const member functions ---------------------
   TGFrame* widget() const;
   
   // ---------- static member functions --------------------

   // ---------- member functions ---------------------------
   void requestForInfo(FWEventItem*);
   void requestForFilter(FWEventItem*);
   void requestForError(FWEventItem*);

private:
   FWSummaryManager(const FWSummaryManager&);    // stop default

   const FWSummaryManager& operator=(const FWSummaryManager&);    // stop default

   void selectionChanged(const FWSelectionManager&);
   void newItem(FWEventItem* iItem);
   void itemChanged(const FWEventItem*);
   void removeAllItems();
   void changesDone();

   void itemDestroyed(const FWEventItem*);
   
   // ---------- member data --------------------------------
   //TGPack* m_pack;
   TGCompositeFrame* m_pack;
   std::vector<FWCollectionSummaryWidget*> m_collectionWidgets;
   FWGUIManager* m_guiManager;
   bool m_itemChanged;
};


#endif
