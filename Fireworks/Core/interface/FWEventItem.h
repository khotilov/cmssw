#ifndef Fireworks_Core_FWEventItem_h
#define Fireworks_Core_FWEventItem_h
// -*- C++ -*-
//
// Package:     Core
// Class  :     FWEventItem
// 
/**\class FWEventItem FWEventItem.h Fireworks/Core/interface/FWEventItem.h

 Description: Stand in for a top level item in an Event

 Usage:
    <usage>

*/
//
// Original Author:  Chris Jones
//         Created:  Thu Jan  3 14:02:21 EST 2008
// $Id: FWEventItem.h,v 1.10 2008/01/25 04:05:35 chrjones Exp $
//

// system include files
#include <string>
#include <vector>
#include "Reflex/Type.h"
#include <boost/shared_ptr.hpp>

// user include files
#include "Fireworks/Core/interface/FWDisplayProperties.h"
#include "Fireworks/Core/interface/FWPhysicsObjectDesc.h"
#include "Fireworks/Core/interface/FWModelChangeSignal.h"
#include "Fireworks/Core/interface/FWItemChangeSignal.h"

// forward declarations
class TClass;
class FWModelChangeManager;
class FWSelectionManager;
class DetIdToMatrix;
class TVirtualCollectionProxy;

namespace fwlite {
  class Event;
}

class FWEventItem
{

   public:
      struct ModelInfo {
         FWDisplayProperties m_displayProperties;
         bool m_isSelected;
         ModelInfo(const FWDisplayProperties& iProps, bool iIsSelected):
         m_displayProperties(iProps),
         m_isSelected(iIsSelected) {}
         
         const FWDisplayProperties& displayProperties() const {
            return m_displayProperties;
         }
         bool isSelected() const {
            return m_isSelected;
         }
      };
   
      FWEventItem(FWModelChangeManager* iCM,
                  FWSelectionManager* iSM,
                  unsigned int iItemId,
                  const std::string& iName,
		  const TClass* iClass,
		  const FWDisplayProperties& iProperties =
		  FWDisplayProperties(),
		  const std::string& iModuleLabel = std::string(),
		  const std::string& iProductInstanceLabel = std::string(),
		  const std::string& iProcessName = std::string());
   
      FWEventItem(FWModelChangeManager* iCM,
                  FWSelectionManager* iSM,
                  unsigned int iItemId,
                  const FWPhysicsObjectDesc& iDesc);
      //virtual ~FWEventItem();

      // ---------- const member functions ---------------------
#if !defined(__CINT__) && !defined(__MAKECINT__)
      template<class T>
	void get(const T*& oData) const {
	oData=reinterpret_cast<const T*>(data(typeid(T)));
      }
#endif
      const void* data(const std::type_info&) const;
      const FWDisplayProperties& defaultDisplayProperties() const;
   
      /**Unique ID for the item. This number starts at 0 and increments by one for each
       new item.*/
      unsigned int id() const;
      const std::string& name() const;
      const TClass* type() const;

      const std::string& moduleLabel() const;
      const std::string& productInstanceLabel() const;
      const std::string& processName() const;
   
      const TClass* modelType() const;
      const ModelInfo& modelInfo(int iIndex) const;
      size_t size() const;
      const void* modelData(int iIndex) const;
   
      //convenience methods
      FWModelChangeManager* changeManager() const {
         return m_changeManager;
      }
      FWSelectionManager* selectionManager() const {
         return m_selectionManager;
      }
   
      // ---------- static member functions --------------------

      // ---------- member functions ---------------------------
      void setEvent(const fwlite::Event* iEvent);
      void setGeom(const DetIdToMatrix* geom){ m_detIdToGeo = geom; }
      const DetIdToMatrix* getGeom() const { return m_detIdToGeo; }

      void setLabels(const std::string& iModule,
		     const std::string& iProductInstance,
		     const std::string& iProcess);
      void setName(const std::string& iName);

      void unselect(int iIndex) const;
      void select(int iIndex) const;
      void toggleSelect(int iIndex) const;
      void setDisplayProperties(int iIndex, const FWDisplayProperties&) const;
   
      /** connect to this signal if you want to know when models held by the item change */
      mutable FWModelChangeSignal changed_;
   
      /** connect to this signal if you want to know when the data underlying the item changes */
      mutable FWItemChangeSignal itemChanged_;
   
      /** connect to this signal if you want to know immediately when the data underlying the item changes 
       only intended to be used by the FWSelectionManager
       */
      mutable FWItemChangeSignal preItemChanged_;
   private:
      //FWEventItem(const FWEventItem&); // stop default

      //const FWEventItem& operator=(const FWEventItem&); // stop default
      void setData(const void* ) const;
   
      // ---------- member data --------------------------------
      FWModelChangeManager* m_changeManager;
      FWSelectionManager* m_selectionManager;
      unsigned int m_id;
      std::string m_name;
      const TClass* m_type;
      boost::shared_ptr<TVirtualCollectionProxy> m_colProxy; //should be something other than shared_ptr 
      mutable const void * m_data;
      size_t m_collectionOffset;
      FWDisplayProperties m_displayProperties;
      mutable std::vector<ModelInfo> m_itemInfos;

      //This will probably moved to a FWEventItemRetriever class
      std::string m_moduleLabel;
      std::string m_productInstanceLabel;
      std::string m_processName;
      const fwlite::Event* m_event;
      ROOT::Reflex::Type m_wrapperType;
      const DetIdToMatrix* m_detIdToGeo;
};


#endif
