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
// $Id: FWEventItem.h,v 1.27 2008/10/21 19:25:06 chrjones Exp $
//

// system include files
#include <string>
#include <vector>
#include "Reflex/Type.h"
#include <boost/shared_ptr.hpp>
#include <sigc++/connection.h>

// user include files
#include "Fireworks/Core/interface/FWDisplayProperties.h"
#include "Fireworks/Core/interface/FWPhysicsObjectDesc.h"
#include "Fireworks/Core/interface/FWModelChangeSignal.h"
#include "Fireworks/Core/interface/FWItemChangeSignal.h"

#include "Fireworks/Core/interface/FWModelFilter.h"

#include "Fireworks/Core/interface/Context.h"

// forward declarations
class TClass;
class FWModelChangeManager;
class FWSelectionManager;
class DetIdToMatrix;
class TVirtualCollectionProxy;
class FWItemAccessorBase;

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

      FWEventItem(fireworks::Context* iContext,
                  unsigned int iItemId,
                  boost::shared_ptr<FWItemAccessorBase> iAccessor,
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

      unsigned int layer() const;

      const std::string& filterExpression() const;
      /**Unique ID for the item. This number starts at 0 and increments by one for each
       new item.*/
      unsigned int id() const;
      const std::string& name() const;
      const TClass* type() const;
      /** Since the same C++ type can be used for multiple purposes, this string disambiguates them.*/
      const std::string& purpose() const;

      const std::string& moduleLabel() const;
      const std::string& productInstanceLabel() const;
      const std::string& processName() const;

      const TClass* modelType() const;
      ModelInfo modelInfo(int iIndex) const; //return copy for now since want to be able to change visibility
      size_t size() const;
      const void* modelData(int iIndex) const;
      std::string modelName(int iIndex) const;

      bool isCollection() const;

      //convenience methods
      FWModelChangeManager* changeManager() const {
         return m_context->modelChangeManager();
      }
      FWSelectionManager* selectionManager() const {
         return m_context->selectionManager();
      }

      bool hasEvent() const {return 0 != m_event; }

      // hackery methods
      const fwlite::Event *getEvent () const { return  m_event; }


      // ---------- static member functions --------------------

      // ---------- member functions ---------------------------
      void setEvent(const fwlite::Event* iEvent);
      void setGeom(const DetIdToMatrix* geom){ m_detIdToGeo = geom; }
      const DetIdToMatrix* getGeom() const { return m_detIdToGeo; }

      void setLabels(const std::string& iModule,
		     const std::string& iProductInstance,
		     const std::string& iProcess);
      void setName(const std::string& iName);
      void setDefaultDisplayProperties(const FWDisplayProperties&);
      /**Throws an FWExpresionException if there is a problem with the expression */
      void setFilterExpression(const std::string& );

      void unselect(int iIndex) const;
      void select(int iIndex) const;
      void toggleSelect(int iIndex) const;
      void setDisplayProperties(int iIndex, const FWDisplayProperties&) const;

      void destroy() const;
      /** connect to this signal if you want to know when models held by the item change */
      mutable FWModelChangeSignal changed_;

      /** connect to this signal if you want to know when the data underlying the item changes */
      mutable FWItemChangeSignal itemChanged_;

      /** connect to this signal if you want to know immediately when the data underlying the item changes
       only intended to be used by the FWSelectionManager
       */
      mutable FWItemChangeSignal preItemChanged_;

      /** connect to this signal if you want to know that the default display properties of the item have changed.
       This is only useful if you are displaying these properties and not just the underlying models.*/
      mutable FWItemChangeSignal defaultDisplayPropertiesChanged_;

      /** connect to this signal if you need to know that this item is going to be destroyed.
       */
      mutable FWItemChangeSignal goingToBeDestroyed_;
   private:
      //FWEventItem(const FWEventItem&); // stop default

      //const FWEventItem& operator=(const FWEventItem&); // stop default
      void setData(const ROOT::Reflex::Object& ) const;

      void getPrimaryData() const;
      void runFilter();
      // ---------- member data --------------------------------
      fireworks::Context* m_context;
      unsigned int m_id;
      std::string m_name;
      const TClass* m_type;
      std::string m_purpose;
      boost::shared_ptr<FWItemAccessorBase> m_accessor;
      FWDisplayProperties m_displayProperties;
      unsigned int m_layer;
      mutable std::vector<ModelInfo> m_itemInfos;

      //This will probably moved to a FWEventItemRetriever class
      std::string m_moduleLabel;
      std::string m_productInstanceLabel;
      std::string m_processName;
      const fwlite::Event* m_event;
      ROOT::Reflex::Type m_wrapperType;
      const DetIdToMatrix* m_detIdToGeo;

      FWModelFilter m_filter;
      sigc::connection m_shouldFilterConnection;
      mutable bool m_printedNoDataError;
      mutable std::string m_fullBranchName;
};


#endif
