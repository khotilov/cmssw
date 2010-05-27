#ifndef Fireworks_Calo_FWCaloTowerProxyBuilder_h
#define Fireworks_Calo_FWCaloTowerProxyBuilder_h
// -*- C++ -*-
//
// Package:     Calo
// Class  :     FWCaloTowerProxyBuilderBase
//
/**\class FWCaloTowerProxyBuilderBase FWCaloTowerProxyBuilderBase.h Fireworks/Calo/interface/FWCaloTowerProxyBuilderBase.h

   Description: <one line class summary>

   Usage:
    <usage>

 */
//
// Original Author:  Chris Jones
//         Created:  Wed Dec  3 11:28:08 EST 2008
// $Id: FWCaloTowerProxyBuilder.h,v 1.8 2010/05/10 11:49:40 amraktad Exp $
//

#include "Rtypes.h"
#include <string>

#include "Fireworks/Core/interface/FWProxyBuilderBase.h"
#include "DataFormats/CaloTowers/interface/CaloTower.h"
#include "DataFormats/CaloTowers/interface/CaloTowerFwd.h"

class TH2F;

class FWCaloTowerProxyBuilderBase : public FWProxyBuilderBase {

public:
   FWCaloTowerProxyBuilderBase();
   virtual ~FWCaloTowerProxyBuilderBase();

   // ---------- const member functions ---------------------
   virtual const std::string histName() const = 0;
   virtual double getEt(const CaloTower&) const = 0;

   virtual bool willHandleInteraction() const { return true; }
   // ---------- static member functions --------------------

   // ---------- member functions ---------------------------

private:
   FWCaloTowerProxyBuilderBase(const FWCaloTowerProxyBuilderBase&); // stop default

   const FWCaloTowerProxyBuilderBase& operator=(const FWCaloTowerProxyBuilderBase&); // stop default


   virtual void build(const FWEventItem* iItem,
                      TEveElementList* product, const FWViewContext*);


   virtual void modelChanges(const FWModelIds&, Product*);
   virtual void applyChangesToAllModels(Product*);
   virtual void itemBeingDestroyed(const FWEventItem*);

   void clearCaloDataSelection();

   // ---------- member data --------------------------------
   TEveCaloDataHist* m_caloData;
   TH2F* m_hist;
   Int_t m_sliceIndex;
   const CaloTowerCollection* m_towers;
};

//
// Ecal
//

class FWECalCaloTowerProxyBuilder : public FWCaloTowerProxyBuilderBase {
public:
   FWECalCaloTowerProxyBuilder() {
   }
   virtual ~FWECalCaloTowerProxyBuilder() {
   }

   // ---------- const member functions ---------------------
   virtual const std::string histName() const {
      return "ECal";
   }

   virtual double getEt(const CaloTower& iTower) const {
      return iTower.emEt();
   }

   REGISTER_PROXYBUILDER_METHODS();
private:
   FWECalCaloTowerProxyBuilder(const FWECalCaloTowerProxyBuilder&); // stop default
   const FWECalCaloTowerProxyBuilder& operator=(const FWECalCaloTowerProxyBuilder&); // stop default
};


//
// Ecal
//

class FWHCalCaloTowerProxyBuilder : public FWCaloTowerProxyBuilderBase {
public:
   FWHCalCaloTowerProxyBuilder() {
   }
   virtual ~FWHCalCaloTowerProxyBuilder(){
   }

   // ---------- const member functions ---------------------
   virtual const std::string histName() const {
      return "HCal";
   }

   virtual double getEt(const CaloTower& iTower) const {
      return iTower.hadEt()+iTower.outerEt();
   }

   REGISTER_PROXYBUILDER_METHODS();
private:
   FWHCalCaloTowerProxyBuilder(const FWHCalCaloTowerProxyBuilder&); // stop default

   const FWHCalCaloTowerProxyBuilder& operator=(const FWHCalCaloTowerProxyBuilder&); // stop default
};


#endif
