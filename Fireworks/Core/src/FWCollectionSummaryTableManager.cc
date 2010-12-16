// -*- C++ -*-
//
// Package:     Core
// Class  :     FWCollectionSummaryTableManager
// 
// Implementation:
//     <Notes on implementation>
//
// Original Author:  Chris Jones
//         Created:  Sun Feb 22 10:13:39 CST 2009
// $Id: FWCollectionSummaryTableManager.cc,v 1.7 2010/06/18 10:17:14 yana Exp $
//

// system include files
#include <sstream>
#include <boost/bind.hpp>
#include "TClass.h"

// user include files
#include "Fireworks/Core/src/FWCollectionSummaryTableManager.h"
#include "Fireworks/Core/interface/FWEventItem.h"
#include "Fireworks/Core/interface/FWItemValueGetter.h"
#include "Fireworks/Core/src/FWCollectionSummaryWidget.h"

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
FWCollectionSummaryTableManager::FWCollectionSummaryTableManager(FWEventItem* iItem, const TGGC* iContext, const TGGC* iHighlightContext,
FWCollectionSummaryWidget* iWidget):
m_collection(iItem),
m_renderer(iContext,iHighlightContext),
m_bodyRenderer(iContext, iHighlightContext, FWTextTableCellRenderer::kJustifyRight),
m_widget(iWidget)
{
   m_collection->changed_.connect(boost::bind(&FWTableManagerBase::dataChanged,this));
   m_collection->itemChanged_.connect(boost::bind(&FWCollectionSummaryTableManager::dataChanged,this));
   
   //try to find the default columns
   std::vector<std::pair<std::string,std::string> > s_names;
   ROOT::Reflex::Type type = ROOT::Reflex::Type::ByTypeInfo(*(m_collection->modelType()->GetTypeInfo()));

   if ( type.Name() == "CaloTower" ){
     if ( m_collection->purpose() == "ECal" ){
       s_names.push_back(std::pair<std::string,std::string>("emEt","GeV"));
       boost::shared_ptr<FWItemValueGetter> trans( new FWItemValueGetter(type,s_names));
       if(trans->isValid()) m_valueGetters.push_back(trans);
     }
     else if ( m_collection->purpose() == "HCal" ){
       s_names.push_back(std::pair<std::string,std::string>("hadEt","GeV"));
       boost::shared_ptr<FWItemValueGetter> hadEt( new FWItemValueGetter(type,s_names));
       if(hadEt->isValid()) m_valueGetters.push_back(hadEt);
     }
     else if (m_collection->purpose() == "HCal Outer"){
        s_names.push_back(std::pair<std::string,std::string>("outerEt","GeV"));
        boost::shared_ptr<FWItemValueGetter> outerEt( new FWItemValueGetter(type,s_names));
        if(outerEt->isValid()) m_valueGetters.push_back(outerEt);
     }
   } else {
     s_names.push_back(std::pair<std::string,std::string>("pt","GeV"));
     s_names.push_back(std::pair<std::string,std::string>("et","GeV"));
     s_names.push_back(std::pair<std::string,std::string>("energy","GeV"));
     boost::shared_ptr<FWItemValueGetter> trans( new FWItemValueGetter(type,s_names));
     if(trans->isValid()) m_valueGetters.push_back(trans);
   }

   
   s_names.clear();
   s_names.push_back(std::pair<std::string,std::string>("eta",""));
   boost::shared_ptr<FWItemValueGetter> eta( new FWItemValueGetter(type,s_names));
   if(eta->isValid()) {
      s_names.clear();
      s_names.push_back(std::pair<std::string,std::string>("phi",""));
      boost::shared_ptr<FWItemValueGetter> phi( new FWItemValueGetter(type,s_names));
      if(phi->isValid()) {
         m_valueGetters.push_back(eta);
         m_valueGetters.push_back(phi);
      }
   }
   
   dataChanged();
}

// FWCollectionSummaryTableManager::FWCollectionSummaryTableManager(const FWCollectionSummaryTableManager& rhs)
// {
//    // do actual copying here;
// }

FWCollectionSummaryTableManager::~FWCollectionSummaryTableManager()
{
}

//
// assignment operators
//
// const FWCollectionSummaryTableManager& FWCollectionSummaryTableManager::operator=(const FWCollectionSummaryTableManager& rhs)
// {
//   //An exception safe implementation is
//   FWCollectionSummaryTableManager temp(rhs);
//   swap(rhs);
//
//   return *this;
// }

//
// member functions
//
namespace {
   template<typename S>
   void doSort(const FWEventItem& iItem,
               FWItemValueGetter& iGetter,
               std::multimap<double,int,S>& iMap,
               std::vector<int>& oNewSort) {
      int size = iItem.size();
      for(int index = 0; index < size; ++index) {
         iMap.insert(std::make_pair(iGetter.valueFor(iItem.modelData(index)),
                                       index));
      }
      std::vector<int>::iterator itVec = oNewSort.begin();
      for(typename std::map<double,int,S>::iterator it = iMap.begin(), itEnd = iMap.end();
          it != itEnd;
          ++it,++itVec) {
         *itVec = it->second;
      }
   }
}

void 
FWCollectionSummaryTableManager::implSort(int iCol, bool iSortOrder)
{
   if(iSortOrder) {
      std::multimap<double,int, std::greater<double> > s;
      doSort(*m_collection, *(m_valueGetters[iCol]), s, m_sortedToUnsortedIndicies);
   } else {
      std::multimap<double,int, std::less<double> > s;
      doSort(*m_collection, *(m_valueGetters[iCol]), s, m_sortedToUnsortedIndicies);
   }
}

void 
FWCollectionSummaryTableManager::buttonReleasedInRowHeader(Int_t row, Event_t* event, Int_t relX, Int_t relY)
{
   Int_t realRow = unsortedRowNumber(row);
   int hit = m_renderer.clickHit(relX,relY);
   if(hit == FWCollectionSummaryModelCellRenderer::kMiss) {
      return;
   }
   if(hit == FWCollectionSummaryModelCellRenderer::kHitColor) {
      m_widget->itemColorClicked(realRow,event->fXRoot, event->fYRoot+12-relY);
      return;
   }
   FWEventItem::ModelInfo mi = m_collection->modelInfo(realRow);
   FWDisplayProperties dp = mi.displayProperties();
   if(hit == FWCollectionSummaryModelCellRenderer::kHitCheck) {
      dp.setIsVisible(!dp.isVisible());
   }
   m_collection->setDisplayProperties(realRow,dp);
}

//
// const member functions
//
int 
FWCollectionSummaryTableManager::numberOfRows() const
{
   return m_collection->size();
}

int 
FWCollectionSummaryTableManager::numberOfColumns() const {
   return m_valueGetters.size();
}

std::vector<std::string> 
FWCollectionSummaryTableManager::getTitles() const {
   std::vector<std::string> titles;
   titles.reserve(m_valueGetters.size());
   for(std::vector<boost::shared_ptr<FWItemValueGetter> >::const_iterator it = m_valueGetters.begin(), itEnd=m_valueGetters.end();
       it != itEnd;
       ++it) {
      titles.push_back((*it)->valueName());
   }
   return titles;
}

int 
FWCollectionSummaryTableManager::unsortedRowNumber(int iSortedRowNumber) const
{
   return m_sortedToUnsortedIndicies[iSortedRowNumber];
}

FWTableCellRendererBase* 
FWCollectionSummaryTableManager::cellRenderer(int iSortedRowNumber, int iCol) const
{
   if(iCol >= static_cast<int>(m_valueGetters.size())) {
      return 0;
   }
   if(iSortedRowNumber >= static_cast<int>(m_collection->size())) {
      m_bodyRenderer.setData("",false);
      return &m_bodyRenderer;
   }
   int index = m_sortedToUnsortedIndicies[iSortedRowNumber];
   std::stringstream s;
   s.setf(std::ios_base::fixed,std::ios_base::floatfield);
   s.precision(1);
   double v = m_valueGetters[iCol]->valueFor(m_collection->modelData(index));
   s <<v;
   m_bodyRenderer.setData(s.str(),
                          m_collection->modelInfo(index).isSelected());
   return &m_bodyRenderer;
}

bool 
FWCollectionSummaryTableManager::hasRowHeaders() const
{
   return true;
}

FWTableCellRendererBase* 
FWCollectionSummaryTableManager::rowHeader(int iSortedRowNumber) const
{
   if(iSortedRowNumber >= static_cast<int>(m_collection->size())) {
      return 0;
   }
   int index = m_sortedToUnsortedIndicies[iSortedRowNumber];
   m_renderer.setData(m_collection,
                      index);
   return &m_renderer;
}

void
FWCollectionSummaryTableManager::dataChanged() 
{
   m_sortedToUnsortedIndicies.clear();
   m_sortedToUnsortedIndicies.reserve(m_collection->size());
   for(int i=0; i< static_cast<int>(m_collection->size());++i) {
      m_sortedToUnsortedIndicies.push_back(i);
   }
   FWTableManagerBase::dataChanged();
}
//
// static member functions
//
