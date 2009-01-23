// -*- C++ -*-
//
// Package:     Core
// Class  :     FWSimpleProxyHelper
//
// Implementation:
//     <Notes on implementation>
//
// Original Author:  Chris Jones
//         Created:  Tue Dec  2 15:13:22 EST 2008
// $Id: FWSimpleProxyHelper.cc,v 1.1 2008/12/02 21:11:53 chrjones Exp $
//

// system include files
#include <sstream>

#include "Reflex/Object.h"
#include "Reflex/Type.h"
#include "TClass.h"

// user include files
#include "Fireworks/Core/interface/FWSimpleProxyHelper.h"
#include "Fireworks/Core/interface/FWEventItem.h"


//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
FWSimpleProxyHelper::FWSimpleProxyHelper(const std::type_info& iType) :
   m_itemType(&iType),
   m_objectOffset(0)
{
}

// FWSimpleProxyHelper::FWSimpleProxyHelper(const FWSimpleProxyHelper& rhs)
// {
//    // do actual copying here;
// }

//FWSimpleProxyHelper::~FWSimpleProxyHelper()
//{
//}

//
// assignment operators
//
// const FWSimpleProxyHelper& FWSimpleProxyHelper::operator=(const FWSimpleProxyHelper& rhs)
// {
//   //An exception safe implementation is
//   FWSimpleProxyHelper temp(rhs);
//   swap(rhs);
//
//   return *this;
// }

//
// member functions
//
void
FWSimpleProxyHelper::itemChanged(const FWEventItem* iItem)
{
   if(0!=iItem) {
      using namespace ROOT::Reflex;
      Type myType = Type::ByTypeInfo(*m_itemType);
      Object dummy(Type::ByTypeInfo(*(iItem->modelType()->GetTypeInfo())),
                   reinterpret_cast<void*>(0xFFFF));
      Object castTo = dummy.CastObject(myType);
      assert(0!=castTo.Address());
      m_objectOffset=static_cast<char*>(dummy.Address())-static_cast<char*>(castTo.Address());
   }
}

//
// const member functions
//
void
FWSimpleProxyHelper::fillTitle(const FWEventItem& iItem, int iIndex, std::string& oTitle) const
{
   oTitle = iItem.modelName(iIndex);
   std::stringstream s;
   if(iItem.haveInterestingValue()) {
      s<<oTitle<<", "<<iItem.modelInterestingValueAsString(iIndex);
      oTitle=s.str();
   }
}


//
// static member functions
//
