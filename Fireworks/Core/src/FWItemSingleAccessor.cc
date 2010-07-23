// -*- C++ -*-
//
// Package:     Core
// Class  :     FWItemSingleAccessor
//
// Implementation:
//     <Notes on implementation>
//
// Original Author:  Chris Jones
//         Created:  Sat Oct 18 11:36:44 EDT 2008
// $Id: FWItemSingleAccessor.cc,v 1.4 2009/01/23 21:35:43 amraktad Exp $
//

// system include files
#include <assert.h>
#include "Reflex/Object.h"

// user include files
#include "Fireworks/Core/src/FWItemSingleAccessor.h"


//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
FWItemSingleAccessor::FWItemSingleAccessor(const TClass* iClass) :
   m_type(iClass), m_data(0)
{
}

// FWItemSingleAccessor::FWItemSingleAccessor(const FWItemSingleAccessor& rhs)
// {
//    // do actual copying here;
// }

FWItemSingleAccessor::~FWItemSingleAccessor()
{
}

//
// assignment operators
//
// const FWItemSingleAccessor& FWItemSingleAccessor::operator=(const FWItemSingleAccessor& rhs)
// {
//   //An exception safe implementation is
//   FWItemSingleAccessor temp(rhs);
//   swap(rhs);
//
//   return *this;
// }

//
// member functions
//
void
FWItemSingleAccessor::setData(const ROOT::Reflex::Object& product)
{
   if (product.Address() == 0)
   {
      reset();
      return;
   }
   
   using ROOT::Reflex::Object;
   
   if(product.TypeOf().IsTypedef())
      m_data = Object(product.TypeOf().ToType(), product.Address()).Address();
   else
      m_data = product.Address();
   assert(0!=m_data);
}

void
FWItemSingleAccessor::reset()
{
   m_data = 0;
}

//
// const member functions
//
const void*
FWItemSingleAccessor::modelData(int iIndex) const
{
   if(0==iIndex) {
      return m_data;
   }
   return 0;
}

const void*
FWItemSingleAccessor::data() const
{
   return m_data;
}

unsigned int
FWItemSingleAccessor::size() const
{
   return 0 == m_data ? 0 : 1;
}

const TClass*
FWItemSingleAccessor::modelType() const
{
   return m_type;
}

const TClass*
FWItemSingleAccessor::type() const
{
   return m_type;
}

bool
FWItemSingleAccessor::isCollection() const
{
   return false;
}

//
// static member functions
//
