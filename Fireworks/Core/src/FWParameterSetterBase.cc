// -*- C++ -*-
//
// Package:     Core
// Class  :     FWParameterSetterBase
// 
// Implementation:
//     <Notes on implementation>
//
// Original Author:  Chris Jones
//         Created:  Fri Mar  7 14:16:20 EST 2008
// $Id: FWParameterSetterBase.cc,v 1.3 2008/06/29 13:10:03 chrjones Exp $
//

// system include files
#include "Reflex/Type.h"
#include "Reflex/Object.h"
#include "TGedFrame.h"

#include <assert.h>
#include <iostream>
#include <boost/bind.hpp>

// user include files
#include "FWCore/Utilities/interface/TypeID.h"

#include "Fireworks/Core/interface/FWParameterSetterBase.h"
#include "Fireworks/Core/interface/FWParameterBase.h"
#include "Fireworks/Core/interface/FWParameterSetterEditorBase.h"

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
FWParameterSetterBase::FWParameterSetterBase():
m_frame(0)
{
}

// FWParameterSetterBase::FWParameterSetterBase(const FWParameterSetterBase& rhs)
// {
//    // do actual copying here;
// }

FWParameterSetterBase::~FWParameterSetterBase()
{
}

//
// assignment operators
//
// const FWParameterSetterBase& FWParameterSetterBase::operator=(const FWParameterSetterBase& rhs)
// {
//   //An exception safe implementation is
//   FWParameterSetterBase temp(rhs);
//   swap(rhs);
//
//   return *this;
// }

//
// member functions
//

void 
FWParameterSetterBase::attach(FWParameterBase* iBase, FWParameterSetterEditorBase* iFrame)
{
   m_frame=iFrame;
   attach(iBase);
}


//
// const member functions
//
void 
FWParameterSetterBase::update() const
{
  if (m_frame != 0)
    m_frame->updateEditor();
}

//
// static member functions
//
boost::shared_ptr<FWParameterSetterBase>
FWParameterSetterBase::makeSetterFor(FWParameterBase* iParam)
{
   static std::map<edm::TypeID,ROOT::Reflex::Type> s_paramToSetterMap;
   edm::TypeID paramType( typeid(*iParam) );
   std::map<edm::TypeID,ROOT::Reflex::Type>::iterator itFind = s_paramToSetterMap.find(paramType);
   if( itFind == s_paramToSetterMap.end() ) {
      ROOT::Reflex::Type paramClass( ROOT::Reflex::Type::ByTypeInfo(typeid(*iParam)) );
      if(paramClass == ROOT::Reflex::Type() ) {
         std::cout << "PROGRAMMING ERROR: the type "<<typeid(*iParam).name()<< " is not known to REFLEX" <<std::endl;
      }
      assert(paramClass != ROOT::Reflex::Type() );
      
      //the corresponding setter has the same name but with 'Setter' at the end
      std::string name = paramClass.Name(ROOT::Reflex::SCOPED);
      name += "Setter";
      
      ROOT::Reflex::Type setterClass( ROOT::Reflex::Type::ByName( name ) );
      if(setterClass == ROOT::Reflex::Type() ) {
         std::cout << "PROGRAMMING ERROR: the type "<<name<< " is not known to REFLEX" <<std::endl;
      }
      assert(setterClass != ROOT::Reflex::Type());
      
      s_paramToSetterMap[paramType]=setterClass;
      itFind = s_paramToSetterMap.find(paramType);
   }
   //create the instance we want
   //NOTE: for some odd reason Reflex 'Construct' uses 'malloc' to allocate the memory.  This means the object
   // can not be deleted using 'delete'!  So we must call Type::Destruct on the object
   ROOT::Reflex::Object setterObj = itFind->second.Construct();
   
   //make it into the base class
   ROOT::Reflex::Type s_setterBaseType( ROOT::Reflex::Type::ByTypeInfo( typeid(FWParameterSetterBase) ) );
   assert(s_setterBaseType != ROOT::Reflex::Type());
   ROOT::Reflex::Object castSetterObj = setterObj.CastObject(s_setterBaseType);
   boost::shared_ptr<FWParameterSetterBase> ptr(reinterpret_cast<FWParameterSetterBase*>( castSetterObj.Address() ),
                                                boost::bind(&ROOT::Reflex::Type::Destruct,itFind->second,setterObj.Address(),true));
   return ptr;
}
