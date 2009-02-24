// -*- C++ -*-
//
// Package:     Utilities
// Class  :     findDataMember
// 
// Implementation:
//     <Notes on implementation>
//
// Original Author:  Chris Jones
//         Created:  Wed Aug 13 10:07:46 EDT 2008
// $Id: findDataMember.cc,v 1.2 2009/01/11 23:37:33 hegner Exp $
//

// system include files
#include "Reflex/Base.h"

// user include files
#include "CommonTools/Utils/src/findDataMember.h"
#include "CommonTools/Utils/src/ErrorCodes.h"

//
// constants, enums and typedefs
//

namespace reco {
   Reflex::Member findDataMember(const Reflex::Type& iType, const std::string& iName, int& oError) {
      using namespace Reflex;
      Member returnValue;
      oError = parser::kNameDoesNotExist;
      Type type = iType;
      if(type) {
         if(type.IsPointer()) {
            type = type.ToType();
         }
         returnValue = type.DataMemberByName(iName);
         if(!returnValue) {
            //check inheriting classes
            for(Base_Iterator b = type.Base_Begin(); b != type.Base_End(); ++ b) {
               returnValue = findDataMember(b->ToType(), iName, oError);
               //only stop if we found it or some other error happened
               if(returnValue || parser::kNameDoesNotExist != oError) {
                  break;
               }
            }
         }
         if(returnValue && !returnValue.IsPublic()) {
            returnValue = Member();
            oError = parser::kIsNotPublic;
         }
      }
      if(returnValue) {
         oError = parser::kNoError;
      }
      return returnValue;
   }
}
