#ifndef Physics_PhysicsListMakerBase_h
#define Physics_PhysicsListMakerBase_h
// -*- C++ -*-
//
// Package:     Physics
// Class  :     PhysicsListMakerBase
// 
/**\class PhysicsListMakerBase PhysicsListMakerBase.h SimG4Core/Physics/interface/PhysicsListMakerBase.h

 Description: Base class for the 'maker' which creates PhysicsLists

 Usage:
    This class is the interface for creating a physics list and for connnecting
 the appropriate OSCAR signals to that physics list

*/
//
// Original Author:  Chris D Jones
//         Created:  Tue Nov 22 13:03:39 EST 2005
// $Id: PhysicsListMakerBase.h,v 1.1 2005/11/22 20:05:22 chrjones Exp $
//

// system include files
#include <memory>

// user include files

// forward declarations
class SimActivityRegistry;
namespace edm{
  class ParameterSet;
}

class PhysicsListMakerBase
{

   public:
      PhysicsListMakerBase() {}
      virtual ~PhysicsListMakerBase() {}

      // ---------- const member functions ---------------------
      virtual std::auto_ptr<PhysicsList> make(G4LogicalVolumeToDDLogicalPartMap&,
					      const edm::ParameterSet&,
					      SimActivityRegistry&) const = 0;

      // ---------- static member functions --------------------

      // ---------- member functions ---------------------------

   private:
      //PhysicsListMakerBase(const PhysicsListMakerBase&); // stop default

      //const PhysicsListMakerBase& operator=(const PhysicsListMakerBase&); // stop default

      // ---------- member data --------------------------------

};


#endif
