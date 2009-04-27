#ifndef Fireworks_Core_FWDetailViewFactory_h
#define Fireworks_Core_FWDetailViewFactory_h
// -*- C++ -*-
//
// Package:     Core
// Class  :     FWDetailViewFactory
//
/**\class FWDetailViewFactory FWDetailViewFactory.h Fireworks/Core/interface/FWDetailViewFactory.h

   Description: <one line class summary>

   Usage:
    <usage>

 */
//
// Original Author:  Chris Jones
//         Created:  Mon Jan 12 09:48:04 EST 2009
// $Id: FWDetailViewFactory.h,v 1.2 2009/01/23 21:35:41 amraktad Exp $
//

// system include files
#include "FWCore/PluginManager/interface/PluginFactory.h"

// user include files

// forward declarations

class FWDetailViewBase;

typedef edmplugin::PluginFactory<FWDetailViewBase*()> FWDetailViewFactory;

#define REGISTER_FWDETAILVIEW(_name_) \
   DEFINE_EDM_PLUGIN(FWDetailViewFactory,_name_,_name_::classRegisterTypeName()+"@" # _name_)


#endif
