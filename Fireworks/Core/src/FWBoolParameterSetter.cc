// -*- C++ -*-
//
// Package:     Core
// Class  :     FWBoolParameterSetter
// 
// Implementation:
//     <Notes on implementation>
//
// Original Author:  Chris Jones
//         Created:  Mon Mar 10 11:22:32 CDT 2008
// $Id: FWBoolParameterSetter.cc,v 1.2 2008/05/18 09:42:48 jmuelmen Exp $
//

// system include files
#include "TGLabel.h"
#include "TGButton.h"
#include <assert.h>
#include <iostream>

// user include files
#include "Fireworks/Core/src/FWBoolParameterSetter.h"
#include "Fireworks/Core/interface/FWBoolParameter.h"


//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
FWBoolParameterSetter::FWBoolParameterSetter():
m_param(0),
m_widget(0)
{
}

// FWBoolParameterSetter::FWBoolParameterSetter(const FWBoolParameterSetter& rhs)
// {
//    // do actual copying here;
// }

FWBoolParameterSetter::~FWBoolParameterSetter()
{
}

//
// assignment operators
//
// const FWBoolParameterSetter& FWBoolParameterSetter::operator=(const FWBoolParameterSetter& rhs)
// {
//   //An exception safe implementation is
//   FWBoolParameterSetter temp(rhs);
//   swap(rhs);
//
//   return *this;
// }

//
// member functions
//

void 
FWBoolParameterSetter::attach(FWParameterBase* iParam)
{
   m_param = dynamic_cast<FWBoolParameter*>(iParam);
   assert(0!=m_param);
}

TGFrame* 
FWBoolParameterSetter::build(TGFrame* iParent)
{
   TGCompositeFrame* frame = new TGHorizontalFrame(iParent);
   
   m_widget = new TGCheckButton(frame, m_param->name().c_str(), 0);
   m_widget->Connect("Clicked()", "FWBoolParameterSetter", this, "doUpdate()");
   frame->AddFrame(m_widget, new TGLayoutHints(kLHintsLeft|kLHintsCenterY,2,0,1,1));
   return frame;
}

void
FWBoolParameterSetter::doUpdate()
{
   assert(0!=m_param);
   assert(0!=m_widget);
   m_param->set(m_widget->IsOn());
   update();
}
//
// const member functions
//

//
// static member functions
//
