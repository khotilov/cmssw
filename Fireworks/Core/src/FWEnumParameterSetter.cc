// -*- C++ -*-
//
// Package:     Core
// Class  :     FWEnumParameterSetter
// 
// Implementation:
//     [Notes on implementation]
//
// Original Author:  matevz
//         Created:  Fri Apr 30 15:17:33 CEST 2010
// $Id: FWEnumParameterSetter.cc,v 1.1 2010/04/30 15:29:44 matevz Exp $
//

// system include files

// user include files
#include "Fireworks/Core/interface/FWEnumParameterSetter.h"
#include "TGComboBox.h"
#include "TGLabel.h"

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
FWEnumParameterSetter::FWEnumParameterSetter() :
   m_param(0),
   m_widget(0)
{}

// FWEnumParameterSetter::FWEnumParameterSetter(const FWEnumParameterSetter& rhs)
// {
//    // do actual copying here;
// }

FWEnumParameterSetter::~FWEnumParameterSetter()
{}

//
// assignment operators
//
// const FWEnumParameterSetter& FWEnumParameterSetter::operator=(const FWEnumParameterSetter& rhs)
// {
//   //An exception safe implementation is
//   FWEnumParameterSetter temp(rhs);
//   swap(rhs);
//
//   return *this;
// }

//
// member functions
//

void
FWEnumParameterSetter::attach(FWParameterBase* iParam)
{
   m_param = dynamic_cast<FWEnumParameter*>(iParam);
   assert(0!=m_param);
}

TGFrame*
FWEnumParameterSetter::build(TGFrame* iParent)
{
   TGCompositeFrame *frame = new TGHorizontalFrame(iParent);

   m_widget = new TGComboBox(iParent);
   std::map<Long_t, std::string>::const_iterator me = m_param->entryMap().begin();
   UInt_t max_len = 0;
   while (me != m_param->entryMap().end())
   {
      m_widget->AddEntry(me->second.c_str(), static_cast<Int_t>(me->first));
      if (me->second.length() > max_len) max_len = me->second.length();
      ++me;
   }
   frame->AddFrame(m_widget, new TGLayoutHints(kLHintsLeft|kLHintsCenterY, 2,8,2,2));
   m_widget->Resize(8*max_len, 20);
   m_widget->Select(static_cast<Int_t>(m_param->value()), kFALSE);

   m_widget->Connect("Selected(Int_t)", "FWEnumParameterSetter", this, "doUpdate(Int_t)");

   // label
   frame->AddFrame(new TGLabel(frame, m_param->name().c_str()),
                   new TGLayoutHints(kLHintsLeft|kLHintsCenterY) );
   return frame;
}

void
FWEnumParameterSetter::doUpdate(Int_t id)
{
   assert(0!=m_param);
   assert(0!=m_widget);
   m_param->set((Long_t) id);
   update();
}

//
// const member functions
//

//
// static member functions
//
