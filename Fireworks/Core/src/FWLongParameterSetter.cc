// -*- C++ -*-
//
// Package:     Core
// Class  :     FWLongParameterSetter
// 
// Implementation:
//     <Notes on implementation>
//
// Original Author:  Chris Jones
//         Created:  Mon Mar 10 11:22:32 CDT 2008
// $Id: FWLongParameterSetter.cc,v 1.1 2008/03/11 02:43:55 chrjones Exp $
//

// system include files
#include "TGLabel.h"
#include "TGNumberEntry.h"

#include <iostream>

// user include files
#include "Fireworks/Core/src/FWLongParameterSetter.h"
#include "Fireworks/Core/interface/FWLongParameter.h"


//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
FWLongParameterSetter::FWLongParameterSetter():
m_param(0),
m_widget(0)
{
}

// FWLongParameterSetter::FWLongParameterSetter(const FWLongParameterSetter& rhs)
// {
//    // do actual copying here;
// }

FWLongParameterSetter::~FWLongParameterSetter()
{
}

//
// assignment operators
//
// const FWLongParameterSetter& FWLongParameterSetter::operator=(const FWLongParameterSetter& rhs)
// {
//   //An exception safe implementation is
//   FWLongParameterSetter temp(rhs);
//   swap(rhs);
//
//   return *this;
// }

//
// member functions
//

void 
FWLongParameterSetter::attach(FWParameterBase* iParam)
{
   m_param = dynamic_cast<FWLongParameter*>(iParam);
   assert(0!=m_param);
}

TGFrame* 
FWLongParameterSetter::build(TGFrame* iParent)
{
   TGCompositeFrame* frame = new TGHorizontalFrame(iParent);
   
   // number entry widget
   TGNumberFormat::ELimit limits = m_param->min()==m_param->max() ? 
   TGNumberFormat::kNELNoLimits : 
   ( m_param->min() > m_param->max()? TGNumberFormat::kNELLimitMin : TGNumberFormat::kNELLimitMinMax);
   double min = 0;
   double max = 1;
   if(m_param->min()!=m_param->max()) {
      min=m_param->min();
      max=m_param->max();
   }
   m_widget = new TGNumberEntry
   (frame, m_param->value(),
    5,                                // number of digits
    0,                                // widget ID
    TGNumberFormat::kNESInteger,      // style
    TGNumberFormat::kNEAAnyNumber,    // input value filter
    limits,                           // specify limits
    min,                              // min value
    max);                             // max value
    
    frame->AddFrame(m_widget, new TGLayoutHints(kLHintsLeft|kLHintsCenterY, 2,8,2,2));
    m_widget->Connect("ValueSet(Long_t)", "FWLongParameterSetter", this, "doUpdate(Long_t)");
   
   // label
   frame->AddFrame(new TGLabel(frame,m_param->name().c_str()),
                   new TGLayoutHints(kLHintsLeft|kLHintsCenterY) );
   return frame;
}

void
FWLongParameterSetter::doUpdate(Long_t)
{
   //std::cout <<"doUpdate called"<<std::endl;
   assert(0!=m_param);
   assert(0!=m_widget);
   //std::cout <<m_widget->GetNumberEntry()->GetNumber()<<std::endl;
   m_param->set(m_widget->GetNumberEntry()->GetIntNumber());
   update();
}
//
// const member functions
//

//
// static member functions
//
