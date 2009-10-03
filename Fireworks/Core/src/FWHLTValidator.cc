// -*- C++ -*-
//
// Package:     Core
// Class  :     FWHLTValidator
// $Id: FWHLTValidator.cc,v 1.1 2009/10/02 17:54:41 dmytro Exp $
//

// system include files
#include <algorithm>
#include <cstring>
#include <boost/regex.hpp>

// user include files
#include "Fireworks/Core/interface/FWHLTValidator.h"
#include "CommonTools/Utils/src/returnType.h"
#include "DataFormats/FWLite/interface/Handle.h"
#include "DataFormats/FWLite/interface/TriggerNames.h"
#include "DataFormats/Common/interface/TriggerResults.h"

void
FWHLTValidator::fillOptions(const char* iBegin, const char* iEnd,
			    std::vector<std::pair<boost::shared_ptr<std::string>, std::string> >& oOptions) const
{
   oOptions.clear();
   std::string part(iBegin,iEnd);
   part = boost::regex_replace(part,boost::regex(".*?(\\&\\&|\\|\\||\\s)+"),"");

   if (m_triggerNames.empty()){
     fwlite::Handle<edm::TriggerResults> hTriggerResults;
     fwlite::TriggerNames const* triggerNames(0);
     try{
       hTriggerResults.getByLabel(m_event,"TriggerResults","","HLT");
       triggerNames = &m_event.triggerNames(*hTriggerResults);
     } catch (...){
       std::cout << "Warning: no trigger results with process name HLT is available" << std::endl;
       return;
     }
     for(unsigned int i=0; i<triggerNames->size(); ++i)
       m_triggerNames.push_back(triggerNames->triggerName(i));
     std::sort(m_triggerNames.begin(),m_triggerNames.end());
   }

   //only use add items which begin with the part of the member we are trying to match
   unsigned int part_size = part.size();
   for(std::vector<std::string>::const_iterator trigger = m_triggerNames.begin();
       trigger != m_triggerNames.end(); ++trigger)
     if(part == trigger->substr(0,part_size) ) {
       oOptions.push_back(std::make_pair(boost::shared_ptr<std::string>(new std::string(*trigger)),
					 trigger->substr(part_size,trigger->size()-part_size)));
     }
}

//
// static member functions
//
