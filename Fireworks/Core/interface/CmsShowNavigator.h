// -*- C++ -*-
#ifndef Fireworks_Core_CmsShowNavigator_h
#define Fireworks_Core_CmsShowNavigator_h
//
// Package:     newVersion
// Class  :     CmsShowNavigator
// $Id: CmsShowNavigator.h,v 1.25 2009/10/27 01:55:28 dmytro Exp $
//

// system include files
#include <string>
#include <sigc++/sigc++.h>

// user include files
#include "DataFormats/FWLite/interface/Event.h"
#include "Fireworks/Core/interface/FWEventSelector.h"
#include "Fireworks/Core/interface/FWConfigurable.h"
#include "Fireworks/Core/interface/FWFileEntry.h"
#include "TEventList.h"

// forward declarations
class TEventList;
class CSGAction;
class CmsShowMain;
class TFile;
class TGWindow;

namespace edm {
   class EventID;
}

class CmsShowNavigator : public FWConfigurable
{
public:

   CmsShowNavigator(const CmsShowMain &);
   virtual ~CmsShowNavigator();

   //configuration management interface
   void addTo(FWConfiguration&) const;
   void setFrom(const FWConfiguration&);

   Int_t realEntry(Int_t rawEntry);
   std::pair<std::deque<FWFileEntry>::iterator,Int_t> realEntry(Int_t run, Int_t event);

   bool loadFile(const std::string& fileName);
   void newRemoteFile(const std::string& fileName);
   void checkPosition();
   void nextEvent();
   void previousEvent();
   void firstEvent();
   void firstEventInTheCurrentFile();
   void lastEvent();
   void lastEventInTheCurrentFile();
   void enableEventFiltering( Bool_t);
   void filterEvents();
   void filterEventsAndReset();
   void goToEvent(Int_t,Int_t);
   void checkForNewFiles();
   void setMaxNumberOfFilesToChain( unsigned int i ){
      m_maxNumberOfFilesToChain = i;
   }
   void showEventFilter(const TGWindow*);

   sigc::signal<void, const fwlite::Event&> newEvent_;
   sigc::signal<void, const fwlite::Event&> oldEvent_;
   sigc::signal<void, const TFile*> fileChanged_;
   sigc::signal<void, bool> atBeginning_;
   sigc::signal<void, bool> atEnd_;

   sigc::signal<void> preFiltering_;
   sigc::signal<void> postFiltering_;

   sigc::signal<void, const std::string&> eventSelectionChanged_;

private:
   CmsShowNavigator(const CmsShowNavigator&);    // stop default

   const CmsShowNavigator& operator=(const CmsShowNavigator&);    // stop default

   void filterEvents(FWFileEntry&, int, std::string);
   bool filterEventsWithCustomParser(FWFileEntry& file, int, std::string);

   // ---------- member data --------------------------------
   unsigned int m_maxNumberOfFilesToChain;
   std::deque<FWFileEntry> m_files;
   std::deque<FWFileEntry>::iterator m_currentFile;
   std::vector<FWEventSelector*> m_selectors;
   std::deque<FWFileEntry>::iterator m_firstSelectedFile;
   std::deque<FWFileEntry>::iterator m_lastSelectedFile;
   std::deque<FWFileEntry>::iterator m_lastFile;
   bool m_filterEvents;
   bool m_globalOR;
   // entry is an event index nubmer which runs from 0 to
   // #events or #selected_events depending on if we filter
   // events or not
   int m_currentEntry;
   int m_lastEntry;
   const CmsShowMain &m_main;
   std::string m_nextFile;
};

#endif
