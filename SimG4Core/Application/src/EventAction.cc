#include "SimG4Core/Application/interface/RunManager.h"
#include "SimG4Core/Application/interface/EventAction.h"
#include "SimG4Core/Application/interface/G4SimTrack.h"
#include "SimG4Core/Application/interface/G4SimVertex.h"
#include "SimG4Core/Notification/interface/BeginOfEvent.h"
#include "SimG4Core/Notification/interface/EndOfEvent.h"

#include <fstream>

using std::cout;
using std::endl;

EventAction::EventAction(const edm::ParameterSet & p) 
    : m_trackManager(p.getParameter<bool>("CollapsePrimaryVertices")),
      m_stopFile(p.getParameter<std::string>("StopFile")),
      m_debug(p.getParameter<bool>("debug"))
{}

EventAction::~EventAction() {}
    
void EventAction::BeginOfEventAction(const G4Event * anEvent)
{
    if (std::ifstream(m_stopFile.c_str()))
    {
        cout << "BeginOfEventAction: termination signal received at event "
             << anEvent->GetEventID() << endl;
        RunManager::instance()->abortRun(true);
    }

    m_trackManager.reset();
    BeginOfEvent e(anEvent);
    m_beginOfEventSignal(&e);
}

void EventAction::EndOfEventAction(const G4Event * anEvent)
{
    if (std::ifstream(m_stopFile.c_str()))
    {
        cout << "EndOfEventAction: termination signal received at event "
             << anEvent->GetEventID() << endl;
        RunManager::instance()->abortRun(true);
    }
    if (anEvent->GetNumberOfPrimaryVertex()==0)
    {
        cout << " EndOfEventAction: event " << anEvent->GetEventID()
             << " must have failed (no G4PrimaryVertices found) and will be skipped " << endl;
        return;
    }

    m_trackManager.storeTracks(RunManager::instance()->simEvent());
    // dispatch now end of event, and only then delete tracks...
    EndOfEvent e(anEvent);
    m_endOfEventSignal(&e);

    m_trackManager.deleteTracks();
}

void EventAction::addTrack(TrackWithHistory* iTrack)
{
  m_trackManager.addTrack(iTrack);
}
