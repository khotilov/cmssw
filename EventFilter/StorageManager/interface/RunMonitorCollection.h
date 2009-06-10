// $Id$

#ifndef StorageManager_RunMonitorCollection_h
#define StorageManager_RunMonitorCollection_h

#include "xdata/UnsignedInteger32.h"

#include "EventFilter/StorageManager/interface/MonitorCollection.h"


namespace stor {

  /**
   * A collection of MonitoredQuantities related to events received
   * in the current run
   *
   * $Author$
   * $Revision$
   * $Date$
   */
  
  class RunMonitorCollection : public MonitorCollection
  {
  private:

    MonitoredQuantity _eventIDsReceived;
    MonitoredQuantity _errorEventIDsReceived;
    MonitoredQuantity _runNumbersSeen;  // Does this make sense?
    MonitoredQuantity _lumiSectionsSeen;


  public:

    explicit RunMonitorCollection(xdaq::Application*);

    const MonitoredQuantity& getEventIDsReceivedMQ() const {
      return _eventIDsReceived;
    }
    MonitoredQuantity& getEventIDsReceivedMQ() {
      return _eventIDsReceived;
    }

    const MonitoredQuantity& getErrorEventIDsReceivedMQ() const {
      return _errorEventIDsReceived;
    }
    MonitoredQuantity& getErrorEventIDsReceivedMQ() {
      return _errorEventIDsReceived;
    }

    const MonitoredQuantity& getRunNumbersSeenMQ() const {
      return _runNumbersSeen;
    }
    MonitoredQuantity& getRunNumbersSeenMQ() {
      return _runNumbersSeen;
    }

    const MonitoredQuantity& getLumiSectionsSeenMQ() const {
      return _lumiSectionsSeen;
    }
    MonitoredQuantity& getLumiSectionsSeenMQ() {
      return _lumiSectionsSeen;
    }


  private:

    //Prevent copying of the RunMonitorCollection
    RunMonitorCollection(RunMonitorCollection const&);
    RunMonitorCollection& operator=(RunMonitorCollection const&);

    virtual void do_calculateStatistics();
    
    virtual void do_updateInfoSpace();
    
    virtual void do_reset();

    xdata::UnsignedInteger32 _runNumber;           // The current run number

    // InfoSpace items which were defined in the old SM
    // xdata::UnsignedInteger32 _receivedEvents;      // Total number of received events
    // xdata::UnsignedInteger32 _receivedErrorEvents; // Total number of received error events

  };
  
} // namespace stor

#endif // StorageManager_RunMonitorCollection_h 


/// emacs configuration
/// Local Variables: -
/// mode: c++ -
/// c-basic-offset: 2 -
/// indent-tabs-mode: nil -
/// End: -
