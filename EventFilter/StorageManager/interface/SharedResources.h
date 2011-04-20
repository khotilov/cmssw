// $Id: SharedResources.h,v 1.6.10.2 2011/02/28 17:56:15 mommsen Exp $
/// @file: SharedResources.h 

#ifndef EventFilter_StorageManager_SharedResources_h
#define EventFilter_StorageManager_SharedResources_h

#include <string>

#include "boost/shared_ptr.hpp"

#include "EventFilter/StorageManager/interface/CommandQueue.h"
#include "EventFilter/StorageManager/interface/DQMEventQueue.h"
#include "EventFilter/StorageManager/interface/DQMEventQueueCollection.h"
#include "EventFilter/StorageManager/interface/EventQueueCollection.h"
#include "EventFilter/StorageManager/interface/FragmentQueue.h"
#include "EventFilter/StorageManager/interface/RegistrationQueue.h"
#include "EventFilter/StorageManager/interface/StreamQueue.h"


namespace stor {

  class Configuration;
  class DiscardManager;
  class DiskWriterResources;
  class DQMEventProcessorResources;
  class InitMsgCollection;
  class RegistrationCollection;
  class SharedResources;
  class StatisticsReporter;


  /**
   * Container for shared resources.
   *
   * $Author: mommsen $
   * $Revision: 1.6.10.2 $
   * $Date: 2011/02/28 17:56:15 $
   */

  struct SharedResources
  {

    // queues
    CommandQueuePtr commandQueue_;
    DQMEventQueuePtr dqmEventQueue_;
    FragmentQueuePtr fragmentQueue_;
    StreamQueuePtr streamQueue_;
    RegistrationQueuePtr registrationQueue_;
    EventQueueCollectionPtr eventQueueCollection_;
    DQMEventQueueCollectionPtr dqmEventQueueCollection_;

    // other
    boost::shared_ptr<Configuration> configuration_;
    boost::shared_ptr<DiscardManager> discardManager_;
    boost::shared_ptr<DiskWriterResources> diskWriterResources_;
    boost::shared_ptr<DQMEventProcessorResources> dqmEventProcessorResources_;
    boost::shared_ptr<InitMsgCollection> initMsgCollection_;
    boost::shared_ptr<StatisticsReporter> statisticsReporter_;
    boost::shared_ptr<RegistrationCollection> registrationCollection_;

    /**
     * Add a Failed state-machine event to the command queue
     */
    void moveToFailedState( xcept::Exception& );

    /**
       Write message to a file in /tmp
       (last resort when everything else fails)
    */
    void localDebug( const std::string& message ) const;

  };

  typedef boost::shared_ptr<SharedResources> SharedResourcesPtr;
  
} // namespace stor

#endif // EventFilter_StorageManager_SharedResources_h 


/// emacs configuration
/// Local Variables: -
/// mode: c++ -
/// c-basic-offset: 2 -
/// indent-tabs-mode: nil -
/// End: -
