#ifndef _FUShmOutputModule_h
#define _FUShmOutputModule_h 

/*
   Description:
     Header file shared memory to be used with FUShmOutputModule.
     See CMS EvF Storage Manager wiki page for further notes.

   $Id: FUShmOutputModule.h,v 1.4 2008/10/14 13:19:50 biery Exp $
*/

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "IOPool/Streamer/interface/InitMsgBuilder.h"
#include "IOPool/Streamer/interface/EventMsgBuilder.h"

#include "EventFilter/ShmBuffer/interface/FUShmBuffer.h"

// Data structure to be shared by all output modules for event serialization
struct SM_SharedMemoryHandle
{
  SM_SharedMemoryHandle():
    shmBuffer_(0)
  { }

  evf::FUShmBuffer* getShmBuffer() {
   if(!shmBuffer_) {
     shmBuffer_ = evf::FUShmBuffer::getShmBuffer();
     return shmBuffer_;
   } else {
     return shmBuffer_;
   }
  }
  void detachShmBuffer() {
   if(!shmBuffer_) {
     // no shared memory was attached to!
   } else {
     shmdt(shmBuffer_);
     shmBuffer_ = 0;
   }
  }

  evf::FUShmBuffer* shmBuffer_;
};

namespace edm
{
  class ParameterSetDescription;
  class FUShmOutputModule
  {
  public:

    FUShmOutputModule(edm::ParameterSet const& ps);
    ~FUShmOutputModule();

    void doOutputHeader(InitMsgBuilder const& initMessage);
    void doOutputEvent(EventMsgBuilder const& eventMessage);
    void start();
    void stop();
    // No parameters.
    static void fillDescription(ParameterSetDescription&) {}

  private:

    evf::FUShmBuffer* shmBuffer_;

    static bool fuIdsInitialized_;
    static uint32 fuProcId_;
    static uint32 fuGuidValue_;

  };
}

#endif
