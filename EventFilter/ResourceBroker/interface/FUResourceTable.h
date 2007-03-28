#ifndef FURESOURCETABLE_H
#define FURESOURCETABLE_H 1


#include "EventFilter/ResourceBroker/interface/FUResource.h"
#include "EventFilter/ResourceBroker/interface/BUProxy.h"
#include "EventFilter/ResourceBroker/interface/SMProxy.h"
#include "EventFilter/ResourceBroker/interface/FUTypes.h"
#include "EventFilter/ShmBuffer/interface/FUShmBuffer.h"
#include "EventFilter/Utilities/interface/Exception.h"

#include "extern/log4cplus/linuxx86/include/log4cplus/logger.h"
#include "toolbox/include/toolbox/lang/Class.h"
#include "toolbox/include/toolbox/task/Action.h"
#include "toolbox/include/toolbox/task/WorkLoop.h"
#include "toolbox/include/BSem.h"

#include <vector>


namespace evf {
  
  class FUResourceTable : public toolbox::lang::Class
  {
  public:
    //
    // construction/destruction
    //
    FUResourceTable(bool   segmentationMode,
		    UInt_t nbRawCells, UInt_t nbRecoCells, UInt_t nbDqmCells,
		    UInt_t rawCellSize,UInt_t recoCellSize,UInt_t dqmCellSize,
		    BUProxy* bu,SMProxy* sm,
		    log4cplus::Logger logger);
    virtual ~FUResourceTable();
    
    
    //
    // member functions
    //
    
    // initialization of the resource queue
    void   initialize(bool   segmentationMode,
		      UInt_t nbRawCells, UInt_t nbRecoCells, UInt_t nbDqmCells,
		      UInt_t rawCellSize,UInt_t recoCellSize,UInt_t dqmCellSize);
    
    // work loop to send data events to storage manager
    void   startSendDataWorkLoop() throw (evf::Exception);
    bool   sendData(toolbox::task::WorkLoop* workLoop);
    
    // work loop to send dqm events to storage manager
    void   startSendDqmWorkLoop() throw (evf::Exception);
    bool   sendDqm(toolbox::task::WorkLoop* workLoop);
    
    // work loop to discard events to builder unit
    void   startDiscardWorkLoop() throw (evf::Exception);
    bool   discard(toolbox::task::WorkLoop* workLoop);
    
    // returns the fuResourceId of the allocated resource
    UInt_t allocateResource();
    
    // process buffer received via I2O_FU_TAKE message
    bool   buildResource(MemRef_t* bufRef);
    
    // process buffer received via I2O_SM_DATA_DISCARD message
    bool   discardDataEvent(MemRef_t* bufRef);
    
    // process buffer received via I2O_SM_DQM_DISCARD message
    bool   discardDqmEvent(MemRef_t* bufRef);
    
    // drop next available event
    void   dropEvent();
    
    // dump event to ascii file
    void   dumpEvent(evf::FUShmRawCell* cell);
    
    // send empty events to notify clients to shutdown
    void   shutDownClients();
    
    // emtpy all containers (resources & ids)
    void   clear();

    // reset the resource table to start over (in its current configuration)
    void   reset();

    // reset event & error counters
    void   resetCounters();

    // tell resources wether to check the crc
    void   setDoCrcCheck(UInt_t doCrcCheck) { doCrcCheck_=doCrcCheck; }

    // tell resources wether to dump events to an ascii file
    void   setDoDumpEvents(UInt_t doDumpEvents) { doDumpEvents_=doDumpEvents; }

    // check if resource table can be savely destroyed
    bool   isReadyToShutDown() const { return isReadyToShutDown_; }
    
    
    // various counters
    UInt_t   nbResources()        const { return resources_.size(); }
    UInt_t   nbFreeSlots()        const { return shmBuffer_->nbRawCellsToWrite(); }
    UInt_t   nbShmClients()       const;
    UInt_t   nbAllocated()        const { return nbAllocated_; }
    UInt_t   nbPending()          const { return nbPending_; }
    UInt_t   nbCompleted()        const { return nbCompleted_; }
    UInt_t   nbProcessed()        const { return nbProcessed_; }
    UInt_t   nbAccepted()         const { return nbAccepted_; }
    UInt_t   nbSent()             const { return nbSent_; }
    UInt_t   nbSentDqm()          const { return nbSentDqm_; }
    UInt_t   nbDiscarded()        const { return nbDiscarded_; }
    UInt_t   nbLost()             const { return nbLost_; }
    
    UInt_t   nbErrors()           const { return nbErrors_; }
    UInt_t   nbCrcErrors()        const { return nbCrcErrors_; }
    UInt_t   nbAllocSent()        const { return nbAllocSent_; }
    
    uint64_t inputSumOfSquares()  const { return inputSumOfSquares_; }
    uint64_t outputSumOfSquares() const { return outputSumOfSquares_; }
    UInt_t   inputSumOfSizes()    const { return inputSumOfSizes_; }
    UInt_t   outputSumOfSizes()   const { return outputSumOfSizes_; }
    
    
    //
    // helpers
    //
    void   sendAllocate();
    void   sendDiscard(UInt_t buResourceId);
    
    void   sendInitMessage(UInt_t  fuResourceId,
			   UChar_t*data,
			   UInt_t  dataSize);

    void   sendDataEvent(UInt_t  fuResourceId,
			 UInt_t  runNumber,
			 UInt_t  evtNumber,
			 UChar_t*data,
			 UInt_t  dataSize);

    void   sendDqmEvent(UInt_t  fuDqmId,
			UInt_t  runNumber,
			UInt_t  evtAtUpdate,
			UInt_t  folderId,
			UChar_t*data,
			UInt_t  dataSize);
    
    bool   isLastMessageOfEvent(MemRef_t* bufRef);
    
    void   lock()   { lock_.take(); }
    void   unlock() { lock_.give(); }
    

  private:
    //
    // member data
    //
    typedef toolbox::task::WorkLoop        WorkLoop_t;
    typedef toolbox::task::ActionSignature ActionSignature_t;

    BUProxy           *bu_;
    SMProxy           *sm_;

    log4cplus::Logger  log_;
    
    WorkLoop_t        *wlSendData_;
    ActionSignature_t *asSendData_;

    WorkLoop_t        *wlSendDqm_;
    ActionSignature_t *asSendDqm_;

    WorkLoop_t        *wlDiscard_;
    ActionSignature_t *asDiscard_;

    FUShmBuffer       *shmBuffer_;
    FUResourceVec_t    resources_;
    
    UInt_t             doCrcCheck_;
    UInt_t             doDumpEvents_;

    UInt_t             nbAllocated_;
    UInt_t             nbPending_;
    UInt_t             nbCompleted_;
    UInt_t             nbProcessed_;
    UInt_t             nbAccepted_;
    UInt_t             nbSent_;
    UInt_t             nbSentDqm_;
    UInt_t             nbDiscarded_;
    UInt_t             nbLost_;
    
    UInt_t             nbClientsToShutDown_;
    bool               isReadyToShutDown_;
    
    UInt_t             nbErrors_;
    UInt_t             nbCrcErrors_;
    UInt_t             nbAllocSent_;
    
    uint64_t           inputSumOfSquares_;
    uint64_t           outputSumOfSquares_;
    UInt_t             inputSumOfSizes_;
    UInt_t             outputSumOfSizes_;
    
    BSem               lock_;
    
  };
  
} // namespace evf


#endif
