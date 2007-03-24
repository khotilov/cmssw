#ifndef FUSHMBUFFER_H
#define FUSHMBUFFER_H 1


#include "EventFilter/ShmBuffer/interface/FUShmRawCell.h"
#include "EventFilter/ShmBuffer/interface/FUShmRecoCell.h"
#include "EventFilter/ShmBuffer/interface/FUShmDqmCell.h"


#include <sys/ipc.h>
#include <sys/types.h>
#include <sys/shm.h>
#include <sys/sem.h>
#include <errno.h>


namespace evf {
  
  // define event data states
  namespace evt {
    enum State_t { EMPTY,
		   RAWWRITING, RAWWRITTEN,
		   RAWREADING, RAWREAD,
		   PROCESSING, PROCESSED,
		   RECOWRITING,RECOWRITTEN,
		   SENDING,    SENT,
		   DISCARDING };
  }
  
  // define dqm data states
  namespace dqm {
    enum State_t { EMPTY,
		   WRITING, WRITTEN,
		   SENDING, SENT,
		   DISCARDING };
  }
  
  
  class FUShmBuffer
  {
    //
    // construction/destruction [-> static 'createShmBuffer'/'getShmBuffer']
    //
  private:
    FUShmBuffer(bool         segmentationMode,
		unsigned int nRawCells,
		unsigned int nRecoCells,
		unsigned int nDqmCells,
		unsigned int rawCellSize,
		unsigned int recoCellSize,
		unsigned int dqmCellSize);
  public:
    ~FUShmBuffer();
    
    
  public:
    //
    // public member functions
    //
    void           initialize(unsigned int shmid,unsigned int semid);
    void           reset();
    
    int            shmid()       const { return shmid_; }
    int            semid()       const { return semid_; }
    unsigned int   nbClients();

    evt::State_t   evtState(unsigned int index);
    dqm::State_t   dqmState(unsigned int index);
    
    int            nbRawCellsToWrite()  const;
    int            nbRawCellsToRead()   const;

    FUShmRawCell*  rawCellToWrite();
    FUShmRawCell*  rawCellToRead();
    FUShmRecoCell* recoCellToRead();
    FUShmDqmCell*  dqmCellToRead();
    FUShmRawCell*  rawCellToDiscard();
    
    void           finishWritingRawCell(FUShmRawCell* cell);
    void           finishReadingRawCell(FUShmRawCell* cell);
    void           finishReadingRecoCell(FUShmRecoCell* cell);
    void           finishReadingDqmCell(FUShmDqmCell* cell);
    
    void           scheduleRawCellForDiscard(unsigned int iCell);
    
    void           discardRawCell(FUShmRawCell* cell);
    void           discardRecoCell(unsigned int iCell);
    void           discardDqmCell(unsigned int iCell);
    
    void           releaseRawCell(FUShmRawCell* cell);
    
    void           writeRawEmptyEvent();
    void           writeRecoEmptyEvent();
    void           writeDqmEmptyEvent();
    
    void           scheduleRawEmptyCellForDiscard();
    void           scheduleRawEmptyCellForDiscard(FUShmRawCell* cell);
    
    bool           writeRecoInitMsg(unsigned char *data,
				    unsigned int   dataSize);

    bool           writeRecoEventData(unsigned int   runNumber,
				      unsigned int   evtNumber,
				      unsigned char *data,
				      unsigned int   dataSize);
    
    bool           writeDqmEventData(unsigned int   runNumber,
				     unsigned int   evtAtUpdate,
				     unsigned int   folderId,
				     unsigned char *data,
				     unsigned int   dataSize);
				     
    void           sem_print();
    void           printEvtState(unsigned int index);
    void           printDqmState(unsigned int index);
    
    
    //
    // static member functions
    //
    static FUShmBuffer* createShmBuffer(bool         semgmentationMode,
					unsigned int nRawCells,
					unsigned int nRecoCells,
					unsigned int nDqmCells,
					unsigned int rawCellSize =0x400000,  //4MB
					unsigned int recoCellSize=0x400000,  //4MB
					unsigned int dqmCellSize =0x400000); //4MB

    static FUShmBuffer* getShmBuffer();
    
    static bool         releaseSharedMemory();
    
    static unsigned int size(bool         segmentationMode,
			     unsigned int nRawCells,
			     unsigned int nRecoCells,
			     unsigned int nDqmCells,
			     unsigned int rawCellSize,
			     unsigned int recoCellSize,
			     unsigned int dqmCellSize);
    
    
    static key_t getShmDescriptorKey();
    static key_t getShmKey();
    static key_t getSemKey();
    
    static int   shm_create(key_t key,int size);
    static int   shm_get(key_t key,int size);
    static void* shm_attach(int shmid);
    static int   shm_nattch(int shmid);
    static int   shm_destroy(int shmid);
    
    static int   sem_create(key_t key,int nsem);
    static int   sem_get(key_t key,int nsem);
    static int   sem_destroy(int semid);
    
    
  private:
    //
    // private member functions
    //
    unsigned int   nextIndex(unsigned int  offset,
			     unsigned int  nCells,
			     unsigned int& iNext);
    void           postIndex(unsigned int  index,
			     unsigned int  offset,
			     unsigned int  nCells,
			     unsigned int& iLast);

    unsigned int   nextRawWriteIndex();
    unsigned int   nextRawReadIndex();
    void           postRawIndexToWrite(unsigned int index);
    void           postRawIndexToRead(unsigned int index);
    
    unsigned int   nextRecoWriteIndex();
    unsigned int   nextRecoReadIndex();
    void           postRecoIndexToWrite(unsigned int index);
    void           postRecoIndexToRead(unsigned int index);
    
    unsigned int   nextDqmWriteIndex();
    unsigned int   nextDqmReadIndex();
    void           postDqmIndexToWrite(unsigned int index);
    void           postDqmIndexToRead(unsigned int index);
    
    unsigned int   indexForEvtNumber(unsigned int evtNumber);
    unsigned int   indexForRecoCellId(unsigned int recoCellId);
    unsigned int   evtNumber(unsigned int index);
    unsigned int   evtRecoId(unsigned int index);
    
    bool           setEvtState(unsigned int index,evt::State_t state);
    bool           setEvtDiscard(unsigned int index,unsigned int discard);
    bool           setEvtNumber(unsigned int index,unsigned int evtNumber);
    bool           setRecoCellId(unsigned int index,unsigned int recoCellId);
    bool           setDqmState(unsigned int index,dqm::State_t state);
    
    FUShmRawCell*  rawCell(unsigned int iCell);
    FUShmRecoCell* recoCell(unsigned int iCell);
    FUShmDqmCell*  dqmCell(unsigned int iCell);
    
    bool           rawCellReadyForDiscard(unsigned int index);

    key_t          shmKey(unsigned int iCell,unsigned int offset);
    key_t          rawCellShmKey(unsigned int iCell);
    key_t          recoCellShmKey(unsigned int iCell);
    key_t          dqmCellShmKey(unsigned int iCell);

    void           sem_init(int isem,int value);
    void           sem_wait(int isem);
    void           sem_post(int isem);

    void           lock()             { sem_wait(0); }
    void           unlock()           { sem_post(0); }
    void           waitRawWrite()     { sem_wait(1); }
    void           postRawWrite()     { sem_post(1); }
    void           waitRawRead()      { sem_wait(2); }
    void           postRawRead()      { sem_post(2); }
    void           waitRawDiscard()   { sem_wait(3); }
    void           postRawDiscard()   { sem_post(3); }
    void           waitRawDiscarded() { sem_wait(4); }
    void           postRawDiscarded() { sem_post(4); }
    void           waitRecoWrite()    { sem_wait(5); }
    void           postRecoWrite()    { sem_post(5); }
    void           waitRecoRead()     { sem_wait(6); }
    void           postRecoRead()     { sem_post(6); }
    void           waitDqmWrite()     { sem_wait(7); }
    void           postDqmWrite()     { sem_post(7); }
    void           waitDqmRead()      { sem_wait(8); }
    void           postDqmRead()      { sem_post(8); }

    
  private:
    //
    // member data
    //
    bool         segmentationMode_;
    int          shmid_;
    int          semid_;

    unsigned int rawWriteNext_;
    unsigned int rawWriteLast_;
    unsigned int rawWriteOffset_;
    unsigned int rawReadNext_;
    unsigned int rawReadLast_;
    unsigned int rawReadOffset_;
    unsigned int rawDiscardIndex_;
    
    unsigned int recoWriteNext_;
    unsigned int recoWriteLast_;
    unsigned int recoWriteOffset_;
    unsigned int recoReadNext_;
    unsigned int recoReadLast_;
    unsigned int recoReadOffset_;
    
    unsigned int dqmWriteNext_;
    unsigned int dqmWriteLast_;
    unsigned int dqmWriteOffset_;
    unsigned int dqmReadNext_;
    unsigned int dqmReadLast_;
    unsigned int dqmReadOffset_;
    
    unsigned int evtStateOffset_;
    unsigned int evtDiscardOffset_;
    unsigned int evtNumberOffset_;
    unsigned int evtRecoIdOffset_;
    unsigned int dqmStateOffset_;
    
    unsigned int nRawCells_;
    unsigned int rawCellPayloadSize_;
    unsigned int rawCellTotalSize_;
    unsigned int rawCellOffset_;
    
    unsigned int recoWriteIndex_;
    unsigned int recoReadIndex_;
    unsigned int nRecoCells_;
    unsigned int recoCellPayloadSize_;
    unsigned int recoCellTotalSize_;
    unsigned int recoCellOffset_;
    
    unsigned int dqmWriteIndex_;
    unsigned int dqmReadIndex_;
    unsigned int nDqmCells_;
    unsigned int dqmCellPayloadSize_;
    unsigned int dqmCellTotalSize_;
    unsigned int dqmCellOffset_;
    
  };

  
} // namespace evf


#endif
