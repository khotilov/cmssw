#ifndef FURAWEVENTFACTORY_H
#define FURAWEVENTFACTORY_H

#include <vector>
#include <stack>
#include <iostream>
#include <pthread.h>
#include <semaphore.h>

class FURawEvent;
class FUAdapter;

/** class FURawEventFactory contains event resources which are distributed to 
    FU and recycled at the end of their lifetime */
class FURawEventFactory
{
  
 public:

  /** virtual destructor */
  virtual ~FURawEventFactory(){}

  /** remove this event */
  void removeEvent(FURawEvent *);

  /** get a built event */
  FURawEvent *getBuiltEvent();

  /** spy a built event */
  unsigned int spyBuiltEvent(unsigned int, unsigned char *);

  /** set minimum events in queue */
  void setMinQueueSize(int min)
    {
      minBuiltInQueue_ = min;
    } 

  /** get the size of the built queue */
  unsigned int queueSize() const;

 private:  

  //these methods are only used by friend class FUAdapter

  /** get an event with a given initial fset and bu instance*/
  virtual FURawEvent *getEvent(unsigned int, unsigned long);

  /** check out event corresponding to this transaction */
  FURawEvent *checkEvent(unsigned int, unsigned int);

  /** take this event as built (after a transaction is completed)*/
  void take(FURawEvent *);

  static int maxHandles;

  std::stack<int> freeRes_;
  std::stack<int> builtRes_;
  unsigned int minBuiltInQueue_;
  std::vector<FURawEvent *> resources_;

  /** constructor: the raw event factory is constructed 
   only by the FUAdapter*/
  FURawEventFactory();
  //@@EM TODO make this a real singleton anyway...

  sem_t empty_; // a POSIX counting semaphore
  pthread_mutex_t mutex_;

  static FUAdapter *adapter_;

  friend class FUAdapter;
};
#endif
