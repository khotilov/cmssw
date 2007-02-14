#include "EventFilter/Unit/interface/EventSink.h"
#include "EventFilter/Unit/interface/FilterUnitFramework.h"
#include "EventFilter/Unit/interface/FURawEvent.h"

FUAdapter *EventSink::fwk_ = 0;
bool EventSink::sinking_ = false;


void EventSink::sink(){
  FURawEvent *event = fwk_->rqstEvent();
  event->reset(true);
}

