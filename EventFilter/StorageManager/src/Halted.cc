// $Id: Halted.cc,v 1.8 2009/08/28 16:41:26 mommsen Exp $
/// @file: Halted.cc

#include "EventFilter/StorageManager/interface/Notifier.h"
#include "EventFilter/StorageManager/interface/StateMachine.h"
#include "EventFilter/StorageManager/interface/TransitionRecord.h"

#include <iostream>

#include "xcept/tools.h"

using namespace std;
using namespace stor;

void Halted::do_entryActionWork()
{
  TransitionRecord tr( stateName(), true );
  outermost_context().updateHistory( tr );
  outermost_context().setExternallyVisibleState( "Halted" );
  outermost_context().getNotifier()->reportNewState( "Halted" );
}

Halted::Halted( my_context c ): my_base(c)
{
  safeEntryAction();
}

void Halted::do_exitActionWork()
{
  TransitionRecord tr( stateName(), false );
  outermost_context().updateHistory( tr );
}

Halted::~Halted()
{
  safeExitAction();
}

string Halted::do_stateName() const
{
  return string( "Halted" );
}

void Halted::do_moveToFailedState( xcept::Exception& exception ) const
{
  outermost_context().getSharedResources()->moveToFailedState( exception );
}



/// emacs configuration
/// Local Variables: -
/// mode: c++ -
/// c-basic-offset: 2 -
/// indent-tabs-mode: nil -
/// End: -
