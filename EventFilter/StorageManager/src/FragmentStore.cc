// $Id$

#include "EventFilter/StorageManager/interface/FragmentStore.h"
#include "EventFilter/StorageManager/interface/Utils.h"

using namespace stor;


const bool FragmentStore::addFragment(I2OChain &chain)
{
  // the trivial case that the chain is already complete
  if ( chain.complete() ) return true;

  const FragKey newKey = chain.fragmentKey();

  // Use efficientAddOrUpdates pattern suggested by Item 24 of 
  // 'Effective STL' by Scott Meyers
  fragmentMap::iterator pos = _store.lower_bound(newKey);

  if(pos != _store.end() && !(_store.key_comp()(newKey, pos->first)))
  {
    // key already exists
    pos->second.addToChain(chain);

    if ( pos->second.complete() )
    {
      chain = pos->second;
      _store.erase(pos);
      return true;
    }
  }
  else
  {
    // The key does not exist in the map, add it to the map
    // Use pos as a hint to insert, so it can avoid another lookup
    _store.insert(pos, fragmentMap::value_type(newKey, chain));
    chain.release();

    // We already handled the trivial case that the chain is complete.
    // Thus, _store will not have a complete event.
  }

  return false;
}

const bool FragmentStore::getStaleEvent(I2OChain &chain, double timeout)
{
  const double cutOffTime = utils::getCurrentTime() - timeout;
  
  fragmentMap::iterator pos = _store.begin();
  fragmentMap::iterator end = _store.end();

  while ( (pos != end) && (pos->second.lastFragmentTime() > cutOffTime ) )
  {
    ++pos;
  }

  if ( pos == end )
  {
    chain.release();
    return false;
  }
  else
  {
    chain = pos->second;
    _store.erase(pos);
    chain.markFaulty();
    return true;
  }
}



/// emacs configuration
/// Local Variables: -
/// mode: c++ -
/// c-basic-offset: 2 -
/// indent-tabs-mode: nil -
/// End: -
