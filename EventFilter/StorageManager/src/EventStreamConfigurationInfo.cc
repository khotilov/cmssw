// $Id: EventStreamConfigurationInfo.cc,v 1.4 2009/09/11 21:07:06 elmer Exp $
/// @file: EventStreamConfigurationInfo.cc

#include "EventFilter/StorageManager/interface/EventStreamConfigurationInfo.h"
#include <ostream>

using stor::EventStreamConfigurationInfo;
using namespace std;

ostream&
stor::operator << ( ostream& os,
		    const EventStreamConfigurationInfo& ci )
{

  os << "EventStreamConfigurationInfo:" << endl
     << " Stream label: " << ci.streamLabel() << endl
     << " Maximum file size, MB: " << ci.maxFileSizeMB() << endl
     << " HLT output: " << ci.outputModuleLabel() << endl
     << " Maximum event size, Bytes: " << ci.maxEventSize() << endl
     << " Fraction of events written to disk: " << ci.fractionToDisk() << endl
     << " Stream Id: " << ci.streamId() << endl;

  if( ci.useCompression() )
    {
      os << " Compression applied at level " << ci.compressionLevel() << endl;
    }
  else
    {
      os << " Compression not applied" << endl;
    }

  os << " Event filters:";
  for( unsigned int i = 0; i < ci.selEvents().size(); ++i )
    {
      os << endl << "  " << ci.selEvents()[i];
    }

    return os;

}
