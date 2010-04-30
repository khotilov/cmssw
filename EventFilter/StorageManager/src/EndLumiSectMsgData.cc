// $Id: EndLumiSectMsgData.cc,v 1.3.2.1 2010/04/22 14:09:04 mommsen Exp $
/// @file: EndLumiSectMsgData.cc

#include "EventFilter/StorageManager/src/ChainData.h"

#include "interface/evb/version.h"
#include "interface/evb/i2oEVBMsgs.h"
#include "interface/shared/i2oXFunctionCodes.h"
#include "interface/shared/version.h"


namespace stor
{

  namespace detail
  {

    EndLumiSectMsgData::EndLumiSectMsgData(toolbox::mem::Reference* pRef):
      #if (INTERFACESHARED_VERSION_MAJOR*1000 + INTERFACESHARED_VERSION_MINOR)>1010
      ChainData(I2O_EVM_LUMISECTION),
      #else
      ChainData(),
      #endif
      _runNumber(0),
      _lumiSection(0)
    {
      addFirstFragment(pRef);

      #if (INTERFACEEVB_VERSION_MAJOR*1000 + INTERFACEEVB_VERSION_MINOR)>1008

      if (validateDataLocation(pRef, INVALID_INITIAL_REFERENCE) &&
          validateMessageCode(pRef, _i2oMessageCode))
      {
        I2O_EVM_END_OF_LUMISECTION_MESSAGE_FRAME* msg_frame =
          (I2O_EVM_END_OF_LUMISECTION_MESSAGE_FRAME*)( pRef->getDataLocation() );
        if (msg_frame)
        {
          _runNumber = msg_frame->runNumber;
          _lumiSection = msg_frame->lumiSection;
        }
      }

      #endif
    }

    uint32 EndLumiSectMsgData::do_runNumber() const
    {
      return _runNumber;
    }

    uint32 EndLumiSectMsgData::do_lumiSection() const
    {
      return _lumiSection;
    }

  } // namespace detail

} // namespace stor



/// emacs configuration
/// Local Variables: -
/// mode: c++ -
/// c-basic-offset: 2 -
/// indent-tabs-mode: nil -
/// End: -
