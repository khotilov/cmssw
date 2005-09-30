
/* \file FEDTrailer.cc
 *
 *  $Date: $
 *  $Revision: $
 *  \author N. Amapane - CERN
 */

#include <DataFormats/FEDRawData/interface/FEDTrailer.h>
#include "fed_trailer.h"

FEDTrailer::FEDTrailer(unsigned char* trailer) :
  theTrailer((fedt_t*)(trailer)) {}


FEDTrailer::~FEDTrailer(){}



int FEDTrailer::lenght(){
  return (theTrailer->eventsize & FED_EVSZ_MASK) >> FED_EVSZ_SHIFT;
}


int FEDTrailer::crc(){
  return ((theTrailer->conscheck & FED_CRCS_MASK) >> FED_CRCS_SHIFT);
}


int FEDTrailer::evtStatus(){
  return ((theTrailer->conscheck & FED_STAT_MASK) >> FED_STAT_SHIFT);
}


int FEDTrailer::ttsBits(){
  return ((theTrailer->conscheck & FED_TTSI_MASK) >> FED_TTSI_SHIFT);  
}

 
bool FEDTrailer::moreTrailers(){
  return ((theTrailer->conscheck & FED_MORE_TRAILERS)!=0);
}


void FEDTrailer::set(int evt_lgth,
		     int crc,  
		     int evt_stat,
		     int tts,
		     bool T){
  // FIXME: should check that input ranges are OK!!!
  theTrailer->eventsize = 
    FED_TCTRLID |
    evt_lgth << FED_EVSZ_SHIFT;
 
  theTrailer->conscheck = 
    crc       << FED_CRCS_SHIFT |
    evt_stat  << FED_STAT_SHIFT |
    tts       << FED_TTSI_SHIFT;
  
  if (T) theTrailer->conscheck |= FED_MORE_TRAILERS;
}


bool FEDTrailer::check(){
  // ...may report with finer detail
  bool result = true;
  result &= ((theTrailer->eventsize & FED_TCTRLID_MASK) == FED_TCTRLID);

  return result;
}
