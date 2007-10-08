/** \file
 *
 *  $Date: 2007/10/08 22:01:06 $
 *  $Revision: 1.11 $
 *  \author A. Tumanov - Rice
 */

#include "DataFormats/MuonDetId/interface/CSCDetId.h"
#include "EventFilter/CSCRawToDigi/src/CSCDigiToRaw.h"
#include "EventFilter/CSCRawToDigi/interface/CSCEventData.h"
#include "EventFilter/CSCRawToDigi/interface/CSCDCCEventData.h"
#include "DataFormats/CSCDigi/interface/CSCStripDigiCollection.h"
#include "DataFormats/CSCDigi/interface/CSCWireDigiCollection.h"
#include "DataFormats/FEDRawData/interface/FEDRawDataCollection.h"
#include "DataFormats/FEDRawData/interface/FEDRawData.h"
#include "DataFormats/FEDRawData/interface/FEDNumbering.h"
#include <boost/dynamic_bitset.hpp>
#include "EventFilter/CSCRawToDigi/src/bitset_append.h"
#include <FWCore/Framework/interface/Event.h>
#include <DataFormats/FEDRawData/interface/FEDHeader.h>
#include <DataFormats/FEDRawData/interface/FEDTrailer.h>
#include "EventFilter/Utilities/interface/Crc.h"
#include "CondFormats/CSCObjects/interface/CSCChamberMap.h"
#include "CondFormats/DataRecord/interface/CSCChamberMapRcd.h"



using namespace edm;
using namespace std;

CSCDigiToRaw::CSCDigiToRaw(){}


CSCDigiToRaw::~CSCDigiToRaw(){}

map<CSCDetId, CSCEventData> 
CSCDigiToRaw::fillChamberDataMap(const CSCStripDigiCollection & stripDigis, 
				 const CSCWireDigiCollection & wireDigis, 
				 const CSCChamberMap* mapping) 
{
  map<CSCDetId, CSCEventData> chamberMap;
  ///iterate over chambers with strip digis in them
  for (CSCStripDigiCollection::DigiRangeIterator j=stripDigis.begin(); j!=stripDigis.end(); ++j)
    {
      CSCDetId const cscDetId=(*j).first;
      CSCDetId chamberID =cscDetId.chamberId();
      /// find the entry into the map
      map<CSCDetId, CSCEventData>::iterator chamberMapItr = chamberMap.find(chamberID);
      if(chamberMapItr == chamberMap.end())
	{
	  /// make an entry, telling it the correct chamberType
	  int istation = cscDetId.station();
	  int iring = cscDetId.ring();
	  int chamberType = 2 * istation + iring; // i=2S+R
	  if ( istation == 1 )
	    {
	      --chamberType;             // ring 1R -> i=1+R (2S+R-1=1+R for S=1)
	      if ( chamberType > 4 ) 
		{
		  chamberType = 1;       // But ring 1A (R=4) -> i=1
		}
	    }   
	  chamberMapItr = chamberMap.insert(pair<CSCDetId, CSCEventData>(chamberID, CSCEventData(chamberType))).first;
	}    
      CSCEventData & cscData = chamberMapItr->second;
      cscData.dmbHeader().setCrateAddress(mapping->crate(cscDetId), mapping->dmb(cscDetId));
      ///add strip digis to that chamber
      std::vector<CSCStripDigi>::const_iterator digiItr = (*j).second.first;
      std::vector<CSCStripDigi>::const_iterator last = (*j).second.second;
      for( ; digiItr != last; ++digiItr) 
	{
	  cscData.add(*digiItr, cscDetId.layer() );    
	}
    }
  ///repeat the same for wire digis
  for (CSCWireDigiCollection::DigiRangeIterator j=wireDigis.begin(); j!=wireDigis.end(); ++j) 
    {
      CSCDetId const cscDetId=(*j).first;
      CSCDetId chamberID =cscDetId.chamberId();
      /// find the entry into the map
      map<CSCDetId, CSCEventData>::iterator chamberMapItr = chamberMap.find(chamberID);
      if(chamberMapItr == chamberMap.end())
	{
	  /// make an entry, telling it the correct chamberType
	  int istation = cscDetId.station();
	  int iring = cscDetId.ring();
	  int chamberType = 2 * istation + iring; // i=2S+R
	  if ( istation == 1 ) 
	    {
	      --chamberType;             // ring 1R -> i=1+R (2S+R-1=1+R for S=1)
	      if ( chamberType > 4 ) 
		{
		  chamberType = 1;       // But ring 1A (R=4) -> i=1
		}
	    }
	  chamberMapItr = chamberMap.insert(pair<CSCDetId, CSCEventData>(chamberID, CSCEventData(chamberType))).first;
	}
      CSCEventData & cscData = chamberMapItr->second;
      cscData.dmbHeader().setCrateAddress(mapping->crate(cscDetId), mapping->dmb(cscDetId));
      ///add strip digis to that chamber
      std::vector<CSCWireDigi>::const_iterator digiItr = (*j).second.first;
      std::vector<CSCWireDigi>::const_iterator last = (*j).second.second;
      for( ; digiItr != last; ++digiItr) 
	{
	  cscData.add(*digiItr, cscDetId.layer() );
	}
    }

  return chamberMap;
}


void CSCDigiToRaw::createFedBuffers(const CSCStripDigiCollection& stripDigis,
				    const CSCWireDigiCollection& wireDigis,
				    FEDRawDataCollection& fed_buffers,
				    const CSCChamberMap* mapping,
				    Event & e)
{

  ///bits of code from ORCA/Muon/METBFormatter - thanks, Rick:)!
  
  ///get fed object from fed_buffers
  /// make a map from the index of a chamber to the event data from it
  map<CSCDetId, CSCEventData> chamberDataMap 
    = fillChamberDataMap(stripDigis, wireDigis, mapping);
  

  int l1a=1; ///need to add increments or get it from lct digis 
  int bx = 0;///same as above
  int startingFED = FEDNumbering::getCSCFEDIds().first;

  for (int idcc=FEDNumbering::getCSCFEDIds().first-startingFED;
       idcc<=FEDNumbering::getCSCFEDIds().second-startingFED;++idcc) 
    {
      ///idcc goes from 0 to 7
      /// @@ if ReadoutMapping changes, this'll have to change
      /// DCCs 1, 2,4,5have 5 DDUs.  Otherwise, 4
      ///int nDDUs = (idcc < 2) || (idcc ==4) || (idcc ==5)
      ///          ? 5 : 4; 
      ///@@ WARNING some DCCs only have 4 DDUs, but I'm giving them all 5, for now
      int nDDUs = 5;
      int indexDDU=0;
      int oldDDUNumber = mapping->ddu(chamberDataMap.begin()->first);///initialize to first DDU

      CSCDCCEventData dccEvent(idcc, nDDUs, bx, l1a);
      /// for every chamber with data, add to a DDU in this DCC Event
      for(map<CSCDetId, CSCEventData>::iterator chamberItr = chamberDataMap.begin();
	chamberItr != chamberDataMap.end(); ++chamberItr)
	{
	  int indexDCC = mapping->slink(chamberItr->first);
	  if (idcc==indexDCC) 
	    { ///fill the right dcc 
	      dccEvent.dduData()[indexDDU].add(chamberItr->second);
	      boost::dynamic_bitset<> dccBits = dccEvent.pack();
	      FEDRawData & fedRawData = fed_buffers.FEDData(startingFED+idcc);
	      fedRawData.resize(dccBits.size());
	      ///fill data with dccEvent
	      bitset_utilities::bitsetToChar(dccBits, fedRawData.data());
	      FEDHeader cscFEDHeader(fedRawData.data());
	      cscFEDHeader.set(fedRawData.data(), 0, e.id().event(), 0, startingFED+idcc);
	      FEDTrailer cscFEDTrailer(fedRawData.data()+(fedRawData.size()-8));
	      cscFEDTrailer.set(fedRawData.data()+(fedRawData.size()-8), 
				fedRawData.size()/8, 
				evf::compute_crc(fedRawData.data(),fedRawData.size()), 0, 0);
	      int dduId = mapping->ddu(chamberItr->first); ///get ddu id based on ChamberId form mapping
	      if (oldDDUNumber!=dduId) 
		{ //if new ddu increment indexDDU counter
		  ++indexDDU;
		  oldDDUNumber = dduId;
		}
	    }
	}
    }
}



