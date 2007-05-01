#include "EventFilter/CSCRawToDigi/interface/CSCAnodeData.h"
#include "EventFilter/CSCRawToDigi/interface/CSCALCTHeader.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include <iostream>
#include <string.h> // for bzero

bool CSCAnodeData::debug = false;

CSCAnodeData::CSCAnodeData() 
: nAFEBs_(0), nTimeBins_(0) 
{
  bzero(this, sizeInWords()*2);
}


CSCAnodeData::CSCAnodeData(const CSCALCTHeader & header)
: nAFEBs_(header.nLCTChipRead()), nTimeBins_(header.NTBins()) 
{
  
  if(debug) edm::LogInfo ("CSCAnodeData") << "Making Anode data " 
					  << nFrames() << " AFEB " << nAFEBs_ 
					  << " TBINS " << nTimeBins_;
  bzero(this, sizeInWords()*2);
  for(int afeb = 0; afeb < nAFEBs_; ++afeb) 
    {
      for(int tbin = 0; tbin < nTimeBins_; ++tbin) 
	{
	  for(int layer = 1; layer <= 6; ++layer) 
	    {
	      for(int halfLayer = 0; halfLayer < 2; ++halfLayer) 
		{
		  rawHit(afeb, tbin, layer, halfLayer)
		    = CSCAnodeDataFrame(afeb, tbin, 0);
		}
	    }
	}
    }
}


// initialize
CSCAnodeData::CSCAnodeData(const CSCALCTHeader & header ,
                               const unsigned short *buf) 
: nAFEBs_(header.nLCTChipRead()), nTimeBins_(header.NTBins()) 
{
  if (debug ) 
    {
      edm::LogInfo ("CSCAnodeData") << "nAFEBs = " << nAFEBs_ << "  nTimeBins = " 
				    << nTimeBins_ << " nFrames = " << nFrames();  
      edm::LogInfo ("CSCAnodeData") << header << " HEADER CHECK " << header.check();
    }  
  memcpy(theDataFrames, buf, sizeInWords()*2);

  if (debug) 
    {
      for (int i = 0; i < 4; ++i)
	{
	  edm::LogInfo ("CSCAnodeData") << std::ios::hex << buf[i+nFrames()] << " ";
	}
    }
}

std::vector<CSCWireDigi> CSCAnodeData::wireDigis(int layer) const 
{
  std::vector<CSCWireDigi> digis;
  uint32_t tbinbits=0;
  uint16_t wireGroup=0;
  for(int afeb = 0; afeb < nAFEBs_; ++afeb) 
    {
      for(int halfLayer = 0; halfLayer <2; ++halfLayer) 
	{
	  for (int j=0;j<8;++j) 
	    {
	      for(int tbin = 0; tbin < nTimeBins_; ++tbin) 
		{
		  const CSCAnodeDataFrame & frame = rawHit(afeb,tbin,layer, halfLayer);
		  // see if there's anything in 1st 8 bits.  Usually zero
		  if(frame.data() != 0) 
		    {
		      if(frame.isHit(j)) 
			{
			  // wireGroup; 5th e time; FD time; beamCrossingTag; adcCounts
			  wireGroup = (afeb*16+halfLayer*8+j)+1;
			  tbinbits=tbinbits + (1<<tbin);      
			}
		    }
		}//end of tbin loop
	      if (tbinbits !=0 ) 
		{
		  CSCWireDigi digi(wireGroup, tbinbits);
		  if (debug)
		    edm::LogInfo ("CSCAnodeData") << "Layer " << layer << " " << digi;
		  digis.push_back(digi);
		  tbinbits=0;
		}
	    }
	}
    }
  return digis;
}

std::vector < std::vector<CSCWireDigi> > CSCAnodeData::wireDigis() const 
{
  std::vector < std::vector<CSCWireDigi> > result;
  for (int layer = 1; layer <= 6; ++layer) 
    {
      result.push_back(wireDigis(layer));
    }
  return result;
}

void CSCAnodeData::add(const CSCWireDigi & digi, int layer) 
{

  int wireGroup = digi.getWireGroup();
  int bxn=digi.getBeamCrossingTag(); 
  int alctBoard  = (wireGroup-1) / 16;
  int localGroup = (wireGroup-1) % 16;

  // crash if there's a bad strip number, but don't freak out
  // if a time bin is out of range 
  //  assert(alctBoard < nAFEBs_);
  if(alctBoard > nAFEBs_)
    {
      edm::LogError("CSCAnodeData") << "Bad Wire Number for this digi.";
      return;
    }

  if(bxn >= 0 && bxn < nTimeBins_) 
    {
      // 12 16-bit words per time bin, two per layer
      // wiregroups 0-7 go on the first line, 8-15 go on the 2nd.
      unsigned halfLayer = (localGroup > 7);
      unsigned bitNumber = localGroup % 8;
      // and pack it in the 8 bits allocated
      rawHit(alctBoard, bxn, layer, halfLayer).addHit(bitNumber);
    } 
  else 
    {
      edm::LogWarning("CSCAnodeData")<< "warning: not saving anode data in bx " << bxn 
				     << ": out of range ";
    }
}

