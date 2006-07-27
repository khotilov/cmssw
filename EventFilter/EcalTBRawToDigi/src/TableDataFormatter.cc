#include "TableDataFormatter.h"

using namespace edm;
using namespace std;

#include <iostream>

TableDataFormatter::TableDataFormatter () {
}

void TableDataFormatter::interpretRawData( const FEDRawData & fedData, 
					   EcalTBEventHeader& tbEventHeader)
{
  const ulong * buffer = ( reinterpret_cast<ulong*>(const_cast<unsigned char*> ( fedData.data())));
  int fedLenght                        = fedData.size(); // in Bytes
  
  // check ultimate fed size and strip off fed-header and -trailer
  if (fedLenght != (nWordsPerEvent *4) )
    {
      LogError("TableDataFormatter") << "TableData has size "  <<  fedLenght
				       <<" Bytes as opposed to expected " 
				       << (nWordsPerEvent *4)
				       << ". Returning."<< endl;
      return;
    }

  ulong a=1; // used to extract an 8 Bytes word from fed 
  ulong b=1; // used to manipulate the 8 Bytes word and get what needed

  int wordCounter =0;
  wordCounter +=4;

  a = buffer[wordCounter];wordCounter++;
  b = (a& 0xffffffff);
  tbEventHeader.setThetaTableIndex(b);
  LogDebug("TableDataFormatter") << "Table theta position:\t" << b << endl;
  a = buffer[wordCounter];wordCounter++;
  b = (a& 0xffffffff);
  tbEventHeader.setPhiTableIndex(b);
  LogDebug("TableDataFormatter") << "Table phi position:\t" << b << endl;
  a = buffer[wordCounter];wordCounter++;
  b = (a& 0xffff);
  tbEventHeader.setCrystalInBeam(EBDetId(1,b,EBDetId::SMCRYSTALMODE));
  LogDebug("TableDataFormatter") << "Actual Current crystal in beam:\t" << b << endl;
  b = (a& 0xffff0000);
  b = b >> 16;
  tbEventHeader.setNominalCrystalInBeam(EBDetId(1,b,EBDetId::SMCRYSTALMODE));
  LogDebug("TableDataFormatter") << "Nominal Current crystal in beam:\t" << b << endl;
  a = buffer[wordCounter];wordCounter++;
  b = (a& 0xffff);
  tbEventHeader.setNextCrystalInBeam(EBDetId(1,b,EBDetId::SMCRYSTALMODE));
  LogDebug("TableDataFormatter") << "Next crystal in beam:\t" << b << endl;
}
