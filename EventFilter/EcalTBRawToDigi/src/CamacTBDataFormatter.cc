/*  
 *
 *  \author G. Franzoni
 *
 */

#include "CamacTBDataFormatter.h"

using namespace edm;
using namespace std;

#include <iostream>


// pro-memo:
// "ff" = 1 Byte
// 64 bits = 8 Bytes = 16 hex carachters
// for now: event is  ( 114 words x 32 bits ) = 448 Bytes
// size of ulong = 4 Bytes. Thus, 1 (32 bit) word = 1 (ulong) 

struct hodo_fibre_index 
{
  int nfiber;
  int ndet;
};

// nHodoscopes = 2; nFibres = 64 
const static struct hodo_fibre_index hodoFiberMap[2][64] = {
  { // Hodo 0
    // unit 1A
    {23,44}, {29,47}, {31,48}, {21,43},
    { 5,35}, {15,40}, { 7,36}, {13,39},
    { 1,33}, {11,38}, { 3,34}, { 9,37},
    { 6, 3}, {16, 8}, { 8, 4}, {14, 7},
    // unit 1C
    {17,41}, {19,42}, {27,46}, {25,45},
    {32,16}, {22,11}, {24,12}, {30,15},
    {12, 6}, { 2, 1}, { 4, 2}, {10, 5},
    {28,14}, {18, 9}, {20,10}, {26,13},
    // unit 2A
    {54,27}, {56,28}, {64,32}, {62,31},
    {49,57}, {59,62}, {51,58}, {57,61},
    {53,59}, {63,64}, {55,60}, {61,63},
    {45,55}, {39,52}, {37,51}, {47,56},
    // unit 2C
    {34,17}, {42,21}, {44,22}, {36,18},
    {50,25}, {60,30}, {58,29}, {52,26},
    {38,19}, {40,20}, {48,24}, {46,23},
    {41,53}, {35,50}, {33,49}, {43,54}
  },
  { // Hodo 1
    // unit 1A
    {31,48}, {29,47}, {23,44}, {21,43},
    { 5,35}, { 7,36}, {15,40}, {13,39},
    { 1,33}, { 3,34}, {11,38}, { 9,37},
    { 6, 3}, { 8, 4}, {16, 8}, {14, 7},
    // unit 1C
    {17,41}, {27,46}, {19,42}, {25,45},
    {24,12}, {22,11}, {32,16}, {30,15},
    { 4, 2}, { 2, 1}, {12, 6}, {10, 5},
    {20,10}, {18, 9}, {28,14}, {26,13},
    // unit 2A
    {54,27}, {64,32}, {56,28}, {62,31},
    {49,57}, {51,58}, {59,62}, {57,61},
    {53,59}, {55,60}, {63,64}, {61,63},
    {45,55}, {47,56}, {37,51}, {39,52},
    // unit 2C
    {34,17}, {42,21}, {36,18}, {44,22},
    {50,25}, {52,26}, {58,29}, {60,30},
    {38,19}, {48,24}, {40,20}, {46,23},
    {41,53}, {43,54}, {33,49}, {35,50}
  }
};



CamacTBDataFormatter::CamacTBDataFormatter () {
  nWordsPerEvent = 148;
}



void CamacTBDataFormatter::interpretRawData( const FEDRawData & fedData, 
					     EcalTBEventHeader& tbEventHeader,
					     EcalTBHodoscopeRawInfo& hodoRaw,
					     EcalTBTDCRawInfo& tdcRawInfo )
{
  

  const ulong * buffer = ( reinterpret_cast<ulong*>(const_cast<unsigned char*> ( fedData.data())));
  int fedLenght                        = fedData.size(); // in Bytes
  
  // check ultimate fed size and strip off fed-header and -trailer
  if (fedLenght != (nWordsPerEvent *4) )
    {
      LogError("CamacTBDataFormatter") << "CamacTBData has size "  <<  fedLenght
				       <<" Bytes as opposed to expected " 
				       << (nWordsPerEvent *4)
				       << ". Returning."<< endl;
      return;
    }

  
  
  ulong a=1; // used to extract an 8 Bytes word from fed 
  ulong b=1; // used to manipulate the 8 Bytes word and get what needed

  // initializing array of statuses
  for (int wordNumber=0; wordNumber<nWordsPerEvent; wordNumber++)
    { statusWords[wordNumber -1] = true;}

  //  for (int wordNumber=0; wordNumber<nWordsPerEvent; wordNumber++)
  //    { checkStatus( buffer[wordNumber],  wordNumber);}
  
  //   for (int wordNumber=0; wordNumber<nWordsPerEvent; wordNumber++)
  //     {
  //       if (! statusWords[wordNumber])
  // 	{
  // 	  LogError("CamacTBDataFormatter") << "bad status in some of the event words; returning;" << endl;	  
  // 	}
  //     }
  
  

  int wordCounter =0;
  wordCounter +=4;


  // read first word
  a = buffer[wordCounter];wordCounter++;
  LogDebug("CamacTBDataFormatter") << "\n\nword:\t" << a << endl;
  
  b = (a& 0xff000000);
  b = b >> 24;
  LogDebug("CamacTBDataFormatter") << "format  ver:\t" << b << endl;

  b = (a& 0xff0000);
  b = b >> 16;
  LogDebug("CamacTBDataFormatter") << "major:\t" << b << endl;

  b = (a& 0xff00);
  b = b >> 8;
  LogDebug("CamacTBDataFormatter") << "minor:\t" << b << endl;

  a = buffer[wordCounter];wordCounter++;
  LogDebug("CamacTBDataFormatter") << "\n\n word:\t" << a << endl;
  LogDebug("CamacTBDataFormatter") << "time stamp secs: "<<a << endl;

  a = buffer[wordCounter];wordCounter++;
  LogDebug("CamacTBDataFormatter") << "\n\n word:\t" << a << endl;
  LogDebug("CamacTBDataFormatter") << "time stamp musecs: " <<a << endl;


  a = buffer[wordCounter];wordCounter++;
  LogDebug("CamacTBDataFormatter") << "\n\n word:\t" << a << endl;
  b = (a& 0xffffff);
  LogDebug("CamacTBDataFormatter") << "LV1A: "<< b << endl;
  int lv1 = b;

  a = buffer[wordCounter];wordCounter++;
  LogDebug("CamacTBDataFormatter") << "\n\n word:\t" << a << endl;
  b = (a& 0xffff0000);
  b = b >> 16;
  LogDebug("CamacTBDataFormatter") << "run number: "<< b << endl;
  int run = b;
  b = (a& 0xffff);
  LogDebug("CamacTBDataFormatter") << "spill number: "<< b << endl;
  int spill = b;

  a = buffer[wordCounter];wordCounter++;
  b = (a& 0xffff);
  LogDebug("CamacTBDataFormatter") << "event number in spill: "<< b << endl;

  a = buffer[wordCounter];wordCounter++;
  b = (a& 0xffffff);
  LogDebug("CamacTBDataFormatter") << "internal event number: "<< b << endl;

  a = buffer[wordCounter];wordCounter++;
  LogDebug("CamacTBDataFormatter") << "\n\n word:\t" << a << endl;
  b = (a& 0xffff0000);
  b = b >> 16;
  LogDebug("CamacTBDataFormatter") << "vme errors: "<< b << endl;
  b = (a& 0xffff);
  LogDebug("CamacTBDataFormatter") << "camac errors: "<< b << endl;
  // if any of these numbers >~50% raise alarm


  // skip 2 reserved words
  wordCounter +=2;

  /**********************************
  // acessing the hodoscope block
  **********************************/

  // getting 16 words buffer and checking words statuses
  ulong bufferHodo[16]; 
  bool hodoAreGood = true;
  for (int hodo=0; hodo<16; hodo++)
    {
      hodoAreGood = hodoAreGood && checkStatus(buffer[wordCounter], wordCounter);

      a                 = buffer[wordCounter];
      bufferHodo[hodo]  = buffer[wordCounter];
      wordCounter++;
            
      b = (a& 0xffffff);
      LogDebug("CamacTBDataFormatter") << "hodo: " << hodo << "\t: " << b << endl;
    }

  // unpacking the hodo data
  if (hodoAreGood){
  for (int iplane=0; iplane<nHodoPlanes; iplane++) 
    {         
      int detType = 1;       // new mapping for electronics channels  
               
      for (int fiber=0; fiber<nHodoFibers; fiber++) { hodoHits[iplane][fiber] = 0; }            
               
      int ch=0;
      
      // loop on [4-24bits words] = 1 plane 
      for(int j=0; j<hodoRawLen; j++) 
	{
	  int word=  bufferHodo[  j+iplane*hodoRawLen  ]  &0xffff;
	  for(int i=1; i<0x10000; i<<=1) 
	    {
	      if ( word & i ) 
		{
		  // map electronics channel to No of fibre
		  hodoHits[iplane][ hodoFiberMap[detType][ch].nfiber - 1]++;
		}
	      ch ++;
	    }
	} 
    }

  
  // building the hodo infos (returning decoded hodoscope hits information)
  hodoRaw.setPlanes((unsigned int)nHodoPlanes);
  for (int ipl = 0; ipl < nHodoPlanes; ipl++) 
    {             
      EcalTBHodoscopePlaneRawHits theHodoPlane;
      theHodoPlane.setChannels((unsigned int)nHodoFibers);
      for (int fib = 0; fib < nHodoFibers; fib++){ theHodoPlane.setHit((unsigned int)fib, (bool)hodoHits[ipl][fib]); }
      hodoRaw.setPlane((unsigned int)ipl, theHodoPlane);
    }
  }
  else
    {
      LogWarning("CamacTBDataFormatter") << "hodoscope block has hardware problems or is partly unused at LV1: "
					 << lv1 << " spill: " << spill 
					 << "run: " << run 
					 << ". Skipping digi." << endl;
    }
  
  



  /**********************************
  // acessing the scalers block
  **********************************/

  // getting 72 words buffer and checking words statuses

  bool scalersAreGood = true;
  for (int scaler=0; scaler<72; scaler++)
    {
      scalersAreGood = scalersAreGood && checkStatus(buffer[wordCounter], wordCounter);

      a = buffer[wordCounter];      wordCounter++;
      b = (a& 0xffffff);
      LogDebug("CamacTBDataFormatter") << "scaler: " << scaler << "\t: " << b << endl;
    }
  if (scalersAreGood){
    ;  }
  else
    {
      LogWarning("CamacTBDataFormatter") << "scalers block has hardware problems  or is partly unused at LV1: "
					 << lv1 << " spill: " << spill 
					 << "run: " << run << endl;
    }
  




  /**********************************
  // acessing the fingers block
  **********************************/

  LogDebug("CamacTBDataFormatter") <<"\n";
  bool fingersAreGood = true;
  for (int finger=0; finger<2; finger++)
    {
      fingersAreGood = fingersAreGood && checkStatus(buffer[wordCounter], wordCounter);

      a = buffer[wordCounter];      wordCounter++;
      b = (a& 0xffffff);
      LogDebug("CamacTBDataFormatter") << "finger: " << finger << "\t: " << b << endl;
    }
  if (fingersAreGood){
    ;  }
  else
    {
      LogWarning("CamacTBDataFormatter") << "fingers block has hardware problems  or is partly unused at LV1: "
					 << lv1 << " spill: " << spill 
					 << "run: " << run << endl;
    }
  



  /**********************************
  // acessing the multi stop TDC block
  **********************************/

  a = buffer[wordCounter];      wordCounter++;
  LogDebug("CamacTBDataFormatter") << "\n\n word:\t" << a << endl;
  b = (a& 0x000000ff);
  LogDebug("CamacTBDataFormatter") << "number of words used in multi stop TDC words: "<< b << endl;
  
  int numberTDCwords = b;
  numberTDCwords = 16;
  bool multiStopTDCIsGood = true;
  for (int tdc=0; tdc< numberTDCwords ; tdc++)
    {
      multiStopTDCIsGood =  multiStopTDCIsGood && checkStatus(buffer[wordCounter], wordCounter);

      a = buffer[wordCounter];      wordCounter++;
      b =a;
      LogDebug("CamacTBDataFormatter") << "tdc: " << tdc << "\t: " << b << endl;
    }
  if ( multiStopTDCIsGood ){
    ;  }
  else
    {
      LogWarning("CamacTBDataFormatter") << "multi stop TDC block has hardware problems or is partly unused at LV1: "
					 << lv1 << " spill: " << spill 
					 << "run: " << run << endl;
    }
  
  // skip the unused words in multi stop TDC block
  wordCounter += (16 - numberTDCwords);

  

  
  /**********************************
  // acessing table in position bit
  **********************************/
  a = buffer[wordCounter];      wordCounter++;
  b = (a & 0x00000001);  //1= table is moving; 0=table is still
  bool tableIsMoving;
  if ( b ){
    LogWarning("CamacTBDataFormatter") << " table is not in position."  << endl;
    tableIsMoving = false;
  }
  else
    {
    LogWarning("CamacTBDataFormatter") << " table is in position."  << endl;
    tableIsMoving = true;
    }
  tbEventHeader.setTableIsMoving( tableIsMoving );


  wordCounter += 3;

  
  
  /**********************************
   // acessing ADC block
   **********************************/
  // skip 10 reserved words
  wordCounter += 10;
  bool ADCIsGood = true;
  ADCIsGood =  ADCIsGood && checkStatus(buffer[wordCounter], wordCounter);
  a = buffer[wordCounter];      wordCounter++;
  b = (a&0x00ffffff);
  LogDebug("CamacTBDataFormatter") << "ADC word1: " << a << "\t ADC2: " << b << " word is: " << (wordCounter-1) << endl;
  ADCIsGood =  ADCIsGood && checkStatus(buffer[wordCounter], wordCounter);
  a = buffer[wordCounter];      wordCounter++;
  b = (a&0xffffff);
  LogDebug("CamacTBDataFormatter") << "ADC word2: " << a << "\t ADC2: " << b << endl;


  
  /**********************************
   // acessing TDC block
   **********************************/
  // skip 6 reserved words
  wordCounter += 6;
  bool TDCIsGood = true;
  TDCIsGood =  ADCIsGood && checkStatus(buffer[wordCounter], wordCounter);
  a = buffer[wordCounter];      wordCounter++;
  b = (a & 0xfffff);
  LogDebug("CamacTBDataFormatter") << "TDC word1: " << a << "\t TDC2: " << b << endl;
  TDCIsGood =  ADCIsGood && checkStatus(buffer[wordCounter], wordCounter);
  a = buffer[wordCounter];      wordCounter++;
  b = (a & 0xfffff);
  LogDebug("CamacTBDataFormatter") << "TDC word2: (ext_val_trig - LHC_clock) " 
				   << a << "\t (ext_val_trig - LHC_clock): "
				   << b << endl;
  
  tdcRawInfo.setSize(1);
  int sampleNumber =1;
  EcalTBTDCSample theTdc(sampleNumber, b);
  tdcRawInfo.setSample(0, theTdc);


  a = buffer[wordCounter];      wordCounter++;
  LogDebug("CamacTBDataFormatter") << "\n\n word:\t" << a << endl;
  b = a;
  LogDebug("CamacTBDataFormatter") << "last word of event: "<< b << endl;


}








// given a data word with 8 msb as status, checks status

bool CamacTBDataFormatter::checkStatus(ulong word, int wordNumber){
  

  if ( wordNumber > nWordsPerEvent)
    { 
      LogWarning("CamacTBDataFormatter::checkStatus") << "checking word number: "
						    <<  wordNumber << " which is out of allowed range (" 
						    << nWordsPerEvent << ")" << endl;
    }

  bool isOk = true;

  if  (word & 0x80000000) // daq item not used
    { 
      LogWarning("CamacTBDataFormatter::checkStatus") << "daq item not used at word: "<<  wordNumber << endl;
      statusWords[wordNumber -1] = false;      
      isOk = false;
    }
  
  if (word & 0x40000000) // vme error on data
    { 
      LogWarning("CamacTBDataFormatter::checkStatus") << "vme error on word: "<<  wordNumber << endl;
      statusWords[wordNumber -1] = false;      
      isOk = false;
    }
    
  if (word & 0x20000000) // vme error on status
    { 
      LogWarning("CamacTBDataFormatter::checkStatus") << "vme status error at word: "<<  wordNumber << endl;
      statusWords[wordNumber -1] = false;      
      isOk = false;
    }
    
  if (word & 0x10000000) // camac error (no X)
    { 
      LogWarning("CamacTBDataFormatter::checkStatus") << "camac error (no X) at word: "<<  wordNumber << endl;
      statusWords[wordNumber -1] = false;      
      isOk = false;
    }
    
  if (word & 0x08000000) // camac error (no Q)
    { 
      LogWarning("CamacTBDataFormatter::checkStatus") << "camac error (no Q) at word: "<<  wordNumber << endl;
      statusWords[wordNumber -1] = false;      
      isOk = false;
    }
 
  if (word & 0x04000000) // no camac check error
    { 
      LogWarning("CamacTBDataFormatter::checkStatus") << "no camac check error at word: "<<  wordNumber << endl;
      statusWords[wordNumber -1] = false;      
      isOk = false;
    }

  return isOk;

}
