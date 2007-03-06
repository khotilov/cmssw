//-------------------------------------------------
//
//   Class: DTTFFEDReader
//
//   L1 DT Track Finder Raw-to-Digi
//
//
//   $Date: 2006/06/01 00:00:00 $
//   $Revision: 1.1 $
//
//   Author :
//   J. Troconiz  UAM Madrid
//   E. Delmeire  UAM Madrid
//
//--------------------------------------------------

#include "EventFilter/DTTFRawToDigi/interface/DTTFFEDReader.h"

#include "FWCore/Framework/interface/Handle.h"
#include "FWCore/Framework/interface/Event.h"

#include <DataFormats/L1DTTrackFinder/interface/L1MuDTChambPhContainer.h>
#include <DataFormats/L1DTTrackFinder/interface/L1MuDTChambThContainer.h>
#include <DataFormats/FEDRawData/interface/FEDRawData.h>
#include <DataFormats/FEDRawData/interface/FEDRawDataCollection.h>

#include <iostream>
#include <string>

using namespace std;


DTTFFEDReader::DTTFFEDReader(const edm::ParameterSet& pset) {

  produces<L1MuDTChambPhContainer>();
  produces<L1MuDTChambThContainer>();

}

DTTFFEDReader::~DTTFFEDReader(){}

void DTTFFEDReader::produce(edm::Event& e, const edm::EventSetup& c) {

  auto_ptr<L1MuDTChambPhContainer> phi_product(new L1MuDTChambPhContainer);
  auto_ptr<L1MuDTChambThContainer> the_product(new L1MuDTChambThContainer);

  L1MuDTChambPhContainer::Phi_Container phi_data;
  L1MuDTChambThContainer::The_Container the_data;

  if (!fillRawData(e, phi_data, the_data)) return;

  phi_product->setContainer(phi_data);
  the_product->setContainer(the_data);

  e.put(phi_product);
  e.put(the_product);

}

bool DTTFFEDReader::fillRawData(edm::Event& e,
                                L1MuDTChambPhContainer::Phi_Container& phi_data,
                                L1MuDTChambThContainer::The_Container& the_data) {

  analyse(e);

  phi_data = p_data();
  the_data = t_data();

  return true;
}

//--------------
// Operations --
//--------------
void DTTFFEDReader::analyse(edm::Event& e) {
  clear();
  process(e);
  return;
}

// process data
void DTTFFEDReader::process(edm::Event& e) {

  // Container
  vector<long long> DTTFWordContainer; 
  vector<long long>::iterator DTTFiterator;

  // Header constituents
  int BOEevTy, DTTFId;

  // DTTF Payload constituents 
  long long DTTFWord;
  int DTTFChan, bitsID;

  // Trailer constituents
  int evtLgth , CRC;

  bool goOn = true;

  //--> Header

  edm::Handle<FEDRawDataCollection> data;
  e.getByType(data);
  FEDRawData dttfdata = data->FEDData(0x030C);

  long long* dataWord = new long long;
  unsigned char* LineFED=dttfdata.data();
  *dataWord=*((long long*)LineFED);
  int lines  = 1; // already counting header

  BOEevTy = ((*dataWord)&0xFF00000000000000)>>56; // positions 57 -> 64
  DTTFId  = ((*dataWord)&0x00000000000FFF00)>>8;  // positions 9 ->20

  if( (BOEevTy != 0x50) || ( DTTFId != 0x030C) ){
    cout << "Not a DTTF header " << hex << *dataWord << endl;
    goOn = false;
  }

  int newCRC =  0xFFFF;
  calcCRC(*dataWord, newCRC);  


  //--> DTTF data 

  LineFED+=8;
  *dataWord=*((long long*)LineFED);
  int chkEOE = ((*dataWord)&0xFFF0000000000000)>>52; 
  lines++;

  while(chkEOE != 0xA00){

    calcCRC(*dataWord, newCRC);

    DTTFWord = *dataWord;
    DTTFWordContainer.push_back(DTTFWord);

    LineFED+=8;
    *dataWord=*((long long*)LineFED);
    chkEOE     = ((*dataWord)&0xFFF0000000000000)>>52; 
    lines++;

    if(lines > 3026){
      cout << "Warning : number of DTTF lines > 3026 " << endl; // 3026 = 1(header) + 3024(max # PHTF-ETTF 64 bits words) + 1(trailer)
      goOn = false;
    }

  } // end while-Data loop


  //--> Trailer

  evtLgth   = ((*dataWord)&0x00FFFFFF00000000)>>32; // positions 33 ->56
  CRC       = ((*dataWord)&0x00000000FFFF0000)>>16; // positions 17->32

  calcCRC((*dataWord)&0xFFFFFFFF0000FFFF, newCRC);
  if( newCRC != CRC){
    cout << "Calculated CRC " ;
    cout << hex << newCRC << " differs from CRC in trailer " << hex << CRC << endl;
    goOn = false;
  }

  if( lines != evtLgth){
    cout << "Number of words read != event length " << dec << lines << " " << evtLgth << endl;
    goOn = false;
  }


  // --> analyse event    

  if( !goOn ) {
    delete dataWord;
    return;
  }

  for( DTTFiterator =  DTTFWordContainer.begin();
       DTTFiterator != DTTFWordContainer.end();
       DTTFiterator++ ){

    DTTFChan = ((*DTTFiterator)&0xFF00000000000000)>>56;
    bitsID   = ((*DTTFiterator)&0x00000000F0000000)>>28;

    int bxID     = bxNr(DTTFChan);
    if(bxID     == -999) continue;
    int wheelID  = wheel(DTTFChan);    
    if(wheelID  == -999) continue;
    int sectorID = sector(DTTFChan);    
    if(sectorID == -999) continue;

    //Input
    if(wheelID!=0 && bitsID<=0x9){   

      int wheelPh   = (abs(wheelID)-1)*wheelID/abs(wheelID); 
      int stationID = 0;
      int     ra    = 0;
      int     ba    = 0;
      int tsqual    = 0;
      int ts2tag    = 0;

      if ( ( bitsID >> 1 ) == 0 ){ stationID = 1;}
      if ( ( bitsID >> 1 ) == 1 ){ stationID = 2;}
      if ( ( bitsID >> 1 ) == 4 ){ stationID = 3;}
      if ( ( bitsID >> 1 ) == 2 ){ stationID = 4;}

      if(stationID != 3){ 
	  
        ts2tag = (bitsID)&0x1;
	tsqual = (~(*DTTFiterator)&0x07)-1;
	ba     = (~(*DTTFiterator)&0x1FF8)>>3;
        if( ba>0x1FF) ba-=0x400;
	ra     = (~(*DTTFiterator)&0x1FFE000)>>13;
        if( ra>0x7FF) ra-=0x1000;
      }   
      else{ 

        ts2tag = (bitsID)&0x1;
	tsqual = (~(*DTTFiterator)&0x07)-1;
	ra     = (~(*DTTFiterator)&0x7FF8)>>3;
        if( ra>0x7FF) ra-=0x1000;
      }

      if(tsqual!=7 && wheelID!=-1){
        phiSegments.push_back(
		    L1MuDTChambPhDigi( bxID+ts2tag, wheelPh, sectorID, stationID,
		    ra, ba, tsqual, ts2tag, 0) );
      }
    }
    //Input

  } // end for-loop container content

  delete dataWord;
  return;
}

// access data
const L1MuDTChambPhContainer::Phi_Container& DTTFFEDReader::p_data() {
  return phiSegments;
}

const L1MuDTChambThContainer::The_Container& DTTFFEDReader::t_data() {
  return theSegments;
}

void DTTFFEDReader::clear() {
  phiSegments.clear();
  theSegments.clear();
  return;
}

int DTTFFEDReader::channel( int wheel, int sector,  int bx ){

  // wheel  :  -3 -2 -1 +1 +2 +3 <=> PHTF's : N2, N1, N0, P0, P1, P2
  //                           0 <=> ETTF  
  // sector :  0 -> 11
  // bx     : -1 -> +1

  int myChannel = 255;

  if ( abs(bx) > 1)               { return myChannel; }
  if ( sector < 0 || sector > 11) { return myChannel; }
  if ( abs(wheel) > 3)            { return myChannel; }

   myChannel = sector*21 + wheel*3 + bx + 10 ; 

  return myChannel;
}

int DTTFFEDReader::bxNr( int channel ){

  if (channel < 0 || channel > 252 ){ return -999; }

  int myBx = (channel%3)-1;

  return myBx;
}

int DTTFFEDReader::sector( int channel ){

  if (channel < 0 || channel > 252 ){ return -999; }

  return channel/21;
}

int DTTFFEDReader::wheel( int channel ){

  if (channel < 0 || channel > 252 ){ return -999; }    

  int myWheel = ((channel%21)/3)-3;

  return myWheel;
}

void DTTFFEDReader::calcCRC(long long myD, int &myC){

  int myCRC[16],D[64],C[16];

  for( int i=0; i < 64; i++ ){ D[i]=(myD>>i)&0x1; }
  for( int i=0; i < 16; i++ ){ C[i]=(myC>>i)&0x1; }

  myCRC[0] = ( D[63] + D[62] + D[61] + D[60] + D[55] + D[54] +
               D[53] + D[52] + D[51] + D[50] + D[49] + D[48] +
               D[47] + D[46] + D[45] + D[43] + D[41] + D[40] +
               D[39] + D[38] + D[37] + D[36] + D[35] + D[34] +
               D[33] + D[32] + D[31] + D[30] + D[27] + D[26] +
               D[25] + D[24] + D[23] + D[22] + D[21] + D[20] +
               D[19] + D[18] + D[17] + D[16] + D[15] + D[13] +
               D[12] + D[11] + D[10] + D[9]  + D[8]  + D[7]  +
               D[6]  + D[5]  + D[4]  + D[3]  + D[2]  + D[1]  +
               D[0]  + C[0]  + C[1]  + C[2]  + C[3]  + C[4]  +
               C[5]  + C[6]  + C[7]  + C[12] + C[13] + C[14] +
               C[15] )%2;

  myCRC[1] = ( D[63] + D[62] + D[61] + D[56] + D[55] + D[54] +
               D[53] + D[52] + D[51] + D[50] + D[49] + D[48] +
               D[47] + D[46] + D[44] + D[42] + D[41] + D[40] +
               D[39] + D[38] + D[37] + D[36] + D[35] + D[34] +
               D[33] + D[32] + D[31] + D[28] + D[27] + D[26] +
	       D[25] + D[24] + D[23] + D[22] + D[21] + D[20] +
	       D[19] + D[18] + D[17] + D[16] + D[14] + D[13] +
	       D[12] + D[11] + D[10] + D[9]  + D[8]  + D[7]  +
	       D[6]  + D[5]  + D[4]  + D[3]  + D[2]  + D[1]  +
	       C[0]  + C[1]  + C[2]  + C[3]  + C[4]  + C[5]  +
	       C[6]  + C[7]  + C[8]  + C[13] + C[14] + C[15] )%2;

  myCRC[2] = ( D[61] + D[60] + D[57] + D[56] + D[46] + D[42] +
	       D[31] + D[30] + D[29] + D[28] + D[16] + D[14] +
	       D[1]  + D[0]  + C[8]  + C[9]  + C[12] + C[13] )%2;

  myCRC[3] = ( D[62] + D[61] + D[58] + D[57] + D[47] + D[43] +
	       D[32] + D[31] + D[30] + D[29] + D[17] + D[15] +
	       D[2]  + D[1]  + C[9]  + C[10] + C[13] + C[14] )%2;

  myCRC[4] = ( D[63] + D[62] + D[59] + D[58] + D[48] + D[44] +
	       D[33] + D[32] + D[31] + D[30] + D[18] + D[16] + 
	       D[3]  + D[2]  + C[0]  + C[10] + C[11] + C[14] +
	       C[15] )%2;

  myCRC[5] = ( D[63] + D[60] + D[59] + D[49] + D[45] + D[34] +
	       D[33] + D[32] + D[31] + D[19] + D[17] + D[4]  +
	       D[3]  + C[1]  + C[11] + C[12] + C[15] )%2;

  myCRC[6] = ( D[61] + D[60] + D[50] + D[46] + D[35] + D[34] +
	       D[33] + D[32] + D[20] + D[18] + D[5]  + D[4]   +
	       C[2]  + C[12] + C[13] )%2;

  myCRC[7] = ( D[62] + D[61] + D[51] + D[47] + D[36] + D[35] +
	       D[34] + D[33] + D[21] + D[19] + D[6]  + D[5]   +
	       C[3]  + C[13] + C[14] )%2;

  myCRC[8] = ( D[63] + D[62] + D[52] + D[48] + D[37] + D[36] +
	       D[35] + D[34] + D[22] + D[20] + D[7]  + D[6]   +
	       C[0]  + C[4]  + C[14] + C[15] )%2;

  myCRC[9] = ( D[63] + D[53] + D[49] + D[38] + D[37] + D[36] +
	       D[35] + D[23] + D[21] + D[8]  + D[7]  + C[1]  +
	       C[5]  + C[15] )%2;

  myCRC[10] = ( D[54] + D[50] + D[39] + D[38] + D[37] + D[36] + 
		D[24] + D[22] + D[9]  + D[8]  + C[2]  + C[6] )%2;

  myCRC[11] = ( D[55] + D[51] + D[40] + D[39] + D[38] + D[37] +
		D[25] + D[23] + D[10] + D[9]  + C[3]  + C[7] )%2;

  myCRC[12] = ( D[56] + D[52] + D[41] + D[40] + D[39] + D[38] +
		D[26] + D[24] + D[11] + D[10] + C[4] + C[8] )%2;

  myCRC[13] = ( D[57] + D[53] + D[42] + D[41] + D[40] + D[39] +
		D[27] + D[25] + D[12] + D[11] + C[5]  + C[9] )%2;

  myCRC[14] = ( D[58] + D[54] + D[43] + D[42] + D[41] + D[40] +
		D[28] + D[26] + D[13] + D[12] + C[6]  + C[10] )%2;

  myCRC[15] = ( D[63] + D[62] + D[61] + D[60] + D[59] + D[54] +
		D[53] + D[52] + D[51] + D[50] + D[49] + D[48] + 
		D[47] + D[46] + D[45] + D[44] + D[42] + D[40] +
		D[39] + D[38] + D[37] + D[36] + D[35] + D[34] + 
		D[33] + D[32] + D[31] + D[30] + D[29] + D[26] +
		D[25] + D[24] + D[23] + D[22] + D[21] + D[20] + 
		D[19] + D[18] + D[17] + D[16] + D[15] + D[14] +
		D[12] + D[11] + D[10] + D[9]  + D[8]  + D[7]  + 
		D[6]  + D[5]  + D[4]  + D[3]  + D[2]  + D[1]  +
		D[0]  + C[0]  + C[1]  + C[2]  + C[3]  + C[4]  + 
		C[5]  + C[6]  + C[11] + C[12] + C[13] + C[14] +
		C[15] )%2;

  int tempC = 0x0;  
  for(int i=0; i<16 ; i++){ tempC = tempC + (myCRC[i]<<i); }
  myC = tempC ;
  return;
}
