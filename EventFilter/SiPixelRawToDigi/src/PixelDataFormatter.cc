#include "EventFilter/SiPixelRawToDigi/interface/PixelDataFormatter.h"

#include "CondFormats/SiPixelObjects/interface/PixelFEDCabling.h"
#include "CondFormats/SiPixelObjects/interface/PixelFEDLink.h"
#include "CondFormats/SiPixelObjects/interface/PixelROC.h"
#include "DataFormats/SiPixelDigi/interface/PixelDigi.h"
#include "DataFormats/FEDRawData/interface/FEDRawData.h"

#include "FWCore/Utilities/interface/Exception.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include <bitset>
#include <sstream>

using namespace std;

const int PixelDataFormatter::LINK_bits = 6;
const int PixelDataFormatter::ROC_bits  = 5;
const int PixelDataFormatter::DCOL_bits = 5;
const int PixelDataFormatter::PXID_bits = 8;
const int PixelDataFormatter::ADC_bits  = 8;

const int PixelDataFormatter::ADC_shift  = 0;
const int PixelDataFormatter::PXID_shift = ADC_shift + ADC_bits;
const int PixelDataFormatter::DCOL_shift = PXID_shift + PXID_bits;
const int PixelDataFormatter::ROC_shift  = DCOL_shift + DCOL_bits;
const int PixelDataFormatter::LINK_shift = ROC_shift + ROC_bits;


PixelDataFormatter::PixelDataFormatter()
{
  int s32 = sizeof(Word32);
  int s64 = sizeof(Word64);
  int s8  = sizeof(char);
  if ( s8 != 1 || s32 != 4*s8 || s64 != 2*s32) {
     edm::LogError("**PixelDataFormatter**")
          <<" unexpected sizes: "
          <<"  size of char is: " << s8
          <<", size of Word32 is: " << s32
          <<", size of Word64 is: " << s64
          <<", send exception" ;
  }
}

void PixelDataFormatter::interpretRawData(const PixelFEDCabling& fed, const FEDRawData& rawData, Digis& digis)
{

  try {
    int nWords = rawData.size()/8;
    LogDebug(" PixelDataFormatter ") <<"input size: "<<nWords<<" (8-bytes words)";
    if (nWords !=0) {
      const Word64 * word = reinterpret_cast<const Word64* >(rawData.data());
      for (int i=0; i<nWords; i++) {
        LogDebug("PixelDataFormatter") << "word64: " 
            << *reinterpret_cast<const bitset<64>*>(word);
        static const Word64 WORD32_mask  = 0xffffffff;
        Word32 w1 =  *word >> 32 & WORD32_mask;
        Word32 w2 =  *word       & WORD32_mask;
        word2digi(fed, w1, digis);
        word2digi(fed, w2, digis );
        word++;
      }
    }
  }
  catch ( cms::Exception & err) {
    edm::LogError("PixelDataFormatter, exception") <<err.what();
  }
}

FEDRawData * PixelDataFormatter::formatData( 
    PixelFEDCabling & fed, 
    const Digis & digis)
{

  vector<Word32> words;
  for (int idxLink = 0; idxLink < fed.numberOfLinks(); idxLink++) {
    PixelFEDLink * link = fed.link(idxLink);
    int numberOfRocs = link->numberOfROCs();
    for(int idxRoc = 0; idxRoc < numberOfRocs; idxRoc++) {
      PixelROC * roc = link->roc(idxRoc);
      Digis::const_iterator im= digis.find(roc->rawId());
      if (im == digis.end() ) continue;
      Range range(im->second.begin(), im->second.end());
      roc2words(*roc, range, words);
    }
  }

  //
  // debug only
  //
//{
//  if (numDigi != words.size() ) {
//    cout << " ** HERE PixelDataFormatter** PROBLEM !!!!"
//               <<" numDigi: "<< numDigi
//               <<" words.size(): " << words.size()
//               <<", send exception" << endl;
//    throw Genexception("** PixelDataFormatter, numDigi != words.size()");
//  }
//}
 
  //
  // since digis are writted in the form og 64-bit packets
  // add extra 32-bit word to make number of digis even
  //
  if (words.size() %2 != 0) words.push_back( Word32(0) );


  //
  // size in Bytes; create output structure
  //
  int dataSize = words.size() * sizeof(Word32);
  if (dataSize == 0) return new FEDRawData(0);
  FEDRawData * rawData = new FEDRawData(dataSize);
/*
  edm::LogDebug("kuku")<< " ** PixelDataFormatter dataSize is: " << dataSize 
                  <<"("<<words.size()<<"*"<<sizeof(Word32)<<")"
                  <<" size of raw data: "<< sizeof(*rawData) ;
*/

  //
  // write data
  //
  Word64 * word = reinterpret_cast<Word64* >(rawData->data());
  for (unsigned int i=0; i < words.size(); i+=2) {
    *word = (Word64(words[i]) << 32 ) | words[i+1];
//                <<*reinterpret_cast<bitset<32>*>(&(words[i]))
//                <<*reinterpret_cast<bitset<32>*>(&(words[i+1]))<<endl;
      LogDebug("PixelDataFormatter")  <<"word64:  "
       << *reinterpret_cast<bitset<64>*> (word);
    word++;
  }

  //
  // check memory
  //
  if (word != reinterpret_cast<Word64* >(rawData->data()+dataSize)) {
    string s = "** PROBLEM in PixelDataFormatter !!!";
    throw cms::Exception(s);
  }

  return rawData;
}
 

void PixelDataFormatter::roc2words(
    PixelROC &roc, 
    const Range & range,
    vector<Word32> &words) const
{
  for (DetDigis::const_iterator it = range.first; it != range.second; it++) {
    const PixelDigi & pd = (*it);
    PixelROC::GlobalPixel glo = { pd.row(), pd.column() };
    PixelROC::LocalPixel  loc = roc.toLocal(glo);
    if (! roc.inside(loc) ) continue;
    LogDebug("PixelDataFormatter")<< "DIGI: row: " << pd.row()
               <<", col: " << pd.column()
               <<", adc: " << pd.adc() ;
    Word32 word;
    word =   (roc.link()->id() << LINK_shift)
           | (roc.idInLink() << ROC_shift)
           | (loc.dcol << DCOL_shift)
           | (loc.pxid << PXID_shift)
           | (pd.adc() << ADC_shift);
    words.push_back(word);
  }
}

void PixelDataFormatter::word2digi(const PixelFEDCabling & fed, 
    const Word32 & word, Digis & digis) const
{
  // do not interpret false digis
  if (word == 0 ) return;

  static const Word32 LINK_mask = ~(~Word32(0) << LINK_bits);
  static const Word32 ROC_mask  = ~(~Word32(0) << ROC_bits);
  static const Word32 DCOL_mask = ~(~Word32(0) << DCOL_bits);
  static const Word32 PXID_mask = ~(~Word32(0) << PXID_bits);
  static const Word32 ADC_mask  = ~(~Word32(0) << ADC_bits);

  PixelROC::LocalPixel loc;
  loc.dcol = (word >> DCOL_shift) & DCOL_mask;
  loc.pxid = (word >> PXID_shift) & PXID_mask;
  int ilink = (word >> LINK_shift) & LINK_mask;
  int iroc   = (word >> ROC_shift) & ROC_mask;
  int adc   = (word >> ADC_shift) & ADC_mask;


  PixelFEDLink * link = fed.link(ilink);
  if (!link) {
    stringstream stm;
    stm << "FED shows no link of id= " << ilink;
    throw cms::Exception(stm.str());
  }
  PixelROC * roc = link->roc(iroc);
  if (!roc) {
    stringstream stm;
    stm << "Link=" << ilink << " shows no ROC with id=" << iroc;
    throw cms::Exception(stm.str());
  }

  PixelROC::GlobalPixel glo = roc->toGlobal(loc);
  PixelDigi pd( glo.row, glo.col, adc);
//
//static bool testOut = uv.testOut;
  LogDebug("PixelDataFormatter")<< "DIGI: row: " << pd.row()
                          <<", col: " << pd.column()
                          <<", adc: " << pd.adc();
  //
  //
  uint32_t detid = roc->rawId();
  digis[detid].push_back(pd);
}
