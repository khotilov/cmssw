// -*- C++ -*-
//
// Package:     CSCTFTBRawFormat
// Module:      CSCTFTBDTData
// 
// Description: Source file for DT Event Data sent to SP,
//              rewritten to use bitfields.
// Implementation:
//     <Notes on implementation>
//
// Author:      Lindsey Gray
// Created:     13.1.2005
//
// $Id: CSCTFTBDTData.cc,v 1.2 2005/03/03 18:14:48 lgray Exp $
//
// Revision History
// $Log: CSCTFTBDTData.cc,v $
// Revision 1.2  2005/03/03 18:14:48  lgray
// Added ability to pack data back into raw form. Added test program for this as well.
//
// Revision 1.1  2005/02/14 21:01:32  lgray
// First Commit from UF
//
//
//

// System include files
#include <string.h>
// Package include files
#include "TBDataFormats/CSCTFTBRawData/interface/CSCTFTBDTData.h"
// External package include files
//#include "Utilities/GenUtil/interface/BitVector.h"

CSCTFTBDTData::CSCTFTBDTData()
{}

CSCTFTBDTData::CSCTFTBDTData(unsigned quality, unsigned bend, unsigned flag,
			     unsigned calib, unsigned phi, unsigned bx, 
			     unsigned bc0):
  quality_(quality),bend_(bend),flag_(flag),calib_(calib),
  phi_(phi),bx_(bx),bc0_(bc0)
{
  zero1_ = zero2_ = zero3_ = zero4_ = 0;
}

CSCTFTBDTData::CSCTFTBDTData(const CSCTFTBDTData &parent)
{
  memcpy(this,&parent,size()*sizeof(unsigned short));
}

CSCTFTBDTData::~CSCTFTBDTData()
{}

/* put back in later if BitVector Comes back to life. -LG
BitVector CSCTFTBDTData::packVector() const
{
  return BitVector(reinterpret_cast<const unsigned*>(this),size()*sizeof(short)*8);
}
*/

std::ostream & operator<<(std::ostream & stream, const CSCTFTBDTData & dtData)
{
  stream << " DT Track Stub:\n\t  " <<"Replace me with pretty print!\n";//dtData.dtDataStub();
  return stream;
}
