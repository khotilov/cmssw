#if !defined(CSCTFTBRAWFORMAT_CSCTFTBFRONTBLOCK_H)
#define CSCTFTBRAWFORMAT_CSCTFTBFRONTBLOCK_H
// -*- C++ -*-
//
// Package:     CSCTFTBRawFormat
// Module:      CSCTFTBFrontBlock
// 
// Description: Header file for Front Event Data Block
//
// Implementation:
//     <Notes on implementation>
//
// Author:      Lindsey Gray
// Created:     13.5.2003
//
// $Id: CSCTFTBFrontBlock.h,v 1.2 2005/03/03 18:14:48 lgray Exp $
//
// Revision History
// $Log: CSCTFTBFrontBlock.h,v $
// Revision 1.2  2005/03/03 18:14:48  lgray
// Added ability to pack data back into raw form. Added test program for this as well.
//
// Revision 1.1  2005/02/14 20:59:46  lgray
// First Commit from UF
//
// Revision 1.2  2004/05/18 08:00:10  tfcvs
// DEA: touch base
//
// Revision 1.1  2004/05/17 08:25:52  tfcvs
// DEA: switch to SR BX data
//
// Revision 1.2  2003/08/27 22:07:51  tfcvs
// Added pretty-print - Rick
//
// Revision 1.1  2003/05/25 10:13:02  tfcvs
// first working version -DEA
//
// Revision 1.7  2003/05/20 22:13:06  tfcvs
// HS - Added Darin's changes
//
// Revision 1.6  2003/05/19 23:23:12  tfcvs
// HS - Commit after some changes
//
// Revision 1.4  2003/05/19 15:47:18  tfcvs
// HS - Some cleanup
//
// Revision 1.3  2003/05/19 00:25:56  tfcvs
// DEA: committed, but may not compile
//
// Revision 1.2  2003/05/15 23:58:40  tfcvs
// HS - Some cosmetics
//
//
//

// System include files
#include <vector>
#include <iostream>

// Package include files
#include "TBDataFormats/CSCTFTBRawData/interface/CSCTFTBFrontData.h"
#include "TBDataFormats/CSCTFTBRawData/interface/CSCTFTBFrontHeader.h"

// External package include files

// STL classes

// Forward declarations
//class BitVector;
class CSCTFTBEventData;
class CSCTFTBEventHeader;

class CSCTFTBFrontBlock 
{

// Friend classses and functions

// Public part
   public:
        // Constants, enums and typedefs

        // Constructors and destructor
        CSCTFTBFrontBlock();

        CSCTFTBFrontBlock(unsigned short *buf, int bx,
			  const CSCTFTBEventHeader& hdr);

	~CSCTFTBFrontBlock();

        // Member functions
  
        // Const member functions

	/// return this SR block's header
	CSCTFTBFrontHeader frontHeader() const {return frontHeader_;}

	/// return the SR data, in vector of vectors
        std::vector<std::vector<CSCTFTBFrontData> > frontData() const {return srdata_;}
	
	/// return one mpc's data
        std::vector<CSCTFTBFrontData> frontData(unsigned mpc) const;

	/// return one link's data
	CSCTFTBFrontData frontData(unsigned mpc, unsigned link) const;

	/// size of data bank in 16-bit words
        int size() const {return size_;};

	/// (relative) BX assigned to this bank
	int BX() const {return myBX_;};

	/// make a bit vector
	//BitVector pack();

        /// pretty-print
        friend std::ostream & operator<<(std::ostream & stream, const CSCTFTBFrontBlock &);
// Protected part
   protected:

        // Protected member functions
        int unpackData(unsigned short * buf, const CSCTFTBEventHeader&);


        // Protected const member functions
	
        /// SRData is unpacked and stored in this vector
	std::vector<std::vector<CSCTFTBFrontData> > srdata_;
	CSCTFTBFrontHeader frontHeader_;

	int size_;
	int myBX_;

// Private part
   private:

        // Constructors and destructor

        // Assignment operator(s)

        // Private member functions

        // Private const member functions

        // Data members

        // Static data members

        // Inline function definitions

};

#endif
