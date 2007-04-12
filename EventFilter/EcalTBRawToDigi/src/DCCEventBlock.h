// Date   : 30/05/2005
// Author : N.Almeida (LIP)
// falta fazer o update dos block sizes

#ifndef DCCEVENTBLOCK_HH
#define DCCEVENTBLOCK_HH


#include "DCCBlockPrototype.h"

class DCCTowerBlock;
class DCCDataParser;
class DCCTrailerBlock;
class DCCTCCBlock;
class DCCSRPBlock;

class DCCEventBlock : public DCCBlockPrototype {
	
	public :
		
		DCCEventBlock(
			DCCDataParser * parser, 
			ulong * buffer, 
			ulong numbBytes, 
			ulong wordsToEnd, 
			ulong wordBufferOffset = 0 , 
			ulong wordEventOffset = 0 
		);
		
		~DCCEventBlock();
		
		void dataCheck(); 
		
		std::vector< DCCTowerBlock * > & towerBlocks();
		std::vector< DCCTCCBlock *   > & tccBlocks();
		DCCSRPBlock               * srpBlock();
		DCCTrailerBlock           * trailerBlock();
		std::vector< DCCTowerBlock * >   towerBlocksById(ulong towerId);
		std::pair<bool,std::string> compare(DCCEventBlock * );

		bool eventHasErrors();
		std::string eventErrorString();
		void displayEvent(std::ostream & os=std::cout);
	
		
	protected :
		enum dccFields{ 
			
			PHYSICTRIGGER        = 1,
			CALIBRATIONTRIGGER   = 2,
			TESTTRIGGER          = 3,
			TECHNICALTRIGGER     = 4,
			
			CH_ENABLED           = 0,
			CH_DISABLED          = 1,
			CH_TIMEOUT           = 2,
			CH_SUPPRESS          = 7,
			
			SR_NREAD              = 0,
			
		
			BOE                  = 0x5, 
			
			DCCERROR_EMPTYEVENT  = 0x1, 
			
			TOWERHEADER_SIZE     = 8, 
			TRAILER_SIZE         = 8
	
		
		};		

		std::vector< DCCTowerBlock * > towerBlocks_      ;
		std::vector< DCCTCCBlock   * > tccBlocks_        ;
		DCCTrailerBlock       *   dccTrailerBlock_  ;
		DCCSRPBlock           *   srpBlock_;
		ulong wordBufferOffset_;
		bool emptyEvent;
};


inline std::vector< DCCTowerBlock * > & DCCEventBlock::towerBlocks()  { return towerBlocks_;     }
inline std::vector< DCCTCCBlock * >   & DCCEventBlock::tccBlocks()    { return tccBlocks_;       }
inline DCCSRPBlock               * DCCEventBlock::srpBlock()     { return srpBlock_;        }
inline DCCTrailerBlock           * DCCEventBlock::trailerBlock() { return dccTrailerBlock_; }

#endif
