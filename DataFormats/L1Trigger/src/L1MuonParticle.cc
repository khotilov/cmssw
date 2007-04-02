// -*- C++ -*-
//
// Package:     L1Trigger
// Class  :     L1MuonParticle
// 
/**\class L1MuonParticle \file L1MuonParticle.cc DataFormats/L1Trigger/src/L1MuonParticle.cc \author Werner Sun
*/
//
// Original Author:  Werner Sun
//         Created:  Tue Jul 25 17:51:21 EDT 2006
// $Id: L1MuonParticle.cc,v 1.5 2006/11/13 00:18:09 wsun Exp $
//

// system include files

// user include files
#include "DataFormats/L1Trigger/interface/L1MuonParticle.h"

using namespace l1extra ;

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
L1MuonParticle::L1MuonParticle()
{
}

L1MuonParticle::L1MuonParticle(
   Charge q,
   const LorentzVector& p4,
   const L1MuGMTExtendedCand& aCand )
  : LeafCandidate( q, p4 ),
     cand_( aCand )
{
   isolated_ = cand_.isol() ;
   mip_ = cand_.mip() ;
   forward_ = cand_.isFwd() ;
   rpc_ = cand_.isRPC() ;
}

L1MuonParticle::L1MuonParticle( Charge q,
				const LorentzVector& p4,
				bool isolated,
				bool mip,
				bool forward,
				bool rpc,
				unsigned int detector )
   : LeafCandidate( q, p4 ),
     isolated_( isolated ),
     mip_( mip ),
     forward_( forward ),
     rpc_( rpc ),
     cand_( L1MuGMTExtendedCand() )
{
}

// L1MuonParticle::L1MuonParticle(const L1MuonParticle& rhs)
// {
//    // do actual copying here;
// }

L1MuonParticle::~L1MuonParticle()
{
}

//
// assignment operators
//
// const L1MuonParticle& L1MuonParticle::operator=(const L1MuonParticle& rhs)
// {
//   //An exception safe implementation is
//   L1MuonParticle temp(rhs);
//   swap(rhs);
//
//   return *this;
// }

//
// member functions
//

//
// const member functions
//

//
// static member functions
//
