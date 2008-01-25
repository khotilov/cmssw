#ifndef StarterKit_HistoJet_h
#define StarterKit_HistoJet_h

//------------------------------------------------------------
// Title: HistoJet.h
// Purpose: To histogram Jets
//
// Authors:
// Liz Sexton-Kennedy <sexton@fnal.gov>
// Eric Vaandering <ewv@fnal.gov >
// Petar Maksimovic <petar@jhu.edu>
// Sal Rappoccio <rappocc@fnal.gov>
//------------------------------------------------------------
//
// Interface:
//
//   HistoJet ( TFile * file );
//   Description: Constructor.
//
//   void fill( TK::Jet * );
//   Description: Fill object. Will fill relevant jet variables
//
//   void write();
//   Description: Write object to file in question.
//
//   ~HistoJet
//    Description: Destructor. Deallocates memory.
//
//------------------------------------------------------------
//
// Modification History:
//
//   -29Nov07: Sal Rappoccio: Creation of the object
//------------------------------------------------------------


// CMSSW include files
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "PhysicsTools/StarterKit/interface/HistoGroup.h"

// STL include files
#include <string>

// ROOT include files
#include <TH1D.h>
#include <TFile.h>

namespace pat {

  class HistoJet : public HistoGroup<Jet> {

  public:
    HistoJet( std::string dir = "jet",
	      double pt1=0, double pt2=200, double m1=0, double m2=200 );
    virtual ~HistoJet();

    void fill( const Jet * jet, uint iPart = 0 );
    void fill( const Jet & jet, uint iPart = 0 ) { fill(&jet, iPart); }


    // Clear ntuple cache
    void clearVec();
  protected:
    PhysVarHisto *    h_jetFlavour_;
    PhysVarHisto *    h_BDiscriminant_;
    PhysVarHisto *    h_jetCharge_;
    PhysVarHisto *    h_nTrk_;

  };

}
#endif
