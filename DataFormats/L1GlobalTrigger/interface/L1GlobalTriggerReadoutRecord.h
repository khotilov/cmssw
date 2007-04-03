#ifndef L1GlobalTrigger_L1GlobalTriggerReadoutRecord_h
#define L1GlobalTrigger_L1GlobalTriggerReadoutRecord_h

/**
 * \class L1GlobalTriggerReadoutRecord
 * 
 * 
 * Description: readout record for L1 Global Trigger.  
 *
 * Implementation:
 *    <TODO: enter implementation details>
 *   
 * \author: N. Neumeister        - HEPHY Vienna - ORCA version 
 * \author: Vasile Mihai Ghete   - HEPHY Vienna - CMSSW version 
 * 
 * $Date$
 * $Revision$
 *
 */

// system include files
#include <string>
#include <vector>
#include <iostream>

// user include files
#include "DataFormats/L1GlobalTrigger/interface/L1GlobalTriggerReadoutSetupFwd.h"
#include "DataFormats/L1GlobalTrigger/interface/L1GlobalTriggerReadoutSetup.h"

#include "DataFormats/L1GlobalMuonTrigger/interface/L1MuGMTReadoutCollection.h"
#include "DataFormats/L1GlobalMuonTrigger/interface/L1MuGMTCand.h"
#include "DataFormats/L1GlobalMuonTrigger/interface/L1MuGMTExtendedCand.h"

#include "DataFormats/L1GlobalCaloTrigger/interface/L1GctEmCand.h"
#include "DataFormats/L1GlobalCaloTrigger/interface/L1GctJetCand.h"
#include "DataFormats/L1GlobalCaloTrigger/interface/L1GctEtSums.h"
#include "DataFormats/L1GlobalCaloTrigger/interface/L1GctJetCounts.h"

#include "DataFormats/L1GlobalTrigger/interface/L1GtfeWord.h"
#include "DataFormats/L1GlobalTrigger/interface/L1GtFdlWord.h"
#include "DataFormats/L1GlobalTrigger/interface/L1GtPsbWord.h"

#include "DataFormats/Common/interface/RefProd.h"

// forward declarations
namespace edm {
    template <typename T> class Handle;
}

// class interface

class L1GlobalTriggerReadoutRecord
{

public:

    /// constructors
    L1GlobalTriggerReadoutRecord();

    L1GlobalTriggerReadoutRecord(int NumberBxInEvent);

    /// copy constructor
    L1GlobalTriggerReadoutRecord(const L1GlobalTriggerReadoutRecord&);

    /// destructor
    virtual ~L1GlobalTriggerReadoutRecord();
      
    /// assignment operator
    L1GlobalTriggerReadoutRecord& operator=(const L1GlobalTriggerReadoutRecord&);
  
    /// equal operator
    bool operator==(const L1GlobalTriggerReadoutRecord&) const;

    /// unequal operator
    bool operator!=(const L1GlobalTriggerReadoutRecord&) const;

public:

    /// get Global Trigger decision and the decision word
    ///   overloaded w.r.t. bxInEvent argument
    ///   bxInEvent not given: for bunch cross with L1Accept
    const bool decision(unsigned int bxInEvent) const;
    const bool decision() const;

    const DecisionWord decisionWord(unsigned int bxInEvent) const;
    const DecisionWord decisionWord() const;
  
    /// set global decision and the decision word
    void setDecision(const bool& t, unsigned int bxInEvent);
    void setDecision(const bool& t);

    void setDecisionWord(const DecisionWord& decisionWordValue, unsigned int bxInEvent);
    void setDecisionWord(const DecisionWord& decisionWordValue);

    /// print global decision and algorithm decision word
    void printGtDecision(std::ostream& myCout, unsigned int bxInEventValue) const;
    void printGtDecision(std::ostream& myCout) const;
    
    /// print technical triggers
    void printTechnicalTrigger(std::ostream& myCout, unsigned int bxInEventValue) const;
    void printTechnicalTrigger(std::ostream& myCout) const;
    
// TODO after inserting the configuration
//    ///  set decision for given configuration
//    ///  disable trigger item(s)

    //**************************************************************************
    // get/set physical candidates 
    //     indexCand    index of candidate
    //     bxInEvent    bunch cross in the event (3bx records, 5bx records) 
    //**************************************************************************

    /// get candidate 
    ///     two arguments:   candidate index, bunch cross index
    ///     one argument:    candidate index, for bunch cross with L1Accept

    /// get all non-empty candidates
    ///     one argument:  bunch cross index
    ///     no argument:   for bunch cross with L1Accept

    /// muon 
    
    /// get / set reference to L1MuGMTReadoutCollection
    const edm::RefProd<L1MuGMTReadoutCollection> muCollectionRefProd() const;
    void setMuCollectionRefProd(edm::Handle<L1MuGMTReadoutCollection>&);
    
    
    const L1MuGMTExtendedCand muonCand(unsigned int indexCand, unsigned int bxInEvent) const;
    const L1MuGMTExtendedCand muonCand(unsigned int indexCand) const;

    std::vector<L1MuGMTExtendedCand> muonCands(unsigned int bxInEvent) const;
    std::vector<L1MuGMTExtendedCand> muonCands() const;

    /// electron

    const L1GctEmCand electronCand(unsigned int indexCand, unsigned int bxInEvent) const;
    const L1GctEmCand electronCand(unsigned int indexCand) const;

    std::vector<L1GctEmCand> electronCands(unsigned int bxInEvent) const;
    std::vector<L1GctEmCand> electronCands() const;

    /// isolated electron
    const L1GctEmCand isolatedElectronCand(unsigned int indexCand, unsigned int bxInEvent) const;
    const L1GctEmCand isolatedElectronCand(unsigned int indexCand) const;

    std::vector<L1GctEmCand> isolatedElectronCands(unsigned int bxInEvent) const;
    std::vector<L1GctEmCand> isolatedElectronCands() const;

    /// central jet
    const L1GctJetCand centralJetCand(unsigned int indexCand, unsigned int bxInEvent) const;
    const L1GctJetCand centralJetCand(unsigned int indexCand) const;

    std::vector<L1GctJetCand> centralJetCands(unsigned int bxInEvent) const;
    std::vector<L1GctJetCand> centralJetCands() const;

    /// forward jet
    const L1GctJetCand forwardJetCand(unsigned int indexCand, unsigned int bxInEvent) const;
    const L1GctJetCand forwardJetCand(unsigned int indexCand) const;

    std::vector<L1GctJetCand> forwardJetCands(unsigned int bxInEvent) const;
    std::vector<L1GctJetCand> forwardJetCands() const;

    /// tau jet
    const L1GctJetCand tauJetCand(unsigned int indexCand, unsigned int bxInEvent) const;
    const L1GctJetCand tauJetCand(unsigned int indexCand) const;

    std::vector<L1GctJetCand> tauJetCands(unsigned int bxInEvent) const;
    std::vector<L1GctJetCand> tauJetCands() const;

    /// missing Et
    const L1GctEtMiss missingEt(unsigned int bxInEvent) const;
    const L1GctEtMiss missingEt() const;

    /// total Et
    const L1GctEtTotal totalEt(unsigned int bxInEvent) const;
    const L1GctEtTotal totalEt() const;

    /// total calibrated Et in jets
    const L1GctEtHad totalHt(unsigned int bxInEvent) const;
    const L1GctEtHad totalHt() const;

    /// jet counts
    const L1GctJetCounts jetCounts(unsigned int bxInEvent) const;
    const L1GctJetCounts jetCounts() const;
  
    /// set candidate data words (all non-empty candidates)
    ///     two arguments:   candidate index, bunch cross index
    ///     one argument:    candidate index, for bunch cross with L1Accept
  
    /// muon

    void setMuons(const std::vector<MuonDataWord>&, unsigned int bxInEvent);
    void setMuons(const std::vector<MuonDataWord>&);

    /// electron
    void setElectrons(const std::vector<CaloDataWord>&, unsigned int bxInEvent);
    void setElectrons(const std::vector<CaloDataWord>&);
  
    /// isolated electron
    void setIsolatedElectrons(const std::vector<CaloDataWord>&, unsigned int bxInEvent); 
    void setIsolatedElectrons(const std::vector<CaloDataWord>&); 
  
    /// central jets
    void setCentralJets(const std::vector<CaloDataWord>&, unsigned int bxInEvent); 
    void setCentralJets(const std::vector<CaloDataWord>&); 
  
    /// forward jets
    void setForwardJets(const std::vector<CaloDataWord>&, unsigned int bxInEvent);
    void setForwardJets(const std::vector<CaloDataWord>&);
  
    /// tau jets
    void setTauJets(const std::vector<CaloDataWord>&, unsigned int bxInEvent);
    void setTauJets(const std::vector<CaloDataWord>&);
  
    /// missing Et
    void setMissingEt(const CaloMissingEtWord&, unsigned int bxInEvent);
    void setMissingEt(const CaloMissingEtWord&);
  
    /// total Et
    void setTotalEt(const CaloDataWord&, unsigned int bxInEvent);
    void setTotalEt(const CaloDataWord&);

    /// total calibrated Et
    void setTotalHt(const CaloDataWord&, unsigned int bxInEvent);
    void setTotalHt(const CaloDataWord&);
  
    /// jet count
    void setJetCounts(const CaloJetCountsWord&, unsigned int bxInEvent);
    void setJetCounts(const CaloJetCountsWord&);

    
    /// print all L1 Trigger Objects
    void printL1Objects(std::ostream& myCout, unsigned int bxInEventValue) const;
    void printL1Objects(std::ostream& myCout) const;


    //**************************************************************************
    // get/set hardware-related words
    // 
    //**************************************************************************


    /// get / set GTFE word (record) in the GT readout record
    const L1GtfeWord gtfeWord() const;
    void setGtfeWord(const L1GtfeWord&);

    /// get / set FDL word (record) in the GT readout record
    const L1GtFdlWord gtFdlWord(unsigned int bxInEvent) const;
    const L1GtFdlWord gtFdlWord() const;

    void setGtFdlWord(const L1GtFdlWord&, unsigned int bxInEvent);
    void setGtFdlWord(const L1GtFdlWord&);

    /// get / set PSB word (record) in the GT readout record
    const L1GtPsbWord gtPsbWord(unsigned int bxInEvent) const;
    const L1GtPsbWord gtPsbWord() const;

    void setGtPsbWord(const L1GtPsbWord&, unsigned int bxInEvent);
    void setGtPsbWord(const L1GtPsbWord&);

    // other methods
    
    /// clear the record
    void reset(); 
  
    /// output stream operator
    friend std::ostream& operator<<(std::ostream&, const L1GlobalTriggerReadoutRecord&);


private:

    L1GtfeWord m_gtfeWord;
    
    std::vector<L1GtFdlWord> m_gtFdlWord;
    
    std::vector<L1GtPsbWord> m_gtPsbWord;

    edm::RefProd<L1MuGMTReadoutCollection> m_muCollRefProd;        
       
};


#endif
