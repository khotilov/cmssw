#ifndef L1RCTCrate_h
#define L1RCTCRate_h

#include <vector>
#include "L1Trigger/RegionalCaloTrigger/interface/L1RCTReceiverCard.h"
#include "L1Trigger/RegionalCaloTrigger/interface/L1RCTElectronIsolationCard.h"
#include "L1Trigger/RegionalCaloTrigger/interface/L1RCTJetSummaryCard.h"

class L1RCTCrate {

 public:
  L1RCTCrate(int crtNo);

  int crateNumber(){return crtNo;}

  //For sharing information between crates.
  //It passes the pointers to the cards rather than the copies of the cards
  //because we need to modify the actual regions when setting their
  //neighbors rather than just copies.
  //Working in non garbage collected languages can really suck sometimes.
  L1RCTReceiverCard* getReceiverCard(int i) { return &receiverCards.at(i);}
  
  //This method receives the input from the L1RCT class and distributes
  //the RCInput to the 7 receiver cards and sends the HFInput straight
  //to the JSC for this crate.  The RCs never see the HF data.  Instead
  //the JSC acts like a primitive RC for these regions.
  void input(vector<vector<unsigned short> > RCInput,
	     vector<unsigned short> HFInput);
  //The two following are methods for running the actual data processing
  //in the RCs and the EICs.  They're to be called for each card
  //from the L1RCT process method
  void processReceiverCards();
  void fillElectronIsolationCards();
  void processElectronIsolationCards();
  //Pulls the information from the RCs and EICs and sends it to the 
  //JSC.
  void fillJetSummaryCard();
  void processJetSummaryCard();
  void print();
  void printJSC(){
    jetSummaryCard.print();
  }
  void printRC(int i){
    receiverCards.at(i).print();
  }
  void printEIC(int i){
    electronCards.at(i).print();
  }
  void printEICEdges(int i){
    electronCards.at(i).printEdges();
  }

  vector<unsigned short> getJetRegions(){
    jetSummaryCard.getJetRegions();
  }
  
  vector<unsigned short> getIsolatedEGObjects(){
    jetSummaryCard.getIsolatedEGObjects();
  }
  vector<unsigned short> getNonisolatedEGObjects(){
    jetSummaryCard.getNonisolatedEGObjects();
  }

 private:
  //The seven RCs and EICs
  //They are laid out according to CMS IN 2004/008
  //Increasing number towards higher absolute eta
  //The seventh card is always sideways with respect to the
  //other six.
  vector<L1RCTReceiverCard> receiverCards;
  vector<L1RCTElectronIsolationCard> electronCards;
  //The JSC receives the jet and electron information from the 
  //RCs and EICs.  There is only one per crate.
  L1RCTJetSummaryCard jetSummaryCard;
  
  int crtNo;

  L1RCTCrate();

  //L1RCTJetCaptureCard jetCaptureCard;
};
#endif
