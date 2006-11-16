#ifndef L1GCTJETFINALSTAGE_H_
#define L1GCTJETFINALSTAGE_H_

#include "L1Trigger/GlobalCaloTrigger/interface/L1GctJet.h"
#include "L1Trigger/GlobalCaloTrigger/interface/L1GctProcessor.h"
#include "L1Trigger/GlobalCaloTrigger/interface/L1GctWheelJetFpga.h"
//#include "L1Trigger/GlobalCaloTrigger/src/L1Gct.h"

#include <vector>

/*!
* \class L1GctJetFinalStage
* \brief Represents the final stage of L1 jet processing.
*
*  Takes as input the jet data from the two Wheel Jet FPGAs
*  and outputs the top four of each type of jet - central, forward,
*  and tau - for the whole of the CMS detector.
* 
* \author Jim Brooke & Robert Frazier
* \date June 2006
*/ 

class L1GctJetFinalStage : public L1GctProcessor
{
public:
  typedef std::vector<L1GctJet> JetVector;
  static const unsigned int MAX_WHEEL_FPGAS; ///< Max number of wheel FPGA pointers

  /// Takes a vector of 2 wheel jet FPGA pointers, with which to get input data from
	L1GctJetFinalStage(std::vector<L1GctWheelJetFpga*> m_wheelFpgas);
	~L1GctJetFinalStage();

  /// Overload << operator
  friend std::ostream& operator << (std::ostream& os, const L1GctJetFinalStage& fpga);

  /// clear internal buffers
  virtual void reset();

  /// get input data from sources
  virtual void fetchInput();

  /// process the data, fill output buffers
  virtual void process();
    	
  void setInputCentralJet(int i, L1GctJet jet);  ///< set the central jets input data
  void setInputForwardJet(int i, L1GctJet jet);  ///< set the forward jets input data
  void setInputTauJet(int i, L1GctJet jet);      ///< set the tau jets input data

  JetVector getInputCentralJets() const { return m_inputCentralJets; } ///< get the central jets input data
  JetVector getInputForwardJets() const { return m_inputForwardJets; } ///< get the forward jets input data
  JetVector getInputTauJets() const { return m_inputTauJets; }         ///< get the tau jets input data

  JetVector getCentralJets() const { return m_centralJets; } ///< get the central jets output data
  JetVector getForwardJets() const { return m_forwardJets; } ///< get the forward jets output data
  JetVector getTauJets() const { return m_tauJets; }         ///< get the tau jets output data
private:

  static const int MAX_JETS_IN;  ///< Max number of jets of each type coming in
  static const int MAX_JETS_OUT; ///< Max number of jets of each type going out
  
  /// wheel jet FPGAs
  std::vector<L1GctWheelJetFpga*> m_wheelFpgas;

  // input data
  JetVector m_inputCentralJets;
  JetVector m_inputForwardJets;
  JetVector m_inputTauJets;

  // output data
  JetVector m_centralJets;
  JetVector m_forwardJets;
  JetVector m_tauJets;

  //PRIVATE MEMBER FUNCTIONS
  ///Enters jets into the specified storageVector, according to which wheel card we are taking them from.
  void storeJets(JetVector& storageVector, JetVector jets, unsigned short iWheel);
  
};

std::ostream& operator << (std::ostream& os, const L1GctJetFinalStage& fpga);

#endif /*L1GCTJETFINALSTAGE_H_*/
