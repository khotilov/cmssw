#ifndef L1GCTELECTRONSORTER_H_
#define L1GCTELECTRONSORTER_H_

#include "DataFormats/L1GlobalCaloTrigger/interface/L1GctDigis.h"
#include "L1Trigger/GlobalCaloTrigger/interface/L1GctProcessor.h"
#include <vector>
#include <functional>
#include <ostream>


/*
 * \Class L1GctElectronSorter
 * \Class that sorts electron candidates
 *
 * This class can be constructed to sort iso or non-iso
 * electron candidates from 4 or 5 source cards.
 * The electrons are sorted in ascending order and the 4
 * highest in rank will be returned.
 * It represents the 1st stage sorter FPGA's on the electron
 * leaf cards.
 *
 * \author  Maria Hansen
 * \date    21/04/06
 */

using std::binary_function;

class L1GctSourceCard;

class L1GctElectronSorter : public L1GctProcessor
{
 public:
  
  /// constructor; set type (isolated or non-isolated)
  L1GctElectronSorter(int nInputs, bool iso,  std::vector<L1GctSourceCard*> sCards);
  /// constructor for test purposes, takes no of inputs and type
  L1GctElectronSorter(int nInputs, bool iso);
  ///   
  ~L1GctElectronSorter();
  ///
  /// clear internal buffers
  virtual void reset();
  ///
  /// get input data from sources
  virtual void fetchInput();
  ///
  /// process the data, fill output buffers
  virtual void process();
  ///
  /// set input candidate
  void setInputEmCand(int i, L1GctEmCand cand);
  ///	
  /// get input candidates
  inline std::vector<L1GctEmCand> getInputCands() { return m_inputCands; }
  ///
  /// get output candidates
  inline std::vector<L1GctEmCand> getOutputCands() { return m_outputCands; }
  ///
  /// overload of cout operator
  friend std::ostream& operator<<(std::ostream& s,const L1GctElectronSorter& card);  
 private:
  /// comparison operator for sort
  struct rank_gt : public binary_function<L1GctEmCand, L1GctEmCand, bool> {
    bool operator()(const L1GctEmCand& x, const L1GctEmCand& y) { return x.rank() > y.rank(); }
  };
  
 private:
  ///
  /// algo ID (is it FPGA 1 or 2 processing)
  int m_id;
  ///
  /// type of electrons to sort (isolated or non isolated)
  bool m_emCandType;
  ///
  /// source card input
  std::vector<L1GctSourceCard*> m_theSCs;
  ///
  /// input data
  std::vector<L1GctEmCand> m_inputCands;
  ///
  /// output data
  std::vector<L1GctEmCand> m_outputCands;
  
};

std::ostream& operator<<(std::ostream& s,const L1GctElectronSorter& card); 
#endif /*L1GCTELECTRONSORTER_H_*/
