#ifndef JetMETCorrections_METPUSubtraction_GBRForestWriter_h
#define JetMETCorrections_METPUSubtraction_GBRForestWriter_h

/** \class GBRForestWriter
 *
 * Read GBRForest objects from ROOT file input
 * and store it in SQL-lite output file
 *
 * \authors Christian Veelken, LLR
 *
 * \version $Revision: 1.3 $
 *
 * $Id: GBRForestWriter.h,v 1.3 2013/02/25 15:08:39 veelken Exp $
 *
 */

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include <vector>
#include <string>

class GBRForestWriter : public edm::EDAnalyzer 
{
 public:
  GBRForestWriter(const edm::ParameterSet&);
  ~GBRForestWriter();
  
 private:
  virtual void analyze(const edm::Event&, const edm::EventSetup&);

  std::string moduleLabel_;

  bool hasRun_;

  struct jobEntryType
  {
    edm::FileInPath inputFileName_;
    std::string gbrForestName_;
    std::string outputRecord_;
  };
  std::vector<jobEntryType*> jobs_;
};

#endif
