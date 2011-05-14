#ifndef TauAnalysis_RecoTools_VertexMultiplicityReweightExtractor_h
#define TauAnalysis_RecoTools_VertexMultiplicityReweightExtractor_h

/** \class VertexMultiplicityReweightProducer
 *
 * Reweight Monte Carlo events simulated with pile-up to match the
 * vertex multiplicity distribution observed in the analyzed data sample
 *
 * NOTE:
 *      o weight > 1: fraction of events which given vertex multiplicity higher in data than in (pile-up) MC
 *      o weight = 1: fraction of events which given vertex multiplicity same   in data than in (pile-up) MC
 *      o weight < 1: fraction of events which given vertex multiplicity lower  in data than in (pile-up) MC
 *
 * \authors Christian Veelken
 *
 * \version $Revision: 1.1 $
 *
 * $Id: VertexMultiplicityReweightExtractor.h,v 1.1 2011/04/13 17:05:13 veelken Exp $
 *
 */

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"

#include "CommonTools/Utils/interface/StringObjectFunction.h"

#include "TauAnalysis/BgEstimationTools/interface/ObjValExtractorBase.h"

#include <TFile.h>
#include <TH1.h>

class VertexMultiplicityReweightExtractor : public ObjValExtractorBase
{
 public:
  explicit VertexMultiplicityReweightExtractor(const edm::ParameterSet&);
  ~VertexMultiplicityReweightExtractor();

  double operator()(const edm::Event&) const;

 private:
  edm::InputTag src_;
  
  TFile* inputFile_;
  TH1* lut_;

  enum { kUndefined, kGenLevel, kRecLevel };
  int type_;
};

#endif

