#include <utility>
#include <vector>

#include "DataFormats/Common/interface/Wrapper.h"
#include "SimDataFormats/JetMatching/interface/JetFlavour.h"
#include "SimDataFormats/JetMatching/interface/JetFlavourMatching.h"
#include "SimDataFormats/JetMatching/interface/MatchedPartons.h"
#include "SimDataFormats/JetMatching/interface/JetMatchedPartons.h"

namespace {
  namespace {
    reco::JetFlavour                                 jf;
    reco::JetFlavourMatchingCollectionBase           jfmcb;
    reco::JetFlavourMatchingCollection               jfmc;
    reco::JetFlavourMatchingRef                      jfmr;
    reco::JetFlavourMatchingRefProd                  jfmrp;
    reco::JetFlavourMatchingRefVector                jfrv;
    edm::Wrapper<reco::JetFlavourMatchingCollection> wjfmc;
    reco::MatchedPartons                             mp;
    reco::JetMatchedPartonsCollectionBase           jmpcb;
    reco::JetMatchedPartonsCollection               jmpc;
    reco::JetMatchedPartonsRef                      jmpr;
    reco::JetMatchedPartonsRefProd                  jmprp;
    reco::JetMatchedPartonsRefVector                jmpv;
    edm::Wrapper<reco::JetMatchedPartonsCollection> wjmpc;
  }
}
