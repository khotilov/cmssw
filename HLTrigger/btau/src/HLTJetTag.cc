/** \class HLTJetTag
 *
 *  This class is an HLTFilter (a spcialized EDFilter) implementing 
 *  tagged multi-jet trigger for b and tau. 
 *  It should be run after the normal multi-jet trigger.
 *
 *  $Date: 2012/02/06 10:06:49 $
 *  $Revision: 1.12 $
 *
 *  \author Arnaud Gay, Ian Tomalin
 *  \maintainer Andrea Bocci
 *
 */

#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/Common/interface/RefToBase.h"
#include "DataFormats/JetReco/interface/CaloJetCollection.h"
#include "DataFormats/BTauReco/interface/JetTag.h"
#include "DataFormats/HLTReco/interface/TriggerFilterObjectWithRefs.h"
#include "HLTrigger/HLTcore/interface/HLTFilter.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"

#include "HLTJetTag.h"

#include<vector>
#include<string>
#include<typeinfo>

//
// constructors and destructor
//

template<typename T>
HLTJetTag<T>::HLTJetTag(const edm::ParameterSet & config) : HLTFilter(config),
  m_Jets   (config.getParameter<edm::InputTag> ("Jets") ),
  m_JetTags(config.getParameter<edm::InputTag> ("JetTags") ),
  m_MinTag (config.getParameter<double>        ("MinTag") ),
  m_MaxTag (config.getParameter<double>        ("MaxTag") ),
  m_MinJets(config.getParameter<int>           ("MinJets") ),
  m_TriggerType(config.getParameter<int>       ("TriggerType") )
{

  edm::LogInfo("") << " (HLTJetTag) trigger cuts: " << std::endl
                   << "\ttype of        jets used: " << m_Jets.encode() << std::endl
                   << "\ttype of tagged jets used: " << m_JetTags.encode() << std::endl
                   << "\tmin/max tag value: [" << m_MinTag << ".." << m_MaxTag << "]" << std::endl
                   << "\tmin no. tagged jets: " << m_MinJets
		   << "\tTriggerType: " << m_TriggerType << std::endl;
}

template<typename T>
HLTJetTag<T>::~HLTJetTag()
{
}

template<typename T>
void
HLTJetTag<T>::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;
  makeHLTFilterDescription(desc);
  desc.add<edm::InputTag>("Jets",edm::InputTag("hltJetCollection"));
  desc.add<edm::InputTag>("JetTags",edm::InputTag("hltJetTagCollection"));
  desc.add<double>("MinTag",2.0);
  desc.add<double>("MaxTag",999999.0);
  desc.add<int>("MinJets",1);
  desc.add<int>("TriggerType",0);
  descriptions.add(std::string("hlt")+std::string(typeid(HLTJetTag<T>).name()),desc);
}

//
// member functions
//

// ------------ method called to produce the data  ------------
template<typename T>
bool
HLTJetTag<T>::hltFilter(edm::Event& event, const edm::EventSetup& setup, trigger::TriggerFilterObjectWithRefs & filterproduct)
{
  using namespace std;
  using namespace edm;
  using namespace reco;

  typedef vector<T> TCollection;
  typedef Ref<TCollection> TRef;

  edm::Handle<TCollection> h_Jets;
  event.getByLabel(m_Jets, h_Jets);
  if (saveTags()) filterproduct.addCollectionTag(m_Jets);

  edm::Handle<JetTagCollection> h_JetTags;
  event.getByLabel(m_JetTags, h_JetTags);

  TRef jetRef;

  // Look at all jets in decreasing order of Et.
  int nJet = 0;
  int nTag = 0;
  for (JetTagCollection::const_iterator jet = h_JetTags->begin(); jet != h_JetTags->end(); ++jet) {
    jetRef = TRef(h_Jets,jet->first.key());
    LogTrace("") << "Jet " << nJet
                 << " : Et = " << jet->first->et()
                 << " , tag value = " << jet->second;
    ++nJet;
    // Check if jet is tagged.
    if ( (m_MinTag <= jet->second) and (jet->second <= m_MaxTag) ) {
      ++nTag;

      // Store a reference to the jets which passed tagging cuts
      // N.B. this *should* work as we start from a CaloJet in HLT
      filterproduct.addObject(static_cast<trigger::TriggerObjectType>(m_TriggerType),jetRef);
    }
  }

  // filter decision
  bool accept = (nTag >= m_MinJets);

  edm::LogInfo("") << " trigger accept ? = " << accept
                   << " nTag/nJet = " << nTag << "/" << nJet << std::endl;

  return accept;
}
