#include <vector>
#include <string>

#include "FWCore/Utilities/interface/Exception.h"
#include "DataFormats/BTauReco/interface/BaseTagInfo.h"
#include "RecoBTag/MCTools/interface/JetFlavour.h"
#include "RecoBTag/Analysis/interface/BaseTagInfoPlotter.h"

using namespace std;
using namespace reco;

void BaseTagInfoPlotter::analyzeTag(const BaseTagInfo * tagInfo, const BTagMCTools::JetFlavour & jetFlavour)
{
  throw cms::Exception("MissingVirtualMethod")
  	<< "No analyzeTag method overloaded from BaseTagInfoPlotter." << endl;
}

void BaseTagInfoPlotter::analyzeTag(const vector<const BaseTagInfo *> &tagInfos, const BTagMCTools::JetFlavour & jetFlavour)
{
  if (tagInfos.size() != 1)
    throw cms::Exception("MismatchedTagInfos")
    	<< tagInfos.size() << " BaseTagInfos passed, but only one expected." << endl;

  analyzeTag(tagInfos[0], jetFlavour);
}

void BaseTagInfoPlotter::setEventSetup(const edm::EventSetup & setup)
{
}

vector<string> BaseTagInfoPlotter::tagInfoRequirements() const
{
  return vector<string>();
}
