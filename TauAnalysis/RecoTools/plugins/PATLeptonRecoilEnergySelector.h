//
// $Id: PATLeptonRecoilEnergySelector.h,v 1.2 2009/10/25 12:38:23 veelken Exp $
//

#ifndef TauAnalysis_RecoTools_PATLeptonRecoilEnergySelector_h
#define TauAnalysis_RecoTools_PATLeptonRecoilEnergySelector_h

#include "DataFormats/Common/interface/RefVector.h"

#include "CommonTools/UtilAlgos/interface/StringCutObjectSelector.h"
#include "CommonTools/UtilAlgos/interface/SingleObjectSelector.h"
#include "CommonTools/UtilAlgos/interface/ObjectSelector.h"
#include "CommonTools/UtilAlgos/interface/SingleElementCollectionSelector.h"

#include "AnalysisDataFormats/TauAnalysis/interface/PATLeptonRecoilEnergy.h"

#include <vector>

typedef SingleObjectSelector<
            std::vector<PATTauRecoilEnergyFromJets>,
            StringCutObjectSelector<PATTauRecoilEnergyFromJets>
        > PATTauRecoilEnergyFromJetsSelector;
typedef SingleObjectSelector<
            std::vector<PATTauRecoilEnergyFromCaloTowers>,
            StringCutObjectSelector<PATTauRecoilEnergyFromCaloTowers>
        > PATTauRecoilEnergyFromCaloTowersSelector;

typedef SingleObjectSelector<
            std::vector<PATTauRecoilEnergyFromJets>,
            StringCutObjectSelector<PATTauRecoilEnergyFromJets>,
            edm::RefVector<std::vector<PATTauRecoilEnergyFromJets> >
        > PATTauRecoilEnergyFromJetsRefSelector;
typedef SingleObjectSelector<
            std::vector<PATTauRecoilEnergyFromCaloTowers>,
            StringCutObjectSelector<PATTauRecoilEnergyFromCaloTowers>,
            edm::RefVector<std::vector<PATTauRecoilEnergyFromCaloTowers> >
        > PATTauRecoilEnergyFromCaloTowersRefSelector;

#endif
