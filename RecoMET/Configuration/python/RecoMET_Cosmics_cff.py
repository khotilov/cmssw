import FWCore.ParameterSet.Config as cms

# Name:   RecoMET.cff
# Original Author: R.Cavanaugh
# Date:   05.11.2006
# Notes:  CaloMET.cfi assumes that a product with label "caloTowers" is
#         already written into the event.
# Modification by F. Ratnikov and R. Remington
# Date: 10/21/08 
# Addition of MET significance by F.Blekman
# Date: 10/23/08
# Addition of HCAL noise by JP Chou
# Date:  3/26/09

import FWCore.ParameterSet.Config as cms
from RecoMET.Configuration.RecoMET_cff import *

tcMetP5 = tcMet.clone(trackInputTag = 'ctfWithMaterialTracksP5LHCNavigation')
hcalnoise_cosmics = hcalnoise.clone(fillTracks = False)

metreco_cosmics = cms.Sequence(
        met+metNoHF+metHO+metNoHFHO+
            calotoweroptmaker+metOpt+metOptNoHF+calotoweroptmakerWithHO+metOptHO+metOptNoHFHO+
            htMetSC5+htMetSC7+htMetKT4+htMetKT6+htMetIC5+muonMETValueMapProducer+corMetGlobalMuons+muonTCMETValueMapProducer+tcMetP5
            )

metrecoPlusHCALNoise_cosmics = cms.Sequence( metreco_cosmics + hcalnoise )



