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

from RecoJets.Configuration.CaloTowersES_cfi import *
from RecoMET.METProducers.CaloTowersOpt_cfi import *
from RecoMET.METProducers.CaloMET_cfi import *
from RecoMET.METProducers.HTMET_cfi import *
from RecoMET.METProducers.CaloMETSignif_cfi import *
from RecoMET.METProducers.TCMET_cfi import *
#sequence metreco = {met, metsig, htMetIC5, htMetMC5}
metreco = cms.Sequence(
        met+metNoHF+metHO+metNoHFHO+
            calotoweroptmaker+metOpt+metOptNoHF+calotoweroptmakerWithHO+metOptHO+metOptNoHFHO+
            htMetSC5+htMetSC7+htMetKT4+htMetKT6+htMetIC5+tcMet
            )

