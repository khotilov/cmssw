import FWCore.ParameterSet.Config as cms

from DQMOffline.Muon.muonQualityTests_cff import *
from DQMOffline.Ecal.ecal_dqm_client_offline_cosmic_cff import *
from DQM.HcalMonitorModule.hcal_dqm_client_fileT0_cff import *
from DQM.SiPixelCommon.SiPixelOfflineDQM_client_cff import *
from DQMOffline.Trigger.DQMOffline_Trigger_Client_cff import *
from DQM.SiStripMonitorClient.SiStripClientConfig_Tier0_cff import *
from DQM.DTMonitorClient.dtDQMOfflineClients_cff import *
from DQMOffline.Trigger.DQMOffline_HLT_Client_cff import *
from DQM.RPCMonitorClient.RPCTier0Client_cff import *
from DQMOffline.EGamma.photonOfflineDQMClient_cff import *

from DQMServices.Components.DQMFEDIntegrityClient_cff import *

DQMOfflineCosmics_SecondStep = cms.Sequence(ecal_dqm_client_offline*cosmicMuonQualityTests*hcalOfflineDQMClient*sipixelEDAClient*triggerOfflineDQMClient*SiStripOfflineDQMClient*dtClients*hltOfflineDQMClient*rpcTier0Client*photonOfflineDQMClient*dqmFEDIntegrityClient)

