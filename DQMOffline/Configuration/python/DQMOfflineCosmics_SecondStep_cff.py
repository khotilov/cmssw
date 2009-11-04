import FWCore.ParameterSet.Config as cms

from CondTools.DQM.DQMReferenceHistogramRootFileEventSetupAnalyzer_cfi import *
from DQMServices.Components.DQMMessageLoggerClient_cfi import *

from DQMOffline.Ecal.ecal_dqm_client_offline_cosmic_cff import *
from DQM.HcalMonitorModule.hcal_dqm_client_fileT0_cff import *
from DQM.SiStripMonitorClient.SiStripClientConfig_Tier0_cff import *
from DQM.SiPixelCommon.SiPixelOfflineDQM_client_cff import *
from DQM.DTMonitorClient.dtDQMOfflineClients_Cosmics_cff import *
from DQM.RPCMonitorClient.RPCTier0Client_cff import *
from DQM.CSCMonitorModule.csc_dqm_offlineclient_cosmics_cff import *
from DQM.EcalPreshowerMonitorClient.es_dqm_client_offline_cosmic_cff import *
from DQMServices.Components.DQMFEDIntegrityClient_cff import *

DQMOfflineCosmics_SecondStep_PreDPG = cms.Sequence( ecal_dqm_client_offline *
                                                    hcalOfflineDQMClient *
                                                    SiStripOfflineDQMClient *
                                                    sipixelEDAClient *
                                                    dtClientsCosmics *
                                                    rpcTier0Client *
                                                    cscOfflineCosmicsClients *
                                                    es_dqm_client_offline *
                                                    dqmFEDIntegrityClient )


DQMOfflineCosmics_SecondStepDPG = cms.Sequence( dqmRefHistoRootFileGetter *
                                                DQMOfflineCosmics_SecondStep_PreDPG *
                                                DQMMessageLoggerClient )

from DQMOffline.Muon.muonQualityTests_cff import *
from DQMOffline.EGamma.photonOfflineDQMClient_cff import *
from DQMOffline.Trigger.DQMOffline_Trigger_Client_cff import *
from DQMOffline.Trigger.DQMOffline_HLT_Client_cff import *

DQMOfflineCosmics_SecondStep_PrePOG = cms.Sequence( cosmicMuonQualityTests *
                                                    photonOfflineDQMClient *
                                                    triggerOfflineDQMClient *
                                                    hltOfflineDQMClient )
 
DQMOfflineCosmics_SecondStepPOG = cms.Sequence( dqmRefHistoRootFileGetter *
                                                DQMOfflineCosmics_SecondStep_PrePOG *
                                                DQMMessageLoggerClient )

DQMOfflineCosmics_SecondStep = cms.Sequence( dqmRefHistoRootFileGetter *
                                             DQMOfflineCosmics_SecondStep_PreDPG *
                                             DQMOfflineCosmics_SecondStep_PrePOG *
                                             DQMMessageLoggerClient )

