import FWCore.ParameterSet.Config as cms

from DQM.EcalEndcapMonitorTasks.EEOccupancyTask_cfi import *
from DQM.EcalEndcapMonitorTasks.EEIntegrityTask_cfi import *
from DQM.EcalEndcapMonitorTasks.EEStatusFlagsTask_cfi import *
from DQM.EcalEndcapMonitorTasks.EECosmicTask_cfi import *
from DQM.EcalEndcapMonitorTasks.EELaserTask_cfi import *
from DQM.EcalEndcapMonitorTasks.EELedTask_cfi import *
from DQM.EcalEndcapMonitorTasks.EEPedestalOnlineTask_cfi import *
from DQM.EcalEndcapMonitorTasks.EEPedestalTask_cfi import *
from DQM.EcalEndcapMonitorTasks.EETestPulseTask_cfi import *
from DQM.EcalEndcapMonitorTasks.EETriggerTowerTask_cfi import *
from DQM.EcalEndcapMonitorTasks.EETimingTask_cfi import *
from DQM.EcalEndcapMonitorTasks.EEBeamHodoTask_cfi import *
from DQM.EcalEndcapMonitorTasks.EEBeamCaloTask_cfi import *
from DQM.EcalEndcapMonitorTasks.EEClusterTask_cfi import *
from DQM.EcalEndcapMonitorTasks.EESelectiveReadoutTask_cfi import *

ecalEndcapDefaultTasksSequence = cms.Sequence(ecalEndcapOccupancyTask*ecalEndcapIntegrityTask*ecalEndcapStatusFlagsTask*ecalEndcapLaserTask*ecalEndcapLedTask*ecalEndcapPedestalOnlineTask*ecalEndcapPedestalTask*ecalEndcapTestPulseTask*ecalEndcapTriggerTowerTask*ecalEndcapTimingTask)

ecalEndcapCosmicTasksSequence = cms.Sequence(ecalEndcapDefaultTasksSequence*ecalEndcapCosmicTask)

ecalEndcapTestBeamTasksSequence = cms.Sequence(ecalEndcapDefaultTasksSequence*ecalEndcapBeamHodoTask*ecalEndcapBeamCaloTask)

