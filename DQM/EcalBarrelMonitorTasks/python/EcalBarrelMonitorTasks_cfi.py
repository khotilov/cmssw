import FWCore.ParameterSet.Config as cms

from DQM.EcalBarrelMonitorTasks.EBOccupancyTask_cfi import *
from DQM.EcalBarrelMonitorTasks.EBIntegrityTask_cfi import *
from DQM.EcalBarrelMonitorTasks.EBStatusFlagsTask_cfi import *
from DQM.EcalBarrelMonitorTasks.EBCosmicTask_cfi import *
from DQM.EcalBarrelMonitorTasks.EBLaserTask_cfi import *
from DQM.EcalBarrelMonitorTasks.EBPedestalOnlineTask_cfi import *
from DQM.EcalBarrelMonitorTasks.EBPedestalTask_cfi import *
from DQM.EcalBarrelMonitorTasks.EBTestPulseTask_cfi import *
from DQM.EcalBarrelMonitorTasks.EBTriggerTowerTask_cfi import *
from DQM.EcalBarrelMonitorTasks.EBTimingTask_cfi import *
from DQM.EcalBarrelMonitorTasks.EBBeamHodoTask_cfi import *
from DQM.EcalBarrelMonitorTasks.EBBeamCaloTask_cfi import *
from DQM.EcalBarrelMonitorTasks.EBClusterTask_cfi import *
from DQM.EcalBarrelMonitorTasks.EBSelectiveReadoutTask_cfi import *

ecalBarrelDefaultTasksSequence = cms.Sequence(ecalBarrelOccupancyTask*ecalBarrelIntegrityTask*ecalBarrelStatusFlagsTask*ecalBarrelLaserTask*ecalBarrelPedestalOnlineTask*ecalBarrelPedestalTask*ecalBarrelTestPulseTask*ecalBarrelTriggerTowerTask*ecalBarrelTimingTask)

ecalBarrelCosmicTasksSequence = cms.Sequence(ecalBarrelDefaultTasksSequence*ecalBarrelCosmicTask)

ecalBarrelTestBeamTasksSequence = cms.Sequence(ecalBarrelDefaultTasksSequence*ecalBarrelBeamHodoTask*ecalBarrelBeamCaloTask)

