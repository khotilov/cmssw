#include "FWCore/PluginManager/interface/ModuleDef.h"

#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/SourceFactory.h"

#include "CalibMuon/DTCalibration/plugins/DTMapGenerator.h"
#include "CalibMuon/DTCalibration/plugins/DTTTrigCalibration.h"
#include "CalibMuon/DTCalibration/plugins/DTTTrigWriter.h"
#include "CalibMuon/DTCalibration/plugins/DTTTrigCorrection.h"
#include "CalibMuon/DTCalibration/plugins/DTT0Calibration.h"
#include "CalibMuon/DTCalibration/plugins/DTT0CalibrationNew.h"
#include "CalibMuon/DTCalibration/plugins/DTTPDeadWriter.h"
#include "CalibMuon/DTCalibration/plugins/DTVDriftCalibration.h"
#include "CalibMuon/DTCalibration/plugins/DTVDriftWriter.h"
#include "CalibMuon/DTCalibration/plugins/DTNoiseComputation.h"
#include "CalibMuon/DTCalibration/plugins/DTNoiseCalibration.h"
#include "CalibMuon/DTCalibration/plugins/DTFakeTTrigESProducer.h"
#include "CalibMuon/DTCalibration/plugins/DTFakeT0ESProducer.h"
#include "CalibMuon/DTCalibration/plugins/DTFakeVDriftESProducer.h"

DEFINE_SEAL_MODULE();
DEFINE_ANOTHER_FWK_MODULE(DTMapGenerator);
DEFINE_ANOTHER_FWK_MODULE(DTTTrigCalibration);
DEFINE_ANOTHER_FWK_MODULE(DTTTrigWriter);
DEFINE_ANOTHER_FWK_MODULE(DTTTrigCorrection);
DEFINE_ANOTHER_FWK_MODULE(DTT0Calibration);
DEFINE_ANOTHER_FWK_MODULE(DTT0CalibrationNew);
DEFINE_ANOTHER_FWK_MODULE(DTTPDeadWriter);
DEFINE_ANOTHER_FWK_MODULE(DTVDriftCalibration);
DEFINE_ANOTHER_FWK_MODULE(DTVDriftWriter);
DEFINE_ANOTHER_FWK_MODULE(DTNoiseComputation);
DEFINE_ANOTHER_FWK_MODULE(DTNoiseCalibration);
DEFINE_ANOTHER_FWK_EVENTSETUP_SOURCE(DTFakeTTrigESProducer);
DEFINE_ANOTHER_FWK_EVENTSETUP_SOURCE(DTFakeT0ESProducer);
DEFINE_ANOTHER_FWK_EVENTSETUP_SOURCE(DTFakeVDriftESProducer);

