/** \file
 *
 *  $Date: 2009/01/16 16:43:46 $
 *  $Revision: 1.2 $
 */

#include "MagneticField/GeomBuilder/plugins/VolumeBasedMagneticFieldESProducer.h"
#include "MagneticField/VolumeBasedEngine/interface/VolumeBasedMagneticField.h"

#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"

#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/EventSetup.h"

#include "DetectorDescription/Core/interface/DDCompactView.h"
#include "FWCore/Framework/interface/ModuleFactory.h"
#include "MagneticField/GeomBuilder/src/MagGeoBuilderFromDDD.h"

#include <string>
#include <vector>
#include <iostream>

using namespace std;
using namespace magneticfield;

VolumeBasedMagneticFieldESProducer::VolumeBasedMagneticFieldESProducer(const edm::ParameterSet& iConfig) : pset(iConfig)
{
  setWhatProduced(this, pset.getUntrackedParameter<std::string>("label",""));
}




// ------------ method called to produce the data  ------------
std::auto_ptr<MagneticField> VolumeBasedMagneticFieldESProducer::produce(const IdealMagneticFieldRecord & iRecord)
{
  edm::ESHandle<DDCompactView> cpv;
  iRecord.get("magfield",cpv );
  MagGeoBuilderFromDDD builder(pset.getParameter<std::string>("version"),
			       pset.getUntrackedParameter<bool>("debugBuilder", false));

  // Get scaling factors
  vector<int> keys = pset.getParameter<vector<int> >("scalingVolumes");
  vector<double> values = pset.getParameter<vector<double> >("scalingFactors");

  if (keys.size() != 0) {
    builder.setScaling(keys, values);
  }

  builder.build(*cpv);


  // Get slave field
  edm::ESHandle<MagneticField> paramField;
  if (pset.getParameter<bool>("useParametrizedTrackerField")) {;
    iRecord.get("parametrizedField",paramField);
  }
  std::auto_ptr<MagneticField> s(new VolumeBasedMagneticField(pset,builder.barrelLayers(), builder.endcapSectors(), builder.barrelVolumes(), builder.endcapVolumes(), builder.maxR(), builder.maxZ(), paramField.product(), false));
  return s;
}

DEFINE_FWK_EVENTSETUP_MODULE(VolumeBasedMagneticFieldESProducer);
