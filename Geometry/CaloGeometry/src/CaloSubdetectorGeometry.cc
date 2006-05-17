#include "Geometry/CaloGeometry/interface/CaloCellGeometry.h"
#include "Geometry/CaloGeometry/interface/CaloSubdetectorGeometry.h"
#include "Geometry/CaloGeometry/interface/CaloGeometry.h"


void CaloSubdetectorGeometry::addCell(const DetId& id, const CaloCellGeometry* ccg) {
  cellGeometries_.insert(std::pair<DetId,const CaloCellGeometry*>(id,ccg));
}

const CaloCellGeometry* CaloSubdetectorGeometry::getGeometry(const DetId& id) const {
  std::map<DetId, const CaloCellGeometry*>::const_iterator i=cellGeometries_.find(id);
  return i==(cellGeometries_.end())?(0):(i->second);
}

bool CaloSubdetectorGeometry::present(const DetId& id) const {
  std::map<DetId, const CaloCellGeometry*>::const_iterator i=cellGeometries_.find(id);
  return i!=cellGeometries_.end();
}


std::vector<DetId> CaloSubdetectorGeometry::getValidDetIds(DetId::Detector det, int subdet) const {
  if (validIds_.empty()) {
    std::map<DetId, const CaloCellGeometry*>::const_iterator i;
    for (i=cellGeometries_.begin(); i!=cellGeometries_.end(); i++)
      validIds_.push_back(i->first);
  }

  return validIds_;    
}

const DetId CaloSubdetectorGeometry::getClosestCell(const GlobalPoint& r) const 
{
  return DetId(0);
}

