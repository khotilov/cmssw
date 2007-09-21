#include "Geometry/HcalTowerAlgo/interface/HcalDDDGeometryLoader.h"
#include "DataFormats/HcalDetId/interface/HcalDetId.h"
#include "Geometry/CaloGeometry/interface/CaloCellGeometry.h"
#include "Geometry/CaloGeometry/interface/IdealObliquePrism.h"
#include "Geometry/CaloGeometry/interface/IdealZPrism.h"
#include "Geometry/CaloGeometry/interface/CaloGeometryLoader.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include<string>

HcalDDDGeometryLoader::HcalDDDGeometryLoader(const DDCompactView & cpv) {

  std::string name = "HcalHits";
  numberingFromDDD = new HcalNumberingFromDDD(name, cpv);
}

HcalDDDGeometryLoader::~HcalDDDGeometryLoader() {
  delete numberingFromDDD;
}


std::auto_ptr<CaloSubdetectorGeometry> HcalDDDGeometryLoader::load(DetId::Detector det, int subdet) {

  HcalSubdetector  hsub        = static_cast<HcalSubdetector>(subdet);
  HcalDDDGeometry* geometryDDD = new HcalDDDGeometry();
  std::auto_ptr<CaloSubdetectorGeometry> hg(geometryDDD);
  fill (hsub, geometryDDD, hg.get());
  return hg;
}

std::auto_ptr<CaloSubdetectorGeometry> HcalDDDGeometryLoader::load() {

  HcalDDDGeometry* geometryDDD = new HcalDDDGeometry();
  std::auto_ptr<CaloSubdetectorGeometry> hg(geometryDDD);
  fill(HcalBarrel,  geometryDDD, hg.get()); 
  fill(HcalEndcap,  geometryDDD, hg.get()); 
  fill(HcalForward, geometryDDD, hg.get()); 
  fill(HcalOuter,   geometryDDD, hg.get());
  return hg;
}

void HcalDDDGeometryLoader::fill(HcalSubdetector subdet, 
				 HcalDDDGeometry* geometryDDD,
				 CaloSubdetectorGeometry* geom) {

  // start by making the new HcalDetIds
  std::vector<HcalCellType::HcalCellType> hcalCells = numberingFromDDD->HcalCellTypes(subdet);
  geometryDDD->insertCell(hcalCells);
  LogDebug("HCalGeom") << "HcalDDDGeometryLoader::fill gets " 
		       << hcalCells.size() << " cells for subdetector " 
		       << subdet;
  if( geom->cornersMgr() == 0 ) geom->allocateCorners( 11000 ) ;
  if( geom->parMgr()     == 0 ) geom->allocatePar( 500, 3 ) ;
			 
  // Make the new HcalDetIds and the cells
  double deg = M_PI/180.;
  std::vector<HcalDetId> hcalIds;
  for (unsigned int i=0; i<hcalCells.size(); i++) {
    int etaRing  = hcalCells[i].etaBin();
    int depthBin = hcalCells[i].depthSegment();
    int phiInc   = 4/hcalCells[i].nPhiModule();
    unsigned int iphi = 1;
    double  dphi = (hcalCells[i].phiBinWidth())*deg;
    double   phi =-(hcalCells[i].phiOffset())*deg + 0.5*dphi;
    LogDebug("HCalGeom") << "HcalDDDGeometryLoader: Subdet " << subdet
			 << " eta " << etaRing << " depth " << depthBin
			 << " modules " << hcalCells[i].nPhiModule() << " "
			 << phiInc << " phi " << phi/deg << " " << dphi/deg;
    for (int k = 0; k < hcalCells[i].nPhiBins(); k++) {
      LogDebug("HCalGeom") << "HcalDDDGeometryLoader::fill Cell " << i
			   << " eta " << etaRing << " phi " << iphi << "("
			   << phi/deg << ", " << dphi/deg << ") depth "
			   << depthBin;
      HcalDetId id(subdet, etaRing, iphi, depthBin);
      hcalIds.push_back(id);
      const CaloCellGeometry * geometry = makeCell(id,hcalCells[i],phi,dphi,geom);
      geom->addCell(id, geometry);
      if (hcalCells[i].nHalves() > 1) {
	LogDebug("HCalGeom") << "HcalDDDGeometryLoader::fill Cell " << i
			     << " eta " << -etaRing << " phi " << iphi << " ("
			     << phi/deg << ", " << dphi/deg << ") depth "
			     << depthBin;
	HcalDetId id(subdet, -etaRing, iphi, depthBin);
	hcalIds.push_back(id);
	const CaloCellGeometry * geometry = makeCell(id,hcalCells[i],phi,dphi,geom);
	geom->addCell(id, geometry);
      }
      iphi += phiInc;
      phi  += dphi;
    }
  }
  
  edm::LogInfo("HCalGeom") << "Number of HCAL DetIds made for " << subdet
			   << " is " << hcalIds.size();
}

const CaloCellGeometry* 
HcalDDDGeometryLoader::makeCell(const HcalDetId& detId,
				HcalCellType::HcalCellType hcalCell,
				double phi, 
				double dphi,
				CaloSubdetectorGeometry* geom) const 
{
  // the two eta boundaries of the cell
  double          eta1   = hcalCell.etaMin();
  double          eta2   = hcalCell.etaMax();
  HcalSubdetector subdet = detId.subdet();
  double          eta    = 0.5*(eta1+eta2) * detId.zside();
  double          deta   = (eta2-eta1);
  double          theta  = 2.0*atan(exp(-eta));

  // barrel vs forward
  bool isBarrel = (subdet == HcalBarrel || subdet == HcalOuter);

  double          z, r, thickness;

  if (isBarrel) {
    r          = hcalCell.depthMin();
    thickness  = hcalCell.depthMax() - r;
    LogDebug("HCalGeom") << "HcalDDDGeometryLoader::makeCell SubDet " << subdet
			 << " eta = " << eta << " theta = " << theta
			 << " r = " << r << " thickness = " << thickness;
    z          = r * sinh(eta); // sinh(eta) == 1/tan(theta)
    thickness *= cosh(eta);     // cosh(eta) == 1/sin(theta)
  } else {
    z          = hcalCell.depthMin();
    thickness  = hcalCell.depthMax() - z;
    LogDebug("HCalGeom") << "HcalDDDGeometryLoader::makeCell SubDet " << subdet
			 << " eta = " << eta << " theta = " << theta
			 << " z = " << z << " thickness = " << thickness;
    z         *= detId.zside(); // get the sign right.
    r          = z * tan(theta);
  }

  double x = r * cos(phi);
  double y = r * sin(phi);
  GlobalPoint point(x,y,z);

  LogDebug("HCalGeom") << "HcalDDDGeometryLoader::makeCell for " << detId
		       << " Point " << point << " deta = " << deta 
		       << " dphi = " << dphi << " thickness = " << thickness
		       << " isBarrel = " << isBarrel;

  
  if (subdet==HcalForward) 
  {
     std::vector<double> hf ;
     hf.reserve(3) ;
     hf.push_back( deta ) ;
     hf.push_back( dphi ) ;
     hf.push_back( thickness ) ;
     return new calogeom::IdealZPrism( 
	point, 
	geom->cornersMgr(),
	CaloCellGeometry::getParmPtr( hf, 
				      geom->parMgr(), 
				      geom->parVecVec() ) );
  } 
  else 
  { 
     const double sign ( isBarrel ? 1 : -1 ) ;
     std::vector<double> hh ;
     hh.reserve(3) ;
     hh.push_back( deta ) ;
     hh.push_back( dphi ) ;
     hh.push_back( sign*thickness ) ;
     return new calogeom::IdealObliquePrism(
	point,
	geom->cornersMgr(),
	CaloCellGeometry::getParmPtr( hh, 
				      geom->parMgr(), 
				      geom->parVecVec() ) );
 }
}
