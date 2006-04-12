#include <L1Trigger/CSCTrackFinder/interface/CSCSectorReceiverLUT.h>
#include <L1Trigger/CSCCommonTrigger/interface/CSCTriggerGeometry.h>
#include <L1Trigger/CSCCommonTrigger/interface/CSCTriggerGeomManager.h>
#include <L1Trigger/CSCCommonTrigger/interface/CSCPatternLUT.h>
#include <L1Trigger/CSCCommonTrigger/interface/CSCFrontRearLUT.h>
#include <L1Trigger/CSCCommonTrigger/interface/CSCBitWidths.h>
#include <L1Trigger/CSCCommonTrigger/interface/CSCConstants.h>

#include <Geometry/CSCGeometry/interface/CSCLayerGeometry.h>
#include <Geometry/Vector/interface/LocalPoint.h>
#include <Geometry/Vector/interface/GlobalPoint.h>

#include <DataFormats/MuonDetId/interface/CSCTriggerNumbering.h>

#include <FWCore/MessageLogger/interface/MessageLogger.h>

#include <fstream>

lclphidat* CSCSectorReceiverLUT::me_lcl_phi = NULL;
bool CSCSectorReceiverLUT::me_lcl_phi_loaded = false;

CSCSectorReceiverLUT::CSCSectorReceiverLUT(int endcap, int sector, int subsector, int station,
					   const edm::ParameterSet & pset):_endcap(endcap),_sector(sector),
									   _subsector(subsector),
									   _station(station)
{
  LUTsFromFile = pset.getUntrackedParameter<bool>("ReadLUTs",false);
  isBinary = pset.getUntrackedParameter<bool>("Binary",false);
  lut_path = pset.getUntrackedParameter<std::string>("LUTPath","./");
  me_global_eta = NULL;
  me_global_phi = NULL;
  if(LUTsFromFile) readLUTsFromFile();
}

CSCSectorReceiverLUT::CSCSectorReceiverLUT(const CSCSectorReceiverLUT& lut):_endcap(lut._endcap),
									    _sector(lut._sector),
									    _subsector(lut._subsector),
									    _station(lut._station),
									    lut_path(lut.lut_path),
									    LUTsFromFile(lut.LUTsFromFile),
									    isBinary(lut.isBinary)
{
  if(lut.me_global_phi)
    {
      me_global_phi = new gblphidat[1<<CSCBitWidths::kGlobalPhiAddressWidth];
      memcpy(me_global_phi, lut.me_global_phi, (1<<CSCBitWidths::kGlobalPhiAddressWidth)*sizeof(gblphidat));
    }

  if(lut.me_global_eta)
    {
      me_global_eta = new gbletadat[1<<CSCBitWidths::kGlobalEtaAddressWidth];
      memcpy(me_global_eta, lut.me_global_eta, (1<<CSCBitWidths::kGlobalEtaAddressWidth)*sizeof(gbletadat));
    }
}
										       
CSCSectorReceiverLUT& CSCSectorReceiverLUT::operator=(const CSCSectorReceiverLUT& lut)
{
  if(this != &lut)
    {
      _endcap = lut._endcap;
      _sector = lut._sector;
      _subsector = lut._subsector;
      _station = lut._station;
      lut_path = lut.lut_path;
      LUTsFromFile = lut.LUTsFromFile;
      isBinary = lut.isBinary;

      
      if(lut.me_global_eta)
	{
	  me_global_eta = new gbletadat[1<<CSCBitWidths::kGlobalEtaAddressWidth];
	  memcpy(me_global_eta, lut.me_global_eta, (1<<CSCBitWidths::kGlobalEtaAddressWidth)*sizeof(gbletadat));
	}
      else me_global_eta = NULL;
    }
  return *this;
}

CSCSectorReceiverLUT::~CSCSectorReceiverLUT()
{
  if(me_lcl_phi_loaded)
    {
      delete me_lcl_phi;
      me_lcl_phi = NULL;
      me_lcl_phi_loaded = false;
    }
  if(me_global_eta)
    {
      delete me_global_eta;
      me_global_eta = NULL;
    }
}

lclphidat CSCSectorReceiverLUT::calcLocalPhi(const lclphiadd& theadd) const
{
  lclphidat data;

  int maxPhiL = 1<<CSCBitWidths::kLocalPhiDataBitWidth;
  double binPhiL = static_cast<double>(maxPhiL)/(2.*CSCConstants::MAX_NUM_STRIPS);

  memset(&data,0,sizeof(lclphidat));
  double patternOffset = CSCPatternLUT::getPosition(theadd.clct_pattern);
  
  if(theadd.strip < 2*CSCConstants::MAX_NUM_STRIPS)
    if(theadd.pattern_type == 1) // if halfstrip
      data.phi_local = static_cast<unsigned>((0.5 + theadd.strip + patternOffset)*binPhiL);
    else // if distrip
      data.phi_local = static_cast<unsigned>((2 + theadd.strip + 4.*patternOffset)*binPhiL);
  else // set out of bounds values
    if(theadd.pattern_type == 1)
      data.phi_local = static_cast<unsigned>((0.5 + (2*CSCConstants::MAX_NUM_STRIPS-1) + patternOffset)*binPhiL);
    else
      data.phi_local = static_cast<unsigned>((2 + (2*CSCConstants::MAX_NUM_STRIPS-1) + 4.*patternOffset)*binPhiL);
  
  /// Local Phi Bend is always zero. Until we start using it.
  data.phi_bend_local = 0;

  return data; //return LUT result
}


void CSCSectorReceiverLUT::fillLocalPhiLUT()
{ 
  // read data in from a file... Add this later.
}

lclphidat CSCSectorReceiverLUT::localPhi(int strip, int pattern, int quality, int lr) const
{
  lclphiadd theadd;

  theadd.strip = strip;
  theadd.clct_pattern = pattern;
  theadd.quality = quality;
  theadd.lr = lr;
  theadd.spare = 0;

  return localPhi(theadd);
}

lclphidat CSCSectorReceiverLUT::localPhi(unsigned address) const
{
  lclphidat result;
  lclphiadd theadd(address);

  if(LUTsFromFile) result = me_lcl_phi[address];
  else result = calcLocalPhi(theadd);

  return result;
}

lclphidat CSCSectorReceiverLUT::localPhi(lclphiadd address) const
{
  lclphidat result;
  
  if(LUTsFromFile) result = me_lcl_phi[address.toint()];
  else result = calcLocalPhi(address);
  
  return result;
}

double CSCSectorReceiverLUT::getGlobalPhiValue(const CSCLayer* thelayer, const unsigned& strip, const unsigned& wire_group) const
{
  double result = 0.0;
  CSCLayerGeometry* thegeom;
  LocalPoint lp;
  GlobalPoint gp;

  try
    {
      thegeom = const_cast<CSCLayerGeometry*>(thelayer->geometry());
      lp = thegeom->stripWireGroupIntersection(strip+1, wire_group+1);
      gp = thelayer->surface().toGlobal(lp);
      result = gp.phi();
      if (result < 0.) result += 2*M_PI;
    }
  catch(edm::Exception& e)
    {
      LogDebug("CSCSectorReceiverLUT|getGlobalPhiValue") << e.what();
    }

  return result;
}

gblphidat CSCSectorReceiverLUT::calcGlobalPhiME(const gblphiadd& address) const
{
  gblphidat result(0);
  CSCTriggerGeomManager* thegeom = CSCTriggerGeometry::get();
  CSCChamber* thechamber = NULL;
  CSCLayer* thelayer = NULL;
  CSCLayerGeometry* layergeom = NULL;
  unsigned cscid = address.cscid;
  unsigned wire_group = address.wire_group;
  unsigned local_phi = address.phi_local;
  const double sectorOffset = (CSCConstants::SECTOR1_CENT_RAD-CSCConstants::SECTOR_RAD/2.) + (_sector-1)*M_PI/3.;
  
  //Number of global phi units per radian.
  int    maxPhiG = 1<<CSCBitWidths::kGlobalPhiDataBitWidth;
  double binPhiG = static_cast<double>(maxPhiG)/CSCConstants::SECTOR_RAD;

  // We will use chamberWidth to convert the local phi into radians.
  int maxPhiL = 1<<CSCBitWidths::kLocalPhiDataBitWidth;
  const double binPhiL = static_cast<double>(maxPhiL)/(2.*CSCConstants::MAX_NUM_STRIPS);

  if(cscid < CSCTriggerNumbering::minTriggerCscId() || cscid > CSCTriggerNumbering::maxTriggerCscId())
    {
      LogDebug("CSCSectorReceiverLUT|getGlobalPhiValue") << " warning: cscId " << cscid
							 << " is out of bounds (1-" << CSCTriggerNumbering::maxTriggerCscId();
      cscid = CSCTriggerNumbering::minTriggerCscId();
    }
  if(wire_group >= 1<<5)
    {
      LogDebug("CSCSectorReceiverLUT|getGlobalPhiValue") << "warning: wire_group" << wire_group
							 << " is out of bounds (1-" << ((1<<5)-1);
      wire_group = (1<<5) - 1;
    }
  if(local_phi >= 1<<CSCBitWidths::kLocalPhiDataBitWidth)
    {
      LogDebug("CSCSectorReceiverLUT|getGlobalPhiValue") << "warning: local_phi" << local_phi
                                                         << " is out of bounds (1-" << ((1<<CSCBitWidths::kLocalPhiDataBitWidth)-1);
      local_phi = (1<<CSCBitWidths::kLocalPhiDataBitWidth) - 1;
    }
  
  try
    {
      thechamber = thegeom->chamber(_endcap,_station,_sector,_subsector,cscid);
      if(thechamber)
	{
	  layergeom = const_cast<CSCLayerGeometry*>(thechamber->layer(3)->geometry());
	  thelayer = const_cast<CSCLayer*>(thechamber->layer(3));
	  const int nStrips = layergeom->numberOfStrips();
	  // PhiL is the strip number converted into some units between 0 and
	  // 1023.  When we did the conversion in fillLocalPhiTable(), we did
	  // not know for which chamber we do it (and, therefore, how many strips
	  // it has), and always used the maximum possible number of strips
	  // per chamber, MAX_NUM_STRIPS=80.  Now, since we know the chamber id
	  // and how many strips the chamber has, we can re-adjust the scale.
	  const double scale = static_cast<double>(CSCConstants::MAX_NUM_STRIPS)/nStrips;

	  int strip = 0, halfstrip = 0;

          halfstrip = static_cast<int>(local_phi/binPhiL);
          strip     = halfstrip/2;
	  
	  // Find the phi width of the chamber and the position of its "left"
	  // (lower phi) edge (both in radians).
	  // Phi positions of the centers of the first and of the last strips
	  // in the chamber.
	  const double phi_f = getGlobalPhiValue(thelayer, 0, wire_group);
	  const double phi_l = getGlobalPhiValue(thelayer, nStrips - 1, wire_group);
	  // Phi widths of the half-strips at both ends of the chamber;
	  // surprisingly, they are not the same.
	  const double hsWidth_f = fabs(getGlobalPhiValue(thelayer, 1, wire_group) - phi_f)/2.;
	  const double hsWidth_l = fabs(phi_l - getGlobalPhiValue(thelayer, nStrips-2, wire_group))/2.;
	  // The "natural" match between the strips and phi values -- when
	  // a larger strip number corresponds to a larger phi value, i.e. strips
	  // are counted clockwise if we look at them from the inside of the
	  // detector -- is reversed for some stations.  At the moment, these
	  // are stations 3 and 4 of the 1st endcap, and stations 1 and 2 of
	  // the 2nd endcap.  Instead of using
	  // if ((theEndcap == 1 && theStation <= 2) ||
	  // (theEndcap == 2 && theStation >= 3)),
	  // we get the order from the phi values of the first and the last strip
	  // in a chamber, just in case the counting scheme changes in the future.
	  // Once we know how the strips are counted, we can go from the middle
	  // of the strips to their outer edges.
	  bool   clockwiseOrder;
	  double leftEdge, rightEdge;
	  if (fabs(phi_f - phi_l) < M_PI) 
	    {
	      if (phi_f < phi_l) clockwiseOrder = true;
	      else clockwiseOrder = false;
	    }
	  else 
	    { // the chamber crosses the phi = pi boundary
	      if (phi_f < phi_l) clockwiseOrder = false;
	      else clockwiseOrder = true;
	    }
	  if (clockwiseOrder) 
	    {
	      leftEdge  = phi_f - hsWidth_f;
	      rightEdge = phi_l + hsWidth_l;
	    }
	  else 
	    {
	      leftEdge  = phi_l - hsWidth_l;
	      rightEdge = phi_f + hsWidth_f;
	    }
	  if (fabs(phi_f - phi_l) >= M_PI) {rightEdge += 2.*M_PI;}
	  double chamberWidth = (rightEdge - leftEdge);

	  // Chamber offset, relative to the edge of the sector.
	  double chamberOffset = leftEdge - sectorOffset;
	  if (chamberOffset < -M_PI) chamberOffset += 2*M_PI;

	  double temp_phi = 0.0, strip_phi = 0.0, delta_phi = 0.0;
	  double distFromHalfStripCenter = 0.0, halfstripWidth = 0.0;
	  
	  if (strip < nStrips) 
	    {
	      // Approximate distance from the center of the half-strip to the center
	      // of this phil bin, in units of half-strip width.
	      distFromHalfStripCenter = (local_phi+0.5)/binPhiL - halfstrip - 0.5;
	      // Half-strip width (in rad), calculated as the half-distance between
	      // the adjacent strips.  Since in the current ORCA implementation
	      // the half-strip width changes from strip to strip, base the choice
	      // of the adjacent strip on the half-strip number.
	      if ((halfstrip%2 == 0 && halfstrip != 0) || halfstrip == 2*nStrips-1) {
		halfstripWidth = 
		  fabs(getGlobalPhiValue(thelayer, strip, wire_group) - getGlobalPhiValue(thelayer, strip - 1, wire_group)) / 2.;
	      }
	      else 
		{
		  halfstripWidth = 
		    fabs(getGlobalPhiValue(thelayer, strip, wire_group) - getGlobalPhiValue(thelayer, strip+1, wire_group)) / 2.;
		}
	      // Correction for the strips crossing the 180 degree boundary.
	      if (halfstripWidth > M_PI/2.) halfstripWidth = M_PI - halfstripWidth;
	      // Phi at the center of the strip.
	      strip_phi = getGlobalPhiValue(thelayer, strip, wire_group);
	      // Distance between the center of the strip and the phil position.
	      delta_phi = halfstripWidth*(((halfstrip%2)-0.5)+distFromHalfStripCenter);
	      if (clockwiseOrder)
		temp_phi = strip_phi + delta_phi;
	      else
		temp_phi = strip_phi - delta_phi;
	    }
	  else 
	    {
	      // PhiL values that do not have corresponding strips (the chamber
	      // has less than 80 strips assumed in fillLocalPhi).  It does not
	      // really matter what we do with these values; at the moment, just
	      // set them to the phis of the edges of the chamber.
	      if (clockwiseOrder) temp_phi = rightEdge;
	      else temp_phi = leftEdge;
	    }

	  // Finally, subtract the sector offset and convert to the scale of
	  // the global phi.
	  temp_phi -= sectorOffset;
	  if (temp_phi < 0.) temp_phi += 2.*M_PI;
	  temp_phi *= binPhiG;

	  if (temp_phi < 0) 
	    {
	      result.global_phi = 0;
	    }
	  else if (temp_phi >= maxPhiG) 
	    {
	      result.global_phi = maxPhiG - 1;
	    }
	  else 
	    {
	     result.global_phi = static_cast<unsigned short>(temp_phi);
	    }
	}
    }
  catch(edm::Exception& e)
    {
      edm::LogError("CSCSectorReceiverLUT|getGlobalPhiValue") << e.what();
    }

  return result;
}

gblphidat CSCSectorReceiverLUT::globalPhiME(int phi_local, int wire_group, int cscid) const
{
  gblphidat result;
  gblphiadd theadd;
  theadd.phi_local = phi_local;
  theadd.wire_group = ((1<<5)-1)&(wire_group >> 2); // want 2-7 of wg
  theadd.cscid = cscid;
  
  if(LUTsFromFile) result = me_global_phi[theadd.toint()];
  else result = calcGlobalPhiME(theadd);

  return result;
}

gblphidat CSCSectorReceiverLUT::globalPhiME(unsigned address) const
{
  gblphidat result;

  if(LUTsFromFile) result = me_global_phi[address];
  else result = calcGlobalPhiME(gblphiadd(address));
  
  return result;
}

gblphidat CSCSectorReceiverLUT::globalPhiME(gblphiadd address) const
{
  gblphidat result;

  if(LUTsFromFile) result = me_global_phi[address.toint()];
  else result = calcGlobalPhiME(address);

  return result;
}

double CSCSectorReceiverLUT::getGlobalEtaValue(const unsigned& thecscid, const unsigned& thewire_group, const unsigned& thephi_local) const
{
  double result = 0.0;
  unsigned wire_group = thewire_group;
  int cscid = thecscid;
  unsigned phi_local = thephi_local;

  if(cscid < CSCTriggerNumbering::minTriggerCscId() || cscid > CSCTriggerNumbering::maxTriggerCscId())
    {
      LogDebug("CSCSectorReceiverLUT|getEtaValue") << " warning: cscId " << cscid
						   << " is out of bounds (1-" << CSCTriggerNumbering::maxTriggerCscId();
      cscid = CSCTriggerNumbering::maxTriggerCscId();
    }

  CSCTriggerGeomManager* thegeom = CSCTriggerGeometry::get();
  CSCLayerGeometry* layerGeom = NULL;
  const unsigned numBins = 1 << 2; // 4 local phi bins
  
  if(phi_local > numBins - 1)
    {
      LogDebug("CSCSectorReceiverLUT|getEtaValue") << "warning: phiL " << phi_local
						   << " is out of bounds (0-" << numBins - 1;
      phi_local = numBins - 1;
    }
  try 
    {    
      const CSCChamber* thechamber = thegeom->chamber(_endcap,_station,_sector,_subsector,cscid);     
      if(thechamber) 
	{
	  if(_station != 1 && CSCTriggerNumbering::ringFromTriggerLabels(_station, cscid) != 1)
	    {
	      layerGeom = const_cast<CSCLayerGeometry*>(thechamber->layer(3)->geometry());
	      const unsigned nWireGroups = layerGeom->numberOfWireGroups();

	      if(wire_group > nWireGroups)
		{
		  LogDebug("CSCSectorReceiverLUT|getEtaValue") << "warning: wireGroup "
                                                               << wire_group << " is out of bounds (0-"
                                                               << nWireGroups;
                  wire_group = nWireGroups - 1;
		}

	      result = thechamber->layer(3)->centerOfWireGroup(wire_group+1).eta();
	    }
	  else
	    {
	      layerGeom = const_cast<CSCLayerGeometry*>(thechamber->layer(3)->geometry());
	      
	      const unsigned nStrips = layerGeom->numberOfStrips();
	      const unsigned nWireGroups = layerGeom->numberOfWireGroups();
	      const unsigned nStripsPerBin = CSCConstants::MAX_NUM_STRIPS/numBins;
	      
	      
	      if(wire_group > nWireGroups) // apply maximum limit
		{
		  LogDebug("CSCSectorReceiverLUT|getEtaValue") << "warning: wireGroup "
							       << wire_group << " is out of bounds (0-"
							       << nWireGroups;
		  wire_group = nWireGroups - 1;
		}
	      
	      /**
	       * Calculate Eta correction
	       */
	      
	      // Check that no strips will be left out.
	      
	      if (nStrips%numBins != 0 || CSCConstants::MAX_NUM_STRIPS%numBins != 0)
		LogDebug("CSCSectorReceiverLUT|EtaCorrectionWarning") << "calcEtaCorrection warning: number of strips "
								      << nStrips << " (" << CSCConstants::MAX_NUM_STRIPS
								      << ") is not divisible by numBins " << numBins
								      << " Station " << _station << " sector " << _sector
								      << " subsector " << _subsector << " cscid " << cscid;
	      
	      unsigned    maxStripPrevBin = 0, maxStripThisBin = 0;
	      unsigned    correctionStrip;
	      LocalPoint  lPoint;
	      GlobalPoint gPoint;
	      // Bins phi_local and find the the middle strip for each bin.
	      maxStripPrevBin = nStripsPerBin * (phi_local);
	      maxStripThisBin = nStripsPerBin * (phi_local+1);
	      if (maxStripThisBin <= nStrips) 
		{
		  correctionStrip = nStripsPerBin/2 * (2*phi_local+1);
		  maxStripPrevBin = maxStripThisBin;
		}
	      else 
		{
		  // If the actual number of strips in the chamber is smaller than
		  // the number of strips corresponding to the right edge of this phi
		  // local bin, we take the middle strip between number of strips
		  // at the left edge of the bin and the actual number of strips.
		  correctionStrip = (nStrips+maxStripPrevBin)/2;
		}
	      
	      lPoint = layerGeom->stripWireGroupIntersection(correctionStrip, wire_group+1);
	      if(thechamber) gPoint = thechamber->layer(3)->surface().toGlobal(lPoint);
	      
	      // end calc of eta correction.
	      result = gPoint.eta();
	    }
	}
    }
  catch (cms::Exception &e)
    {
      LogDebug("CSCSectorReceiver:OutofBoundInput") << e.what();
    }
  
  return std::fabs(result);
}


gbletadat CSCSectorReceiverLUT::calcGlobalEtaME(const gbletaadd& address) const
{
  gbletadat result;
  double float_eta = getGlobalEtaValue(address.cscid, address.wire_group, address.phi_local);
  unsigned int_eta = 0;
  unsigned bend_global = 0; // not filled yet... will change when it is.
  const double etaPerBin = (CSCConstants::maxEta - CSCConstants::minEta)/CSCConstants::etaBins;
  const unsigned me12EtaCut = 56;
    
  if ((float_eta < CSCConstants::minEta) || (float_eta >= CSCConstants::maxEta)) 
    {
      
      LogDebug("CSCSectorReceiverLUT:OutOfBounds") << "L1MuCSCSectorReceiverLUT warning: float_eta = " << float_eta
						   << " minEta = " << CSCConstants::minEta << " maxEta = " << CSCConstants::maxEta
						   << "   station " << _station << " sector " << _sector
						   << " chamber "   << address.cscid << " wire group " << address.wire_group;
      
      if (float_eta < CSCConstants::minEta) 
	result.global_eta = 0;
      else if (float_eta >= CSCConstants::maxEta) 
	result.global_eta = CSCConstants::etaBins - 1;
    }
  else
    {  
      float_eta -= CSCConstants::minEta;
      float_eta = float_eta/etaPerBin;
      int_eta = static_cast<unsigned>(float_eta);
      /* Commented until I find out its use.
      // Fine-tune eta boundary between DT and CSC.
      if ((intEta == L1MuCSCSetup::CscEtaStart() && (L1MuCSCSetup::CscEtaStartCorr() > 0.) ) ||
	  (intEta == L1MuCSCSetup::CscEtaStart() - 1 && (L1MuCSCSetup::CscEtaStartCorr() < 0.) ) ) {
	bitEta = (thisEta-minEta-L1MuCSCSetup::CscEtaStartCorr())/EtaPerBin;
	intEta = static_cast<int>(bitEta);
      }
      */
      if (_station == 1 && address.cscid >= static_cast<unsigned>(CSCTriggerNumbering::minTriggerCscId()) 
	  && address.cscid <= static_cast<unsigned>(CSCTriggerNumbering::maxTriggerCscId()) )
	{
	  unsigned ring = CSCTriggerNumbering::ringFromTriggerLabels(_station, address.cscid);
	  
	  if      (ring == 1 && int_eta <  me12EtaCut) {int_eta = me12EtaCut;}
	  else if (ring == 2 && int_eta >= me12EtaCut) {int_eta = me12EtaCut-1;}
	}
      result.global_eta = int_eta;
    }
  result.global_bend = bend_global;

  return result;
}

gbletadat CSCSectorReceiverLUT::globalEtaME(int tphi_bend, int tphi_local, int twire_group, int tcscid) const
{
  gbletadat result;
  gbletaadd theadd;

  theadd.phi_bend = tphi_bend;
  theadd.phi_local = (tphi_local>>(CSCBitWidths::kLocalPhiDataBitWidth - 2)) & 0x3; // want 2 msb of local phi
  theadd.wire_group = twire_group;
  theadd.cscid = tcscid;
  
  if(LUTsFromFile) result = me_global_eta[theadd.toint()];
  else result = calcGlobalEtaME(theadd);

  //  if(address.wire_group == 0 && address.phi_local==0) std::cout << result.global_eta << std::endl;

  return result;
}

gbletadat CSCSectorReceiverLUT::globalEtaME(unsigned address) const
{
  gbletadat result;
  gbletaadd theadd(address);

  if(LUTsFromFile) result = me_global_eta[address];
  else result = calcGlobalEtaME(theadd);
  return result;
}

gbletadat CSCSectorReceiverLUT::globalEtaME(gbletaadd address) const
{
  gbletadat result;

  if(LUTsFromFile) result = me_global_eta[address.toint()];
  else result = calcGlobalEtaME(address);
  return result;
}

std::string CSCSectorReceiverLUT::encodeFileIndex() const {
  std::string fileName = "";
  if (_station == 1) {
    if (_subsector == 1) fileName += "1a";
    if (_subsector == 2) fileName += "1b";
  }
  else if (_station == 2) fileName += "2";
  else if (_station == 3) fileName += "3";
  else if (_station == 4) fileName += "4";
  fileName += "End";
  if (_endcap == 1) fileName += "1";
  else                fileName += "2";
  fileName += "Sec";
  if      (_sector == 1) fileName += "1";
  else if (_sector == 2) fileName += "2";
  else if (_sector == 3) fileName += "3";
  else if (_sector == 4) fileName += "4";
  else if (_sector == 5) fileName += "5";
  else if (_sector == 6) fileName += "6";
  fileName += "LUT";
  return fileName;
}

void CSCSectorReceiverLUT::readLUTsFromFile()
{
  if(!me_lcl_phi_loaded)
    {
      me_lcl_phi = new lclphidat[1<<CSCBitWidths::kLocalPhiAddressWidth];
      memset(me_lcl_phi, 0, (1<<CSCBitWidths::kLocalPhiAddressWidth)*sizeof(short));
      std::string fName("LocalPhiLUT");
      std::ifstream LocalPhiLUT;

      edm::LogInfo("CSCSectorReceiverLUT|loadLUT") << "Loading SR LUT: " << fName; 

      fName += ((isBinary) ? ".bin" : ".dat");

      if(isBinary)
	{
	  LocalPhiLUT.open((lut_path+fName).c_str(),std::ios::binary);
          LocalPhiLUT.seekg(0,std::ios::end);
          int length = LocalPhiLUT.tellg();
          if(length == (1<<CSCBitWidths::kLocalPhiAddressWidth)*sizeof(short))
	    {
	      LocalPhiLUT.seekg(0,std::ios::beg);
	      LocalPhiLUT.read(reinterpret_cast<char*>(me_lcl_phi),length);
	      LocalPhiLUT.close();
	    }
	  else
	    edm::LogError("CSCSectorReceiverLUT|loadLUT") << "File "<< fName << " is incorrect size!";
	  LocalPhiLUT.close();
	}
      else
        {
          LocalPhiLUT.open((lut_path+fName).c_str());
      	  unsigned i = 0;
	  unsigned short temp = 0;
          while(!LocalPhiLUT.eof() && i < 1<<CSCBitWidths::kLocalPhiAddressWidth)
	    {
	      LocalPhiLUT >> temp;
	      me_lcl_phi[i++] = (*reinterpret_cast<lclphidat*>(&temp));
	    }
	  LocalPhiLUT.close();
	}
    }
  if(!me_global_phi)
    {
      me_global_phi = new gblphidat[1<<CSCBitWidths::kGlobalPhiAddressWidth];
      memset(me_global_phi, 0, (1<<CSCBitWidths::kGlobalPhiAddressWidth)*sizeof(short));
      std::string fName = std::string("GlobalPhiME") + encodeFileIndex();
      std::ifstream GlobalPhiLUT;

      edm::LogInfo("CSCSectorReceiverLUT|loadLUT") << "Loading SR LUT: " << fName;

      fName += ((isBinary) ? ".bin" : ".dat");

      std::cout << "OPENING " << lut_path + fName << std::endl;

      if(isBinary)
        {
          GlobalPhiLUT.open((lut_path + fName).c_str(),std::ios::binary);
          GlobalPhiLUT.seekg(0,std::ios::end);
          int length = GlobalPhiLUT.tellg();
          if(length == (1<<CSCBitWidths::kGlobalPhiAddressWidth)*sizeof(short))
            {
              GlobalPhiLUT.seekg(0,std::ios::beg);
              GlobalPhiLUT.read(reinterpret_cast<char*>(me_global_phi),length);
            }
          else
            edm::LogError("CSCSectorReceiverLUT|loadLUT") << "File "<< fName << " is incorrect size!";
          GlobalPhiLUT.close();
        }
      else
        {
          GlobalPhiLUT.open((lut_path + fName).c_str());
          unsigned short temp = 0;
          unsigned i = 0;
          while(!GlobalPhiLUT.eof() && i < 1<<CSCBitWidths::kGlobalPhiAddressWidth)
	    {
	      GlobalPhiLUT >> temp;
	      me_global_phi[i++] = (*reinterpret_cast<gblphidat*>(&temp));
	    }
          GlobalPhiLUT.close();
        }
    }
  if(!me_global_eta) 
    {
      me_global_eta = new gbletadat[1<<CSCBitWidths::kGlobalEtaAddressWidth];
      memset(me_global_eta, 0, (1<<CSCBitWidths::kGlobalEtaAddressWidth)*sizeof(short));
      std::string fName = std::string("GlobalEtaME") + encodeFileIndex();
      std::ifstream GlobalEtaLUT;
      
      edm::LogInfo("CSCSectorReceiverLUT|loadLUT") << "Loading SR LUT: " << fName;

      fName += ((isBinary) ? ".bin" : ".dat");

      if(isBinary)
	{
	  GlobalEtaLUT.open((lut_path + fName).c_str(),std::ios::binary);
	  GlobalEtaLUT.seekg(0,std::ios::end);
	  int length = GlobalEtaLUT.tellg();
	  if(length == (1<<CSCBitWidths::kGlobalEtaAddressWidth)*sizeof(short))
	    {
	      GlobalEtaLUT.seekg(0,std::ios::beg);
	      GlobalEtaLUT.read(reinterpret_cast<char*>(me_global_eta),length);
	    }
	  else
	    edm::LogError("CSCSectorReceiverLUT|loadLUT") << "File "<< fName << " is incorrect size!";
	  GlobalEtaLUT.close();
	}
      else
	{
	  GlobalEtaLUT.open((lut_path + fName).c_str());
	  unsigned short temp = 0;
	  unsigned i = 0;
	  while(!GlobalEtaLUT.eof() && i < 1<<CSCBitWidths::kGlobalEtaAddressWidth)
	  {
	    GlobalEtaLUT >> temp;
	    me_global_eta[i++] = (*reinterpret_cast<gbletadat*>(&temp));
	  }
	  GlobalEtaLUT.close();
	}
    }
}

