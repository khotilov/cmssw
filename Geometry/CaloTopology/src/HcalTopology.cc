#include "Geometry/CaloTopology/interface/HcalTopology.h"
#include <cmath>
#include <iostream>


static const int IPHI_MAX=72;

HcalTopology::HcalTopology() :
  excludeHB_(false),
  excludeHE_(false),
  excludeHO_(false),
  excludeHF_(false),
  firstHBRing_(1),
  lastHBRing_(16),
  firstHERing_(16),
  lastHERing_(29),
  firstHFRing_(29),
  lastHFRing_(41),
  firstHORing_(1),
  lastHORing_(15),
  firstHEDoublePhiRing_(21),
  firstHFQuadPhiRing_(40),
  firstHETripleDepthRing_(27),
  singlePhiBins_(72),
  doublePhiBins_(36)
{
}


bool HcalTopology::valid(const HcalDetId& id) const {
  // check the raw rules
  bool ok=validRaw(id);

  ok=ok && !isExcluded(id);

  return ok;
}

bool HcalTopology::isExcluded(const HcalDetId& id) const {
  bool exed=false;
  // first, check the full detector exclusions...  (fast)
  switch (id.subdet()) {
  case(HcalBarrel): exed=excludeHB_; break;
  case(HcalEndcap): exed=excludeHE_; break;
  case(HcalOuter): exed=excludeHO_; break;
  case(HcalForward): exed=excludeHF_; break;
  default: exed=false;
  }
  // next, check the list (slower)
  if (!exed && !exclusionList_.empty()) {
    std::vector<HcalDetId>::const_iterator i=std::lower_bound(exclusionList_.begin(),exclusionList_.end(),id);
    if (i!=exclusionList_.end() && *i==id) exed=true;
  }
  return exed;
}

void HcalTopology::exclude(const HcalDetId& id) {
  std::vector<HcalDetId>::iterator i=std::lower_bound(exclusionList_.begin(),exclusionList_.end(),id);
  if (i==exclusionList_.end() || *i!=id) {
    exclusionList_.insert(i,id);
  }
}

void HcalTopology::excludeSubdetector(HcalSubdetector subdet) {
  switch (subdet) {
  case(HcalBarrel): excludeHB_=true; break;
  case(HcalEndcap): excludeHE_=true; break;
  case(HcalOuter): excludeHO_=true; break;
  case(HcalForward): excludeHF_=true; break;
  default: break;
  }
}



std::vector<DetId> HcalTopology::east(const DetId& id) const
{
  std::vector<DetId> vNeighborsDetId;
  HcalDetId neighbors[2];
  for (int i=0;i<decIEta(HcalDetId(id),neighbors);i++)
    vNeighborsDetId.push_back(DetId(neighbors[i].rawId()));
  return vNeighborsDetId;
}

std::vector<DetId> HcalTopology::west(const DetId& id) const
{
  std::vector<DetId> vNeighborsDetId;
  HcalDetId neighbors[2];
  for (int i=0;i<incIEta(HcalDetId(id),neighbors);i++)
    vNeighborsDetId.push_back(DetId(neighbors[i].rawId()));
  return  vNeighborsDetId;
}

std::vector<DetId> HcalTopology::north(const DetId& id) const
{
  std::vector<DetId> vNeighborsDetId;
  HcalDetId neighbor;
  if (incIPhi(HcalDetId(id),neighbor))
    vNeighborsDetId.push_back(DetId(neighbor.rawId()));
  return  vNeighborsDetId;
}

std::vector<DetId> HcalTopology::south(const DetId& id) const
{
  std::vector<DetId> vNeighborsDetId;
  HcalDetId neighbor;
  if (decIPhi(HcalDetId(id),neighbor))
    vNeighborsDetId.push_back(DetId(neighbor.rawId()));
  return  vNeighborsDetId;
}

std::vector<DetId> HcalTopology::up(const DetId& id) const
{
  HcalDetId neighbor = id;
  incrementDepth(neighbor);
  std::vector<DetId> vNeighborsDetId;
  if(incrementDepth(neighbor)) 
  {
    vNeighborsDetId.push_back(neighbor);
  }
  return  vNeighborsDetId;
}

std::vector<DetId> HcalTopology::down(const DetId& id) const
{
  std::cout << "HcalTopology::down() not yet implemented" << std::endl; 
  std::vector<DetId> vNeighborsDetId;
  return  vNeighborsDetId;
}

int HcalTopology::exclude(HcalSubdetector subdet, int ieta1, int ieta2, int iphi1, int iphi2, int depth1, int depth2) {

  bool exed=false;
  // first, check the full detector exclusions...  (fast)
  switch (subdet) {
  case(HcalBarrel): exed=excludeHB_; break;
  case(HcalEndcap): exed=excludeHE_; break;
  case(HcalOuter): exed=excludeHO_; break;
  case(HcalForward): exed=excludeHF_; break;
  default: exed=false;
  }
  if (exed) return 0; // if the whole detector is excluded...

  int ieta_l=std::min(ieta1,ieta2);
  int ieta_h=std::max(ieta1,ieta2);
  int iphi_l=std::min(iphi1,iphi2);
  int iphi_h=std::max(iphi1,iphi2);
  int depth_l=std::min(depth1,depth2);
  int depth_h=std::max(depth1,depth2);

  int n=0;
  for (int ieta=ieta_l; ieta<=ieta_h; ieta++) 
    for (int iphi=iphi_l; iphi<=iphi_h; iphi++) 
      for (int depth=depth_l; depth<=depth_h; depth++) {
	HcalDetId id(subdet,ieta,iphi,depth);
	if (validRaw(id)) { // use 'validRaw' to include check validity in "uncut" detector
	  exclude(id);  
	  n++;
	}
      }
  return n;
}

  /** Basic rules used to derive this code:
      
  HB has 72 towers in iphi.  Ieta 1-14 have depth=1, Ieta 15-16 have depth=1 or 2.

  HE ieta=16-20 have 72 towers in iphi
     ieta=21-29 have 36 towers in iphi
     ieta=16 is depth 3 only
     ieta=17 is depth 1 only
     ieta=18-26 & 29 have depth 1 and 2
     ieta=27-28 has depth 1-3

  HF ieta=29-39 have 36 in iphi
     ieta=40-41 have 18 in iphi
     all have two depths

  HO has 15 towers in ieta and 72 in iphi and depth = 4 (one value)
  */

  /** Is this a valid cell id? */
  bool HcalTopology::validRaw(const HcalDetId& id) const {
    bool ok=true;
    int ieta=id.ieta();
    int aieta=id.ietaAbs();
    int depth=id.depth();
    int iphi=id.iphi();

    if ((ieta==0 || iphi<=0 || iphi>IPHI_MAX) || aieta>41) return false; // outer limits
    
    if (ok) {
      HcalSubdetector subdet=id.subdet();
      if (subdet==HcalBarrel) {
	if (aieta>16 || depth>2 || (aieta<=14 && depth>1)) ok=false;	    
      } else if (subdet==HcalEndcap) {
	if (aieta<16 || aieta>29 ||
	    (aieta==16 && depth!=3) ||
	    (aieta==17 && depth!=1) ||
	    (((aieta>=18 && aieta<=26) || aieta==29) && depth>2) ||
	    ((aieta==27 || aieta==28) && depth>3) ||
	    (aieta>=21 && (iphi%2)==0)) ok=false;
      } else if (subdet==HcalOuter) {
	if (aieta>15 || iphi>IPHI_MAX || depth!=4) ok=false;
      } else if (subdet==HcalForward) {
	if (aieta<29 || aieta>41 ||
	    ((iphi%2)==0) ||
	    (depth>2) ||
	    (aieta>=40 && ((iphi+1)%4)==0)) ok=false;
      } else ok=false;
    }
    
    return ok;
  }

  
  bool HcalTopology::incIPhi(const HcalDetId& id, HcalDetId &neighbor) const {
    bool ok=valid(id);
    if (ok) {
      switch (id.subdet()) {
      case (HcalBarrel):
      case (HcalOuter):
	if (id.iphi()==IPHI_MAX) neighbor=HcalDetId(id.subdet(),id.ieta(),1,id.depth()); 
	else neighbor=HcalDetId(id.subdet(),id.ieta(),id.iphi()+1,id.depth()); 
	break;
      case (HcalEndcap):
	if (id.ietaAbs()>=21) {
	  if (id.iphi()==IPHI_MAX-1) neighbor=HcalDetId(HcalEndcap,id.ieta(),1,id.depth()); 
	  else neighbor=HcalDetId(HcalEndcap,id.ieta(),id.iphi()+2,id.depth()); 
	} else {
	  if (id.iphi()==IPHI_MAX) neighbor=HcalDetId(HcalEndcap,id.ieta(),1,id.depth()); 
	  else neighbor=HcalDetId(HcalEndcap,id.ieta(),id.iphi()+1,id.depth()); 
	}	
	break;
      case (HcalForward):
	if (id.ietaAbs()>=40) {
	  if (id.iphi()==IPHI_MAX-3) neighbor=HcalDetId(HcalEndcap,id.ieta(),1,id.depth()); 
	  else neighbor=HcalDetId(HcalEndcap,id.ieta(),id.iphi()+4,id.depth()); 
	} else {
	  if (id.iphi()==IPHI_MAX-1) neighbor=HcalDetId(HcalEndcap,id.ieta(),1,id.depth()); 
	  else neighbor=HcalDetId(HcalEndcap,id.ieta(),id.iphi()+2,id.depth()); 
	}
	break;
      default: ok=false;
      }
    } 
    return ok;
  }
  /** Get the neighbor (if present) of the given cell with lower iphi */
  bool HcalTopology::decIPhi(const HcalDetId& id, HcalDetId &neighbor) const {
    bool ok=valid(id);
    if (ok) {
      switch (id.subdet()) {
      case (HcalBarrel):
      case (HcalOuter):
	if (id.iphi()==1) neighbor=HcalDetId(id.subdet(),id.ieta(),IPHI_MAX,id.depth()); 
	else neighbor=HcalDetId(id.subdet(),id.ieta(),id.iphi()-1,id.depth()); 
	break;
      case (HcalEndcap):
	if (id.ietaAbs()>=21) {
	  if (id.iphi()==1) neighbor=HcalDetId(HcalEndcap,id.ieta(),IPHI_MAX-1,id.depth()); 
	  else neighbor=HcalDetId(HcalEndcap,id.ieta(),id.iphi()-2,id.depth()); 
	} else {
	  if (id.iphi()==1) neighbor=HcalDetId(HcalEndcap,id.ieta(),IPHI_MAX,id.depth()); 
	  else neighbor=HcalDetId(HcalEndcap,id.ieta(),id.iphi()-1,id.depth()); 
	}
	break;
      case (HcalForward):
	if (id.ietaAbs()>=40) {
	  if (id.iphi()==1) neighbor=HcalDetId(HcalEndcap,id.ieta(),IPHI_MAX-3,id.depth()); 
	  else neighbor=HcalDetId(HcalEndcap,id.ieta(),id.iphi()-4,id.depth()); 
	} else {
	  if (id.iphi()==1) neighbor=HcalDetId(HcalEndcap,id.ieta(),IPHI_MAX-1,id.depth()); 
	  else neighbor=HcalDetId(HcalEndcap,id.ieta(),id.iphi()-2,id.depth()); 
	}
	break;
      default: ok=false;
      }
    } 
    return ok;
  }

  int HcalTopology::incIEta(const HcalDetId& id, HcalDetId neighbors[2]) const {
    if (id.zside()==1) return incAIEta(id,neighbors);
    else return decAIEta(id,neighbors);
  }

  int HcalTopology::decIEta(const HcalDetId& id, HcalDetId neighbors[2]) const {
    if (id.zside()==1) return decAIEta(id,neighbors);
    else return incAIEta(id,neighbors);
  }

  /** Increasing in |ieta|, there is always at most one neighbor */
  int HcalTopology::incAIEta(const HcalDetId& id, HcalDetId neighbors[2]) const {
    int n=1;
    int aieta=id.ietaAbs();

    if (aieta==20 && (id.iphi()%2)==0) 
      neighbors[0]=HcalDetId(id.subdet(),(aieta+1)*id.zside(),id.iphi()-1,id.depth());
    else if (aieta==39 && ((id.iphi()+1)%4)==0) 
      neighbors[0]=HcalDetId(id.subdet(),(aieta+1)*id.zside(),id.iphi()-2,id.depth());
    else
      neighbors[0]=HcalDetId(id.subdet(),(aieta+1)*id.zside(),id.iphi(),id.depth());
    
    if (!valid(neighbors[0])) n=0;
    return n;
  }

  /** Decreasing in |ieta|, there are be two neighbors of 40 and 21*/
  int HcalTopology::decAIEta(const HcalDetId& id, HcalDetId neighbors[2]) const {
    int n=1;
    int aieta=id.ietaAbs();

    if (aieta==21) { 
      n=2;
      neighbors[0]=HcalDetId(id.subdet(),(aieta-1)*id.zside(),id.iphi(),id.depth());
      neighbors[1]=HcalDetId(id.subdet(),(aieta-1)*id.zside(),id.iphi()+1,id.depth());
    } else if (aieta==40) {
      n=2;
      neighbors[0]=HcalDetId(id.subdet(),(aieta-1)*id.zside(),id.iphi(),id.depth());
      neighbors[1]=HcalDetId(id.subdet(),(aieta-1)*id.zside(),id.iphi()+2,id.depth());
    } else
      neighbors[0]=HcalDetId(id.subdet(),(aieta-1)*id.zside(),id.iphi(),id.depth());
    
    if (!valid(neighbors[0]) && n==2) {
      if (!valid(neighbors[1])) n=0;
      else {
	n=1;
	neighbors[0]=neighbors[1];
      }
    }
    if (n==2 && !valid(neighbors[1])) n=1;

    return n;
  }


void HcalTopology::depthBinInformation(HcalSubdetector subdet, int etaRing,
                                       int & nDepthBins, int & startingBin) const {
  if(subdet == HcalBarrel) {
    if (etaRing<=14) {
      nDepthBins = 1;
      startingBin = 1;
    } else {
      nDepthBins = 2;
      startingBin = 1;
    }
  } else if(subdet == HcalEndcap) {
    if (etaRing==16) {
      nDepthBins = 1;
      startingBin = 3;
    } else if (etaRing==17) {
      nDepthBins = 1;
      startingBin = 1;
    } else if (etaRing==lastHERing()) {
      nDepthBins = 2;
      startingBin = 1;
    }
    else {
      nDepthBins = (etaRing >= firstHETripleDepthRing_) ? 3 : 2;
      startingBin = 1;
    }
  }

  else if(subdet == HcalForward) {
    nDepthBins = 2;
    startingBin = 1;
  }

  else if(subdet == HcalOuter) {
    nDepthBins = 1;
    startingBin = 4;
  }

  else {
    std::cerr << "Bad HCAL subdetector " << subdet << std::endl;
  }
}


bool HcalTopology::incrementDepth(HcalDetId & detId) const
{
  HcalSubdetector subdet = detId.subdet();
  int ieta = detId.ieta();
  int etaRing = detId.ietaAbs();
  int depth = detId.depth();
  int nDepthBins, startingBin;
  depthBinInformation(subdet, etaRing, nDepthBins, startingBin);

  // see if the new depth bin exists
  ++depth;
  if(depth > nDepthBins)
  {
    // handle on a case-by-case basis
    if(subdet == HcalBarrel && etaRing < lastHORing())
    {
      // HO
      subdet = HcalOuter;
      depth = 4;
    }
    else if(subdet == HcalBarrel && etaRing == lastHBRing())
    {
      // overlap
      subdet = HcalEndcap;
    }
    else if(subdet == HcalEndcap && etaRing ==  lastHERing()-1)
    {
      // guard ring HF29 is behind HE 28
      subdet = HcalForward;
      (ieta > 0) ? ++ieta : --ieta;
      depth = 1;
    }
    else if(subdet == HcalEndcap && etaRing ==  lastHERing())
    {
      // split cells go to bigger granularity.  Ring 29 -> 28
      (ieta > 0) ? --ieta : ++ieta;
    }
    else 
    {
      // no more chances
      detId = HcalDetId();
      return false;
    }
  }
  detId = HcalDetId(subdet, ieta, detId.iphi(), depth);
  assert(validRaw(detId));
  return true;
}


int HcalTopology::nPhiBins(int etaRing) const {
  int lastPhiBin=singlePhiBins_;
  if (etaRing>= firstHFQuadPhiRing_) lastPhiBin=doublePhiBins_/2;
  else if (etaRing>= firstHEDoublePhiRing_) lastPhiBin=doublePhiBins_;
  return lastPhiBin;
}



