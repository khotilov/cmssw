#include "SimG4CMS/Muon81/interface/MuonSlaveSD.h"
#include "Geometry/MuonNumbering/interface/MuonSubDetector.h"

#include "SimG4Core/Application81/interface/SimTrackManager.h"

#include <iostream>

//#define DEBUG


MuonSlaveSD::MuonSlaveSD(MuonSubDetector* d,
			 const SimTrackManager* manager): 
  TrackingSlaveSD(d->name() ), m_trackManager(manager)
{
  detector=d;

}

MuonSlaveSD::~MuonSlaveSD() { 
}

void MuonSlaveSD::clearHits()
{
#ifdef DEBUG
    std::cout << " MuonSlaveSD::clearHits "<< detector->name() << std::endl;
#endif
    hits_.clear();
}

bool MuonSlaveSD::format()
{
#ifdef DEBUG
  std::cout << " MuonSlaveSD "<<detector->name()<<" formatting " << hits_.size() <<" hits."<< std::endl;
#endif
  if (detector->isBarrel()) {
    sort(hits_.begin(),hits_.end(), FormatBarrelHits());
  } else if (detector->isEndcap()) {
    sort(hits_.begin(),hits_.end(), FormatEndcapHits());
  } else if (detector->isRpc()) {
    sort(hits_.begin(),hits_.end(), FormatRpcHits());
  } 
  
  return true;
}


bool FormatBarrelHits::operator() (const PSimHit & a, const PSimHit & b)
{
  return (sortId(a)<sortId(b));
}

int FormatBarrelHits::sortId(const PSimHit & a)  const 
{
  return a.detUnitId();
}

bool FormatEndcapHits::operator() (const PSimHit & a, const PSimHit & b)
{
  return (sortId(a)<sortId(b));
}

int FormatEndcapHits::sortId(const PSimHit & a)  const 
{
  return a.detUnitId();
}

bool FormatRpcHits::operator() (const PSimHit & a, const PSimHit & b)
{
  return (sortId(a)<sortId(b));
}

int FormatRpcHits::sortId(const PSimHit & a)  const 
{
  return a.detUnitId();
}

