// $Id: EcalMonitorPrescaler.h,v 1.5 2009/11/06 10:43:12 dellaric Exp $

/*!
  \file EcalMonitorPrescaler.h
  \brief Ecal specific Prescaler 
  \author G. Della Ricca
  \version $Revision: 1.5 $
  \date $Date: 2009/11/06 10:43:12 $
*/

#ifndef EcalMonitorPrescaler_H
#define EcalMonitorPrescaler_H

#include "FWCore/Framework/interface/EDFilter.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

class EcalMonitorPrescaler: public edm::EDFilter {

public:

explicit EcalMonitorPrescaler(edm::ParameterSet const& ps);
virtual ~EcalMonitorPrescaler();

virtual bool filter(edm::Event& e, edm::EventSetup const& c);
void endJob(void);

private:

edm::InputTag EcalRawDataCollection_;

int count_;

// accept one in n

int occupancyPrescaleFactor_;
int integrityPrescaleFactor_;
int statusflagsPrescaleFactor_;

int pedestalonlinePrescaleFactor_;

int laserPrescaleFactor_;
int ledPrescaleFactor_;
int pedestalPrescaleFactor_;
int testpulsePrescaleFactor_;

int pedestaloffsetPrescaleFactor_;

int triggertowerPrescaleFactor_;
int timingPrescaleFactor_;

int cosmicPrescaleFactor_;

int physicsPrescaleFactor_;

int clusterPrescaleFactor_;

};

#endif // EcalMonitorPrescaler_H
