#ifndef CSCTFPacker_h
#define CSCTFPacker_h

#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ParameterSet/interface/InputTag.h"

#include "DataFormats/CSCDigi/interface/CSCCorrelatedLCTDigi.h"
#include "DataFormats/CSCDigi/interface/CSCCorrelatedLCTDigiCollection.h"
#include "DataFormats/FEDRawData/interface/FEDRawData.h"
#include "DataFormats/FEDRawData/interface/FEDNumbering.h"
#include "DataFormats/FEDRawData/interface/FEDRawDataCollection.h"
#include "DataFormats/MuonDetId/interface/CSCDetId.h"

#include <string>

class CSCTFPacker : public edm::EDProducer {
private:
	edm::InputTag lctProducer, trackProducer;

	bool zeroSuppression;
	unsigned short nTBINs;
	unsigned short activeSectors;
	bool putBufferToEvent;

	bool swapME1strips;

	FILE *file;

	int m_minBX, m_maxBX, central_lct_bx, central_sp_bx;

public:
	virtual void produce(edm::Event& e, const edm::EventSetup& c);

	explicit CSCTFPacker(const edm::ParameterSet &conf);
	~CSCTFPacker(void);
};

#endif
