#include "DQM/HcalMonitorTasks/interface/HcalDataIntegrityTask.h"
#include "DataFormats/FEDRawData/interface/FEDNumbering.h"
#include "DataFormats/FEDRawData/interface/FEDTrailer.h"
#include "DQMServices/Core/interface/DQMStore.h"
#include "DQMServices/Core/interface/MonitorElement.h"

HcalDataIntegrityTask::HcalDataIntegrityTask() 
{
  //Initialize phatmap to a vector of vectors of uint64_t 0
  static size_t iphirange = IPHIMAX - IPHIMIN;
  static size_t ietarange = IETAMAX - IETAMIN;
 
  std::vector<uint64_t> phatv (iphirange + 1, 0);
  // ... nothing goes at ieta=0, so an extra bin goes there.
  phatmap = vector< vector < uint64_t> > ( ietarange + 1, phatv);
  HBmap   = vector< vector < uint64_t> > ( ietarange + 1, phatv);
  HEmap   = vector< vector < uint64_t> > ( ietarange + 1, phatv);
  HFmap   = vector< vector < uint64_t> > ( ietarange + 1, phatv);
  HOmap   = vector< vector < uint64_t> > ( ietarange + 1, phatv);
  std::vector<bool> probvect (iphirange + 1, 0);
  // ... nothing goes at ieta=0, so an extra bin goes there.
  problemhere = vector< vector <bool> > ( ietarange + 1, probvect);
  problemHB   = vector< vector <bool> > ( ietarange + 1, probvect);
  problemHE   = vector< vector <bool> > ( ietarange + 1, probvect);
  problemHF   = vector< vector <bool> > ( ietarange + 1, probvect);
  problemHO   = vector< vector <bool> > ( ietarange + 1, probvect);

} // HcalDataIntegrityTask::HcalDataIntegrityTask()

HcalDataIntegrityTask::~HcalDataIntegrityTask() {}

void HcalDataIntegrityTask::reset(){}

void HcalDataIntegrityTask::clearME()
{
  if(m_dbe){
    m_dbe->setCurrentFolder(baseFolder_);
    m_dbe->removeContents();
  }
}


void HcalDataIntegrityTask::setup(const edm::ParameterSet& ps,
				  DQMStore* dbe)
{
  HcalBaseMonitor::setup(ps,dbe);
  
  ievt_=0;
  baseFolder_ = rootFolder_+"HcalDataIntegrityTask";

  if(fVerbosity) 
    cout << "About to pushback fedUnpackList_" << endl;
  
  firstFED_ = FEDNumbering::getHcalFEDIds().first;
  for (int i=FEDNumbering::getHcalFEDIds().first; 
       i<=FEDNumbering::getHcalFEDIds().second;
       ++i) 
    {
      if(fVerbosity) cout << "[DFMon:]Pushback for fedUnpackList_: " << i <<endl;
      fedUnpackList_.push_back(i);
    }

  prtlvl_ = ps.getUntrackedParameter<int>("dfPrtLvl");

  if ( m_dbe ) 
    {
      char* type;
      
      if (fVerbosity)
	cout <<"SET TO HCAL/FEDIntegrity"<<endl;
      m_dbe->setCurrentFolder("Hcal/FEDIntegrity/");
      type="FEDEntries";
      fedEntries_ = m_dbe->book1D(type,"# entries per HCAL FED",32,699.5,731.5);
      type="FEDFatal";
      fedFatal_ = m_dbe->book1D(type,"# fatal errors HCAL FED",32,699.5,731.5);
    } // if (m_dbe)


  return;
}

void HcalDataIntegrityTask::processEvent(const FEDRawDataCollection& rawraw, 
					 const HcalUnpackerReport& report, 
					 const HcalElectronicsMap& emap){
  
  if(!m_dbe) 
    { 
      cout<<"HcalDataIntegrityTask::processEvent DQMStore not instantiated!!!"<<endl;  
      return;
    }
  
  ievt_++;

  // Loop over all FEDs reporting the event, unpacking if good.
  for (vector<int>::const_iterator i=fedUnpackList_.begin();i!=fedUnpackList_.end(); i++) 
    {
      const FEDRawData& fed = rawraw.FEDData(*i);
      if (fed.size()<12) continue; // Was 16. How do such tiny events even get here?
      unpack(fed,emap);
    }

  return;
} //void HcalDataIntegrityTask::processEvent()


// Process one FED's worth (one DCC's worth) of the event data.
void HcalDataIntegrityTask::unpack(const FEDRawData& raw, 
				   const HcalElectronicsMap& emap){
  // get the DCC header
  const HcalDCCHeader* dccHeader=(const HcalDCCHeader*)(raw.data());
  if(!dccHeader) return;

  // get the DCC trailer 
  unsigned char* trailer_ptr = (unsigned char*) (raw.data()+raw.size()-sizeof(uint64_t));
  FEDTrailer trailer = FEDTrailer(trailer_ptr);

  int dccid=dccHeader->getSourceId();

  ////////// Histogram problems with the Common Data Format compliance;////////////
  bool CDFProbThisDCC = false; 
  /* 1 */ //There should always be a second CDF header word indicated.
  if (!dccHeader->thereIsASecondCDFHeaderWord()) 
    {
      CDFProbThisDCC = true; 
    }

  /* 2 */ //Make sure a reference CDF Version value has been recorded for this dccid
  CDFvers_it = CDFversionNumber_list.find(dccid);
  if (CDFvers_it  == CDFversionNumber_list.end()) 
    {
      CDFversionNumber_list.insert(pair<int,short>
				   (dccid,dccHeader->getCDFversionNumber() ) );
      CDFvers_it = CDFversionNumber_list.find(dccid);
    } // then check against it.

  if (dccHeader->getCDFversionNumber()!= CDFvers_it->second) 
    {
      CDFProbThisDCC = true; 
    }
  
  /* 3 */ //Make sure a reference CDF EventType value has been recorded for this dccid
  CDFEvT_it = CDFEventType_list.find(dccid);
  if (CDFEvT_it  == CDFEventType_list.end()) 
    {
      CDFEventType_list.insert(pair<int,short>
			       (dccid,dccHeader->getCDFEventType() ) );
      CDFEvT_it = CDFEventType_list.find(dccid);
    } // then check against it.
  
  if (dccHeader->getCDFEventType()!= CDFEvT_it->second) 
    {
      CDFProbThisDCC = true; 
    }

  /* 4 */ //There should always be a '5' in CDF Header word 0, bits [63:60]
  if (dccHeader->BOEshouldBe5Always()!=5) 
    {
      CDFProbThisDCC = true; 
    }

  /* 5 */ //There should never be a third CDF Header word indicated.
  if (dccHeader->thereIsAThirdCDFHeaderWord()) 
    {
      CDFProbThisDCC = true; 
    }

  /* 6 */ //Make sure a reference value of Reserved Bits has been recorded for this dccid

  CDFReservedBits_it = CDFReservedBits_list.find(dccid);
  if (CDFReservedBits_it  == CDFReservedBits_list.end()) {
    CDFReservedBits_list.insert(pair<int,short>
				(dccid,dccHeader->getSlink64ReservedBits() ) );
    CDFReservedBits_it = CDFReservedBits_list.find(dccid);
  } // then check against it.
  
  if ((int) dccHeader->getSlink64ReservedBits()!= CDFReservedBits_it->second) 
    {
      CDFProbThisDCC = true; 
    }

  /* 7 */ //There should always be 0x0 in CDF Header word 1, bits [63:60]
  if (dccHeader->BOEshouldBeZeroAlways() !=0) 
    {
      CDFProbThisDCC = true; 
    }
  
  /* 8 */ //There should only be one trailer
  if (trailer.moreTrailers()) 
    {
      CDFProbThisDCC = true; 
    }
  //  if trailer.

  /* 9 */ //CDF Trailer [55:30] should be the # 64-bit words in the EvFragment
  if ((uint64_t) raw.size() != ( (uint64_t) trailer.lenght()*sizeof(uint64_t)) )  //The function name is a typo! Awesome.
    {
      CDFProbThisDCC = true; 
    }
  /*10 */ //There is a rudimentary sanity check built into the FEDTrailer class
  if (!trailer.check()) 
    {
      CDFProbThisDCC = true; 
    }

  if (CDFProbThisDCC)
    fedFatal_->Fill(dccid);
  fedEntries_->Fill(dccid);

  return;
} // void HcalDataIntegrityTask::unpack(


