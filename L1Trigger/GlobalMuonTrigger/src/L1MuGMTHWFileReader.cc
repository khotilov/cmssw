//-------------------------------------------------
//
//   \class L1MuGMTHWFileReader
//
//   Description: Puts the GMT input information from 
//                a GMT ascii HW testfile into the Event
//
//
//   $Date$
//   $Revision$
//
//   Author :
//   Tobias Noebauer                 HEPHY Vienna
//   Ivan Mikulec                    HEPHY Vienna
//
//--------------------------------------------------

//-----------------------
// This Class's Header --
//-----------------------
#include "L1Trigger/GlobalMuonTrigger/src/L1MuGMTHWFileReader.h"

//---------------
// C++ Headers --
//---------------
#include <stdexcept>

//-------------------------------
// Collaborating Class Headers --
//-------------------------------
#include "DataFormats/L1GlobalMuonTrigger/interface/L1MuRegionalCand.h"

//----------------
// Constructors --
//----------------
L1MuGMTHWFileReader::L1MuGMTHWFileReader(edm::ParameterSet const& ps,
                                         edm::InputSourceDescription const& desc) :
                                         ExternalInputSource(ps, desc) {

  produces<std::vector<L1MuRegionalCand> >("DT");
  produces<std::vector<L1MuRegionalCand> >("CSC");
  produces<std::vector<L1MuRegionalCand> >("RPCb");
  produces<std::vector<L1MuRegionalCand> >("RPCf");

  if(!fileNames().size()) {
    throw runtime_error("L1MuGMTHWFileReader: no input file");
  }
  cout << "opening file " << fileNames()[0] << endl;
  m_in.open((fileNames()[0].substr(fileNames()[0].find(":")+1)).c_str());
  if(!m_in) {
    throw runtime_error("L1MuGMTHWFileReader: file " + fileNames()[0]
			+ " could not be openned");
  }

}

//--------------
// Destructor --
//--------------
L1MuGMTHWFileReader::~L1MuGMTHWFileReader() {
  m_in.close();
}

//--------------
// Operations --
//--------------
void L1MuGMTHWFileReader::setRunAndEventInfo() {
  readNextEvent();
  setRunNumber(m_evt.getRunNumber());
  setEventNumber(m_evt.getEventNumber());

  cout << "run: " << m_evt.getRunNumber() << 
          "   evt: " << m_evt.getEventNumber() << endl;
}

bool L1MuGMTHWFileReader::produce(edm::Event& e) {
  L1MuRegionalCand empty_mu;

  if(!m_evt.getRunNumber() && !m_evt.getEventNumber()) return false;

  std::auto_ptr<std::vector<L1MuRegionalCand> > DTCands(new std::vector<L1MuRegionalCand>);
  for (unsigned i = 0; i < 4; i++) {
    const L1MuRegionalCand *mu = m_evt.getInputMuon("IND", i);
    if (!mu) mu = &empty_mu;
    DTCands->push_back(*mu);
  }
  e.put(DTCands,"DT");

  std::auto_ptr<std::vector<L1MuRegionalCand> > CSCCands(new std::vector<L1MuRegionalCand>);
  for (unsigned i = 0; i < 4; i++) {
    const L1MuRegionalCand *mu = m_evt.getInputMuon("INC", i);
    if (!mu) mu = &empty_mu;
    CSCCands->push_back(*mu);
  }
  e.put(CSCCands,"CSC");

  std::auto_ptr<std::vector<L1MuRegionalCand> > RPCbCands(new std::vector<L1MuRegionalCand>);
  for (unsigned i = 0; i < 4; i++) {
    const L1MuRegionalCand *mu = m_evt.getInputMuon("INB", i);
    if (!mu) mu = &empty_mu;
    RPCbCands->push_back(*mu);
  }
  e.put(RPCbCands,"RPCb");

  std::auto_ptr<std::vector<L1MuRegionalCand> > RPCfCands(new std::vector<L1MuRegionalCand>);
  for (unsigned i = 0; i < 4; i++) {
    const L1MuRegionalCand *mu = m_evt.getInputMuon("INF", i);
    if (!mu) mu = &empty_mu;
    RPCfCands->push_back(*mu);
  }
  e.put(RPCfCands,"RPCf");

  return true;
}

void L1MuGMTHWFileReader::readNextEvent() {
  m_evt.reset();

  string line_id;
  do {
    int bx = 0;

    m_in >> line_id;
    if (line_id == "--") continue;

    if (line_id == "RUN") {
      unsigned long runnr;
      m_in >> runnr;
      m_evt.setRunNumber(runnr);
    }

    if (line_id == "EVT") {
      unsigned long evtnr;
      m_in >> evtnr;
      m_evt.setEventNumber(evtnr);
    }


    if (line_id == "DT" || line_id == "CSC" || line_id == "BRPC" || line_id == "FRPC")
    {

      // decode input muon

      unsigned inpmu = 0;
      unsigned val;
      m_in >> val; inpmu |= (val & 0x01) << 24; // valid charge
      m_in >> val; inpmu |= (val & 0x01) << 23; // charge
      m_in >> val; inpmu |= (val & 0x01) << 22; // halo / fine
      m_in >> val; inpmu |= (val & 0x3f) << 16; // eta
      m_in >> val; inpmu |= (val & 0x07) << 13; // quality
      m_in >> val; inpmu |= (val & 0x1f) <<  8; // pt
      m_in >> val; inpmu |= (val & 0xff)      ; // phi

      string chipid("IN");
      chipid += line_id[0];

      int type=0;
      if (line_id == "DT") type = 0;
      if (line_id == "CSC") type = 2;
      if (line_id == "BRPC") type = 1;
      if (line_id == "FRPC") type = 3;


      L1MuRegionalCand cand(inpmu);
      cand.setType(type);
      cand.setBx(bx);
      m_evt.addInputMuon(chipid, cand);
    }

    if (line_id == "MIP") {
      int nPairs;
      m_in >> nPairs;
      for (int i=0; i<nPairs; i++) {
        unsigned eta;
        unsigned phi;
        m_in >> eta;
        m_in >> phi;
        if (phi >= 9) phi-=9;
        else phi+=9;
        m_evt.setMipBit(eta, phi, true);
      }
    }

    if (line_id == "NQ") {
      int nPairs;
      m_in >> nPairs;
      for (int i=0; i<nPairs; i++) {
        unsigned eta;
        unsigned phi;
        m_in >> eta;
        m_in >> phi;
        if (phi >= 9) phi-=9;
        else phi+=9;
        m_evt.setIsoBit(eta, phi, false);
      }
    }

    //read the rest of the line
    const int sz=4000; char buf[sz];
    m_in.getline(buf, sz);

  } while (line_id != "NQ" && !m_in.eof());

}
