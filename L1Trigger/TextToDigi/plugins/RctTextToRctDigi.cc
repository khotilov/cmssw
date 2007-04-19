// -*- C++ -*-
//
// Package:    RctTextToRctDigi
// Class:      RctTextToRctDigi
// 
/**\class RctTextToRctDigi RctTextToRctDigi.cc L1Trigger/TextToDigi/src/RctTextToRctDigi.cc

Description: Makes RCT digis from the file format specified by Pam Klabbers

*/
//
// Original Author:  Alex Tapper
//         Created:  Fri Mar  9 19:11:51 CET 2007
// $Id: RctTextToRctDigi.cc,v 1.1 2007/03/21 18:48:23 tapper Exp $

// Rct Input File Format 
// Line 1: Crossing no as "Crossing x" (2)     
// Line 2: isoe0 isoe1 isoe2 isoe3 nonIsoe0 nonIsoe1 nonIso2 nonIso3 (8) 
// Line 3: RC0mip0 RC0mip1 RC1mip0 RC1mip1 RC2mip0 RC2mip1 RC3mip0 RC3mip1 RC4mip0 RC4mip1 
//         RC5mip0 RC5mip1 RC6mip0 RC6mip1 (14)
// Line 4: RC0qt0 RCqt1 RC1qt0 RC1qt1 RC2qt0 RC2qt1 RC3qt0 RC3qt1 RC4qt0 RC4qt1 
//         RC5qt0 RC5qt1 RC6qt0 RC6qt1 (14)
// Line 5: RC0reg0 RC0reg1 RC1reg0 RC1reg1 RC2reg0 RC2reg1 RC3reg0 RC3reg1 RC4reg0 RC4reg1
//         RC5reg0 RC5reg1 RC6reg0 RC6reg1 (14)
// Line 6: HF0eta0 HF0eta1 HF0eta2 HF0eta3 HF1eta0 HF1eta1 HF1eta2 HF1eta3 (8)
//
// NOTE:  CMS IN 2004/009 specifies that cable four provides 8 Quiet (fineGrain) bits for the HF.  These are not
//        detailed in the fileformat above, and are not currently dealt with in any way. Set to true.
// 

#include "RctTextToRctDigi.h"
#include "FWCore/ServiceRegistry/interface/Service.h" // Framework services
#include "FWCore/MessageLogger/interface/MessageLogger.h" // Logger

using namespace edm;
using namespace std;

// Set constant
const static unsigned NUM_RCT_CRATES = 18;

RctTextToRctDigi::RctTextToRctDigi(const edm::ParameterSet& iConfig):
  m_textFileName(iConfig.getParameter<std::string>("TextFileName")),
  m_skipEvents(iConfig.getParameter<int>("SkipEvents")),
  m_nevt(0)
{
  // Produces collections
  produces<L1CaloEmCollection>();
  produces<L1CaloRegionCollection>();

  // Open the input files
  for (unsigned i=0; i<NUM_RCT_CRATES; i++){
    stringstream fileStream;
    fileStream << m_textFileName << i;
    string fileName(fileStream.str());
    m_file[i].open(fileName.c_str(),ios::in);

    if(!m_file[i].good())
      {
        throw cms::Exception("RctTextToRctDigiTextFileOpenError")
          << "RctTextToRctDigi::RctTextToRctDigi : "
          << " couldn't open the file " << fileName << " for reading" << endl;
      }
  }
}

RctTextToRctDigi::~RctTextToRctDigi()
{
  // Close the input files
  for (unsigned i=0; i<NUM_RCT_CRATES; i++){  
    m_file[i].close();
  }
}

// ------------ method called to produce the data  ------------
void RctTextToRctDigi::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{

  // Skip event if required
  if (m_nevt < m_skipEvents){ 
    string tmp;
    for (int i=0; i<6; i++){
      getline(m_file[i],tmp);
    } 
    m_nevt++;
    return;
  }

  // New collections
  auto_ptr<L1CaloEmCollection> em (new L1CaloEmCollection);
  auto_ptr<L1CaloRegionCollection> rgn (new L1CaloRegionCollection);

  // Loop over RCT crates
  for (unsigned i=0; i<NUM_RCT_CRATES; i++){  

    // Check we're not at the end of the file
    if(m_file[i].eof())
      {
        throw cms::Exception("RctTextToRctDigiTextFileReadError")
          << "RctTextToRctDigi::produce : "
          << " unexpected end of file " << m_textFileName << i << endl;
      }      
    
    // Check we're at the start of an event
    string tmp;
    m_file[i]>> tmp;
    if(tmp!="Crossing")
      {
        throw cms::Exception("RctTextToRctDigiTextFileReadError")
          << "RctTextToRctDigi::produce : "
          << " something screwy happened Crossing!=" << tmp << endl;
      }      

    // Read BX number
    dec(m_file[i]);
    int BXNum;
    m_file[i]>>BXNum;
    
    // Buffers
    unsigned long int uLongBuffer;
    bool mipBitBuffer[14],qBitBuffer[14];

    // All in hex from now on
    hex(m_file[i]); 

    // Isolated electrons
    for (unsigned j=0; j<4; j++){
      m_file[i] >> uLongBuffer;
      em->push_back(L1CaloEmCand(uLongBuffer, i, true));
    }

    // Non-isolated electrons
    for (unsigned j=0; j<4; j++){
      m_file[i] >> uLongBuffer;
      em->push_back(L1CaloEmCand(uLongBuffer, i, false));
    }      
    
    // MIP bits 
    for (unsigned j=0; j<14; j++){
      m_file[i] >> mipBitBuffer[j];
    }   

    // Quiet bits 
    for (unsigned j=0; j<14; j++){
      m_file[i] >> qBitBuffer[j];
    }     

    // Barrel and endcap regions
    for (unsigned j=0; j<14; j++){
      m_file[i] >> uLongBuffer;

      unsigned et = uLongBuffer & 0x3ff;  // put the first 10 bits of rawData into the Et
      uLongBuffer >>= 10;  // shift the remaining bits down to remove the 10 bits of Et
      
      bool overFlow = ((uLongBuffer & 0x1)        != 0); //LSB is now overflow bit
      bool tauVeto  = (((uLongBuffer & 0x2) >> 1) != 0); //2nd bit is tauveto      

      rgn->push_back(L1CaloRegion(et,overFlow,tauVeto,mipBitBuffer[j],qBitBuffer[j],i,j/2,j%2));
    }      
    
    // HF
    for (unsigned j=0; j<8; j++){
      m_file[i] >> uLongBuffer;

      unsigned et = uLongBuffer & 0xff;  // put the first 8 bits into the Et

      rgn->push_back(L1CaloRegion(et,true,i,j));
    }        

  }
  
  iEvent.put(em);
  iEvent.put(rgn);

  m_nevt++;
}



