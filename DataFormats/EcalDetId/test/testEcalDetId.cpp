/**
   \file
   Test suit for EcalDetId

   \version $Id: testEcalDetId.cc,v 1.7 2006/03/14 15:52:55 meridian Exp $

   \note This test is not exaustive     
*/

static const char CVSId[] = "$Id: testEcalDetId.cc,v 1.7 2006/03/14 15:52:55 meridian Exp $";

#include <Utilities/Testing/interface/CppUnit_testdriver.icpp>
#include <cppunit/extensions/HelperMacros.h>
#include <DataFormats/EcalDetId/interface/EBDetId.h>
#include <DataFormats/EcalDetId/interface/EEDetId.h>
#include <DataFormats/EcalDetId/interface/ESDetId.h>
#include <DataFormats/EcalDetId/interface/EcalTrigTowerDetId.h>
#include <DataFormats/EcalDetId/interface/EcalElectronicsId.h>
#include <DataFormats/EcalDetId/interface/EcalPnDiodeDetId.h>


#include <iostream>

class testEcalDetId: public CppUnit::TestFixture {

  CPPUNIT_TEST_SUITE(testEcalDetId);

  CPPUNIT_TEST(testEBDetId);
  CPPUNIT_TEST(testEEDetId);
  CPPUNIT_TEST(testESDetId);
  CPPUNIT_TEST(testEcalTrigTowerDetId);
  CPPUNIT_TEST(testEcalElectronicsId);
  CPPUNIT_TEST(testPnDiodeDetId);

  CPPUNIT_TEST_SUITE_END();

public:


  void setUp(){}
  void tearDown(){}  
  void testEBDetId();
  void testEEDetId();
  void testESDetId();
  void testEcalTrigTowerDetId();
  void testEcalElectronicsId();
  void testPnDiodeDetId();
 
}; 

///registration of the test so that the runner can find it
CPPUNIT_TEST_SUITE_REGISTRATION(testEcalDetId);

void testEcalDetId::testEBDetId(){

  EBDetId smId;

  for (int ieta=EBDetId::MIN_IETA;ieta<=EBDetId::MAX_IETA;ieta++)
    for (int iphi=EBDetId::MIN_IPHI;iphi<=EBDetId::MAX_IPHI;iphi++)
      {
	//EBDetId Zside 1 
	try
	  {
	    {
	      EBDetId aPositiveId(ieta,iphi);
	      CPPUNIT_ASSERT(aPositiveId.ieta()==ieta);
	      CPPUNIT_ASSERT(aPositiveId.iphi()==iphi);
	      CPPUNIT_ASSERT(aPositiveId.zside()==1);
	      CPPUNIT_ASSERT(aPositiveId.ietaAbs()==ieta);

	      smId = EBDetId(aPositiveId.ism(), aPositiveId.ic(),
			     EBDetId::SMCRYSTALMODE);
	      CPPUNIT_ASSERT(aPositiveId==smId);
	      CPPUNIT_ASSERT(aPositiveId.ism()==smId.ism());
	      CPPUNIT_ASSERT(aPositiveId.ic()==smId.ic());
	      CPPUNIT_ASSERT(smId.ism()>=EBDetId::MIN_SM);
	      CPPUNIT_ASSERT(smId.ism()<=EBDetId::MAX_SM);
	      CPPUNIT_ASSERT(smId.ic()>=EBDetId::MIN_C);
	      CPPUNIT_ASSERT(smId.ic()<=EBDetId::MAX_C);
	    }

	    //EBDetId Zside -1 
	    {
	      EBDetId aNegativeId(-1*ieta,iphi);
	      CPPUNIT_ASSERT(aNegativeId.ieta()==-1*ieta);
	      CPPUNIT_ASSERT(aNegativeId.iphi()==iphi);
	      CPPUNIT_ASSERT(aNegativeId.zside()==-1);
	      CPPUNIT_ASSERT(aNegativeId.ietaAbs()==ieta);

	      smId = EBDetId(aNegativeId.ism(), aNegativeId.ic(),
			     EBDetId::SMCRYSTALMODE);
	      CPPUNIT_ASSERT(aNegativeId==smId);
	      CPPUNIT_ASSERT(aNegativeId.ism()==smId.ism());
	      CPPUNIT_ASSERT(aNegativeId.ic()==smId.ic());
	      CPPUNIT_ASSERT(smId.ism()>=EBDetId::MIN_SM);
	      CPPUNIT_ASSERT(smId.ism()<=EBDetId::MAX_SM);
	      CPPUNIT_ASSERT(smId.ic()>=EBDetId::MIN_C);
	      CPPUNIT_ASSERT(smId.ic()<=EBDetId::MAX_C);

	    }
	  }
	catch ( std::runtime_error &e ) 
	  { 
	  }
      }

  
}

void testEcalDetId::testEEDetId(){
  
  
  for (int ix=EEDetId::IX_MIN;ix<=EEDetId::IX_MAX;ix++)
    for (int iy=EEDetId::IY_MIN;iy<=EEDetId::IY_MAX;iy++)
      {
	try
	  {
	    //EEDetId Zside 1 
	    {
	      EEDetId aPositiveId(ix,iy,1);
	      CPPUNIT_ASSERT(aPositiveId.ix()==ix);
	      CPPUNIT_ASSERT(aPositiveId.iy()==iy);
	      CPPUNIT_ASSERT(aPositiveId.zside()==1);
	    }
	    //EEDetId Zside -1 
	    {
	      EEDetId aNegativeId(ix,iy,-1);
	      CPPUNIT_ASSERT(aNegativeId.ix()==ix);
	      CPPUNIT_ASSERT(aNegativeId.iy()==iy);
	      CPPUNIT_ASSERT(aNegativeId.zside()==-1);
	    }
	  }
	catch ( std::runtime_error &e ) 
	  { 
	  }
      }
  
}

void testEcalDetId::testESDetId(){
  
  
  for (int istrip=ESDetId::ISTRIP_MIN;istrip<=ESDetId::ISTRIP_MAX;istrip++)
    for (int ix=ESDetId::IX_MIN;ix<=ESDetId::IX_MAX;ix++)
      for (int iy=ESDetId::IY_MIN;iy<=ESDetId::IY_MAX;iy++)
	{
	  try 
	    {
	      //ESDetId Plane 0 Zside 1 
	      {
		ESDetId aPositiveId(istrip,ix,iy,1,1);
		CPPUNIT_ASSERT(aPositiveId.strip()==istrip);
		CPPUNIT_ASSERT(aPositiveId.six()==ix);
		CPPUNIT_ASSERT(aPositiveId.siy()==iy);
		CPPUNIT_ASSERT(aPositiveId.plane()==1);
		CPPUNIT_ASSERT(aPositiveId.zside()==1);

	      }
	      //ESDetId Plane 1 Zside 1 
	      {
		ESDetId aPositiveId(istrip,ix,iy,2,1);
		CPPUNIT_ASSERT(aPositiveId.strip()==istrip);
		CPPUNIT_ASSERT(aPositiveId.six()==ix);
		CPPUNIT_ASSERT(aPositiveId.siy()==iy);
		CPPUNIT_ASSERT(aPositiveId.plane()==2);
		CPPUNIT_ASSERT(aPositiveId.zside()==1);
	      }
	      //ESDetId Plane 0 Zside -1 
	      {
		ESDetId aNegativeId(istrip,ix,iy,1,-1);
		CPPUNIT_ASSERT(aNegativeId.strip()==istrip);
		CPPUNIT_ASSERT(aNegativeId.six()==ix);
		CPPUNIT_ASSERT(aNegativeId.siy()==iy);
		CPPUNIT_ASSERT(aNegativeId.plane()==1);
		CPPUNIT_ASSERT(aNegativeId.zside()==-1);
	      }
	      //ESDetId Plane 1 Zside -1 
	      {
		ESDetId aNegativeId(istrip,ix,iy,2,-1);
		CPPUNIT_ASSERT(aNegativeId.strip()==istrip);
		CPPUNIT_ASSERT(aNegativeId.six()==ix);
		CPPUNIT_ASSERT(aNegativeId.siy()==iy);
		CPPUNIT_ASSERT(aNegativeId.plane()==2);
		CPPUNIT_ASSERT(aNegativeId.zside()==-1);
	      }
	    }
	  catch ( std::runtime_error &e ) 
	    { 
	    }
	}
}


void testEcalDetId::testEcalTrigTowerDetId() {

  for (int ieta=EcalTrigTowerDetId::MIN_I;ieta<=EcalTrigTowerDetId::MAX_I;ieta++)
    for (int iphi=EcalTrigTowerDetId::MIN_I;iphi<=EcalTrigTowerDetId::MAX_I;iphi++)
      {
	try
	  {
	    //EcalTrigTowerDetId Zside 1 
	    {
	      EcalTrigTowerDetId aPositiveId(1,EcalBarrel,ieta,iphi);
	      CPPUNIT_ASSERT(aPositiveId.subDet()==EcalBarrel);
	      CPPUNIT_ASSERT(aPositiveId.ieta()==ieta);
	      CPPUNIT_ASSERT(aPositiveId.iphi()==iphi);
	      CPPUNIT_ASSERT(aPositiveId.zside()==1);
	      CPPUNIT_ASSERT(aPositiveId.ietaAbs()==ieta);
	    }
	    //EcalTrigTowerDetId Zside -1 
	    {
	      EcalTrigTowerDetId aNegativeId(-1,EcalBarrel,ieta,iphi);
	      CPPUNIT_ASSERT(aNegativeId.subDet()==EcalBarrel);
	      CPPUNIT_ASSERT(aNegativeId.ieta()==-1*ieta);
	      CPPUNIT_ASSERT(aNegativeId.iphi()==iphi);
	      CPPUNIT_ASSERT(aNegativeId.zside()==-1);
	      CPPUNIT_ASSERT(aNegativeId.ietaAbs()==ieta);
	    }
	  }
	catch ( std::runtime_error &e ) 
	  { 
	  }
      }
  for (int ix=EcalTrigTowerDetId::MIN_I;ix<=EcalTrigTowerDetId::MAX_I;ix++)
    for (int iy=EcalTrigTowerDetId::MIN_I;iy<=EcalTrigTowerDetId::MAX_I;iy++)
      {
	try
	  {
	    //EcalTrigTowerDetId Zside 1 
	    {
	      EcalTrigTowerDetId aPositiveId(1,EcalEndcap,ix,iy);
	      CPPUNIT_ASSERT(aPositiveId.ix()==ix);
	      CPPUNIT_ASSERT(aPositiveId.subDet()==EcalEndcap);
	      CPPUNIT_ASSERT(aPositiveId.iy()==iy);
	      CPPUNIT_ASSERT(aPositiveId.zside()==1);
	    }
	    //EcalTrigTowerDetId Zside -1 
	    {
	      EcalTrigTowerDetId aNegativeId(-1,EcalEndcap,ix,iy);
	      CPPUNIT_ASSERT(aNegativeId.subDet()==EcalEndcap);
	      CPPUNIT_ASSERT(aNegativeId.ix()==ix);
	      CPPUNIT_ASSERT(aNegativeId.iy()==iy);
	      CPPUNIT_ASSERT(aNegativeId.zside()==-1);
	    }
	  }
	catch ( std::runtime_error &e ) 
	  { 
	  }
      }
}

void testEcalDetId::testEcalElectronicsId() {
  for (int idcc=EcalElectronicsId::MIN_DCCID;idcc<=EcalElectronicsId::MAX_DCCID;idcc++)
    for (int itower=EcalElectronicsId::MIN_TOWERID;itower<EcalElectronicsId::MAX_TOWERID;itower++)
      for (int ichannel=EcalElectronicsId::MIN_CHANNELID;ichannel<EcalElectronicsId::MAX_CHANNELID;ichannel++)
	{
	  try
	    {
	      EcalElectronicsId electronicsId(idcc,itower,ichannel);
	      CPPUNIT_ASSERT(electronicsId.dccId()==idcc);
	      CPPUNIT_ASSERT(electronicsId.towerId()==itower);
	      CPPUNIT_ASSERT(electronicsId.channelId()==ichannel);
	    }
	  catch ( std::runtime_error &e ) 
	    { 
	    }
	}
}

void testEcalDetId::testPnDiodeDetId() {

  for (int iEcalSubDet=EcalBarrel ; iEcalSubDet<=EcalEndcap ; iEcalSubDet++)
    for (int iDCC=1;iDCC<=36;iDCC++)
      for (int iPn=1;iPn<=10;iPn++)
	{
	  try
	    {
	      //EcalTrigTowerDetId Zside 1 
	      {
		EcalPnDiodeDetId aId(iEcalSubDet,iDCC,iPn);
		CPPUNIT_ASSERT(aId.iEcalSubDetectorId()==iEcalSubDet);
		CPPUNIT_ASSERT(aId.iDCCId()==iDCC);
		CPPUNIT_ASSERT(aId.iPnId()==iPn);
	      }
	    }
	  catch ( std::runtime_error &e ) 
	    { 
	    }
	}
}
