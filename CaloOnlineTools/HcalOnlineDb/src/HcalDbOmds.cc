// -*- C++ -*-
//
// Original Author:  Gena Kukartsev Mar 11, 2009
// Adapted from HcalDbASCIIIO.cc,v 1.41
// $Id: HcalDbOmds.cc,v 1.9 2009/03/26 18:57:40 kukartse Exp $
//
//
#include <vector>
#include <string>

#include "DataFormats/HcalDetId/interface/HcalGenericDetId.h"
#include "DataFormats/HcalDetId/interface/HcalZDCDetId.h"
#include "DataFormats/HcalDetId/interface/HcalCastorDetId.h"
#include "DataFormats/HcalDetId/interface/HcalElectronicsId.h"
#include "CalibFormats/HcalObjects/interface/HcalText2DetIdConverter.h"

#include "CondFormats/HcalObjects/interface/AllObjects.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "CaloOnlineTools/HcalOnlineDb/interface/HcalDbOmds.h"
#include "CaloOnlineTools/HcalOnlineDb/interface/RooGKCounter.h"

// get the proper detId from the result of the oracle query
// assumed that channel info comes in from of the ResultSet
// in the following order (some may not be filled):
// 1.objectname, ( values: HcalDetId, HcalCalibDetId, HcalTrigTowerDetId, HcalZDCDetId or HcalCastorDetId)
// 2.subdet, 3.ieta, 4.iphi, 5.depth, 6.type, 7.section, 8.ispositiveeta, 9.sector, 10.module, 11.channel 
DetId HcalDbOmds::getId(oracle::occi::ResultSet * rs){
  std::string _name = rs->getString(1);
  //cerr << "DEBUG: name - " << _name << endl;
  if (rs->getString(1).find("HcalDetId")!=std::string::npos){
    //cerr << "DEBUG: HcalDetId" << endl;
    return HcalDetId(get_subdetector(rs->getString(2)),
		     rs->getInt(3),
		     rs->getInt(4),
		     rs->getInt(5));
  }
  else if (rs->getString(1).find("HcalCalibDetId")!=std::string::npos){
    //cerr << "DEBUG: HcalCalibDetId" << endl;
    return HcalCalibDetId(get_subdetector(rs->getString(2)),
			  rs->getInt(3),
			  rs->getInt(4),
			  rs->getInt(6));
  }
  else if (rs->getString(1).find("HcalTrigTowerDetId")!=std::string::npos){
    //cerr << "DEBUG: HcalTrigTowerDetId" << endl;
    return HcalTrigTowerDetId(rs->getInt(3),
			      rs->getInt(4));
  }
  else if (rs->getString(1).find("HcalZDCDetId")!=std::string::npos){
    //cerr << "DEBUG: HcalZDCDetId" << endl;
    return HcalZDCDetId(get_zdc_section(rs->getString(7)),
			rs->getInt(8)>0,
			rs->getInt(11));
  }
  else if (rs->getString(1).find("HcalCastorDetId")!=std::string::npos){
    //cerr << "DEBUG: HcalCastorDetId" << endl;
    return HcalCastorDetId(rs->getInt(8)>0,
			   rs->getInt(9)>0,
			   rs->getInt(10));
  }
  else return 0;
}

bool HcalDbOmds::getObject (oracle::occi::Connection * connection, 
			    const std::string & fTag, 
			    const std::string & fVersion,
			    const int fSubversion,
			    const std::string & fQuery,
			    HcalPedestals* fObject) {
  bool result=true;
  if (!fObject) fObject = new HcalPedestals;
  int _unit=0;
  try {
    oracle::occi::Statement* stmt = connection->createStatement(fQuery);
    stmt->setString(1,fTag);
    stmt->setString(2,fVersion);
    //stmt->setInt(3,fSubversion);

    ResultSet *rs = stmt->executeQuery();

    RooGKCounter _row(1,100);
    _row.setMessage("HCAL channels processed: ");
    _row.setPrintCount(true);
    _row.setNewLine(true);
    while (rs->next()) {
      _row.count();
      _unit = rs->getInt(1);
      float cap0 = rs->getFloat(2);
      float cap1 = rs->getFloat(3);
      float cap2 = rs->getFloat(4);
      float cap3 = rs->getFloat(5);
      float variance0 = rs->getFloat(6);
      float variance1 = rs->getFloat(7);
      float variance2 = rs->getFloat(8);
      float variance3 = rs->getFloat(9);
      int ieta = rs->getInt(10);
      int iphi = rs->getInt(11);
      int depth = rs->getInt(12);
      HcalSubdetector subdetector = get_subdetector(rs->getString(13));
      HcalDetId id(subdetector,ieta,iphi,depth);
      //cout << "DEBUG: " << id << " " << cap0 << " " << cap1 << " " << cap2 << " " << cap3 << endl;
      HcalPedestal * fCondObject = new HcalPedestal(id.rawId(), cap0, cap1, cap2, cap3, variance0, variance1, variance2, variance3);
      fObject->addValues(*fCondObject);
      delete fCondObject;
    }
    //Always terminate statement
    connection->terminateStatement(stmt);
  } catch (SQLException& e) {
    throw cms::Exception("ReadError") << ::toolbox::toString("Oracle  exception : %s",e.getMessage().c_str()) << std::endl;
  }
  bool unit_is_adc = false;
  if (_unit!=0) unit_is_adc = true;
  fObject->setUnitADC(unit_is_adc);
  return result;
}


bool HcalDbOmds::getObject (oracle::occi::Connection * connection, 
			    const std::string & fTag, 
			    const std::string & fVersion,
			    const int fSubversion,
			    const std::string & fQuery,
			    HcalPedestalWidths* fObject) {
  bool result=true;
  if (!fObject) fObject = new HcalPedestalWidths;
  int _unit=0;
  try {
    oracle::occi::Statement* stmt = connection->createStatement(fQuery);
    stmt->setString(1,fTag);
    stmt->setString(2,fVersion);
    //stmt->setInt(3,fSubversion);

    ResultSet *rs = stmt->executeQuery();

    rs->setMaxColumnSize(1,128);
    rs->setMaxColumnSize(2,128);
    rs->setMaxColumnSize(3,128);
    rs->setMaxColumnSize(4,128);
    rs->setMaxColumnSize(5,128);
    rs->setMaxColumnSize(6,128);
    rs->setMaxColumnSize(7,128);
    rs->setMaxColumnSize(8,128);
    rs->setMaxColumnSize(9,128);
    rs->setMaxColumnSize(10,128);
    rs->setMaxColumnSize(11,128);

    RooGKCounter _row(1,100);
    _row.setMessage("HCAL channels processed: ");
    _row.setPrintCount(true);
    _row.setNewLine(true);
    //
    // The query result must be ordered in the following way
    // 1.objectname, ( values: HcalDetId, HcalCalibDetId, HcalTrigTowerDetId, HcalZDCDetId or HcalCastorDetId)
    // 2.subdet, 3.ieta, 4.iphi, 5.depth, 6.type, 7.section, 8.ispositiveeta, 9.sector, 10.module, 11.channel 
    // 12. is_ADC(int), 13-28. covariance00-01-...-10-11-...33
    //
    while (rs->next()) {
      _row.count();
      DetId id = getId(rs);
      _unit = rs->getInt(12);
      float covariance_00 = rs->getFloat(13);
      float covariance_01 = rs->getFloat(14);
      float covariance_02 = rs->getFloat(15);
      float covariance_03 = rs->getFloat(16);
      float covariance_10 = rs->getFloat(17);
      float covariance_11 = rs->getFloat(18);
      float covariance_12 = rs->getFloat(19);
      float covariance_13 = rs->getFloat(20);
      float covariance_20 = rs->getFloat(21);
      float covariance_21 = rs->getFloat(22);
      float covariance_22 = rs->getFloat(23);
      float covariance_23 = rs->getFloat(24);
      float covariance_30 = rs->getFloat(25);
      float covariance_31 = rs->getFloat(26);
      float covariance_32 = rs->getFloat(27);
      float covariance_33 = rs->getFloat(28);
      HcalPedestalWidth * fCondObject = new HcalPedestalWidth(id);
      fCondObject->setSigma(0,0,covariance_00);
      fCondObject->setSigma(0,1,covariance_01);
      fCondObject->setSigma(0,2,covariance_02);
      fCondObject->setSigma(0,3,covariance_03);
      fCondObject->setSigma(1,0,covariance_10);
      fCondObject->setSigma(1,1,covariance_11);
      fCondObject->setSigma(1,2,covariance_12);
      fCondObject->setSigma(1,3,covariance_13);
      fCondObject->setSigma(2,0,covariance_20);
      fCondObject->setSigma(2,1,covariance_21);
      fCondObject->setSigma(2,2,covariance_22);
      fCondObject->setSigma(2,3,covariance_23);
      fCondObject->setSigma(3,0,covariance_30);
      fCondObject->setSigma(3,1,covariance_31);
      fCondObject->setSigma(3,2,covariance_32);
      fCondObject->setSigma(3,3,covariance_33);
      fObject->addValues(*fCondObject);
      delete fCondObject;
    }
    //Always terminate statement
    connection->terminateStatement(stmt);
  } catch (SQLException& e) {
    throw cms::Exception("ReadError") << ::toolbox::toString("Oracle  exception : %s",e.getMessage().c_str()) << std::endl;
  }
  bool unit_is_adc = false;
  if (_unit!=0) unit_is_adc = true;
  fObject->setUnitADC(unit_is_adc);
  return result;
}


bool HcalDbOmds::getObject (oracle::occi::Connection * connection, 
			    const std::string & fTag, 
			    const std::string & fVersion,
			    const int fSubversion,
			    const std::string & fQuery,
			    HcalGains* fObject) {
  bool result=true;
  if (!fObject) fObject = new HcalGains;
  try {
    oracle::occi::Statement* stmt = connection->createStatement(fQuery);
    stmt->setString(1,fTag);
    stmt->setString(2,fVersion);
    //stmt->setInt(3,fSubversion);

    ResultSet *rs = stmt->executeQuery();

    RooGKCounter _row(1,100);
    _row.setMessage("HCAL channels processed: ");
    _row.setPrintCount(true);
    _row.setNewLine(true);
    while (rs->next()) {
      _row.count();
      float cap0 = rs->getFloat(1);
      float cap1 = rs->getFloat(2);
      float cap2 = rs->getFloat(3);
      float cap3 = rs->getFloat(4);
      int ieta = rs->getInt(5);
      int iphi = rs->getInt(6);
      int depth = rs->getInt(7);
      HcalSubdetector subdetector = get_subdetector(rs->getString(8));
      HcalDetId id(subdetector,ieta,iphi,depth);
      //cout << "DEBUG: " << id << " " << cap0 << " " << cap1 << " " << cap2 << " " << cap3 << endl;
      HcalGain * fCondObject = new HcalGain(id, cap0, cap1, cap2, cap3);
      fObject->addValues(*fCondObject);
      delete fCondObject;
    }
    //Always terminate statement
    connection->terminateStatement(stmt);
  } catch (SQLException& e) {
    throw cms::Exception("ReadError") << ::toolbox::toString("Oracle  exception : %s",e.getMessage().c_str()) << std::endl;
  }
  return result;
}


bool HcalDbOmds::getObject (oracle::occi::Connection * connection, 
			    const std::string & fTag, 
			    const std::string & fVersion,
			    const int fSubversion,
			    const std::string & fQuery,
			    HcalGainWidths* fObject) {
  bool result=true;
  if (!fObject) fObject = new HcalGainWidths;
  try {
    oracle::occi::Statement* stmt = connection->createStatement(fQuery);
    stmt->setString(1,fTag);
    stmt->setString(2,fVersion);
    //stmt->setInt(3,fSubversion);

    ResultSet *rs = stmt->executeQuery();

    RooGKCounter _row(1,100);
    _row.setMessage("HCAL channels processed: ");
    _row.setPrintCount(true);
    _row.setNewLine(true);
    while (rs->next()) {
      _row.count();
      float cap0 = rs->getFloat(1);
      float cap1 = rs->getFloat(2);
      float cap2 = rs->getFloat(3);
      float cap3 = rs->getFloat(4);
      int ieta = rs->getInt(5);
      int iphi = rs->getInt(6);
      int depth = rs->getInt(7);
      HcalSubdetector subdetector = get_subdetector(rs->getString(8));
      HcalDetId id(subdetector,ieta,iphi,depth);
      //cout << "DEBUG: " << id << " " << cap0 << " " << cap1 << " " << cap2 << " " << cap3 << endl;
      HcalGainWidth * fCondObject = new HcalGainWidth(id, cap0, cap1, cap2, cap3);
      fObject->addValues(*fCondObject);
      delete fCondObject;
    }
    //Always terminate statement
    connection->terminateStatement(stmt);
  } catch (SQLException& e) {
    throw cms::Exception("ReadError") << ::toolbox::toString("Oracle  exception : %s",e.getMessage().c_str()) << std::endl;
  }
  return result;
}


bool HcalDbOmds::getObject (oracle::occi::Connection * connection, 
			    const std::string & fTag, 
			    const std::string & fVersion,
			    const int fSubversion,
			    const std::string & fQuery,
			    HcalQIEData* fObject) {
  bool result=true;
  if (!fObject) fObject = new HcalQIEData;
  try {
    oracle::occi::Statement* stmt = connection->createStatement(fQuery);
    stmt->setString(1,fTag);
    stmt->setString(2,fVersion);
    //stmt->setInt(3,fSubversion);

    ResultSet *rs = stmt->executeQuery();

    rs->setMaxColumnSize(1,128);
    rs->setMaxColumnSize(2,128);
    rs->setMaxColumnSize(3,128);
    rs->setMaxColumnSize(4,128);
    rs->setMaxColumnSize(5,128);
    rs->setMaxColumnSize(6,128);
    rs->setMaxColumnSize(7,128);
    rs->setMaxColumnSize(8,128);
    rs->setMaxColumnSize(9,128);
    rs->setMaxColumnSize(10,128);
    rs->setMaxColumnSize(11,128);

    RooGKCounter _row(1,100);
    _row.setMessage("HCAL channels processed: ");
    _row.setPrintCount(true);
    _row.setNewLine(true);
    //
    // The query result must be ordered in the following way
    // 1.objectname, ( values: HcalDetId, HcalCalibDetId, HcalTrigTowerDetId, HcalZDCDetId or HcalCastorDetId)
    // 2.subdet, 3.ieta, 4.iphi, 5.depth, 6.type, 7.section, 8.ispositiveeta, 9.sector, 10.module, 11.channel 
    // 13-27. cap0_range0_slope, cap0_range1_slope... 33, 28-43. cap0_range0_offset, cap0_range1_offset...
    //
    while (rs->next()) {
      _row.count();
      DetId id = getId(rs);
      fObject->sort();
      float items[32];
      for (int _i=0; _i!=32; _i++) items[_i] = rs->getFloat(_i+12);
      HcalQIECoder * fCondObject = new HcalQIECoder(id.rawId());
      for (unsigned capid = 0; capid < 4; capid++) {
	for (unsigned range = 0; range < 4; range++) {
	  fCondObject->setSlope (capid, range, items[capid*4+range]);
	}
      }
      for (unsigned capid = 0; capid < 4; capid++) {
	for (unsigned range = 0; range < 4; range++) {
	  fCondObject->setOffset (capid, range, items[16+capid*4+range]);
	}
      }
      fObject->addCoder(*fCondObject);
      delete fCondObject;
    }
    //Always terminate statement
    connection->terminateStatement(stmt);
  } catch (SQLException& e) {
    throw cms::Exception("ReadError") << ::toolbox::toString("Oracle  exception : %s",e.getMessage().c_str()) << std::endl;
  }
  return result;
}


bool HcalDbOmds::getObject (oracle::occi::Connection * connection, 
			    const std::string & fTag, 
			    const std::string & fVersion,
			    const int fSubversion,
			    const std::string & fQuery,
			    HcalCalibrationQIEData* fObject) {
  std::cerr << "NOT IMPLEMENTED!" << std::endl;
  return false;
}


bool HcalDbOmds::getObject (oracle::occi::Connection * connection, 
			    const std::string & fTag, 
			    const std::string & fVersion,
			    const int fSubversion,
			    const std::string & fQuery,
			    HcalElectronicsMap* fObject) {
  bool result=true;
  if (!fObject) fObject = new HcalElectronicsMap;
  try {
    oracle::occi::Statement* stmt = connection->createStatement(fQuery);
    stmt->setString(1,fTag);
    stmt->setString(2,fVersion);
    //stmt->setInt(3,fSubversion);

    ResultSet *rs = stmt->executeQuery();

    // protection agains NULL values
    for (int _i=1; _i!=20; _i++) rs->setMaxColumnSize(_i,128);

    RooGKCounter _row(1,100);
    _row.setMessage("HCAL channels processed: ");
    _row.setPrintCount(true);
    _row.setNewLine(true);
    //
    // The query result must be ordered in the following way
    // 1.objectname, ( values: HcalDetId, HcalCalibDetId, HcalTrigTowerDetId, HcalZDCDetId or HcalCastorDetId)
    // 2.subdet, 3.ieta, 4.iphi, 5.depth, 6.type, 7.section, 8.ispositiveeta, 9.sector, 10.module, 11.channel 
    // 12. i, 13. crate, 14. slot, 15. tb, 16. dcc, 17. spigot, 18. fiber(slb), 19. fiberchan(slbchan)
    //
    while (rs->next()) {
      _row.count();
      DetId id = getId(rs);
      std::string _obj_name = rs->getString(1);
      int crate = rs->getInt(13);
      int slot = rs->getInt(14);
      int dcc = rs->getInt(16);
      int spigot = rs->getInt(17);
      std::string tb = rs->getString(15);
      int top = 1;
      if (tb.find("b")!=std::string::npos) top = 0;
      HcalElectronicsId * fCondObject = 0;
      if (_obj_name.find("HcalTrigTowerDetId")!=std::string::npos){
	int slbCh = rs->getInt(19);
	int slb = rs->getInt(18);
	fCondObject = new HcalElectronicsId(slbCh, slb, spigot, dcc,crate,slot,top);
      }
      else{
	int fiberCh = rs->getInt(19);
	int fiber = rs->getInt(18);
	fCondObject = new HcalElectronicsId(fiberCh, fiber, spigot, dcc);
	fCondObject->setHTR(crate,slot,top);
      }
      if (_obj_name.find("HcalTrigTowerDetId")!=std::string::npos){
	fObject->mapEId2tId(*fCondObject, id);
      }
      else if (_obj_name.find("HcalDetId")!=std::string::npos ||
	       _obj_name.find("HcalCalibDetId")!=std::string::npos ||
	       _obj_name.find("HcalZDCDetId")!=std::string::npos){
	fObject->mapEId2chId(*fCondObject, id);
      }
      else {
	edm::LogWarning("Format Error") << "HcalElectronicsMap-> Unknown subdetector: " 
					<< _obj_name << std::endl; 
      }
      delete fCondObject;
    }
    //Always terminate statement
    connection->terminateStatement(stmt);
  } catch (SQLException& e) {
    throw cms::Exception("ReadError") << ::toolbox::toString("Oracle  exception : %s",e.getMessage().c_str()) << std::endl;
  }
  fObject->sort ();
  return result;
}


bool HcalDbOmds::getObject (oracle::occi::Connection * connection, 
			    const std::string & fTag, 
			    const std::string & fVersion,
			    const int fSubversion,
			    const std::string & fQuery,
			    HcalChannelQuality* fObject) {
  bool result=true;
  if (!fObject) fObject = new HcalChannelQuality;
  try {
    oracle::occi::Statement* stmt = connection->createStatement(fQuery);
    stmt->setString(1,fTag);
    stmt->setString(2,fVersion);
    //stmt->setInt(3,fSubversion);

    ResultSet *rs = stmt->executeQuery();

    RooGKCounter _row(1,100);
    _row.setMessage("HCAL channels processed: ");
    _row.setPrintCount(true);
    _row.setNewLine(true);
    while (rs->next()) {
      _row.count();
      int value = rs->getInt(1);
      int ieta = rs->getInt(2);
      int iphi = rs->getInt(3);
      int depth = rs->getInt(4);
      HcalSubdetector subdetector = get_subdetector(rs->getString(5));
      HcalDetId id(subdetector,ieta,iphi,depth);
      //cout << "DEBUG: " << id << " " << zs << endl;
      HcalChannelStatus * fCondObject = new HcalChannelStatus(id, value);
      fObject->addValues(*fCondObject);
      delete fCondObject;
    }
    //Always terminate statement
    connection->terminateStatement(stmt);
  } catch (SQLException& e) {
    throw cms::Exception("ReadError") << ::toolbox::toString("Oracle  exception : %s",e.getMessage().c_str()) << std::endl;
  }
  return result;
}


bool HcalDbOmds::getObject (oracle::occi::Connection * connection, 
			    const std::string & fTag, 
			    const std::string & fVersion,
			    const int fSubversion,
			    const std::string & fQuery,
			    HcalRespCorrs* fObject) {
  bool result=true;
  if (!fObject) fObject = new HcalRespCorrs;
  try {
    oracle::occi::Statement* stmt = connection->createStatement(fQuery);
    stmt->setString(1,fTag);
    stmt->setString(2,fVersion);
    //stmt->setInt(3,fSubversion);

    ResultSet *rs = stmt->executeQuery();

    RooGKCounter _row(1,100);
    _row.setMessage("HCAL channels processed: ");
    _row.setPrintCount(true);
    _row.setNewLine(true);
    while (rs->next()) {
      _row.count();
      float value = rs->getFloat(1);
      int ieta = rs->getInt(2);
      int iphi = rs->getInt(3);
      int depth = rs->getInt(4);
      HcalSubdetector subdetector = get_subdetector(rs->getString(5));
      HcalDetId id(subdetector,ieta,iphi,depth);
      //cout << "DEBUG: " << id << " " << value << endl;
      HcalRespCorr * fCondObject = new HcalRespCorr(id, value);
      fObject->addValues(*fCondObject);
      delete fCondObject;
    }
    //Always terminate statement
    connection->terminateStatement(stmt);
  } catch (SQLException& e) {
    throw cms::Exception("ReadError") << ::toolbox::toString("Oracle  exception : %s",e.getMessage().c_str()) << std::endl;
  }
  return result;
}

// Oracle database connection ownership is transferred here, DO terminate after use
bool HcalDbOmds::getObject (oracle::occi::Connection * connection, 
			    const std::string & fTag, 
			    const std::string & fVersion,
			    const int fSubversion,
			    const std::string & fQuery,
			    HcalZSThresholds* fObject) {
  bool result=true;
  if (!fObject) fObject = new HcalZSThresholds;
  try {
    oracle::occi::Statement* stmt = connection->createStatement(fQuery);
    stmt->setString(1,fTag);
    stmt->setString(2,fVersion);
    //stmt->setInt(3,fSubversion);

    ResultSet *rs = stmt->executeQuery();

    RooGKCounter _row(1,100);
    _row.setMessage("HCAL channels processed: ");
    _row.setPrintCount(true);
    _row.setNewLine(true);
    while (rs->next()) {
      _row.count();
      int zs = rs->getInt(1);
      int ieta = rs->getInt(2);
      int iphi = rs->getInt(3);
      int depth = rs->getInt(4);
      HcalSubdetector subdetector = get_subdetector(rs->getString(5));
      HcalDetId id(subdetector,ieta,iphi,depth);
      //cout << "DEBUG: " << id << " " << zs << endl;
      HcalZSThreshold * fCondObject = new HcalZSThreshold(id, zs);
      fObject->addValues(*fCondObject);
      delete fCondObject;
    }
    //Always terminate statement
    connection->terminateStatement(stmt);
  } catch (SQLException& e) {
    throw cms::Exception("ReadError") << ::toolbox::toString("Oracle  exception : %s",e.getMessage().c_str()) << std::endl;
  }
  return result;
}


bool HcalDbOmds::getObject (oracle::occi::Connection * connection, 
			    const std::string & fTag, 
			    const std::string & fVersion,
			    const int fSubversion,
			    const std::string & fQuery,
			    HcalL1TriggerObjects* fObject) {
  bool result=true;
  if (!fObject) fObject = new HcalL1TriggerObjects;
  std::string _tag;
  std::string _algo;
  try {
    oracle::occi::Statement* stmt = connection->createStatement(fQuery);
    stmt->setString(1,fTag);
    stmt->setString(2,fVersion);
    //stmt->setInt(3,fSubversion);

    ResultSet *rs = stmt->executeQuery();

    rs->setMaxColumnSize(1,128);
    rs->setMaxColumnSize(2,128);
    rs->setMaxColumnSize(3,128);
    rs->setMaxColumnSize(4,128);
    rs->setMaxColumnSize(5,128);
    rs->setMaxColumnSize(6,128);
    rs->setMaxColumnSize(7,128);
    rs->setMaxColumnSize(8,128);
    rs->setMaxColumnSize(9,128);
    rs->setMaxColumnSize(10,128);
    rs->setMaxColumnSize(11,128);

    RooGKCounter _row(1,100);
    _row.setMessage("HCAL channels processed: ");
    _row.setPrintCount(true);
    _row.setNewLine(true);
    //
    // The query result must be ordered in the following way
    // 1.objectname, ( values: HcalDetId, HcalCalibDetId, HcalTrigTowerDetId, HcalZDCDetId or HcalCastorDetId)
    // 2.subdet, 3.ieta, 4.iphi, 5.depth, 6.type, 7.section, 8.ispositiveeta, 9.sector, 10.module, 11.channel 
    // 12. AVERAGE_PEDESTAL, 13. RESPONSE_CORRECTED_GAIN, 14. FLAG, 15. LUT_TAG_NAME, 16. ALGO_NAME
    //
    bool _is_first_entry=true;
    while (rs->next()) {
      _row.count();
      DetId id = getId(rs);
      float average_pedestal = rs->getFloat(12);
      float response_corrected_gain = rs->getFloat(13);
      int flag = rs->getInt(14);
      if (_is_first_entry){
	_tag = rs->getString(15);
	_algo = rs->getString(16);
	_is_first_entry=false;
      }
      HcalL1TriggerObject * fCondObject = new HcalL1TriggerObject(id, average_pedestal, response_corrected_gain, flag);
      fObject->addValues(*fCondObject);
      delete fCondObject;
    }
    //Always terminate statement
    connection->terminateStatement(stmt);
  } catch (SQLException& e) {
    throw cms::Exception("ReadError") << ::toolbox::toString("Oracle  exception : %s",e.getMessage().c_str()) << std::endl;
  }
  fObject->setTagString(_tag);
  fObject->setAlgoString(_algo);
  return result;
}


bool HcalDbOmds::dumpObject (std::ostream& fOutput, const HcalZSThresholds& fObject) {
  return true;
}


HcalSubdetector HcalDbOmds::get_subdetector( string _det )
{
  HcalSubdetector result;
  if      ( _det.find("HB") != string::npos ) result = HcalBarrel;
  else if ( _det.find("HE") != string::npos ) result = HcalEndcap;
  else if ( _det.find("HF") != string::npos ) result = HcalForward;
  else if ( _det.find("HO") != string::npos ) result = HcalOuter;
  else                                        result = HcalOther;  

  return result;
}

HcalZDCDetId::Section HcalDbOmds::get_zdc_section( string _section )
{
  HcalZDCDetId::Section result;
  if      ( _section.find("Unknown") != string::npos ) result = HcalZDCDetId::Unknown;
  else if ( _section.find("EM") != string::npos )   result = HcalZDCDetId::EM;
  else if ( _section.find("HAD") != string::npos )  result = HcalZDCDetId::HAD;
  else if ( _section.find("LUM") != string::npos ) result = HcalZDCDetId::LUM;
  else                                              result = HcalZDCDetId::Unknown;  
  
  return result;
}
