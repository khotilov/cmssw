// -*- C++ -*-
//
// Package:    LumiProducer
// Class:      LumiCorrectionSource
// 
/**\class LumiCorrectionSource LumiCorrectionSource.cc RecoLuminosity/LumiProducer/src/LumiCorrectionSource.cc
Description: A essource/esproducer for lumi correction factor and run parameters needed to deduce the corrections
      Author: Zhen Xie
*/
// $Id: LumiCorrectionSource.cc,v 1.8 2012/08/21 14:18:03 xiezhen Exp $

//#include <memory>
//#include "boost/shared_ptr.hpp"
#include "FWCore/Framework/interface/SourceFactory.h"
#include "FWCore/Framework/interface/ESProducer.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/EventSetup.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/LuminosityBlock.h"
#include "FWCore/Framework/interface/Run.h"
#include "DataFormats/Provenance/interface/LuminosityBlockID.h"
#include "FWCore/ServiceRegistry/interface/Service.h"

#include "FWCore/Framework/interface/IOVSyncValue.h"

#include "CoralBase/Exception.h"
#include "CoralBase/AttributeList.h"
#include "CoralBase/Attribute.h"
#include "CoralBase/AttributeSpecification.h"
#include "CoralBase/Exception.h"
#include "RelationalAccess/ISessionProxy.h"
#include "RelationalAccess/ITransaction.h"
#include "RelationalAccess/AccessMode.h"
#include "RelationalAccess/ITypeConverter.h"
#include "RelationalAccess/IQuery.h"
#include "RelationalAccess/ICursor.h"
#include "RelationalAccess/ISchema.h"
#include "RelationalAccess/ITable.h"
#include "RecoLuminosity/LumiProducer/interface/DBService.h"
#include "RecoLuminosity/LumiProducer/interface/Exception.h"
#include "RecoLuminosity/LumiProducer/interface/ConstantDef.h"
#include "RecoLuminosity/LumiProducer/interface/LumiCorrectionParam.h"
#include "RecoLuminosity/LumiProducer/interface/LumiCorrectionParamRcd.h"
#include "RecoLuminosity/LumiProducer/interface/RevisionDML.h"
#include "RecoLuminosity/LumiProducer/interface/NormDML.h"
#include "LumiCorrectionSource.h"
#include <iostream>
#include <sstream>
#include <string>
#include <memory>
#include <algorithm>
#include <vector>
#include <cstring>
#include <iterator>
#include <boost/foreach.hpp>
#include <boost/tokenizer.hpp>
#include <xercesc/dom/DOM.hpp>
#include <xercesc/parsers/XercesDOMParser.hpp>
#include <xercesc/util/PlatformUtils.hpp>
#include <xercesc/util/XMLString.hpp>

#include "boost/filesystem/path.hpp"
#include "boost/filesystem/operations.hpp"
std::string 
LumiCorrectionSource::x2s(const XMLCh *toTranscode)const{
  std::string tmp(xercesc::XMLString::transcode(toTranscode));
  return tmp;
}

XMLCh*  
LumiCorrectionSource::s2x( const std::string& temp )const{
  XMLCh* buff = xercesc::XMLString::transcode(temp.c_str());    
  return  buff;
}

std::string
LumiCorrectionSource::toParentString(const xercesc::DOMNode &nodeToConvert)const{
  std::ostringstream oss;
  xercesc::DOMNodeList *childList = nodeToConvert.getChildNodes();

  unsigned int numNodes = childList->getLength ();
  for (unsigned int i = 0; i < numNodes; ++i){
    xercesc::DOMNode *childNode = childList->item(i);
    if (childNode->getNodeType() != xercesc::DOMNode::ELEMENT_NODE)
      continue;
    xercesc::DOMElement *child = static_cast < xercesc::DOMElement *> (childNode);
    xercesc::DOMNamedNodeMap *attributes = child->getAttributes();
    unsigned int numAttributes = attributes->getLength ();
    for (unsigned int j = 0; j < numAttributes; ++j){
      xercesc::DOMNode *attributeNode = attributes->item(j);
      if (attributeNode->getNodeType() != xercesc::DOMNode::ATTRIBUTE_NODE)
	continue;
      xercesc::DOMAttr *attribute = static_cast <xercesc::DOMAttr *> (attributeNode);
      
      oss << "(" << x2s(child->getTagName()) << 
	x2s(attribute->getName()) << "=" << 
	x2s(attribute->getValue()) << ")";
    }
  }
  return oss.str();
}

const std::string
LumiCorrectionSource::servletTranslation(const std::string& servlet) const{
  std::string frontierConnect;
  std::string realconnect;
  xercesc::XMLPlatformUtils::Initialize();  
  std::auto_ptr< xercesc::XercesDOMParser > parser(new xercesc::XercesDOMParser);
  try{
    parser->setValidationScheme(xercesc::XercesDOMParser::Val_Auto);
    parser->setDoNamespaces(false);
    parser->parse(m_siteconfpath.c_str());
    xercesc::DOMDocument* doc=parser->getDocument();
    if(!doc){
      return "";
    }
   
    xercesc::DOMNodeList *frontierConnectList=doc->getElementsByTagName(s2x("frontier-connect"));
    if (frontierConnectList->getLength()>0){
      xercesc::DOMElement *frontierConnectElement=static_cast < xercesc::DOMElement *> (frontierConnectList->item (0));
      frontierConnect = toParentString(*frontierConnectElement);
    }
    // Replace the last component of every "serverurl=" piece (up to the
    //   next close-paren) with the servlet
    std::string::size_type nextparen = 0;
    std::string::size_type serverurl, lastslash;
    std::string complexstr = "";
    while ((serverurl=frontierConnect.find("(serverurl=", nextparen)) != std::string::npos){
      realconnect.append(frontierConnect, nextparen, serverurl - nextparen);
      nextparen=frontierConnect.find(')', serverurl);
      lastslash=frontierConnect.rfind('/', nextparen);
      realconnect.append(frontierConnect,serverurl,lastslash-serverurl+1);
      realconnect.append(servlet);
    }
    realconnect.append(frontierConnect, nextparen,frontierConnect.length()-nextparen);
  }catch(xercesc::DOMException &e){
  }
  return realconnect;
}

LumiCorrectionSource::LumiCorrectionSource(const edm::ParameterSet& iConfig):m_connectStr(""),m_datatag(""),m_normtag(""),m_paramcachedrun(0),m_cachesize(0){
  setWhatProduced(this,&LumiCorrectionSource::produceLumiCorrectionParam);
  findingRecord<LumiCorrectionParamRcd>();
  std::string connectStr=iConfig.getParameter<std::string>("connect");
  m_datatag=iConfig.getUntrackedParameter<std::string>("datatag","");
  m_normtag=iConfig.getUntrackedParameter<std::string>("normtag","");
  m_cachesize=iConfig.getUntrackedParameter<unsigned int>("ncacheEntries",3);
  const std::string fproto("frontier://");
  //test if need frontier servlet site-local translation  
  if(connectStr.substr(0,fproto.length())==fproto){
    std::string::size_type startservlet=fproto.length();
    std::string::size_type endservlet=connectStr.find("(",startservlet);
    if(endservlet==std::string::npos){
      endservlet=connectStr.rfind('/',connectStr.length());
    }
    std::string servlet=connectStr.substr(startservlet,endservlet-startservlet);
    if( (servlet !="")&& (servlet.find_first_of(":/)[]")==std::string::npos)){
      if(servlet=="cms_conditions_data") servlet="";
      
      std::string siteconfpath=iConfig.getUntrackedParameter<std::string>("siteconfpath","");
      if(siteconfpath.length()==0){
	std::string url=(boost::filesystem::path("SITECONF")/boost::filesystem::path("local")/boost::filesystem::path("JobConfig")/boost::filesystem::path("site-local-config.xml")).string();
	char * tmp = getenv ("CMS_PATH");
	if(tmp){
	  m_siteconfpath = (boost::filesystem::path(tmp)/boost::filesystem::path(url)).string();
	}
      }else{
	if(!boost::filesystem::exists(boost::filesystem::path(siteconfpath))){
	  throw cms::Exception("Non existing path ")<<siteconfpath;
	}
	m_siteconfpath = (boost::filesystem::path(siteconfpath)/boost::filesystem::path("site-local-config.xml")).string();
      }
      //std::cout<<"servlet : "<<servlet<<std::endl;
      m_connectStr=fproto+servletTranslation(servlet)+connectStr.substr(endservlet);
    }else{
      m_connectStr=connectStr;
    }
  }else{
    m_connectStr=connectStr;
  }
}

LumiCorrectionSource::ReturnParamType
LumiCorrectionSource::produceLumiCorrectionParam(const LumiCorrectionParamRcd&)  
{ 
  unsigned int currentrun=m_pcurrentTime->eventID().run();
  if(currentrun==0||currentrun==4294967295){ 
    return  boost::shared_ptr<LumiCorrectionParam>(new LumiCorrectionParam());
  }
  if(m_paramcachedrun!=currentrun){//i'm in a new run
    fillparamcache(currentrun);//fill cache
  }else{ //i'm in an old run
    if(m_paramcache.find(currentrun)==m_paramcache.end()){//i'm not cached 
      fillparamcache(currentrun);// 
    }
  }
  if(m_paramcache.empty()){
    return boost::shared_ptr<LumiCorrectionParam>(new LumiCorrectionParam());
  }
  m_paramresult=m_paramcache[currentrun];
  if(m_paramresult.get()==0){
    return boost::shared_ptr<LumiCorrectionParam>(new LumiCorrectionParam());
  }
  return m_paramresult;
}

void 
LumiCorrectionSource::setIntervalFor( 
    const edm::eventsetup::EventSetupRecordKey& iKey, 
    const edm::IOVSyncValue& iTime, 
    edm::ValidityInterval& oValidity ) {
  m_pcurrentTime=&iTime;
  oValidity.setFirst(iTime);
  oValidity.setLast(iTime);
}

void
LumiCorrectionSource::fillparamcache(unsigned int runnumber){
  m_paramcache.clear();
  m_paramcachedrun=runnumber;
  edm::Service<lumi::service::DBService> mydbservice;
  if( !mydbservice.isAvailable() ){
    throw cms::Exception("Non existing service lumi::service::DBService");
  }
  coral::ISessionProxy* session=mydbservice->connectReadOnly(m_connectStr);
  coral::ITypeConverter& tconverter=session->typeConverter();
  tconverter.setCppTypeForSqlType(std::string("float"),std::string("FLOAT(63)"));
  tconverter.setCppTypeForSqlType(std::string("unsigned int"),std::string("NUMBER(10)"));
  tconverter.setCppTypeForSqlType(std::string("unsigned short"),std::string("NUMBER(1)"));
  boost::shared_ptr<LumiCorrectionParam> result(new LumiCorrectionParam(LumiCorrectionParam::HF));
  try{
    session->transaction().start(true);
    coral::ISchema& schema=session->nominalSchema();
    lumi::RevisionDML dml;
    unsigned long long tagid=dml.currentHFDataTagId(schema);//get datatag id
    lumi::RevisionDML::DataID dataid=dml.dataIDForRun(schema,runnumber,tagid);//get data id
    unsigned int lumiid=dataid.lumi_id;
    if(lumiid==0){
      result->setNBX(0);
      m_paramcache.insert(std::make_pair(runnumber,result));
      session->transaction().commit();
      mydbservice->disconnect(session);
      return;
    }

    coral::AttributeList lumidataBindVariables;
    lumidataBindVariables.extend("dataid",typeid(unsigned long long));
    lumidataBindVariables["dataid"].data<unsigned long long>()=lumiid;   
    std::string conditionStr("DATA_ID=:dataid");
    coral::AttributeList lumiparamOutput;
    lumiparamOutput.extend("NCOLLIDINGBUNCHES",typeid(unsigned int));
    coral::IQuery* lumiparamQuery=schema.newQuery();
    lumiparamQuery->addToTableList(std::string("LUMIDATA"));
    lumiparamQuery->setCondition(conditionStr,lumidataBindVariables);
    coral::ICursor& lumiparamcursor=lumiparamQuery->execute();
    unsigned int ncollidingbx=0;
    while( lumiparamcursor.next() ){
      const coral::AttributeList& row=lumiparamcursor.currentRow();
      if(!row["NCOLLIDINGBUNCHES"].isNull()){
	ncollidingbx=row["NCOLLIDINGBUNCHES"].data<unsigned int>();
      }
      result->setNBX(ncollidingbx);
    }
    delete lumiparamQuery;
    lumi::NormDML normdml;
    unsigned long long normid=0;
    std::map<std::string,unsigned long long> normidmap;
    if (m_normtag.empty()){
      normdml.normIdByType(schema,normidmap,lumi::NormDML::HF,true);
      m_normtag=normidmap.begin()->first;
      normid=normidmap.begin()->second;
    }else{
      normid=normdml.normIdByName(schema,m_normtag);
    }

    std::map< unsigned int,lumi::NormDML::normData > normDataMap;
    normdml.normById(schema,normid,normDataMap); 
    session->transaction().commit();
    std::map< unsigned int,lumi::NormDML::normData >::iterator normIt=normDataMap.upper_bound(runnumber);
    --normIt;
    result->setNormtag(m_normtag);
    result->setcorrFunc(normIt->second.corrfunc);
    result->setnonlinearCoeff(normIt->second.coefficientmap);
    result->setafterglows(normIt->second.afterglows);
    result->setdescription(normIt->second.amodetag,normIt->second.beamegev);
    m_paramcache.insert(std::make_pair(runnumber,result));
  }catch(const coral::Exception& er){
    session->transaction().rollback();
    mydbservice->disconnect(session);
    throw cms::Exception("DatabaseError ")<<er.what();
  }
  mydbservice->disconnect(session);
}

LumiCorrectionSource::~LumiCorrectionSource(){}
//define this as a plug-in
DEFINE_FWK_EVENTSETUP_SOURCE(LumiCorrectionSource);
