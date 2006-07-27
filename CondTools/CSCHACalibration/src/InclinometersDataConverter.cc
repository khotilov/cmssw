#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CondCore/DBCommon/interface/Exception.h"
#include "CondCore/DBOutputService/interface/PoolDBOutputService.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "CondCore/DBOutputService/interface/PoolDBOutputService.h"
#include "CondFormats/OptAlignObjects/interface/Inclinometers.h"
#include "CondTools/Utilities/interface/CSVDataLineParser.h"
#include "CondTools/Utilities/interface/CSVHeaderLineParser.h"
#include "CondTools/Utilities/interface/CSVBlankLineParser.h"
#include "InclinometersDataConverter.h"
#include <fstream>
#include <iostream>
InclinometersDataConverter::InclinometersDataConverter(const edm::ParameterSet& iConfig):m_inFileName( iConfig.getUntrackedParameter< std::string >("inputFile") ){}
void InclinometersDataConverter::endJob()
{
  std::ifstream myfile(m_inFileName.c_str());
  if ( !myfile.is_open() ) throw cms::Exception("unable to open file");
  Inclinometers* myobj=new Inclinometers;
  std::vector<std::string> fieldNames,fieldTypes;
  CSVHeaderLineParser headerParser;
  CSVDataLineParser dataParser;
  int counter=0;
  while (! myfile.eof() ){
    std::string line;
    std::getline (myfile,line);
    CSVBlankLineParser blank;
    if(blank.isBlank(line)){
      continue;
    }
    if(counter<2){//two lines of header
      if(counter==0) {
	if(!headerParser.parse(line)) {
	  throw cms::Exception("unable to parse header: ")<<line;
	}
	fieldNames=headerParser.result();
      }
      if(counter==1) {
	if(!headerParser.parse(line)) {
	}
	fieldTypes=headerParser.result();
	int idx=0;
	for(std::vector<std::string>::iterator it=fieldTypes.begin(); it!=fieldTypes.end(); ++it, ++idx){
	  std::cout<<fieldNames[idx]<<":"<<*it<<std::endl;
	  m_fieldMap.push_back(fieldNames[idx],*it);
	}
      }
      ++counter;
      continue;
    }
    if(!dataParser.parse(line)) {
      throw cms::Exception("unable to parse data")<<line;
    }
    std::vector<boost::any> result=dataParser.result();
    Inclinometers::Item data;
    int idx=0;
    for(std::vector<boost::any>::iterator it=result.begin(); 
	it!=result.end(); ++it, ++idx){
      std::string fieldName=m_fieldMap.fieldName(idx);
      if(fieldName=="Sensor_type"){
	//std::cout<<"fieldName "<<fieldName<<" field type "<<m_fieldMap.fieldTypeName(idx)<<std::endl;
	if( m_fieldMap.fieldType(idx)!= typeid(std::string) ) throw cond::Exception("unexpected type");
	data.Sensor_type=boost::any_cast<std::string>(*it);
      }
      if(fieldName=="Sensor_number"){
	if( m_fieldMap.fieldType(idx) != typeid(int) ) throw cond::Exception("unexpected type");
	data.Sensor_number=boost::any_cast<int>(*it);
      }
      if(fieldName=="ME_layer"){
	if( m_fieldMap.fieldType(idx) != typeid(std::string) ) throw cond::Exception("unexpected type");
	data.ME_layer=boost::any_cast<std::string>(*it);
      }
      if(fieldName=="Logical_Alignment_Name"){
	if( m_fieldMap.fieldType(idx) != typeid(std::string) ) throw cond::Exception("unexpected type");
	data.Logical_Alignment_Name=boost::any_cast<std::string>(*it);
      }
      if(fieldName=="CERN_Designator"){
	if( m_fieldMap.fieldType(idx) != typeid(std::string) ) throw cond::Exception("unexpected type");
	data.CERN_Designator=boost::any_cast<std::string>(*it);
      }
      if(fieldName=="CERN_Barcode"){
	if( m_fieldMap.fieldType(idx) != typeid(std::string) ) throw cond::Exception("unexpected type");
	data.CERN_Barcode=boost::any_cast<std::string>(*it);
      }
      if(fieldName=="Inclination_Direction"){
	if( m_fieldMap.fieldType(idx) != typeid(std::string) ) throw cond::Exception("unexpected type");
	data.Inclination_Direction=boost::any_cast<std::string>(*it);
      }
      if(fieldName=="Abs_Slope"){
	if( m_fieldMap.fieldType(idx) != typeid(float) ) throw cond::Exception("unexpected type");
	data.Abs_Slope=(float)boost::any_cast<double>(*it);
      }
      if(fieldName=="Abs_Slope_Error"){
	if( m_fieldMap.fieldType(idx) != typeid(float) ) throw cond::Exception("unexpected type");
	data.Abs_Slope_Error=(float)boost::any_cast<double>(*it);
      }
      if(fieldName=="Norm_Slope"){
	if( m_fieldMap.fieldType(idx) != typeid(float) ) throw cond::Exception("unexpected type");
	data.Norm_Slope=(float)boost::any_cast<double>(*it);
      }
      if(fieldName=="Norm_Slope_Error"){
	if( m_fieldMap.fieldType(idx) != typeid(float) ) throw cond::Exception("unexpected type");
	data.Norm_Slope_Error=(float)boost::any_cast<double>(*it);
      }
      if(fieldName=="Abs_Intercept"){
	if( m_fieldMap.fieldType(idx) != typeid(float) ) throw cond::Exception("unexpected type");
	data.Abs_Intercept=(float)boost::any_cast<double>(*it);
      }
      if(fieldName=="Abs_Intercept_Error"){
	if( m_fieldMap.fieldType(idx) != typeid(float) ) throw cond::Exception("unexpected type");
	data.Abs_Intercept_Error=(float)boost::any_cast<double>(*it);
      }
      if(fieldName=="Norm_Intercept"){
	if( m_fieldMap.fieldType(idx) != typeid(float) ) throw cond::Exception("unexpected type");
	data.Norm_Intercept=(float)boost::any_cast<double>(*it);
      }
      if(fieldName=="Norm_Intercept_Error"){
	if( m_fieldMap.fieldType(idx) != typeid(float) ) throw cond::Exception("unexpected type");
	data.Norm_Intercept_Error=(float)boost::any_cast<double>(*it);
      }
      if(fieldName=="Shifts_due_to_shims_etc."){
	if( m_fieldMap.fieldType(idx) != typeid(float) ) throw cond::Exception("unexpected type");
	data.Shifts_due_to_shims_etc=(float)boost::any_cast<double>(*it);
      }
    }
    myobj->m_inclinometers.push_back(data);
    ++counter;
  }
  myfile.close();
  edm::Service<cond::service::PoolDBOutputService> mydbservice;
  if( mydbservice.isAvailable() ){
    try{
      size_t callbackToken=mydbservice->callbackToken("Inclinometers");
      mydbservice->newValidityForNewPayload<Inclinometers>(myobj,mydbservice->endOfTime(),callbackToken);
    }catch(const cond::Exception& er){
      std::cout<<er.what()<<std::endl;
    }catch(const std::exception& er){
      std::cout<<"caught std::exception "<<er.what()<<std::endl;
    }catch(...){
      std::cout<<"Funny error"<<std::endl;
    }
  }
}
