#ifndef RecoLuminosity_LumiProducer_Lumi2DB_H 
#define RecoLuminosity_LumiProducer_Lumi2DB_H
#include "RelationalAccess/ConnectionService.h"
#include "CoralBase/AttributeList.h"
#include "CoralBase/Attribute.h"
#include "CoralBase/AttributeSpecification.h"
#include "CoralBase/Blob.h"
#include "CoralBase/Exception.h"
#include "RelationalAccess/ISessionProxy.h"
#include "RelationalAccess/ITransaction.h"
#include "RelationalAccess/ITypeConverter.h"
#include "RelationalAccess/ISchema.h"
#include "RelationalAccess/ITable.h"
#include "RelationalAccess/ITableDataEditor.h"
#include "RelationalAccess/IBulkOperation.h"
#include "RecoLuminosity/LumiProducer/interface/LumiRawDataStructures.h"
#include "RecoLuminosity/LumiProducer/interface/DataPipe.h"
#include "RecoLuminosity/LumiProducer/interface/ConstantDef.h"
#include "RecoLuminosity/LumiProducer/interface/LumiNames.h"
#include "RecoLuminosity/LumiProducer/interface/idDealer.h"
#include "RecoLuminosity/LumiProducer/interface/Exception.h"
#include "RecoLuminosity/LumiProducer/interface/DBConfig.h"
#include <iostream>
#include <cstring>
#include <cstdlib>
#include <memory>
#include <map>
#include "TFile.h"
#include "TTree.h"
namespace lumi{
  class Lumi2DB : public DataPipe{
  public:
    const static unsigned int COMMITLSINTERVAL=500; //commit interval in LS,totalrow=nls*(1+nalgo)
    Lumi2DB(const std::string& dest);
    virtual void retrieveData( unsigned int );
    virtual const std::string dataType() const;
    virtual const std::string sourceType() const;
    virtual ~Lumi2DB();
    struct LumiSource{
      unsigned int run;
      unsigned int firstsection;
      char version[8];
      char datestr[9];
    };
    struct PerBXData{
      //int idx;
      float lumivalue;
      float lumierr;
      short lumiquality;
    };
    struct PerLumiData{
      float dtnorm;
      float lhcnorm;
      float instlumi;
      float instlumierror;
      short instlumiquality;
      short lumisectionquality;
      short cmsalive;
      std::string beammode;
      float beamenergy;
      unsigned int cmslsnr;
      unsigned int lumilsnr;
      unsigned int startorbit;
      unsigned int numorbit;
      std::vector<PerBXData> bxET;
      std::vector<PerBXData> bxOCC1;
      std::vector<PerBXData> bxOCC2;
    };
    struct beamData{
      float energy;
      std::string mode;
    };
    typedef std::vector<PerLumiData> LumiResult;
  private:
    void parseSourceString(lumi::Lumi2DB::LumiSource& result)const;
  };//cl Lumi2DB
}//ns lumi

//
//implementation
//
lumi::Lumi2DB::Lumi2DB(const std::string& dest):DataPipe(dest){}

void 
lumi::Lumi2DB::parseSourceString(lumi::Lumi2DB::LumiSource& result)const{
  //parse lumi source file name
  if(m_source.length()==0) throw lumi::Exception("lumi source is not set","parseSourceString","Lumi2DB");
  //std::cout<<"source "<<m_source<<std::endl;
  size_t tempIndex = m_source.find_last_of(".");
  size_t nameLength = m_source.length();
  std::string FileSuffix= m_source.substr(tempIndex+1,nameLength - tempIndex);
  std::string::size_type lastPos=m_source.find_first_not_of("_",0);
  std::string::size_type pos = m_source.find_first_of("_",lastPos);
  std::vector<std::string> pieces;
  while( std::string::npos != pos || std::string::npos != lastPos){
    pieces.push_back(m_source.substr(lastPos,pos-lastPos));
    lastPos=m_source.find_first_not_of("_",pos);
    pos=m_source.find_first_of("_",lastPos);
  }
  if( pieces[1]!="LUMI" || pieces[2]!="RAW" || FileSuffix!="root"){
    throw lumi::Exception("not lumi raw data file CMS_LUMI_RAW","parseSourceString","Lumi2DB");
  }
  std::strcpy(result.datestr,pieces[3].c_str());
  std::strcpy(result.version,pieces[5].c_str());
  //std::cout<<"version : "<<result.version<<std::endl;
  result.run = atoi(pieces[4].c_str());
  //std::cout<<"run : "<<result.run<<std::endl;
  result.firstsection = atoi(pieces[5].c_str());
  //std::cout<<"first section : "<< result.firstsection<<std::endl;
}

void 
lumi::Lumi2DB::retrieveData( unsigned int runnumber){
  lumi::Lumi2DB::LumiResult lumiresult;
  //check filename is in  lumiraw format
  lumi::Lumi2DB::LumiSource filenamecontent;
  try{
    parseSourceString(filenamecontent);
  }catch(const lumi::Exception& er){
    std::cout<<er.what()<<std::endl;
    throw er;
  }
  if( filenamecontent.run!=runnumber ){
    throw lumi::Exception("runnumber in file name does not match requested run number","retrieveData","Lumi2DB");
  }
  TFile* source=TFile::Open(m_source.c_str(),"READ");
  TTree *hlxtree = (TTree*)source->Get("HLXData");
  if(!hlxtree){
    throw lumi::Exception(std::string("non-existing HLXData "),"retrieveData","Lumi2DB");
  }
  //hlxtree->Print();
  std::auto_ptr<HCAL_HLX::LUMI_SECTION> localSection(new HCAL_HLX::LUMI_SECTION);
  HCAL_HLX::LUMI_SECTION_HEADER* lumiheader = &(localSection->hdr);
  HCAL_HLX::LUMI_SUMMARY* lumisummary = &(localSection->lumiSummary);
  HCAL_HLX::LUMI_DETAIL* lumidetail = &(localSection->lumiDetail);
  
  hlxtree->SetBranchAddress("Header.",&lumiheader);
  hlxtree->SetBranchAddress("Summary.",&lumisummary);
  hlxtree->SetBranchAddress("Detail.",&lumidetail);
   
  size_t nentries=hlxtree->GetEntries();
  //source->GetListOfKeys()->Print();
  std::map<unsigned int, Lumi2DB::beamData> dipmap;
  TTree *diptree= (TTree*)source->Get("DIPCombined");
  if(diptree){
    //throw lumi::Exception(std::string("non-existing DIPData "),"retrieveData","Lumi2DB");
    std::auto_ptr<HCAL_HLX::DIP_COMBINED_DATA> dipdata(new HCAL_HLX::DIP_COMBINED_DATA);
    diptree->SetBranchAddress("DIPCombined.",&dipdata);
    size_t ndipentries=diptree->GetEntries();
    for(size_t i=0;i<ndipentries;++i){
      diptree->GetEntry(i);
      beamData b;
      //std::cout<<"Beam Mode : "<<dipdata->beamMode<<"\n";
      //std::cout<<"Beam Energy : "<<dipdata->Energy<<"\n";
      //std::cout<<"sectionUmber : "<<dipdata->sectionNumber<<"\n";
      unsigned int dipls=dipdata->sectionNumber;
      if (std::string(dipdata->beamMode).empty()){
	b.mode="N/A";
      }else{
	b.mode=dipdata->beamMode;
      }
      b.energy=dipdata->Energy;
      dipmap.insert(std::make_pair(dipls,b));
    }
  }else{
    for(size_t i=0;i<nentries;++i){
      beamData b;
      b.mode="N/A";
      b.energy=0.0;
      dipmap.insert(std::make_pair(i,b));
    }
  }
  //diptree->Print();
 
  size_t ncmslumi=0;
  std::cout<<"processing total lumi lumisection "<<nentries<<std::endl;
  //size_t lumisecid=0;
  //unsigned int lumilumisecid=0;
  //runnumber=lumiheader->runNumber;
  //
  //hardcode the first LS is always alive
  //
  for(size_t i=0;i<nentries;++i){
    lumi::Lumi2DB::PerLumiData h;
    h.cmsalive=1;
    hlxtree->GetEntry(i);
    //std::cout<<"live flag "<<lumiheader->bCMSLive <<std::endl;
    if( !lumiheader->bCMSLive && i!=0){
      std::cout<<"\t non-CMS LS "<<lumiheader->sectionNumber<<" ";
      h.cmsalive=0;
    }
    ++ncmslumi;
    
    h.bxET.reserve(lumi::N_BX);
    h.bxOCC1.reserve(lumi::N_BX);
    h.bxOCC2.reserve(lumi::N_BX);
    
    //runnumber=lumiheader->runNumber;
    //if(runnumber!=m_run) throw std::runtime_error(std::string("requested run ")+this->int2str(m_run)+" does not match runnumber in the data header "+this->int2str(runnumber));
    h.lumilsnr=lumiheader->sectionNumber;
    std::map<unsigned int , Lumi2DB::beamData >::iterator beamIt=dipmap.find(h.lumilsnr);
    if ( beamIt!=dipmap.end() ){
      h.beammode=beamIt->second.mode;
      h.beamenergy=beamIt->second.energy;
    }else{
      h.beammode="N/A";
      h.beamenergy=0.0;
    }
    h.startorbit=lumiheader->startOrbit;
    h.numorbit=lumiheader->numOrbits;
    if(h.cmsalive==0){
      h.cmslsnr=0; //the dead ls has cmsls number=0
    }else{
      h.cmslsnr=ncmslumi;//we guess cms lumils
    }
    h.instlumi=lumisummary->InstantLumi;
    //std::cout<<"instant lumi "<<lumisummary->InstantLumi<<std::endl;
    h.instlumierror=lumisummary->InstantLumiErr;
    h.lumisectionquality=lumisummary->InstantLumiQlty;
    h.dtnorm=lumisummary->DeadTimeNormalization;
    h.lhcnorm=lumisummary->LHCNormalization;
    for(size_t i=0;i<lumi::N_BX;++i){
      lumi::Lumi2DB::PerBXData bET;
      lumi::Lumi2DB::PerBXData bOCC1;
      lumi::Lumi2DB::PerBXData bOCC2;
      //bET.idx=i+1;
      bET.lumivalue=lumidetail->ETLumi[i];
      bET.lumierr=lumidetail->ETLumiErr[i];
      bET.lumiquality=lumidetail->ETLumiQlty[i];      
      h.bxET.push_back(bET);

      //bOCC1.idx=i+1;
      bOCC1.lumivalue=lumidetail->OccLumi[0][i];
      bOCC1.lumierr=lumidetail->OccLumiErr[0][i];
      bOCC1.lumiquality=lumidetail->OccLumiQlty[0][i]; 
      h.bxOCC1.push_back(bOCC1);
          
      //bOCC2.idx=i+1;
      bOCC2.lumivalue=lumidetail->OccLumi[1][i];
      bOCC2.lumierr=lumidetail->OccLumiErr[1][i];
      bOCC2.lumiquality=lumidetail->OccLumiQlty[1][i]; 
      h.bxOCC2.push_back(bOCC2);
    }
    lumiresult.push_back(h);
  }
  std::cout<<std::endl;
  coral::ConnectionService* svc=new coral::ConnectionService;
  lumi::DBConfig dbconf(*svc);
  if(!m_authpath.empty()){
    dbconf.setAuthentication(m_authpath);
  }
  coral::ISessionProxy* session=svc->connect(m_dest,coral::Update);
  coral::ITypeConverter& tpc=session->typeConverter();
  tpc.setCppTypeForSqlType("unsigned int","NUMBER(10)");

  unsigned int totallumils=lumiresult.size();
  try{
    coral::AttributeList summaryData;
    summaryData.extend<unsigned long long>("LUMISUMMARY_ID");
    summaryData.extend<unsigned int>("RUNNUM");
    summaryData.extend<unsigned int>("CMSLSNUM");
    summaryData.extend<unsigned int>("LUMILSNUM");
    summaryData.extend<std::string>("LUMIVERSION");
    summaryData.extend<float>("DTNORM");
    summaryData.extend<float>("LHCNORM");
    summaryData.extend<float>("INSTLUMI");
    summaryData.extend<float>("INSTLUMIERROR");
    summaryData.extend<short>("INSTLUMIQUALITY");
    summaryData.extend<short>("LUMISECTIONQUALITY");
    summaryData.extend<short>("CMSALIVE");
    summaryData.extend<unsigned int>("STARTORBIT");
    summaryData.extend<unsigned int>("NUMORBIT");
    summaryData.extend<float>("BEAMENERGY");
    summaryData.extend<std::string>("BEAMSTATUS");

    coral::AttributeList detailData;
    detailData.extend("LUMIDETAIL_ID",typeid(unsigned long long));
    detailData.extend("LUMISUMMARY_ID",typeid(unsigned long long));
    detailData.extend("BXLUMIVALUE",typeid(coral::Blob));
    detailData.extend("BXLUMIERROR",typeid(coral::Blob));
    detailData.extend("BXLUMIQUALITY",typeid(coral::Blob));
    detailData.extend("ALGONAME",typeid(std::string));

    //loop over lumi LS
    unsigned long long& lumisummary_id=summaryData["LUMISUMMARY_ID"].data<unsigned long long>();
    unsigned int& lumirunnum = summaryData["RUNNUM"].data<unsigned int>();
    std::string& lumiversion = summaryData["LUMIVERSION"].data<std::string>();
    float& dtnorm = summaryData["DTNORM"].data<float>();
    float& lhcnorm = summaryData["LHCNORM"].data<float>();
    float& instlumi = summaryData["INSTLUMI"].data<float>();
    float& instlumierror = summaryData["INSTLUMIERROR"].data<float>();
    short& instlumiquality = summaryData["INSTLUMIQUALITY"].data<short>();
    short& lumisectionquality = summaryData["LUMISECTIONQUALITY"].data<short>();
    short& cmsalive = summaryData["CMSALIVE"].data<short>();
    unsigned int& lumilsnr = summaryData["LUMILSNUM"].data<unsigned int>();
    unsigned int& cmslsnr = summaryData["CMSLSNUM"].data<unsigned int>();
    unsigned int& startorbit = summaryData["STARTORBIT"].data<unsigned int>();
    unsigned int& numorbit = summaryData["NUMORBIT"].data<unsigned int>();
    float& beamenergy = summaryData["BEAMENERGY"].data<float>();
    std::string& beamstatus = summaryData["BEAMSTATUS"].data<std::string>();

    unsigned long long& lumidetail_id=detailData["LUMIDETAIL_ID"].data<unsigned long long>();
    unsigned long long& d2lumisummary_id=detailData["LUMISUMMARY_ID"].data<unsigned long long>();
    coral::Blob& bxlumivalue=detailData["BXLUMIVALUE"].data<coral::Blob>();
    coral::Blob& bxlumierror=detailData["BXLUMIERROR"].data<coral::Blob>();
    coral::Blob& bxlumiquality=detailData["BXLUMIQUALITY"].data<coral::Blob>();
    std::string& algoname=detailData["ALGONAME"].data<std::string>();
    
    lumi::Lumi2DB::LumiResult::const_iterator lumiIt;
    lumi::Lumi2DB::LumiResult::const_iterator lumiBeg=lumiresult.begin();
    lumi::Lumi2DB::LumiResult::const_iterator lumiEnd=lumiresult.end();

    coral::IBulkOperation* summaryInserter=0;
    coral::IBulkOperation* detailInserter=0;
    //one loop for ids
    //nested transaction doesn't work with bulk inserter
    std::map< unsigned long long,std::vector<unsigned long long> > idallocationtable;
    std::cout<<"\t allocating total lumisummary ids "<<totallumils<<std::endl; 
    std::cout<<"\t allocating total lumidetail ids "<<totallumils*lumi::N_LUMIALGO<<std::endl; 
    session->transaction().start(false);
    unsigned int lumiindx=0;
    lumi::idDealer idg(session->nominalSchema());
    unsigned long long lumisummaryID = idg.generateNextIDForTable(LumiNames::lumisummaryTableName(),totallumils)-totallumils;
    unsigned long long lumidetailID=idg.generateNextIDForTable(LumiNames::lumidetailTableName(),totallumils*lumi::N_LUMIALGO)-totallumils*lumi::N_LUMIALGO;
    session->transaction().commit();
    for(lumiIt=lumiBeg;lumiIt!=lumiEnd;++lumiIt,++lumiindx,++lumisummaryID){
      std::vector< unsigned long long > allIDs;
      allIDs.reserve(1+lumi::N_LUMIALGO);
      allIDs.push_back(lumisummaryID);
      for( unsigned int j=0; j<lumi::N_LUMIALGO; ++j, ++lumidetailID){
	allIDs.push_back(lumidetailID);
      }
      idallocationtable.insert(std::make_pair(lumiindx,allIDs));
    }
    std::cout<<"\t all ids allocated"<<std::endl; 
    lumiindx=0;
    unsigned int comittedls=0;
    for(lumiIt=lumiBeg;lumiIt!=lumiEnd;++lumiIt,++lumiindx){
      if(!session->transaction().isActive()){ 
	session->transaction().start(false);
	coral::ITable& summarytable=session->nominalSchema().tableHandle(LumiNames::lumisummaryTableName());
	summaryInserter=summarytable.dataEditor().bulkInsert(summaryData,totallumils);
	coral::ITable& detailtable=session->nominalSchema().tableHandle(LumiNames::lumidetailTableName());
	detailInserter=detailtable.dataEditor().bulkInsert(detailData,totallumils*lumi::N_LUMIALGO);    
      }
      lumisummary_id=idallocationtable[lumiindx][0];
      lumirunnum = runnumber;
      lumiversion = std::string(filenamecontent.version);
      dtnorm = lumiIt->dtnorm;
      lhcnorm = lumiIt->lhcnorm;
      instlumi = lumiIt->instlumi;
      instlumierror = lumiIt->instlumierror;
      instlumiquality = lumiIt->instlumiquality;
      lumisectionquality = lumiIt->lumisectionquality;
      cmsalive = lumiIt->cmsalive;
      cmslsnr = lumiIt->cmslsnr;
      
      lumilsnr = lumiIt->lumilsnr;
      startorbit = lumiIt->startorbit;
      numorbit = lumiIt->numorbit;
      beamenergy = lumiIt->beamenergy;
      beamstatus = lumiIt->beammode;
      //fetch a new id value 
      //insert the new row
      summaryInserter->processNextIteration();
      summaryInserter->flush();
      unsigned int algoindx=1;
      for( unsigned int j=0; j<lumi:: N_LUMIALGO; ++j,++algoindx ){
	d2lumisummary_id=idallocationtable[lumiindx].at(0);
	lumidetail_id=idallocationtable[lumiindx].at(algoindx);
	std::vector<PerBXData>::const_iterator bxIt;
	std::vector<PerBXData>::const_iterator bxBeg;
	std::vector<PerBXData>::const_iterator bxEnd;
	if(j==0) {
	  algoname=std::string("ET");
	  bxBeg=lumiIt->bxET.begin();
	  bxEnd=lumiIt->bxET.end();
	}
	if(j==1) {
	  algoname=std::string("OCC1");
	  bxBeg=lumiIt->bxOCC1.begin();
	  bxEnd=lumiIt->bxOCC1.end();
	}
	if(j==2) {
	  algoname=std::string("OCC2");
	  bxBeg=lumiIt->bxOCC2.begin();
	  bxEnd=lumiIt->bxOCC2.end();
	}
	float lumivalue[lumi::N_BX];
	std::memset((void*)&lumivalue,0,sizeof(float)*lumi::N_BX );
	float lumierror[lumi::N_BX];
	std::memset((void*)&lumierror,0,sizeof(float)*lumi::N_BX );
	int lumiquality[lumi::N_BX];
	std::memset((void*)&lumiquality,0,sizeof(short)*lumi::N_BX );
	bxlumivalue.resize(sizeof(float)*lumi::N_BX);
	bxlumierror.resize(sizeof(float)*lumi::N_BX);
	bxlumiquality.resize(sizeof(short)*lumi::N_BX);
	void* bxlumivalueStartAddress=bxlumivalue.startingAddress();
	void* bxlumierrorStartAddress=bxlumierror.startingAddress();
	void* bxlumiqualityStartAddress=bxlumiquality.startingAddress();
	unsigned int k=0;
	for( bxIt=bxBeg;bxIt!=bxEnd;++bxIt,++k  ){	    
	  lumivalue[k]=bxIt->lumivalue;
	  lumierror[k]=bxIt->lumierr;
	  lumiquality[k]=bxIt->lumiquality;
	}
	std::memmove(bxlumivalueStartAddress,lumivalue,sizeof(float)*lumi::N_BX);
	std::memmove(bxlumierrorStartAddress,lumierror,sizeof(float)*lumi::N_BX);
	std::memmove(bxlumiqualityStartAddress,lumiquality,sizeof(short)*lumi::N_BX);
	detailInserter->processNextIteration();
      }
      detailInserter->flush();
      ++comittedls;
      if(comittedls==Lumi2DB::COMMITLSINTERVAL){
	std::cout<<"\t committing in LS chunck "<<comittedls<<std::endl; 
	delete summaryInserter;
	summaryInserter=0;
	delete detailInserter;
	detailInserter=0;
	session->transaction().commit();
	comittedls=0;
	std::cout<<"\t committed "<<std::endl; 
      }else if( lumiindx==(totallumils-1) ){
	std::cout<<"\t committing at the end"<<std::endl; 
	delete summaryInserter; summaryInserter=0;
	delete detailInserter; detailInserter=0;
	session->transaction().commit();
	std::cout<<"\t done"<<std::endl; 
      }
    }
  }catch( const coral::Exception& er){
    session->transaction().rollback();
    delete session;
    delete svc;
    throw er;
  }
  delete session;
  delete svc;
}
const std::string lumi::Lumi2DB::dataType() const{
  return "LUMI";
}
const std::string lumi::Lumi2DB::sourceType() const{
  return "DB";
}
lumi::Lumi2DB::~Lumi2DB(){}

#include "RecoLuminosity/LumiProducer/interface/DataPipeFactory.h"
DEFINE_EDM_PLUGIN(lumi::DataPipeFactory,lumi::Lumi2DB,"Lumi2DB");
#endif
