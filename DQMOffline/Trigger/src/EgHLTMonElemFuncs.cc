#include "DQMOffline/Trigger/interface/EgHLTMonElemFuncs.h"

#include "DQMOffline/Trigger/interface/EgHLTMonElemMgrEBEE.h"
#include "DQMOffline/Trigger/interface/EgHLTEgCutCodes.h"
#include "DQMOffline/Trigger/interface/EgHLTDQMCut.h"
#include "DQMOffline/Trigger/interface/EgHLTMonElemWithCutEBEE.h"
#include "DQMOffline/Trigger/interface/EgHLTCutMasks.h"

#include <boost/algorithm/string.hpp>

using namespace egHLT;

void MonElemFuncs::initStdEleHists(std::vector<MonElemManagerBase<OffEle>*>& histVec,const std::string& baseName,const BinData& bins)
{
  addStdHist<OffEle,float>(histVec,baseName+"_et",baseName+" E_{T};E_{T} (GeV)",bins.et,&OffEle::et); 
  addStdHist<OffEle,float>(histVec,baseName+"_etHigh",baseName+" E_{T};E_{T} (GeV)",bins.etHigh,&OffEle::et);
  addStdHist<OffEle,float>(histVec,baseName+"_eta",baseName+" #eta;#eta",bins.eta,&OffEle::detEta);		
  addStdHist<OffEle,float>(histVec,baseName+"_phi",baseName+" #phi;#phi (rad)",bins.phi,&OffEle::phi);
  addStdHist<OffEle,int>(histVec,baseName+"_charge",baseName+" Charge; charge",bins.charge,&OffEle::charge);
  
  addStdHist<OffEle,float>(histVec,baseName+"_hOverE",baseName+" H/E; H/E",bins.hOverE,&OffEle::hOverE);
  addStdHist<OffEle,float>(histVec,baseName+"_dPhiIn",baseName+" #Delta #phi_{in}; #Delta #phi_{in}",bins.dEtaIn,&OffEle::dPhiIn);
  addStdHist<OffEle,float>(histVec,baseName+"_dEtaIn",baseName+" #Delta #eta_{in}; #Delta #eta_{in}",bins.dPhiIn,&OffEle::dEtaIn);
  addStdHist<OffEle,float>(histVec,baseName+"_sigmaIEtaIEta",baseName+"#sigma_{i#etai#eta}; #sigma_{i#etai#eta}",bins.sigEtaEta,&OffEle::scSigmaIEtaIEta);  
  addStdHist<OffEle,float>(histVec,baseName+"_e2x5Over5x5",baseName+"E^{2x5}/E^{5x5}; E^{2x5}/E^{5x5}",bins.e2x5,&OffEle::scE2x5MaxOver5x5);
  addStdHist<OffEle,float>(histVec,baseName+"_e1x5Over5x5",baseName+"E^{1x5}/E^{5x5}; E^{1x5}/E^{5x5}",bins.e1x5,&OffEle::scE1x5Over5x5);
  addStdHist<OffEle,float>(histVec,baseName+"_isolEM",baseName+"Isol EM; Isol EM (GeV)",bins.isolEm,&OffEle::isolEm); 
  addStdHist<OffEle,float>(histVec,baseName+"_isolHad",baseName+"Isol Had; Isol Had (GeV)",bins.isolHad,&OffEle::isolHad);
  addStdHist<OffEle,float>(histVec,baseName+"_isolPtTrks",baseName+"Isol Pt Trks; Isol Pt Tracks (GeV/c)",bins.isolPtTrks,&OffEle::isolPtTrks); 
  addStdHist<OffEle,float>(histVec,baseName+"_hltIsolTrksEle",baseName+"HLT Ele Isol Trks; HLT Ele Iso Tracks (GeV/c)",bins.isolPtTrks,&OffEle::hltIsolTrksEle);  
  addStdHist<OffEle,float>(histVec,baseName+"_hltIsolTrksPho",baseName+"HLT Pho Isol Trks; HLT Pho Iso Tracks (GeV/c)",bins.isolPtTrks,&OffEle::hltIsolTrksPho); 
  addStdHist<OffEle,float>(histVec,baseName+"_hltIsolHad",baseName+"HLT Isol Had; HLT Isol Had (GeV)",bins.isolHad,&OffEle::hltIsolHad);
  
  histVec.push_back(new MonElemManager2D<OffEle,float,float>(baseName+"_etaVsPhi",
							     baseName+" #eta vs #phi;#eta;#phi (rad)",
							     bins.etaVsPhi.nrX,bins.etaVsPhi.xMin,bins.etaVsPhi.xMax,
							     bins.etaVsPhi.nrY,bins.etaVsPhi.yMin,bins.etaVsPhi.yMax,
							     &OffEle::detEta,&OffEle::phi));
 
}

void MonElemFuncs::initStdPhoHists(std::vector<MonElemManagerBase<OffPho>*>& histVec,const std::string& baseName,const BinData& bins)
{
  addStdHist<OffPho,float>(histVec,baseName+"_et",baseName+" E_{T};E_{T} (GeV)",bins.et,&OffPho::et); 
  addStdHist<OffPho,float>(histVec,baseName+"_etHigh",baseName+" E_{T};E_{T} (GeV)",bins.etHigh,&OffPho::et);
  addStdHist<OffPho,float>(histVec,baseName+"_eta",baseName+" #eta;#eta",bins.eta,&OffPho::detEta);		
  addStdHist<OffPho,float>(histVec,baseName+"_phi",baseName+" #phi;#phi (rad)",bins.phi,&OffPho::phi);
  
  addStdHist<OffPho,float>(histVec,baseName+"_hOverE",baseName+" H/E; H/E",bins.hOverE,&OffPho::hOverE);
  addStdHist<OffPho,float>(histVec,baseName+"_r9",baseName+" R9 ; R9",bins.r9,&OffPho::r9);
  addStdHist<OffPho,float>(histVec,baseName+"_sigmaIEtaIEta",baseName+"#sigma_{i#etai#eta}; #sigma_{i#etai#eta}",bins.sigEtaEta,&OffPho::scSigmaIEtaIEta);  
  addStdHist<OffPho,float>(histVec,baseName+"_e2x5Over5x5",baseName+"E^{2x5}/E^{5x5}; E^{2x5}/E^{5x5}",bins.e2x5,&OffPho::scE2x5MaxOver5x5);
  addStdHist<OffPho,float>(histVec,baseName+"_e1x5Over5x5",baseName+"E^{1x5}/E^{5x5}; E^{1x5}/E^{5x5}",bins.e1x5,&OffPho::scE1x5Over5x5);
  addStdHist<OffPho,float>(histVec,baseName+"_isolEM",baseName+"Isol EM; Isol EM (GeV)",bins.isolEm,&OffPho::isolEm); 
  addStdHist<OffPho,float>(histVec,baseName+"_isolHad",baseName+"Isol Had; Isol Had (GeV)",bins.isolHad,&OffPho::isolHad);
  addStdHist<OffPho,float>(histVec,baseName+"_isolPtTrks",baseName+"Isol Pt Trks; Isol Pt Tracks (GeV/c)",bins.isolPtTrks,&OffPho::isolPtTrks);  
  addStdHist<OffPho,int>(histVec,baseName+"_isolNrTrks",baseName+"Isol Nr Trks; Isol Nr Tracks",bins.isolNrTrks,&OffPho::isolNrTrks); 
  addStdHist<OffPho,float>(histVec,baseName+"_hltIsolTrks",baseName+"HLT Isol Trks; HLT Iso Tracks (GeV/c)",bins.isolPtTrks,&OffPho::hltIsolTrks); 
  addStdHist<OffPho,float>(histVec,baseName+"_hltIsolHad",baseName+"HLT Isol Had; HLT Isol Had (GeV)",bins.isolPtTrks,&OffPho::hltIsolHad);
  
  histVec.push_back(new MonElemManager2D<OffPho,float,float>(baseName+"_etaVsPhi",
							     baseName+" #eta vs #phi;#eta;#phi (rad)",
							     bins.etaVsPhi.nrX,bins.etaVsPhi.xMin,bins.etaVsPhi.xMax,
							     bins.etaVsPhi.nrY,bins.etaVsPhi.yMin,bins.etaVsPhi.yMax,
							     &OffPho::detEta,&OffPho::phi));
}
 
void MonElemFuncs::initStdEffHists(std::vector<MonElemWithCutBase<OffEle>*>& histVec,const std::string& baseName,const BinData::Data1D& bins,float (OffEle::*vsVarFunc)()const,const CutMasks& masks)
{
  initStdEffHists(histVec,baseName,bins.nr,bins.min,bins.max,vsVarFunc,masks);
}
  
void MonElemFuncs::initStdEffHists(std::vector<MonElemWithCutBase<OffPho>*>& histVec,const std::string& baseName,const BinData::Data1D& bins,float (OffPho::*vsVarFunc)()const,const CutMasks& masks)
{
  initStdEffHists(histVec,baseName,bins.nr,bins.min,bins.max,vsVarFunc,masks);
}

void MonElemFuncs::initStdEffHists(std::vector<MonElemWithCutBase<OffEle>*>& histVec,const std::string& baseName,int nrBins,double xMin,double xMax,float (OffEle::*vsVarFunc)()const,const CutMasks& masks)
{
  //some convience typedefs, I hate typedefs but atleast here where they are defined is obvious
  typedef EgHLTDQMVarCut<OffEle> VarCut;
  typedef MonElemWithCutEBEE<OffEle,float> MonElemFloat;
  int stdCutCode = masks.stdEle;

  //first do the zero and all cuts histograms
  histVec.push_back(new MonElemFloat(baseName+"_noCuts",baseName+" NoCuts",nrBins,xMin,xMax,vsVarFunc));
  histVec.push_back(new MonElemFloat(baseName+"_allCuts",baseName+" All Cuts",nrBins,xMin,xMax,vsVarFunc,
				     new VarCut(stdCutCode,&OffEle::cutCode))); 
		       
  //now for the n-1
  histVec.push_back(new MonElemFloat(baseName+"_n1_dEtaIn",baseName+" N1 #Delta#eta_{in}",nrBins,xMin,xMax,vsVarFunc,
				     new VarCut(~EgCutCodes::DETAIN&stdCutCode,&OffEle::cutCode)));
  histVec.push_back(new MonElemFloat(baseName+"_n1_dPhiIn",baseName+" N1 #Delta#phi_{in}",nrBins,xMin,xMax,vsVarFunc,
				     new VarCut(~EgCutCodes::DPHIIN&stdCutCode,&OffEle::cutCode)));
  histVec.push_back(new MonElemFloat(baseName+"_n1_sigmaIEtaIEta",baseName+" N1 #sigma_{#ieta#ieta}",nrBins,xMin,xMax,vsVarFunc,
				     new VarCut(~EgCutCodes::SIGMAIETAIETA&stdCutCode,&OffEle::cutCode)));
  histVec.push_back(new MonElemFloat(baseName+"_n1_hOverE",baseName+" N1 H/E",nrBins,xMin,xMax,vsVarFunc,
				     new VarCut(~EgCutCodes::HADEM &stdCutCode,&OffEle::cutCode)));

  histVec.push_back(new MonElemFloat(baseName+"_n1_isolEm",baseName+" N1 Isol Em",nrBins,xMin,xMax,vsVarFunc,
				     new VarCut(~EgCutCodes::ISOLEM&stdCutCode,&OffEle::cutCode)));
  histVec.push_back(new MonElemFloat(baseName+"_n1_isolHad",baseName+" N1 Isol Had",nrBins,xMin,xMax,vsVarFunc,
				     new VarCut(~EgCutCodes::ISOLHAD&stdCutCode,&OffEle::cutCode)));
  histVec.push_back(new MonElemFloat(baseName+"_n1_isolPtTrks",baseName+" N1 Isol Tracks",nrBins,xMin,xMax,vsVarFunc,
				     new VarCut(~EgCutCodes::ISOLPTTRKS&stdCutCode,&OffEle::cutCode)));
  histVec.push_back(new MonElemFloat(baseName+"_n1_hltIsolHad",baseName+" N1 HLT Isol Had",nrBins,xMin,xMax,vsVarFunc,
				     new VarCut(~EgCutCodes::HLTISOLHAD&stdCutCode,&OffEle::cutCode)));
  histVec.push_back(new MonElemFloat(baseName+"_n1_hltIsolTrksEle",baseName+" N1 HLT Isol Tracks Ele ",nrBins,xMin,xMax,vsVarFunc,
				     new VarCut(~EgCutCodes::HLTISOLTRKSELE&stdCutCode,&OffEle::cutCode)));


  histVec.push_back(new MonElemFloat(baseName+"_single_dEtaIn",baseName+" Single #Delta#eta_{in}",nrBins,xMin,xMax,vsVarFunc,
				     new VarCut(EgCutCodes::DETAIN,&OffEle::cutCode)));
  histVec.push_back(new MonElemFloat(baseName+"_single_dPhiIn",baseName+" Single #Delta#phi_{in}",nrBins,xMin,xMax,vsVarFunc,
				     new VarCut(EgCutCodes::DPHIIN,&OffEle::cutCode)));
  histVec.push_back(new MonElemFloat(baseName+"_single_sigmaIEtaIEta",baseName+" Single #sigma_{#ieta#ieta}",nrBins,xMin,xMax,vsVarFunc,
				     new VarCut(EgCutCodes::SIGMAIETAIETA,&OffEle::cutCode)));
  histVec.push_back(new MonElemFloat(baseName+"_single_hOverE",baseName+" Single H/E",nrBins,xMin,xMax,vsVarFunc,
				     new VarCut(EgCutCodes::HADEM,&OffEle::cutCode)));
 

  histVec.push_back(new MonElemFloat(baseName+"_single_isolEm",baseName+" Single Isol Em",nrBins,xMin,xMax,vsVarFunc,
				     new VarCut(EgCutCodes::ISOLEM,&OffEle::cutCode)));
  histVec.push_back(new MonElemFloat(baseName+"_single_isolHad",baseName+" Single Isol Had",nrBins,xMin,xMax,vsVarFunc,
				     new VarCut(EgCutCodes::ISOLHAD,&OffEle::cutCode)));
  histVec.push_back(new MonElemFloat(baseName+"_single_isolPtTrks",baseName+" Single Isol Tracks",nrBins,xMin,xMax,vsVarFunc,
				     new VarCut(EgCutCodes::ISOLPTTRKS,&OffEle::cutCode)));

  histVec.push_back(new MonElemFloat(baseName+"_single_hltIsolHad",baseName+" Single HLT Isol Had",nrBins,xMin,xMax,vsVarFunc,
				     new VarCut(EgCutCodes::HLTISOLHAD,&OffEle::cutCode)));
  histVec.push_back(new MonElemFloat(baseName+"_single_hltIsolTrksEle",baseName+" Single HLT Isol Tracks Ele ",nrBins,xMin,xMax,vsVarFunc,
				     new VarCut(EgCutCodes::HLTISOLTRKSELE,&OffEle::cutCode))); 
  histVec.push_back(new MonElemFloat(baseName+"_single_hltIsolTrksPho",baseName+" Single HLT Isol Tracks Pho ",nrBins,xMin,xMax,vsVarFunc,
				     new VarCut(EgCutCodes::HLTISOLTRKSPHO,&OffEle::cutCode)));
  
}

void MonElemFuncs::initStdEffHists(std::vector<MonElemWithCutBase<OffPho>*>& histVec,const std::string& baseName,int nrBins,double xMin,double xMax,float (OffPho::*vsVarFunc)()const,const CutMasks& masks)
{
  //some convience typedefs, I hate typedefs but atleast here where they are defined is obvious
  typedef EgHLTDQMVarCut<OffPho> VarCut;
  typedef MonElemWithCutEBEE<OffPho,float> MonElemFloat;
  int stdCutCode = masks.stdPho;

  histVec.push_back(new MonElemFloat(baseName+"_noCuts",baseName+" NoCuts",nrBins,xMin,xMax,vsVarFunc));
  histVec.push_back(new MonElemFloat(baseName+"_allCuts",baseName+" All Cuts",nrBins,xMin,xMax,vsVarFunc,
				     new VarCut(stdCutCode,&OffPho::cutCode))); 

  
  histVec.push_back(new MonElemFloat(baseName+"_n1_sigmaIEtaIEta",baseName+" N1 #sigma_{#ieta#ieta}",nrBins,xMin,xMax,vsVarFunc,
				     new VarCut(~EgCutCodes::SIGMAIETAIETA&stdCutCode,&OffPho::cutCode)));
  histVec.push_back(new MonElemFloat(baseName+"_n1_hOverE",baseName+" N1 H/E",nrBins,xMin,xMax,vsVarFunc,
				     new VarCut(~EgCutCodes::R9&stdCutCode,&OffPho::cutCode)));
  histVec.push_back(new MonElemFloat(baseName+"_n1_r9",baseName+" N1 R9",nrBins,xMin,xMax,vsVarFunc,
				     new VarCut(~EgCutCodes::HADEM &stdCutCode,&OffPho::cutCode)));
  histVec.push_back(new MonElemFloat(baseName+"_n1_isolEm",baseName+" N1 Isol Em",nrBins,xMin,xMax,vsVarFunc,
				     new VarCut(~EgCutCodes::ISOLEM&stdCutCode,&OffPho::cutCode)));
  histVec.push_back(new MonElemFloat(baseName+"_n1_isolHad",baseName+" N1 Isol Had",nrBins,xMin,xMax,vsVarFunc,
				     new VarCut(~EgCutCodes::ISOLHAD&stdCutCode,&OffPho::cutCode)));
  histVec.push_back(new MonElemFloat(baseName+"_n1_isolPtTrks",baseName+" N1 Pt Isol Tracks",nrBins,xMin,xMax,vsVarFunc,
				     new VarCut(~EgCutCodes::ISOLPTTRKS&stdCutCode,&OffPho::cutCode))); 
  histVec.push_back(new MonElemFloat(baseName+"_n1_isolNrTrks",baseName+" N1 Nr Isol Tracks",nrBins,xMin,xMax,vsVarFunc,
				     new VarCut(~EgCutCodes::ISOLNRTRKS&stdCutCode,&OffPho::cutCode)));
 
  histVec.push_back(new MonElemFloat(baseName+"_single_sigmaIEtaIEta",baseName+" Single #sigma_{#ieta#ieta}",nrBins,xMin,xMax,vsVarFunc,
				     new VarCut(EgCutCodes::SIGMAIETAIETA,&OffPho::cutCode)));
  histVec.push_back(new MonElemFloat(baseName+"_single_hltIsolHad",baseName+" N1 HLT Isol Had",nrBins,xMin,xMax,vsVarFunc,
				     new VarCut(EgCutCodes::HLTISOLHAD,&OffPho::cutCode)));
  histVec.push_back(new MonElemFloat(baseName+"_single_hltIsolTrksPho",baseName+" N1 HLT Isol Tracks Pho ",nrBins,xMin,xMax,vsVarFunc,
				     new VarCut(EgCutCodes::HLTISOLTRKSPHO,&OffPho::cutCode)));

  
}


//we own the passed in cut, so we give it to the first mon element and then clone it after that
//only currently used for trigger tag and probe
void MonElemFuncs::initStdEleCutHists(std::vector<MonElemWithCutBase<OffEle>*>& histVec,const std::string& baseName,const BinData& bins,EgHLTDQMCut<OffEle>* cut)
{
  histVec.push_back(new MonElemWithCutEBEE<OffEle,float>(baseName+"_et",
							 baseName+" E_{T};E_{T} (GeV)",
							 bins.et.nr,bins.et.min,bins.et.max,&OffEle::et,cut));
  histVec.push_back(new MonElemWithCutEBEE<OffEle,float>(baseName+"_eta",
							 baseName+" #eta;#eta",
							 bins.eta.nr,bins.eta.min,bins.eta.max,
							 &OffEle::detEta,cut ? cut->clone(): NULL));		
  histVec.push_back(new MonElemWithCutEBEE<OffEle,float>(baseName+"_phi",
							 baseName+" #phi;#phi (rad)",
							 bins.phi.nr,bins.phi.min,bins.phi.max,
							 &OffEle::phi,cut ? cut->clone():NULL));
  histVec.push_back(new MonElemWithCutEBEE<OffEle,int>(baseName+"_charge",
						       baseName+" Charge; charge",
						       bins.charge.nr,bins.charge.min,bins.charge.max,
						       &OffEle::charge,cut ? cut->clone():NULL)); 
}



//we transfer ownership of eleCut to addTrigLooseTrigHist which transfers it to the eleMonElems
void MonElemFuncs::initTightLooseTrigHists(std::vector<MonElemContainer<OffEle>*>& eleMonElems,const std::vector<std::string>& tightLooseTrigs,const BinData& bins,EgHLTDQMCut<OffEle>* eleCut)
{
  for(size_t trigNr=0;trigNr<tightLooseTrigs.size();trigNr++){
    std::vector<std::string> splitString;
    boost::split(splitString,tightLooseTrigs[trigNr],boost::is_any_of(":"));
    if(splitString.size()!=2) continue; //format incorrect
    const std::string& tightTrig = splitString[0];
    const std::string& looseTrig = splitString[1];
    //this step is necessary as we want to transfer ownership of eleCut to the addTrigLooseTrigHist func on the last iteration
    //but clone it before that
    //perhaps my object ownership rules need to be re-evalulated
    if(trigNr!=tightLooseTrigs.size()-2) addTightLooseTrigHist(eleMonElems,tightTrig,looseTrig,eleCut->clone(),"gsfEle",bins);
    else addTightLooseTrigHist(eleMonElems,tightTrig,looseTrig,eleCut,"gsfEle",bins);
  }
} 


void MonElemFuncs::initTightLooseTrigHists(std::vector<MonElemContainer<OffPho>*>& phoMonElems,const std::vector<std::string>& tightLooseTrigs,const BinData& bins,EgHLTDQMCut<OffPho>* phoCut)
{
  for(size_t trigNr=0;trigNr<tightLooseTrigs.size();trigNr++){
    std::vector<std::string> splitString;
    boost::split(splitString,tightLooseTrigs[trigNr],boost::is_any_of(":"));
    if(splitString.size()!=2) continue; //format incorrect
    const std::string& tightTrig = splitString[0];
    const std::string& looseTrig = splitString[1];

    //this step is necessary as we want to transfer ownership of phoCut to the addTrigLooseTrigHist func on the last iteration
    //but clone it before that
    //perhaps my object ownership rules need to be re-evalulated
    if(trigNr!=tightLooseTrigs.size()-2) addTightLooseTrigHist(phoMonElems,tightTrig,looseTrig,phoCut->clone(),"pho",bins);
    else addTightLooseTrigHist(phoMonElems,tightTrig,looseTrig,phoCut,"pho",bins);
  }
}
   
//there is nothing electron specific here, will template at some point
void MonElemFuncs::addTightLooseTrigHist(std::vector<MonElemContainer<OffEle>*>& eleMonElems,
					 const std::string& tightTrig,const std::string& looseTrig,
					 EgHLTDQMCut<OffEle>* eleCut,
					 const std::string& histId,const BinData& bins)
{
  MonElemContainer<OffEle>* passMonElem = NULL;
  passMonElem = new MonElemContainer<OffEle>(tightTrig+"_"+looseTrig+"_"+histId+"_passTrig","",
						  &(*(new EgMultiCut<OffEle>) << 
						    new EgObjTrigCut<OffEle>(TrigCodes::getCode(tightTrig+":"+looseTrig),EgObjTrigCut<OffEle>::AND)  <<
						    eleCut->clone()));
  
  
  MonElemContainer<OffEle>* failMonElem = NULL;
  failMonElem = new MonElemContainer<OffEle>(tightTrig+"_"+looseTrig+"_"+histId+"_failTrig","",
						  &(*(new EgMultiCut<OffEle>) << 
						    new EgObjTrigCut<OffEle>(TrigCodes::getCode(looseTrig),EgObjTrigCut<OffEle>::AND,TrigCodes::getCode(tightTrig))  << 
						    eleCut));
  
  MonElemFuncs::initStdEleHists(passMonElem->monElems(),passMonElem->name(),bins);
  MonElemFuncs::initStdEleHists(failMonElem->monElems(),failMonElem->name(),bins);
  eleMonElems.push_back(passMonElem);
  eleMonElems.push_back(failMonElem);
}

 
//there is nothing photon specific here, will template at some point
void MonElemFuncs::addTightLooseTrigHist(std::vector<MonElemContainer<OffPho>*>& phoMonElems,
					 const std::string& tightTrig,const std::string& looseTrig,
					 EgHLTDQMCut<OffPho>* phoCut,
					 const std::string& histId,const BinData& bins)
{
  MonElemContainer<OffPho>* passMonElem = NULL;
  passMonElem = new MonElemContainer<OffPho>(tightTrig+"_"+looseTrig+"_"+histId+"_passTrig","",
						  &(*(new EgMultiCut<OffPho>) << 
						    new EgObjTrigCut<OffPho>(TrigCodes::getCode(tightTrig+":"+looseTrig),EgObjTrigCut<OffPho>::AND)  <<
						    phoCut->clone()));
  
  
  MonElemContainer<OffPho>* failMonElem = NULL;
  failMonElem = new MonElemContainer<OffPho>(tightTrig+"_"+looseTrig+"_"+histId+"_failTrig","",
						  &(*(new EgMultiCut<OffPho>) << 
						    new EgObjTrigCut<OffPho>(TrigCodes::getCode(looseTrig),EgObjTrigCut<OffPho>::AND,TrigCodes::getCode(tightTrig))  << 
						    phoCut));
  
  MonElemFuncs::initStdPhoHists(passMonElem->monElems(),passMonElem->name(),bins);
  MonElemFuncs::initStdPhoHists(failMonElem->monElems(),failMonElem->name(),bins);
  phoMonElems.push_back(passMonElem);
  phoMonElems.push_back(failMonElem);
}
  

//we transfer ownership of eleCut to the monitor elements
void MonElemFuncs::initTightLooseTrigHistsTrigCuts(std::vector<MonElemContainer<OffEle>*>& eleMonElems,const std::vector<std::string>& tightLooseTrigs,const BinData& bins)
{
  for(size_t trigNr=0;trigNr<tightLooseTrigs.size();trigNr++){
    std::vector<std::string> splitString;
    boost::split(splitString,tightLooseTrigs[trigNr],boost::is_any_of(":"));
    if(splitString.size()!=2) continue; //format incorrect
    const std::string& tightTrig = splitString[0];
    const std::string& looseTrig = splitString[1];
    EgHLTDQMCut<OffEle>* eleCut = new EgHLTDQMUserVarCut<OffEle,TrigCodes::TrigBitSet>(&OffEle::trigCutsCutCode,TrigCodes::getCode(tightTrig));
    addTightLooseTrigHist(eleMonElems,tightTrig,looseTrig,eleCut,"gsfEle_trigCuts",bins);
  }
} 

//we transfer ownership of phoCut to the monitor elements
void MonElemFuncs::initTightLooseTrigHistsTrigCuts(std::vector<MonElemContainer<OffPho>*>& phoMonElems,const std::vector<std::string>& tightLooseTrigs,const BinData& bins)
{
  for(size_t trigNr=0;trigNr<tightLooseTrigs.size();trigNr++){
    std::vector<std::string> splitString;
    boost::split(splitString,tightLooseTrigs[trigNr],boost::is_any_of(":"));
    if(splitString.size()!=2) continue; //format incorrect
    const std::string& tightTrig = splitString[0];
    const std::string& looseTrig = splitString[1];
    EgHLTDQMCut<OffPho>* phoCut = new EgHLTDQMUserVarCut<OffPho,TrigCodes::TrigBitSet>(&OffPho::trigCutsCutCode,TrigCodes::getCode(tightTrig));
    addTightLooseTrigHist(phoMonElems,tightTrig,looseTrig,phoCut,"pho_trigCuts",bins);
  }
} 


//we transfer ownership of eleCut to the monitor elements
void MonElemFuncs::initTightLooseDiObjTrigHistsTrigCuts(std::vector<MonElemContainer<OffEle>*>& eleMonElems,const std::vector<std::string>& tightLooseTrigs,const BinData& bins)
{
  for(size_t trigNr=0;trigNr<tightLooseTrigs.size();trigNr++){
    std::vector<std::string> splitString;
    boost::split(splitString,tightLooseTrigs[trigNr],boost::is_any_of(":"));
    if(splitString.size()!=2) continue; //format incorrect
    const std::string& tightTrig = splitString[0];
    const std::string& looseTrig = splitString[1];
    EgHLTDQMCut<OffEle>* eleCut = new EgDiEleUserCut<TrigCodes::TrigBitSet>(&OffEle::trigCutsCutCode,TrigCodes::getCode(tightTrig));
    addTightLooseTrigHist(eleMonElems,tightTrig,looseTrig,eleCut,"gsfEle_trigCuts",bins);
  }
} 


//we transfer ownership of phoCut to the monitor elements
void MonElemFuncs::initTightLooseDiObjTrigHistsTrigCuts(std::vector<MonElemContainer<OffPho>*>& phoMonElems,const std::vector<std::string>& tightLooseTrigs,const BinData& bins)
{
  for(size_t trigNr=0;trigNr<tightLooseTrigs.size();trigNr++){
    std::vector<std::string> splitString;
    boost::split(splitString,tightLooseTrigs[trigNr],boost::is_any_of(":"));
    if(splitString.size()!=2) continue; //format incorrect
    const std::string& tightTrig = splitString[0];
    const std::string& looseTrig = splitString[1];
    EgHLTDQMCut<OffPho>* phoCut = new EgDiPhoUserCut<TrigCodes::TrigBitSet>(&OffPho::trigCutsCutCode,TrigCodes::getCode(tightTrig));
    addTightLooseTrigHist(phoMonElems,tightTrig,looseTrig,phoCut,"pho_trigCuts",bins);
  }
}



//tag and probe trigger efficiencies
//this is to do measure the trigger efficiency with respect to a fully selected offline electron
//using a tag and probe technique (note: this will be different to the trigger efficiency normally calculated) 
void MonElemFuncs::initTrigTagProbeHists(std::vector<MonElemContainer<OffEle>*>& eleMonElems,const std::vector<std::string> filterNames,const BinData& bins)
{
  for(size_t filterNr=0;filterNr<filterNames.size();filterNr++){ 
    
    std::string trigName(filterNames[filterNr]);
    int stdCutCode = EgCutCodes::getCode("detEta:crack:sigmaIEtaIEta:hadem:dPhiIn:dEtaIn"); //will have it non hardcoded at a latter date
    MonElemContainer<OffEle>* monElemCont = new MonElemContainer<OffEle>("trigTagProbe","Trigger Tag and Probe",new EgTrigTagProbeCut(TrigCodes::getCode(trigName),stdCutCode,&OffEle::cutCode));
    MonElemFuncs::initStdEleCutHists(monElemCont->cutMonElems(),trigName+"_"+monElemCont->name()+"_gsfEle_all",bins,new EgGreaterCut<OffEle,float>(15.,&OffEle::etSC));
    MonElemFuncs::initStdEleCutHists(monElemCont->cutMonElems(),trigName+"_"+monElemCont->name()+"_gsfEle_pass",bins,&(*(new EgMultiCut<OffEle>) << new EgGreaterCut<OffEle,float>(15.,&OffEle::etSC) << new EgObjTrigCut<OffEle>(TrigCodes::getCode(trigName),EgObjTrigCut<OffEle>::AND)));
    eleMonElems.push_back(monElemCont);
  } //end filter names loop
   
}


