/*! \file SiPixelInformationExtractor.cc
 *  \brief This class represents ...
 *  
 *  (Documentation under development)
 *  
 */
#include "DQM/SiPixelMonitorClient/interface/SiPixelInformationExtractor.h"
#include "DQM/SiPixelMonitorClient/interface/SiPixelUtility.h"
#include "DQM/SiPixelMonitorClient/interface/SiPixelEDAClient.h"
#include "DQM/SiPixelMonitorClient/interface/ANSIColors.h"
#include "DQM/SiPixelMonitorClient/interface/SiPixelHistoPlotter.h"
#include "DQM/SiPixelCommon/interface/SiPixelFolderOrganizer.h"
#include "DQMServices/Core/interface/DQMStore.h"
#include "DQMServices/Core/interface/MonitorElement.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "DQMServices/WebComponents/interface/CgiReader.h"
#include "FWCore/ParameterSet/interface/FileInPath.h"
#include "DataFormats/DetId/interface/DetId.h"
#include "DataFormats/GeometrySurface/interface/Surface.h"
#include "Geometry/TrackerNumberingBuilder/interface/GeometricDet.h"
#include "Geometry/TrackerGeometryBuilder/interface/TrackerGeometry.h"
#include "Geometry/Records/interface/TrackerDigiGeometryRecord.h"
#include "Geometry/TrackerGeometryBuilder/interface/PixelGeomDetUnit.h"
#include "Geometry/CommonTopologies/interface/PixelTopology.h"

#include "TClass.h"
#include "TText.h"
#include "TROOT.h"
#include "TPad.h"
#include "TStyle.h"
#include "TSystem.h"
#include "TString.h"
#include "TImage.h"
#include "TPaveText.h"
#include "TImageDump.h"
#include "TRandom.h"
#include "TStopwatch.h"
#include "TAxis.h"
#include "TPaveLabel.h"
#include "Rtypes.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TProfile.h"

#include <qstring.h>
#include <qregexp.h>

#include <iostream>
#include <math.h>

#include <cstdlib> // for free() - Root can allocate with malloc() - sigh...
 
using namespace std;
using namespace edm;

//------------------------------------------------------------------------------
/*! \brief Constructor of the SiPixelInformationExtractor class.
 *  
 */
SiPixelInformationExtractor::SiPixelInformationExtractor(std::string summaryXMLfileName) : summaryXMLfileName_(summaryXMLfileName) {
  edm::LogInfo("SiPixelInformationExtractor") << 
    " Creating SiPixelInformationExtractor " << "\n" ;
  
  canvas_ = new TCanvas("PlotCanvas", "Plot Canvas"); 
  readReference_ = false;
  allMods_=0;
  errorMods_=0;
  qflag_=1.;
  histoPlotter_=0;
  histoPlotter_ = new SiPixelHistoPlotter();
}

//------------------------------------------------------------------------------
/*! \brief Destructor of the SiPixelInformationExtractor class.
 *  
 */
SiPixelInformationExtractor::~SiPixelInformationExtractor() {
  edm::LogInfo("SiPixelInformationExtractor") << 
    " Deleting SiPixelInformationExtractor " << "\n" ;
  
  if (canvas_) delete canvas_;
  if (histoPlotter_) delete histoPlotter_;
}

//------------------------------------------------------------------------------
/*! \brief Read Configuration File
 *
 */
void SiPixelInformationExtractor::readConfiguration() { }

//
// -- Select Histograms for a given module
//
void SiPixelInformationExtractor::getSingleModuleHistos(DQMStore * bei, 
                                                        const multimap<string, string>& req_map, 
							xgi::Output * out){

  vector<string> hlist;
  getItemList(req_map,"histo", hlist);

  uint32_t detId = atoi(getItemValue(req_map,"ModId").c_str());
 
  int width  = atoi(getItemValue(req_map, "width").c_str());
  int height = atoi(getItemValue(req_map, "height").c_str());

  string opt =" ";
  
  SiPixelFolderOrganizer folder_organizer;
  string path;
  folder_organizer.getModuleFolder(detId,path);   

  if((bei->pwd()).find("Module_") == string::npos &&
     (bei->pwd()).find("FED_") == string::npos){
    cout<<"This is not a pixel module or FED!"<<endl;
    return;
  }
 
  vector<MonitorElement*> all_mes = bei->getContents(path);
  setHTMLHeader(out);
  *out << path << " ";

  QRegExp rx("(\\w+)_(siPixel|ctfWithMaterialTracks)") ;
  QString theME ;

  for (vector<string>::const_iterator ih = hlist.begin();
       ih != hlist.end(); ih++) {
    for (vector<MonitorElement *>::const_iterator it = all_mes.begin();
	 it!= all_mes.end(); it++) {
      MonitorElement * me = (*it);
      if (!me) continue;
      theME = me->getName();
      string temp_s ; 
      if( rx.search(theME) != -1 ) { temp_s = rx.cap(1).latin1() ; }
      if (temp_s == (*ih)) {
	string full_path = path + "/" + me->getName();
	histoPlotter_->setNewPlot(full_path, opt, width, height);
	*out << me->getName() << " " ;
      }
    }
  }
}

//
// -- Plot Tracker Map MEs
//
void SiPixelInformationExtractor::getTrackerMapHistos(DQMStore* bei, 
                                                      const std::multimap<std::string, std::string>& req_map, 
						      xgi::Output * out) {

//  cout << __LINE__ << ACYellow << ACBold 
//       << "[SiPixelInformationExtractor::getTrackerMapHistos] " << ACPlain << endl ;
  vector<string> hlist;
  string tkmap_name;
  SiPixelConfigParser config_parser;
//  string localPath = string("DQM/SiPixelMonitorClient/test/sipixel_monitorelement_config.xml");
  string localPath = summaryXMLfileName_;
  config_parser.getDocument(edm::FileInPath(localPath).fullPath());
//  if (!config_parser.getMENamesForTrackerMap(tkmap_name, hlist)) return;
//  if (hlist.size() == 0) return;
  if (!config_parser.getMENamesForTrackerMap(tkmap_name, hlist)) 
  {
   cout << __LINE__ << ACYellow << ACBold 
        << "[SiPixelInformationExtractor::getTrackerMapHistos] " 
	<< ACPlain << ACRed << ACPlain 
	<< "getMENamesForTrackerMap return false " 
        << ACPlain << endl ; assert(0) ;
   return;
  }
  if (hlist.size() == 0) 
  {
   cout << __LINE__ << ACYellow << ACBold 
        << "[SiPixelInformationExtractor::getTrackerMapHistos] " 
	<< ACPlain << ACRed << ACPlain 
	<< "hlist.size() == 0 " 
        << ACPlain << endl ;  assert(0) ;
   return;
  }


  uint32_t detId = atoi(getItemValue(req_map,"ModId").c_str());
 
  int width  = atoi(getItemValue(req_map, "width").c_str());
  int height = atoi(getItemValue(req_map, "height").c_str());

  string opt =" ";
  
  SiPixelFolderOrganizer folder_organizer;
  string path;
  folder_organizer.getModuleFolder(detId,path);   
/*
  if((bei->pwd()).find("Module_") == string::npos &&
     (bei->pwd()).find("FED_") == string::npos){
    cout<<"This is not a pixel module or FED!"<<endl;
   cout << __LINE__ << ACYellow << ACBold 
        << "[SiPixelInformationExtractor::getTrackerMapHistos] " 
	<< ACPlain << ACRed << ACPlain 
	<< "This is not a pixel module or FED!" 
        << ACPlain << endl ; assert(0) ;
    return;
  }
*/ 
  vector<MonitorElement*> all_mes = bei->getContents(path);
  setXMLHeader(out);

//  cout << __LINE__ << ACCyan << ACBold 
//       << " [SiPixelInformationExtractor::getTrackerMapHistos()] path "
//       << ACPlain << path << endl ; 
//  cout << __LINE__ << ACCyan << ACBold 
//       << " [SiPixelInformationExtractor::getTrackerMapHistos()] all_mes.size() "
//       << ACPlain << all_mes.size() << endl ; 

  QRegExp rx("(\\w+)_(siPixel|ctfWithMaterialTracks)") ;
  QString theME ;

  *out << "<pathList>" << endl ;
  for (vector<string>::iterator ih = hlist.begin();
       ih != hlist.end(); ih++) {
    for (vector<MonitorElement *>::const_iterator it = all_mes.begin();
	 it!= all_mes.end(); it++) {
      MonitorElement * me = (*it);
      if (!me) 
      { 
//       cout << __LINE__ << ACCyan << ACBold 
//            << " [SiPixelInformationExtractor::getTrackerMapHistos()] skipping "
//        	       << ACPlain << *ih << endl ; 
       continue;
      }
      theME = me->getName();
      string temp_s ; 
      if( rx.search(theME) != -1 ) { temp_s = rx.cap(1).latin1() ; }
//      cout << __LINE__ << ACCyan << ACBold 
//           << " [SiPixelInformationExtractor::getTrackerMapHistos()] temp_s "
//           << ACPlain << temp_s << " <--> " << *ih << " |" << theME << "|" << endl ; 
      if (temp_s == (*ih)) {
	string full_path = path + "/" + me->getName();
	histoPlotter_->setNewPlot(full_path, opt, width, height);
//cout << __LINE__ << ACRed << ACBold 
//     << " [SiPixelInformationExtractor::getTrackerMapHistos()] fullPath: "
//     << ACPlain << full_path << endl ; 
	*out << " <pathElement path='" << full_path << "' />" << endl ;
      }      
    }
  }   
  *out << "</pathList>" << endl ;
//cout << __LINE__ << " [SiPixelInformationExtractor::getTrackerMapHistos()] endlist: " << endl ;
}

//============================================================================================================
// --  Return type of ME
//
std::string  SiPixelInformationExtractor::getMEType(MonitorElement * theMe)
{
  QString qtype = theMe->getRootObject()->IsA()->GetName() ;
  if(         qtype.contains("TH1") > 0 )
  {
    return "TH1" ;
  } else if ( qtype.contains("TH2") > 0  ) {
    return "TH2" ;
  } else if ( qtype.contains("TH3") > 0 ) {
    return "TH3" ;
  }
  return "TH1" ;
}

//------------------------------------------------------------------------------
/*! \brief (Documentation under construction).
 *
 *  This method 
 */
void SiPixelInformationExtractor::readModuleAndHistoList(DQMStore* bei, 
                                                         xgi::Output * out) {
//cout<<"entering SiPixelInformationExtractor::readModuleAndHistoList"<<endl;
   std::map<std::string,std::string> hnames;
   std::vector<std::string> mod_names;
   fillModuleAndHistoList(bei, mod_names, hnames);
   out->getHTTPResponseHeader().addHeader("Content-Type", "text/xml");
  *out << "<?xml version=\"1.0\" ?>" << std::endl;
  *out << "<ModuleAndHistoList>" << endl;
  *out << "<ModuleList>" << endl;
   for (std::vector<std::string>::iterator im = mod_names.begin();
        im != mod_names.end(); im++) {
     *out << "<ModuleNum>" << *im << "</ModuleNum>" << endl;     
   }
   *out << "</ModuleList>" << endl;
   *out << "<HistoList>" << endl;

   for (std::map<std::string,std::string>::iterator ih = hnames.begin();
        ih != hnames.end(); ih++) {
     *out << "<Histo type=\"" 
          << ih->second
	  << "\">" 
	  << ih->first 
	  << "</Histo>" 
	  << endl;     
   }
   *out << "</HistoList>" << endl;
   *out << "</ModuleAndHistoList>" << endl;
//cout<<"leaving SiPixelInformationExtractor::readModuleAndHistoList"<<endl;
}

//------------------------------------------------------------------------------
/*! \brief (Documentation under construction).
 *
 *  This method 
 */
void SiPixelInformationExtractor::fillModuleAndHistoList(DQMStore * bei, 
                                                         vector<string>        & modules,
							 map<string,string>    & histos) {
//cout<<"entering SiPixelInformationExtractor::fillModuleAndHistoList"<<endl;
  string currDir = bei->pwd();
  if (currDir.find("Module_") != string::npos ||
      currDir.find("FED_") != string::npos)  {
    if (histos.size() == 0) {
      //cout<<"currDir="<<currDir<<endl;

      vector<string> contents = bei->getMEs();
          
      for (vector<string>::const_iterator it = contents.begin();
	   it != contents.end(); it++) {
	string hname          = (*it).substr(0, (*it).find("_siPixel"));
	if (hname==" ") hname = (*it).substr(0, (*it).find("_ctfWithMaterialTracks"));
        string fullpathname   = bei->pwd() + "/" + (*it); 

        MonitorElement * me   = bei->get(fullpathname);
	
        string htype          = "undefined" ;
        if (me) 
	{
	 htype = me->getRootObject()->IsA()->GetName() ;
	}
	//cout<<"hname="<<hname<<endl;
        histos[hname] = htype ;
        string mId=" ";
	if(hname.find("ndigis")                !=string::npos) mId = (*it).substr((*it).find("ndigis_siPixelDigis_")+20, 9);
	if(mId==" " && hname.find("nclusters") !=string::npos) mId = (*it).substr((*it).find("nclusters_siPixelClusters_")+26, 9);
        if(mId==" " && hname.find("residualX") !=string::npos) mId = (*it).substr((*it).find("residualX_ctfWithMaterialTracks_")+32, 9);
        if(mId==" " && hname.find("NErrors") !=string::npos) mId = (*it).substr((*it).find("NErrors_siPixelDigis_")+21, 9);
        if(mId==" " && hname.find("ClustX") !=string::npos) mId = (*it).substr((*it).find("ClustX_siPixelRecHit_")+21, 9);
        if(mId==" " && hname.find("pixelAlive") !=string::npos) mId = (*it).substr((*it).find("pixelAlive_siPixelCalibDigis_")+29, 9);
        if(mId==" " && hname.find("Gain1d") !=string::npos) mId = (*it).substr((*it).find("Gain1d_siPixelCalibDigis_")+25, 9);
        if(mId!=" ") modules.push_back(mId);
        //cout<<"mId="<<mId<<endl;
      }    
    }
  } else {  
    vector<string> subdirs = bei->getSubdirs();
    for (vector<string>::const_iterator it = subdirs.begin();
	 it != subdirs.end(); it++) {
      bei->cd(*it);
      fillModuleAndHistoList(bei, modules, histos);
      bei->goUp();
    }
  }
//  fillBarrelList(bei, modules, histos);
//cout<<"leaving SiPixelInformationExtractor::fillModuleAndHistoList"<<endl;
}

//------------------------------------------------------------------------------
/*! \brief (Documentation under construction).
 *
 *  This method 
 */
void SiPixelInformationExtractor::readModuleHistoTree(DQMStore* bei, 
                                                      string& str_name, 
						      xgi::Output * out) {
//cout<<"entering  SiPixelInformationExtractor::readModuleHistoTree"<<endl;
  ostringstream modtree;
  if (goToDir(bei, str_name)) {
    modtree << "<form name=\"IMGCanvasItemsSelection\" "
            << "action=\"javascript:void%200\">" 
	    << endl ;
    modtree << "<ul id=\"dhtmlgoodies_tree\" class=\"dhtmlgoodies_tree\">" << endl;
    printModuleHistoList(bei,modtree);
    modtree <<"</ul>" << endl;   
    modtree <<"</form>" << endl;   
  } else {
    modtree << "Desired Directory does not exist";
  }
  cout << ACYellow << ACBold
       << "[SiPixelInformationExtractor::readModuleHistoTree()]"
       << ACPlain << endl ;
  //     << "html string follows: " << endl ;
  //cout << modtree.str() << endl ;
  //cout << ACYellow << ACBold
  //     << "[SiPixelInformationExtractor::readModuleHistoTree()]"
  //     << ACPlain
  //     << "String complete " << endl ;
  out->getHTTPResponseHeader().addHeader("Content-Type", "text/plain");
  *out << modtree.str();
   bei->cd();
//cout<<"leaving  SiPixelInformationExtractor::readModuleHistoTree"<<endl;
}

//------------------------------------------------------------------------------
/*! \brief (Documentation under construction).
 *
 *  This method 
 */
void SiPixelInformationExtractor::printModuleHistoList(DQMStore * bei, 
                                                       ostringstream& str_val){
//cout<<"entering SiPixelInformationExtractor::printModuleHistoList"<<endl;
  static string indent_str = "";
  string currDir = bei->pwd();
  string dname = currDir.substr(currDir.find_last_of("/")+1);
  str_val << " <li>\n"
	  << "  <a href=\"#\" id=\"" << currDir << "\">\n   " 
	  <<     dname << "\n"
	  << "  </a>\n"
	  << endl << endl;

  vector<string> meVec     = bei->getMEs(); 
  
  vector<string> subDirVec = bei->getSubdirs();
  if ( meVec.size()== 0  && subDirVec.size() == 0 ) {
    str_val << " </li>" << endl;    
    return;
  }
  str_val << "\n   <ul>" << endl; 
  for (vector<string>::const_iterator it  = meVec.begin();
                                      it != meVec.end(); it++) {
    if ((*it).find("_siPixel")!=string::npos || 
        (*it).find("_ctfWithMaterialTracks")!=string::npos) {
      QString qit = (*it) ;
      QRegExp rx("(\\w+)_(siPixel|ctfWithMaterialTracks)") ;
      if( rx.search(qit) > -1 ) {qit = rx.cap(1);} 
      str_val << "    <li class=\"dhtmlgoodies_sheet.gif\">\n"
	      << "     <input id      = \"selectedME\""
	      << "            folder  = \"" << currDir << "\""
	      << "            type    = \"checkbox\""
	      << "            name    = \"selected\""
	      << "            class   = \"smallCheckBox\""
	      << "            value   = \"" << (*it) << "\""
	      << "            onclick = \"javascript:IMGC.selectedIMGCItems()\" />\n"
//	      << "     <a href=\"javascript:IMGC.updateIMGC('" << currDir << "')\">\n       " 
	      << "     <a href=\"javascript:IMGC.plotFromPath('" << currDir << "')\">\n       " 
//	      <<        qit << "\n"
	      <<        (*it) << "\n"
	      << "     </a>\n"
	      << "    </li>" 
	      << endl;
    }
  }
  for (vector<string>::const_iterator ic  = subDirVec.begin();
                                      ic != subDirVec.end(); ic++) {
    bei->cd(*ic);
    printModuleHistoList(bei, str_val);
    bei->goUp();
  }
  str_val << "   </ul>" << endl;  
  str_val << "  </li>"  << endl;  
//cout<<"leaving SiPixelInformationExtractor::printModuleHistoList"<<endl;
}

//------------------------------------------------------------------------------
/*! \brief (Documentation under construction).
 *
 *  This method 
 */
void SiPixelInformationExtractor::readSummaryHistoTree(DQMStore* bei, 
                                                       string& str_name, 
						       xgi::Output * out) {
//cout<<"entering  SiPixelInformationExtractor::readSummaryHistoTree"<<endl;
  ostringstream sumtree;
  if (goToDir(bei, str_name)) {
    sumtree << "<ul id=\"dhtmlgoodies_tree\" class=\"dhtmlgoodies_tree\">" << endl;
    printSummaryHistoList(bei,sumtree);
    sumtree <<"</ul>" << endl;   
  } else {
    sumtree << "Desired Directory does not exist";
  }
  cout << ACYellow << ACBold
       << "[SiPixelInformationExtractor::readSummaryHistoTree()]"
       << ACPlain << endl ;
  //     << "html string follows: " << endl ;
  //cout << sumtree.str() << endl ;
  //cout << ACYellow << ACBold
  //     << "[SiPixelInformationExtractor::readSummaryHistoTree()]"
  //     << ACPlain
  //     << "String complete " << endl ;
  out->getHTTPResponseHeader().addHeader("Content-Type", "text/plain");
  *out << sumtree.str();
   bei->cd();
//cout<<"leaving  SiPixelInformationExtractor::readSummaryHistoTree"<<endl;
}
//------------------------------------------------------------------------------
/*! \brief (Documentation under construction).
 *
 *  Returns a stringstream containing an HTML-formatted list of ME in the current
 *  directory. 
 *  This is a recursive method.
 */
void SiPixelInformationExtractor::printSummaryHistoList(DQMStore * bei, 
                                                        ostringstream& str_val){
//cout<<"entering SiPixelInformationExtractor::printSummaryHistoList"<<endl;
  static string indent_str = "";
  string currDir = bei->pwd();
  string dname = currDir.substr(currDir.find_last_of("/")+1);
  if (dname.find("Module_") ==0 || dname.find("FED_")==0) return;
  str_val << " <li>\n"
          << "  <a href=\"#\" id=\"" << currDir << "\">\n   " 
	  <<     dname 
	  << "  </a>" 
	  << endl;

  vector<string> meVec     = bei->getMEs(); 
  
  vector<string> subDirVec = bei->getSubdirs();
  if ( meVec.size()== 0  && subDirVec.size() == 0 ) {
    str_val << " </li> "<< endl;    
    return;
  }
  str_val << "\n   <ul>" << endl;      
  for (vector<string>::const_iterator it = meVec.begin();
       it != meVec.end(); it++) {
    if ((*it).find("SUM") == 0) {
      QString qit = (*it) ;
      QRegExp rx("(\\w+)_(siPixel|ctfWithMaterialTracks)");
      if( rx.search(qit) > -1 ) {qit = rx.cap(1);} 
      str_val << "    <li class=\"dhtmlgoodies_sheet.gif\">\n"
	      << "     <input id      = \"selectedME\""
	      << "            folder  = \"" << currDir << "\""
	      << "            type    = \"checkbox\""
	      << "            name    = \"selected\""
	      << "            class   = \"smallCheckBox\""
	      << "            value   = \"" << (*it) << "\""
	      << "            onclick = \"javascript:IMGC.selectedIMGCItems()\" />\n"
//              << "     <a href=\"javascript:IMGC.updateIMGC('" << currDir << "')\">\n       " 
              << "     <a href=\"javascript:IMGC.plotFromPath('" << currDir << "')\">\n       " 
	      <<        qit << "\n"
	      << "     </a>\n"
	      << "    </li>" 
	      << endl;
    }
  }

  for (vector<string>::const_iterator ic = subDirVec.begin();
       ic != subDirVec.end(); ic++) {
    bei->cd(*ic);
    printSummaryHistoList(bei, str_val);
    bei->goUp();
  }
  str_val << "   </ul> "<< endl;  
  str_val << "  </li> "<< endl;  
//cout<<"leaving SiPixelInformationExtractor::printSummaryHistoList"<<endl;
}


//------------------------------------------------------------------------------
/*! \brief (Documentation under construction).
 *
 *  This method 
 */
void SiPixelInformationExtractor::readAlarmTree(DQMStore* bei, 
                                                string& str_name, 
						xgi::Output * out){
//cout<<"entering SiPixelInformationExtractor::readAlarmTree"<<endl;
  ostringstream alarmtree;
  if (goToDir(bei, str_name)) {
    alarmtree << "<ul id=\"dhtmlgoodies_tree\" class=\"dhtmlgoodies_tree\">" << endl;
    alarmCounter_=0;
    printAlarmList(bei,alarmtree);
    if(alarmCounter_==0) alarmtree <<"<li>No problematic modules found, all ok!</li>" << endl;
    alarmtree <<"</ul>" << endl; 
  } else {
    alarmtree << "Desired Directory does not exist";
  }
  cout << ACYellow << ACBold
       << "[SiPixelInformationExtractor::readAlarmTree()]"
       << ACPlain << endl ;
  //     << "html string follows: " << endl ;
  //cout << alarmtree.str() << endl ;
  //cout << ACYellow << ACBold
  //     << "[SiPixelInformationExtractor::readAlarmTree()]"
  //     << ACPlain
  //     << "String complete " << endl ;
  out->getHTTPResponseHeader().addHeader("Content-Type", "text/plain");
 *out << alarmtree.str();
  bei->cd();
  cout << ACYellow << ACBold
       << "[SiPixelInformationExtractor::readAlarmTree()]"
       << ACPlain 
       << " Done!"
       << endl ;
//cout<<"leaving SiPixelInformationExtractor::readAlarmTree"<<endl;
}
//------------------------------------------------------------------------------
/*! \brief (Documentation under construction).
 *  
 *  Returns a stringstream containing an HTML-formatted list of alarms for the current
 *  directory. 
 *  This is a recursive method.
 */
void SiPixelInformationExtractor::printAlarmList(DQMStore * bei, 
                                                 ostringstream& str_val){
//cout<<"entering SiPixelInformationExtractor::printAlarmList"<<endl;
//   cout << ACRed << ACBold
//        << "[SiPixelInformationExtractor::printAlarmList()]"
//        << ACPlain
//        << " Enter" 
//        << endl ;
  static string indent_str = "";
  string currDir = bei->pwd();
  string dname = currDir.substr(currDir.find_last_of("/")+1);
  string image_name;
  selectImage(image_name,bei->getStatus(currDir));
  if(image_name!="images/LI_green.gif")
    str_val << " <li>\n"
            << "  <a href=\"#\" id=\"" << currDir << "\">\n   " 
	    <<     dname 
	    << "  </a>\n"
	    << "  <img src=\"" 
            <<     image_name 
	    << "\">" << endl;
  vector<string> subDirVec = bei->getSubdirs();

  vector<string> meVec = bei->getMEs();
   
  if (subDirVec.size() == 0 && meVec.size() == 0) {
    str_val <<  "</li> "<< endl;    
    return;
  }
  str_val << "<ul>" << endl;
  for (vector<string>::const_iterator it = meVec.begin();
	   it != meVec.end(); it++) {
    string full_path = currDir + "/" + (*it);

    MonitorElement * me = bei->get(full_path);
    
    if (!me) continue;
    std::vector<QReport *> my_map = me->getQReports();
    if (my_map.size() > 0) {
      string image_name1;
      selectImage(image_name1,my_map);
      if(image_name1!="images/LI_green.gif") {
        alarmCounter_++;
        QString qit = (*it) ;
        QRegExp rx("(\\w+)_(siPixel|ctfWithMaterialTracks)") ;
        if( rx.search(qit) > -1 ) {qit = rx.cap(1);} 
        str_val << "	<li class=\"dhtmlgoodies_sheet.gif\">\n"
        	<< "	 <input id	= \"selectedME\""
        	<< "		folder  = \"" << currDir << "\""
        	<< "		type	= \"checkbox\""
        	<< "		name	= \"selected\""
        	<< "		class	= \"smallCheckBox\""
        	<< "		value	= \"" << (*it) << "\""
        	<< "		onclick = \"javascript:IMGC.selectedIMGCItems()\" />\n"
//        	<< "	 <a href=\"javascript:IMGC.updateIMGC('" << currDir << "')\">\n       " 
        	<< "	 <a href=\"javascript:IMGC.plotFromPath('" << currDir << "')\">\n       " 
        	<<	  qit << "\n"
        	<< "	 </a>\n"
		<< "     <img src=\""
		<<        image_name1 
		<< "\">"
        	<< "	</li>" 
        	<< endl;
	}	
    }
  }
  for (vector<string>::const_iterator ic = subDirVec.begin();
       ic != subDirVec.end(); ic++) {
    bei->cd(*ic);
    printAlarmList(bei, str_val);
    bei->goUp();
  }
  str_val << "</ul> "<< endl;  
  str_val << "</li> "<< endl;  
//   cout << ACGreen << ACBold
//        << "[SiPixelInformationExtractor::printAlarmList()]"
//        << ACPlain
//        << " Done" 
//        << endl ;
//cout<<"leaving SiPixelInformationExtractor::printAlarmList"<<endl;
}


//------------------------------------------------------------------------------
/*! \brief (Documentation under construction).
 *
 *  This method 
 */
void SiPixelInformationExtractor::getItemList(const multimap<string, string>& req_map, 
                                              string item_name,
					      vector<string>& items) {
//cout<<"entering SiPixelInformationExtractor::getItemList"<<endl;
  items.clear();
  for (multimap<string, string>::const_iterator it = req_map.begin();
       it != req_map.end(); it++) {
    if (it->first == item_name) {
      items.push_back(it->second);
    }
  }
//cout<<"leaving SiPixelInformationExtractor::getItemList"<<endl;
}

//------------------------------------------------------------------------------
/*! \brief (Documentation under construction).
 *
 *  This method 
 */
bool SiPixelInformationExtractor::hasItem(multimap<string,string>& req_map,
					  string item_name){
//cout<<"entering SiPixelInformationExtractor::hasItem"<<endl;
  multimap<string,string>::iterator pos = req_map.find(item_name);
  if (pos != req_map.end()) return true;
  return false;  
//cout<<"leaving SiPixelInformationExtractor::hasItem"<<endl;
}

//------------------------------------------------------------------------------
/*! \brief (Documentation under construction).
 *
 *  This method 
 */
std::string SiPixelInformationExtractor::getItemValue(const std::multimap<std::string,std::string>& req_map,
						 std::string item_name){
//cout<<"entering SiPixelInformationExtractor::getItemValue for item: "<<item_name<<endl;
  std::multimap<std::string,std::string>::const_iterator pos = req_map.find(item_name);
  std::string value = " ";
  if (pos != req_map.end()) {
    value = pos->second;
  }
//  cout<<"value = "<<value<<endl;
  return value;
//cout<<"leaving SiPixelInformationExtractor::getItemValue"<<endl;
}
std::string SiPixelInformationExtractor::getItemValue(std::multimap<std::string,std::string>& req_map,
						 std::string item_name){
//cout<<"entering SiPixelInformationExtractor::getItemValue for item: "<<item_name<<endl;
  std::multimap<std::string,std::string>::iterator pos = req_map.find(item_name);
  std::string value = " ";
  if (pos != req_map.end()) {
//  cout<<"item found!"<<endl;
    value = pos->second;
  }
//  cout<<"value = "<<value<<endl;
  return value;
//cout<<"leaving SiPixelInformationExtractor::getItemValue"<<endl;
}

//
// -- Get color  name from status
//
void SiPixelInformationExtractor::selectColor(string& col, int status){
  if (status == dqm::qstatus::STATUS_OK)    col = "#00ff00";
  else if (status == dqm::qstatus::WARNING) col = "#ffff00";
  else if (status == dqm::qstatus::ERROR)   col = "#ff0000";
  else if (status == dqm::qstatus::OTHER)   col = "#ffa500";
  else  col = "#0000ff";
}
//
// -- Get Image name from ME
//
void SiPixelInformationExtractor::selectColor(string& col, vector<QReport*>& reports){
  int istat = 999;
  int status = 0;
  for (vector<QReport*>::const_iterator it = reports.begin(); it != reports.end();
       it++) {
    status = (*it)->getStatus();
    if (status > istat) istat = status;
  }
  selectColor(col, status);
}
//
// -- Get Image name from status
//
void SiPixelInformationExtractor::selectImage(string& name, int status){
  if (status == dqm::qstatus::STATUS_OK) name="images/LI_green.gif";
  else if (status == dqm::qstatus::WARNING) name="images/LI_yellow.gif";
  else if (status == dqm::qstatus::ERROR) name="images/LI_red.gif";
  else if (status == dqm::qstatus::OTHER) name="images/LI_orange.gif";
  else  name="images/LI_blue.gif";
}
//
// -- Get Image name from ME
//
void SiPixelInformationExtractor::selectImage(string& name, vector<QReport*>& reports){
  int istat = 999;
  int status = 0;
  for (vector<QReport*>::const_iterator it = reports.begin(); it != reports.end();
       it++) {
    status = (*it)->getStatus();
    if (status > istat) istat = status;
  }
  selectImage(name, status);
}

//
// -- Get a tagged image 
//
void SiPixelInformationExtractor::getIMGCImage(const multimap<string, string>& req_map, 
                                               xgi::Output * out){
  string path = getItemValue(req_map,"Path");
  string image;
  histoPlotter_->getNamedImageBuffer(path, image);

  out->getHTTPResponseHeader().addHeader("Content-Type", "image/png");
  out->getHTTPResponseHeader().addHeader("Pragma", "no-cache");   
  out->getHTTPResponseHeader().addHeader("Cache-Control", "no-store, no-cache, must-revalidate,max-age=0");
  out->getHTTPResponseHeader().addHeader("Expires","Mon, 26 Jul 1997 05:00:00 GMT");
  *out << image;
}

void SiPixelInformationExtractor::getIMGCImage(multimap<string, string>& req_map, 
                                               xgi::Output * out){
  
  string path = getItemValue(req_map,"Path");
  string image;
  histoPlotter_->getNamedImageBuffer(path, image);

  out->getHTTPResponseHeader().addHeader("Content-Type", "image/png");
  out->getHTTPResponseHeader().addHeader("Pragma", "no-cache");   
  out->getHTTPResponseHeader().addHeader("Cache-Control", "no-store, no-cache, must-revalidate,max-age=0");
  out->getHTTPResponseHeader().addHeader("Expires","Mon, 26 Jul 1997 05:00:00 GMT");
  *out << image;

}


//------------------------------------------------------------------------------
/*! \brief (Documentation under construction).
 *
 *  This method 
 */
bool SiPixelInformationExtractor::goToDir(DQMStore* bei, 
                                          string& sname){ 
//cout<<"entering SiPixelInformationExtractor::goToDir"<<endl;
  bei->cd();
  //if(flg) bei->cd("Collector/Collated");
  bei->cd(sname);
  string dirName = bei->pwd();
  if (dirName.find(sname) != string::npos) return true;
  else return false;  
//cout<<"leaving SiPixelInformationExtractor::goToDir"<<endl;
}

//
// -- Get Warning/Error Messages
//
void SiPixelInformationExtractor::readStatusMessage(DQMStore* bei, 
                                                    std::multimap<std::string, std::string>& req_map, 
						    xgi::Output * out){

  string path = getItemValue(req_map,"Path");

  int width  = atoi(getItemValue(req_map, "width").c_str());
  int height = atoi(getItemValue(req_map, "height").c_str());

  string opt =" ";

  ostringstream test_status;
  
  setXMLHeader(out);
  *out << "<StatusAndPath>" << endl;
  *out << "<PathList>" << endl;
  if (path.size() == 0) {
    *out << "<HPath>" << "NONE" << "</HPath>" << endl;     
    test_status << " ME Does not exist ! " << endl;
  } else {
    vector<MonitorElement*> all_mes = bei->getContents(path);
    *out << "<HPath>" << path << "</HPath>" << endl;     
    for(vector<MonitorElement*>::iterator it=all_mes.begin(); it!=all_mes.end(); it++){
      MonitorElement* me = (*it);
      if (!me) continue;
      string name = me->getName();  

      vector<QReport*> q_reports = me->getQReports();
      if (q_reports.size() == 0) continue;
      string full_path = path + "/" + name;
      histoPlotter_->setNewPlot(full_path, opt, width, height);

      if (q_reports.size() != 0) {
        test_status << " QTest Status for " << name << " : " << endl;
        test_status << " ======================================================== " << endl; 
        for (vector<QReport*>::const_iterator it = q_reports.begin(); it != q_reports.end();
	     it++) {
	  int status = (*it)->getStatus();
	  if (status == dqm::qstatus::WARNING) test_status << " Warning ";
	  else if (status == dqm::qstatus::ERROR) test_status << " Error  ";
	  else if (status == dqm::qstatus::STATUS_OK) test_status << " Ok  ";
	  else if (status == dqm::qstatus::OTHER) test_status << " Other(" << status << ") ";
	  string mess_str = (*it)->getMessage();
	  test_status <<  "&lt;br/&gt;";
	  mess_str = mess_str.substr(mess_str.find(" Test")+5);
	  test_status <<  " QTest Name  : " << mess_str.substr(0, mess_str.find(")")+1) << endl;
	  test_status << "&lt;br/&gt;";
	  test_status <<  " QTest Detail  : " << mess_str.substr(mess_str.find(")")+2) << endl;
	} 
	test_status << " ======================================================== " << endl;
      }
      *out << "<HPath>" << name << "</HPath>" << endl;         
    }    
  }
  *out << "</PathList>" << endl;
  *out << "<StatusList>" << endl;
  *out << "<Status>" << test_status.str() << "</Status>" << endl;      
  *out << "</StatusList>" << endl;
  *out << "</StatusAndPath>" << endl;
}

//------------------------------------------------------------------------------
/*! \brief (Documentation under construction).
 *  
 */
void SiPixelInformationExtractor::computeStatus(MonitorElement      * theME,
                                                double              & colorValue,
						pair<double,double> & norm) 
{
  double normalizationX = 1 ;
  double normalizationY = 1 ;
  double meanX          = 0 ;
  double meanY          = 0 ;
  
  colorValue = 0 ;

  pair<double,double> normX ;
  pair<double,double> normY ;

  QString theMEType = getMEType(theME) ;

//   cout << ACRed << ACReverse
//        << "[SiPixelInformationExtractor::computeStatus()]"
//        << ACPlain
//        << " Computing average for "
//        << theME->getName()
//        << endl ;

  if( theMEType.contains("TH1") > 0 )
  {
   meanX = (double)theME->getMean();
   getNormalization(theME, normX, "TH1") ;
   normalizationX = fabs( normX.second - normX.first) ;
   if( normalizationX == 0 ) {normalizationX=1.E-20;}
   colorValue  = meanX / normalizationX ;
   norm.first  = normX.first ;
   norm.second = normX.second ;
  }
  
  if( theMEType.contains("TH2") > 0 )
  {
   meanX = (double)theME->getMean(1);
   meanY = (double)theME->getMean(2);
   getNormalization2D(theME, normX, normY, "TH2") ;
   normalizationX = fabs( normX.second - normX.first) ;
   normalizationY = fabs( normY.second - normY.first) ;
   if( normalizationX == 0 ) {normalizationX=1.E-20;}
   if( normalizationY == 0 ) {normalizationY=1.E-20;}
   double cVX = meanX / normalizationX ;
   double cVY = meanY / normalizationY ;
   colorValue = sqrt(cVX*cVX + cVY*cVY) ;
   if( normalizationX >= normalizationY )
   { 
    norm.first  = normX.first;
    norm.second = normX.second ;
   } else { 
    norm.first  = normY.first;
    norm.second = normY.second ;
   }
//   cout << ACBlue << ACBold << ACReverse
//        << "[SiPixelInformationExtractor::computeStatus()]"
//	<< ACPlain << "    "
//	<< theME->getName()
//	<< " meanX:Y "
//	<< meanX << ":" << meanY
//	<< " normX:Y " 
//	<< norm.first << ":" << norm.second
//	<< endl ;
  } 
 
  return ;
}

//------------------------------------------------------------------------------
/*! \brief (Documentation under construction).
 *  
 */
void SiPixelInformationExtractor::getNormalization(MonitorElement     * theME, 
                                                   pair<double,double>& norm,
						   QString              theMEType) 
{
  double normLow  = 0 ;
  double normHigh = 0 ;

  if( theMEType.contains("TH1") > 0 )
  {
   normHigh    = (double)theME->getNbinsX() ;
   norm.first  = normLow  ;
   norm.second = normHigh ;
  }
}

//------------------------------------------------------------------------------
/*! \brief (Documentation under construction).
 *  
 */
void SiPixelInformationExtractor::getNormalization2D(MonitorElement     * theME, 
                                                     pair<double,double>& normX,
                                                     pair<double,double>& normY,
						     QString              theMEType) 
{
  double normLow  = 0 ;
  double normHigh = 0 ;

  if( theMEType.contains("TH2") > 0 )
  {
   normHigh    = (double)theME->getNbinsX() ;
   normX.first  = normLow  ;
   normX.second = normHigh ;
   normHigh    = (double)theME->getNbinsY() ;
   normY.first  = normLow  ;
   normY.second = normHigh ;
//   cout << ACCyan << ACBold << ACReverse
//        << "[SiPixelInformationExtractor::getNormalization2D()]"
//	<< ACPlain << " "
//	<< theME->getName()
//	<< " normX: " 
//	<< normX.first << ":" << normX.second
//	<< " normY: " 
//	<< normY.first << ":" << normY.second
//	<< endl ;
  }
}

//------------------------------------------------------------------------------
/*! \brief (Documentation under construction).
 *
 *   
 */
void SiPixelInformationExtractor::selectMEList(DQMStore   * bei,  
					       string	               & theMEName,
					       vector<MonitorElement*> & mes) 
{  
  string currDir = bei->pwd();
   
  QRegExp rx("(\\w+)_(siPixel|ctfWithMaterialTracks)") ;
  QString theME ;
   
  // Get ME from Collector/FU0/Tracker/PixelEndcap/HalfCylinder_pX/Disk_X/Blade_XX/Panel_XX/Module_XX
  if (currDir.find("Module_") != string::npos ||
      currDir.find("FED_") != string::npos)  
  {
    vector<string> contents = bei->getMEs(); 
       
    for (vector<string>::const_iterator it = contents.begin(); it != contents.end(); it++) 
    {
      theME = *it ;
      if( rx.search(theME) == -1 ) {continue ;} // If the ME is not a siPixel or ctfWithMaterialTrack one, skip
      if (rx.cap(1).latin1() == theMEName)  
      {
        string full_path = currDir + "/" + (*it);

        MonitorElement * me = bei->get(full_path.c_str());
	
        if (me) {mes.push_back(me);}
      }
    }
    return;
  } else {  // If not yet reached the desired level in the directory tree, recursively go down one level more
    vector<string> subdirs = bei->getSubdirs();
    for (vector<string>::const_iterator it = subdirs.begin(); it != subdirs.end(); it++) 
    {
      bei->cd(*it);
      selectMEList(bei, theMEName, mes);
      bei->goUp();
    }
  }
}

//------------------------------------------------------------------------------
/*! \brief (Documentation under construction).
 *  
 */
void SiPixelInformationExtractor::sendTkUpdatedStatus(DQMStore  * bei, 
                                                      xgi::Output            * out,
						      std::string            & theMEName,
						      std::string            & theTKType) 
{
  int rval, gval, bval;
  vector<string>          colorMap ;
  vector<MonitorElement*> me_list;
  pair<double,double>     norm ;
  double sts ;
    
  bei->cd();
  selectMEList(bei, theMEName, me_list) ;
  bei->cd();

  string detId = "undefined";

  QRegExp rx("\\w+_\\w+_(\\d+)");
  
//   cout << ACYellow << ACBold
//        << "[SiPixelInformationExtractor::sendTkUpdatedStatus()] "
//        << ACPlain
//        << "Preparing color map update for " 
//        << theMEName
//        << " type "
//        << theTKType
//        << " - List size: "
//        << me_list.size() 
//        << endl ;
  
  int maxEntries = 0 ;
  if( theTKType == "Entries") // In this case find the ME with the highest number of entries
  {			      // first and use that as a vertical scale normalization
   for(vector<MonitorElement*>::iterator it=me_list.begin(); it!=me_list.end(); it++)
   {
    int entries = (int)(*it)->getEntries() ;
    if( entries > maxEntries ) maxEntries = entries ;
   }
  }
  
  int entries = 0 ;
  stringstream jsSnippet ;
  for(vector<MonitorElement*>::iterator it=me_list.begin(); it!=me_list.end(); it++)
  {
    QString meName    = (*it)->getName() ;
    QString theMEType = getMEType(*it) ;
    if( rx.search(meName) != -1 ) 
    {
     detId = rx.cap(1).latin1() ;
     entries = (int)(*it)->getEntries() ;
     if( theTKType == "Averages") 
     {
      computeStatus(*it, sts, norm) ;
      SiPixelUtility::getStatusColor(sts, rval, gval, bval);
     } else if( theTKType == "Entries") {
      sts = (double)entries / (double)maxEntries ;
      SiPixelUtility::getStatusColor(sts, rval, gval, bval);
      if( entries > maxEntries ) maxEntries = entries ;
      norm.first  = 0 ;
      norm.second = maxEntries ;
     } else {
      int status  =  SiPixelUtility::getStatus((*it));
      if(        status == dqm::qstatus::ERROR ) 
      {
       rval = 255; gval =   0; bval =   0;
      } else if (status == dqm::qstatus::WARNING )  {
       rval = 255; gval = 255; bval =   0; 
      } else if (status == dqm::qstatus::OTHER)     {
       rval =   0; gval =   0; bval = 255;
      } else if (status == dqm::qstatus::STATUS_OK) {
       rval =   0; gval = 255; bval =   0;
      } else {  
       rval = 255; gval = 255; bval = 255;
      }
     }
     jsSnippet.str("") ;
     jsSnippet << " <DetInfo DetId='"
     	       << detId
     	       << "' red='"
     	       << rval
     	       << "' green='"
     	       << gval
     	       << "' blue='"
     	       << bval
     	       << "' entries='"
     	       << entries
     	       << "'/>" ;
     colorMap.push_back(jsSnippet.str()) ;
//      if( it == me_list.begin()) // The first should be equal to all others...
//      {
//       getNormalization((*it), norm, theMEType.latin1()) ;
//      }
    }
  }

//  delete random ;
  
//   cout << ACYellow << ACBold
//        << "[SiPixelInformationExtractor::sendTkUpdatedStatus()] "
//        << ACPlain
//        << "Color map consists of "
//        << colorMap.size()
//        << " snippets: start shipping back"
//        << endl ;

  out->getHTTPResponseHeader().addHeader("Content-Type", "text/xml");
  *out << "<?xml version=\"1.0\" ?>" << endl;
  *out << "<TrackerMapUpdate>"       << endl;

  for(vector<string>::iterator it=colorMap.begin(); it!=colorMap.end(); it++)
  {
   *out << *it << endl;
  }
 
  *out << " <theLimits id=\"normalizationLimits\" normLow=\"" 
       << norm.first 
       << "\" normHigh=\""
       << norm.second 
       << "\" />"
       << endl;
  *out << "</TrackerMapUpdate>"              
       << endl;

//   cout << ACYellow << ACBold
//        << "[SiPixelInformationExtractor::sendTkUpdatedStatus()] "
//        << ACPlain
//        << "Color map updated within range " 
//        << norm.first
//        << "-"
//        << norm.second
//        << endl ;
}

//------------------------------------------------------------------------------
/*! \brief (Documentation under construction).
 *  
 *  Given a pointer to ME returns the associated detId 
 */
int SiPixelInformationExtractor::getDetId(MonitorElement * mE) 
{
 QRegExp rx("(\\w+)_(\\w+)_(\\d+)") ;
 QString mEName = mE->getName() ;

 int detId = 0;
 
 if( rx.search(mEName) != -1 )
 {
  detId = rx.cap(3).toInt() ;
 } else {
  cout << ACYellow << ACBold
       << "[SiPixelInformationExtractor::getDetId()] "
       << ACPlain
       << "Could not extract detId from "
       << mEName
       << endl ;
 }
      
  return detId ;
  
}

//------------------------------------------------------------------------------
/*! \brief (Documentation under construction).
 *  
 */
void SiPixelInformationExtractor::getMEList(DQMStore    * bei,  
					    map<string, int>         & mEHash)
{
  string currDir = bei->pwd();
   
//   cout << ACRed << ACBold
//        << "[SiPixelInformationExtractor::getMEList()]"
//        << ACPlain
//        << " Requesting ME list in " 
//        << currDir
//        << endl ;
       
  QRegExp rx("(\\w+)_(siPixel|ctfWithMaterialTracks)") ;
  QString theME ;
   
  // Get ME from Collector/FU0/Tracker/PixelEndcap/HalfCylinder_pX/Disk_X/Blade_XX/Panel_XX/Module_XX
  if (currDir.find("Module_") != string::npos ||
      currDir.find("FED_") != string::npos)  
  {
    vector<string> contents = bei->getMEs(); 
       
    for (vector<string>::const_iterator it = contents.begin(); it != contents.end(); it++) 
    {
      theME = *it ;
//       cout << ACRed << ACReverse
//            << "[SiPixelInformationExtractor::getMEList()]"
//            << ACPlain
//            << " ME: " 
//            << theME
//            << endl ;
      if( rx.search(theME) == -1 ) 
      {
       cout << ACRed << ACBold
            << "[SiPixelInformationExtractor::getMEList()]"
	    << ACPlain
	    << " ----> Skipping " 
	    << theME
	    << endl ;
       continue ;
      } // If the ME is not a Pixel one, skip
      string full_path = currDir + "/" + (*it);
      string mEName = rx.cap(1).latin1() ;
      mEHash[mEName]++ ;
    }
    
    return;
  } else {  // If not yet reached the desired level in the directory tree, recursively go down one level more
    vector<string> subdirs = bei->getSubdirs();
    for (vector<string>::const_iterator it = subdirs.begin(); it != subdirs.end(); it++) 
    {
      bei->cd(*it);
      getMEList(bei, mEHash);
      bei->goUp();
    }
  }
}

//
// -- Get All histograms from a Path
//
void SiPixelInformationExtractor::getHistosFromPath(DQMStore * bei, 
                                                    const std::multimap<std::string, std::string>& req_map, 
						    xgi::Output * out){
//cout<<"Entering SiPixelInformationExtractor::getHistosFromPath: "<<endl;
  string path = getItemValue(req_map,"Path");
//cout<<"Path is: "<<path<<endl;
  if (path.size() == 0) return;

  int width  = atoi(getItemValue(req_map, "width").c_str());
  int height = atoi(getItemValue(req_map, "height").c_str());

  string opt =" ";

  setHTMLHeader(out);
  vector<MonitorElement*> all_mes = bei->getContents(path);
  *out << path << " " ;
  for(vector<MonitorElement*>::iterator it=all_mes.begin(); it!=all_mes.end(); it++){
    MonitorElement* me = (*it);
    //cout<<"I'm in the loop now..."<<endl;
    if (!me) continue;
    string name = me->getName();
    string full_path = path + "/" + name;
//cout<<"Calling HP::setNewPlot now for "<<full_path<<endl;
    histoPlotter_->setNewPlot(full_path, opt, width, height);
    *out << name << " ";
  }
//  cout<<"... leaving SiPixelInformationExtractor::getHistosFromPath!"<<endl;
}

void SiPixelInformationExtractor::bookGlobalQualityFlag(DQMStore * bei) {
//std::cout<<"BOOK GLOBAL QUALITY FLAG MEs!"<<std::endl;
  bei->cd();
  bei->setCurrentFolder("Pixel/EventInfo");
  SummaryReport = bei->bookFloat("reportSummary");
  //SummaryReportMap = bei->book2D("reportSummaryMap","Pixel EtaPhi Summary Map",60,-3.,3.,64,-3.2,3.2);
  SummaryReportMap = bei->book2D("reportSummaryMap","Pixel Summary Map",28,0.,28.,48,0.,48.);
  SummaryReportMap->setAxisTitle("Eta",1);
  SummaryReportMap->setAxisTitle("Phi",2);
  SummaryReportMap->setBinLabel(1,"Endcap",1);
  SummaryReportMap->setBinLabel(2,"+z Disk_2",1);
  SummaryReportMap->setBinLabel(6,"Endcap",1);
  SummaryReportMap->setBinLabel(7,"+z Disk_1",1);
  SummaryReportMap->setBinLabel(12,"Barrel +z",1);
  SummaryReportMap->setBinLabel(16,"Barrel -z",1);
  SummaryReportMap->setBinLabel(20,"Endcap",1);
  SummaryReportMap->setBinLabel(21,"-z Disk_1",1);
  SummaryReportMap->setBinLabel(25,"Endcap",1);
  SummaryReportMap->setBinLabel(26,"-z Disk_2",1);
  SummaryReportMap->setBinLabel(13,"Inside",2);
  SummaryReportMap->setBinLabel(12,"LHC",2);
  SummaryReportMap->setBinLabel(11,"Ring",2);
  SummaryReportMap->setBinLabel(37,"Outside",2);
  SummaryReportMap->setBinLabel(36,"LHC",2);
  SummaryReportMap->setBinLabel(35,"Ring",2);
  bei->setCurrentFolder("Pixel/EventInfo/reportSummaryContents");
  SummaryBarrel = bei->bookFloat("Pixel_Barrel");
  SummaryShellmI = bei->bookFloat("Pixel_Shell_mI");
  SummaryShellmO = bei->bookFloat("Pixel_Shell_mO");
  SummaryShellpI = bei->bookFloat("Pixel_Shell_pI");
  SummaryShellpO = bei->bookFloat("Pixel_Shell_pO");
  SummaryEndcap = bei->bookFloat("Pixel_Endcap");
  SummaryHCmI = bei->bookFloat("Pixel_HalfCylinder_mI");
  SummaryHCmO = bei->bookFloat("Pixel_HalfCylinder_mO");
  SummaryHCpI = bei->bookFloat("Pixel_HalfCylinder_pI");
  SummaryHCpO = bei->bookFloat("Pixel_HalfCylinder_pO");
}

void SiPixelInformationExtractor::computeGlobalQualityFlag(DQMStore * bei, 
                                                           bool init)
{
//cout<<"entering SiPixelInformationExtractor::ComputeGlobalQualityFlag"<<endl;
//   cout << ACRed << ACBold
//        << "[SiPixelInformationExtractor::ComputeGlobalQualityFlag]"
//        << ACPlain
//        << " Enter" 
//        << endl ;
  if(init){
    gotDigis=false;
    allMods_=0; errorMods_=0; qflag_=0.; 
    bpix_mods_=0; err_bpix_mods_=0; bpix_flag_=0.;
    shellmI_mods_=0; err_shellmI_mods_=0; shellmI_flag_=0.;
    shellmO_mods_=0; err_shellmO_mods_=0; shellmO_flag_=0.;
    shellpI_mods_=0; err_shellpI_mods_=0; shellpI_flag_=0.;
    shellpO_mods_=0; err_shellpO_mods_=0; shellpO_flag_=0.;
    fpix_mods_=0; err_fpix_mods_=0; fpix_flag_=0.;
    hcylmI_mods_=0; err_hcylmI_mods_=0; hcylmI_flag_=0.;
    hcylmO_mods_=0; err_hcylmO_mods_=0; hcylmO_flag_=0.;
    hcylpI_mods_=0; err_hcylpI_mods_=0; hcylpI_flag_=0.;
    hcylpO_mods_=0; err_hcylpO_mods_=0; hcylpO_flag_=0.;
  }
  
  string currDir = bei->pwd();
  string dname = currDir.substr(currDir.find_last_of("/")+1);
  
  QRegExp rx("Module_");
 
  if(rx.search(dname)!=-1){
    if(currDir.find("Pixel")!=string::npos) allMods_++;
    if(currDir.find("Barrel")!=string::npos) bpix_mods_++;
    if(currDir.find("Shell_mI")!=string::npos) shellmI_mods_++;
    if(currDir.find("Shell_mO")!=string::npos) shellmO_mods_++;
    if(currDir.find("Shell_pI")!=string::npos) shellpI_mods_++;
    if(currDir.find("Shell_pO")!=string::npos) shellpO_mods_++;
    if(currDir.find("Endcap")!=string::npos) fpix_mods_++;
    if(currDir.find("HalfCylinder_mI")!=string::npos) hcylmI_mods_++;
    if(currDir.find("HalfCylinder_mO")!=string::npos) hcylmO_mods_++;
    if(currDir.find("HalfCylinder_pI")!=string::npos) hcylpI_mods_++;
    if(currDir.find("HalfCylinder_pO")!=string::npos) hcylpO_mods_++;
      
    vector<string> meVec = bei->getMEs();
    //checking for any digis anywhere to decide if Pixel detector is in DAQ:  
    for (vector<string>::const_iterator it = meVec.begin(); it != meVec.end(); it++) {
      string full_path = currDir + "/" + (*it);
      if(full_path.find("ndigis")!=string::npos){
        MonitorElement * me = bei->get(full_path);
        if (!me) continue;
        if(me->getEntries()>0) gotDigis = true;
      }
    }
      
    //checking for FED errors only:
    for (vector<string>::const_iterator it = meVec.begin(); it != meVec.end(); it++) {
      string full_path = currDir + "/" + (*it);
      if(full_path.find("NErrors")!=string::npos){
        MonitorElement * me = bei->get(full_path);
        if (!me) continue;
	
	//if(me->getEntries()>0) cout<<"I am here: "<<currDir<<" and have the NErrors-ME. It has this many entries: "<<me->getEntries()<<endl;
	
        if(me->getEntries()>0){
	
	  full_path = full_path.replace(full_path.find("NErrors"),9,"errorType");
	  me = bei->get(full_path);
	  if(!me) continue;
	  bool type30=false; bool othererror=false; bool reset=false;
	  for(int jj=1; jj<16; jj++){
	    if(me->getBinContent(jj)>0.){
	      if(jj!=6) othererror=true;
              else type30=true;
	    }
	  }
	  if(type30){
	    full_path = full_path.replace(full_path.find("errorType"),10,"TBMMessage");
	    me = bei->get(full_path);
	    if(!me) continue;
	    for(int kk=1; kk<9; kk++){
              if(me->getBinContent(kk)>0.){
		if(kk!=6 && kk!=7) othererror=true;
		else reset=true;
	      }
	    }
	  }
	  
    //if you want to check for QTest results instead:
    //vector<string> meVec = bei->getMEs();
    //bool gotcha = false;
    //for (vector<string>::const_iterator it = meVec.begin();
    //	 it != meVec.end(); it++) {
    //  string full_path = currDir + "/" + (*it);
    // MonitorElement * me = bei->get(full_path);
    // if (!me) continue;
    // std::vector<QReport *> my_map = me->getQReports();
    // if (my_map.size() > 0) {
       // string image_name;
       // selectImage(image_name,my_map);
       // if(image_name!="images/LI_green.gif") {
       
       //   if(!gotcha){
	  if(othererror || (type30 && !reset)){
            if(currDir.find("Pixel")!=string::npos) errorMods_++;
            if(currDir.find("Barrel")!=string::npos) err_bpix_mods_++;
            if(currDir.find("Shell_mI")!=string::npos) err_shellmI_mods_++;
            if(currDir.find("Shell_mO")!=string::npos) err_shellmO_mods_++;
            if(currDir.find("Shell_pI")!=string::npos) err_shellpI_mods_++;
            if(currDir.find("Shell_pO")!=string::npos) err_shellpO_mods_++;
            if(currDir.find("Endcap")!=string::npos) err_fpix_mods_++;
            if(currDir.find("HalfCylinder_mI")!=string::npos) err_hcylmI_mods_++;
            if(currDir.find("HalfCylinder_mO")!=string::npos) err_hcylmO_mods_++;
            if(currDir.find("HalfCylinder_pI")!=string::npos) err_hcylpI_mods_++;
            if(currDir.find("HalfCylinder_pO")!=string::npos) err_hcylpO_mods_++;
	  }
	//  gotcha = true;
        }	
      }
    }
  }
  if(allMods_>0) qflag_ = (float(allMods_)-float(errorMods_))/float(allMods_);
  if(bpix_mods_>0) bpix_flag_ = (float(bpix_mods_)-float(err_bpix_mods_))/float(bpix_mods_);
  if(shellmI_mods_>0) shellmI_flag_ = (float(shellmI_mods_)-float(err_shellmI_mods_))/float(shellmI_mods_);
  if(shellmO_mods_>0) shellmO_flag_ = (float(shellmO_mods_)-float(err_shellmO_mods_))/float(shellmO_mods_);
  if(shellpI_mods_>0) shellpI_flag_ = (float(shellpI_mods_)-float(err_shellpI_mods_))/float(shellpI_mods_);
  if(shellpO_mods_>0) shellpO_flag_ = (float(shellpO_mods_)-float(err_shellpO_mods_))/float(shellpO_mods_);
  if(fpix_mods_>0) fpix_flag_ = (float(fpix_mods_)-float(err_fpix_mods_))/float(fpix_mods_);
  if(hcylmI_mods_>0) hcylmI_flag_ = (float(hcylmI_mods_)-float(err_hcylmI_mods_))/float(hcylmI_mods_);
  if(hcylmO_mods_>0) hcylmO_flag_ = (float(hcylmO_mods_)-float(err_hcylmO_mods_))/float(hcylmO_mods_);
  if(hcylpI_mods_>0) hcylpI_flag_ = (float(hcylpI_mods_)-float(err_hcylpI_mods_))/float(hcylpI_mods_);
  if(hcylpO_mods_>0) hcylpO_flag_ = (float(hcylpO_mods_)-float(err_hcylpO_mods_))/float(hcylpO_mods_);
  if(!gotDigis){
    qflag_ = -1.;
    bpix_flag_ = -1.;
    shellmI_flag_ = -1.;
    shellmO_flag_ = -1.;
    shellpI_flag_ = -1.;
    shellpO_flag_ = -1.;
    fpix_flag_ = -1.;
    hcylmI_flag_ = -1.;
    hcylmO_flag_ = -1.;
    hcylpI_flag_ = -1.;
    hcylpO_flag_ = -1.;
  }
  
  vector<string> subDirVec = bei->getSubdirs();  
  for (vector<string>::const_iterator ic = subDirVec.begin();
       ic != subDirVec.end(); ic++) {
    bei->cd(*ic);
    init=false;
    computeGlobalQualityFlag(bei,init);
    bei->goUp();
  }
  SummaryReport = bei->get("Pixel/EventInfo/reportSummary");
  if(SummaryReport) SummaryReport->Fill(qflag_);
  SummaryBarrel = bei->get("Pixel/EventInfo/reportSummaryContents/Pixel_Barrel");
  if(SummaryBarrel) SummaryBarrel->Fill(bpix_flag_);
  SummaryShellmI = bei->get("Pixel/EventInfo/reportSummaryContents/Pixel_Shell_mI");
  if(SummaryShellmI) SummaryShellmI->Fill(shellmI_flag_);
  SummaryShellmO = bei->get("Pixel/EventInfo/reportSummaryContents/Pixel_Shell_mO");
  if(SummaryShellmO)   SummaryShellmO->Fill(shellmO_flag_);
  SummaryShellpI = bei->get("Pixel/EventInfo/reportSummaryContents/Pixel_Shell_pI");
  if(SummaryShellpI)   SummaryShellpI->Fill(shellpI_flag_);
  SummaryShellpO = bei->get("Pixel/EventInfo/reportSummaryContents/Pixel_Shell_pO");
  if(SummaryShellpO)   SummaryShellpO->Fill(shellpO_flag_);
  SummaryEndcap = bei->get("Pixel/EventInfo/reportSummaryContents/Pixel_Endcap");
  if(SummaryEndcap)   SummaryEndcap->Fill(fpix_flag_);
  SummaryHCmI = bei->get("Pixel/EventInfo/reportSummaryContents/Pixel_HalfCylinder_mI");
  if(SummaryHCmI)   SummaryHCmI->Fill(hcylmI_flag_);
  SummaryHCmO = bei->get("Pixel/EventInfo/reportSummaryContents/Pixel_HalfCylinder_mO");
  if(SummaryHCmO)   SummaryHCmO->Fill(hcylmO_flag_);
  SummaryHCpI = bei->get("Pixel/EventInfo/reportSummaryContents/Pixel_HalfCylinder_pI");
  if(SummaryHCpI)   SummaryHCpI->Fill(hcylpI_flag_);
  SummaryHCpO = bei->get("Pixel/EventInfo/reportSummaryContents/Pixel_HalfCylinder_pO");
  if(SummaryHCpO)   SummaryHCpO->Fill(hcylpO_flag_);

}

void SiPixelInformationExtractor::fillGlobalQualityPlot(DQMStore * bei, bool init, edm::EventSetup const& eSetup)
{
  //calculate eta and phi of the modules and fill a 2D plot:
  
  if(init){
    allmodsMap = new TH2F("allmodsMap","allmodsMap",28,0.,28.,48,0.,48.);
    errmodsMap = new TH2F("errmodsMap","errmodsMap",28,0.,28.,48,0.,48.);
    goodmodsMap = new TH2F("goodmodsMap","goodmodsMap",28,0.,28.,48,0.,48.);
    count=0; errcount=0;
    init=false;
  }
  string currDir = bei->pwd();
  string dname = currDir.substr(currDir.find_last_of("/")+1);
  
  QRegExp rx("Module_");
 
  if(rx.search(dname)!=-1){
    vector<string> meVec = bei->getMEs();
    float detEta=-5.; float detPhi=-5.;
    bool first=true; bool once=true;
    int xbin = -1; int ybin =-1; int xoffA = 0; int xoffB = 0; int yoffA = 0; int yoffB = 0;
    bool gotData=false;
    for (vector<string>::const_iterator it = meVec.begin(); it != meVec.end(); it++) {
      //checking for any digis or FED errors to decide if this module is in DAQ:  
      string full_path = currDir + "/" + (*it);
      if(full_path.find("ndigis")!=string::npos || full_path.find("NErrors")!=string::npos){
        MonitorElement * me = bei->get(full_path);
        if (!me) continue;
        if(me->getMean()>0.) gotData = true;
      }
      
      if(!once) continue;
      MonitorElement * me = bei->get(full_path);
      if (!me) continue;
/*      int id=0;
      if(first){ id = getDetId(me); first=false; }
      DetId detid = DetId(id);
      if(detid.det()!=1) continue;
      edm::ESHandle<TrackerGeometry> pDD;
      eSetup.get<TrackerDigiGeometryRecord>().get( pDD );
      for(TrackerGeometry::DetContainer::const_iterator it = 
	  pDD->dets().begin(); it != pDD->dets().end(); it++){
        if(dynamic_cast<PixelGeomDetUnit*>((*it))!=0){
          DetId detId = (*it)->geographicalId();
	  if(detId!=detid) continue;
	  //if(detId.subdetId()!=1) continue;
          const GeomDetUnit * geoUnit = pDD->idToDetUnit( detId );
          const PixelGeomDetUnit * pixDet = dynamic_cast<const PixelGeomDetUnit*>(geoUnit);
	  float detR = pixDet->surface().position().perp();
	  float detZ = pixDet->surface().position().z();
	  detPhi = pixDet->surface().position().phi();
	  detEta = -1.*log(tan(atan2(detR,detZ)/2.));
	  //cout<<"Module: "<<currDir<<" , Eta= "<<detEta<<" , Phi= "<<detPhi<<endl;
	  once=false;
	}
      }*/
      
      if(full_path.find("Endcap")!=string::npos){ // Endcaps
        if(full_path.find("_m")!=string::npos) xoffA = 19; // the -z endcaps are on the right hand side in x
	if(full_path.find("_m")!=string::npos && full_path.find("Disk_2")!=string::npos) xoffB = 5; // on -z Disk_2 is right of Disk_1
	if(full_path.find("_p")!=string::npos && full_path.find("Disk_1")!=string::npos) xoffB = 5; // on +z Disk_2 is left of Disk_1
	if(full_path.find("_p")!=string::npos){
	  if(full_path.find("Module_4")!=string::npos) xbin = 1+xoffA+xoffB;
	  else if(full_path.find("Module_3")!=string::npos) xbin = 2+xoffA+xoffB;
	  else if(full_path.find("Module_2")!=string::npos) xbin = 3+xoffA+xoffB;
	  else if(full_path.find("Module_1")!=string::npos) xbin = 4+xoffA+xoffB;
	}else if(full_path.find("_m")!=string::npos){
	  if(full_path.find("Module_1")!=string::npos) xbin = 1+xoffA+xoffB;
	  else if(full_path.find("Module_2")!=string::npos) xbin = 2+xoffA+xoffB;
	  else if(full_path.find("Module_3")!=string::npos) xbin = 3+xoffA+xoffB;
	  else if(full_path.find("Module_4")!=string::npos) xbin = 4+xoffA+xoffB;
	}
	if(full_path.find("Panel_2")!=string::npos) yoffA = 1; // Panel_2's are shifted by 1 bin in y w.r.t. their Panel_1 partner
	if(full_path.find("HalfCylinder_mO")!=string::npos || full_path.find("HalfCylinder_pO")!=string::npos) yoffB = 24; // outer HC's have offset off 24 bins in y
	if(full_path.find("Blade_01")!=string::npos) ybin = 1+yoffA+yoffB;
	else if(full_path.find("Blade_02")!=string::npos) ybin = 3+yoffA+yoffB;
	else if(full_path.find("Blade_03")!=string::npos) ybin = 5+yoffA+yoffB;
	else if(full_path.find("Blade_04")!=string::npos) ybin = 7+yoffA+yoffB;
	else if(full_path.find("Blade_05")!=string::npos) ybin = 9+yoffA+yoffB;
	else if(full_path.find("Blade_06")!=string::npos) ybin = 11+yoffA+yoffB;
	else if(full_path.find("Blade_07")!=string::npos) ybin = 13+yoffA+yoffB;
	else if(full_path.find("Blade_08")!=string::npos) ybin = 15+yoffA+yoffB;
	else if(full_path.find("Blade_09")!=string::npos) ybin = 17+yoffA+yoffB;
	else if(full_path.find("Blade_10")!=string::npos) ybin = 19+yoffA+yoffB;
	else if(full_path.find("Blade_11")!=string::npos) ybin = 21+yoffA+yoffB;
	else if(full_path.find("Blade_12")!=string::npos) ybin = 23+yoffA+yoffB;
      }else if(full_path.find("Barrel")!=string::npos){ // Barrel
        if(full_path.find("_m")!=string::npos) xoffA = 1; 
	if(full_path.find("Module_4")!=string::npos) xbin = 11+xoffA*7;
	else if(full_path.find("Module_3")!=string::npos) xbin = 12+xoffA*5;
	else if(full_path.find("Module_2")!=string::npos) xbin = 13+xoffA*3;
	else if(full_path.find("Module_1")!=string::npos) xbin = 14+xoffA*1;
	if(full_path.find("Layer_3")!=string::npos) yoffA = 2;
	if(full_path.find("Layer_2")!=string::npos) yoffA = 8;
	if(full_path.find("Layer_1")!=string::npos) yoffA = 14;
	if(full_path.find("_mO")!=string::npos || full_path.find("_pO")!=string::npos) yoffB = 22;
	if(full_path.find("Ladder_01")!=string::npos) ybin = 1+yoffA+yoffB;
	else if(full_path.find("Ladder_02")!=string::npos) ybin = 2+yoffA+yoffB;
	else if(full_path.find("Ladder_03")!=string::npos) ybin = 3+yoffA+yoffB;
	else if(full_path.find("Ladder_04")!=string::npos) ybin = 4+yoffA+yoffB;
	else if(full_path.find("Ladder_05")!=string::npos) ybin = 5+yoffA+yoffB;
	else if(full_path.find("Ladder_06")!=string::npos) ybin = 6+yoffA+yoffB;
	else if(full_path.find("Ladder_07")!=string::npos) ybin = 7+yoffA+yoffB;
	else if(full_path.find("Ladder_08")!=string::npos) ybin = 8+yoffA+yoffB;
	else if(full_path.find("Ladder_09")!=string::npos) ybin = 9+yoffA+yoffB;
	else if(full_path.find("Ladder_10")!=string::npos) ybin = 10+yoffA+yoffB;
	else if(full_path.find("Ladder_11")!=string::npos) ybin = 11+yoffA+yoffB;
	else if(full_path.find("Ladder_12")!=string::npos) ybin = 12+yoffA+yoffB;
	else if(full_path.find("Ladder_13")!=string::npos) ybin = 13+yoffA+yoffB;
	else if(full_path.find("Ladder_14")!=string::npos) ybin = 14+yoffA+yoffB;
	else if(full_path.find("Ladder_15")!=string::npos) ybin = 15+yoffA+yoffB;
	else if(full_path.find("Ladder_16")!=string::npos) ybin = 16+yoffA+yoffB;
	else if(full_path.find("Ladder_17")!=string::npos) ybin = 17+yoffA+yoffB;
	else if(full_path.find("Ladder_18")!=string::npos) ybin = 18+yoffA+yoffB;
	else if(full_path.find("Ladder_19")!=string::npos) ybin = 19+yoffA+yoffB;
	else if(full_path.find("Ladder_20")!=string::npos) ybin = 20+yoffA+yoffB;
	else if(full_path.find("Ladder_21")!=string::npos) ybin = 21+yoffA+yoffB;
	else if(full_path.find("Ladder_22")!=string::npos) ybin = 22+yoffA+yoffB;
      }
      once=false;
		
    } //end of the for loop over all ME's for a given module
    //cout<<"module: "<<bei->pwd()<<" , xbin="<<xbin<<" , ybin="<<ybin<<endl; 
    //got module ID and eta and phi now! time to count:
    
    // only consider those modules that are actually in the readout.
    if(gotData){
      count++;
      //allmodsEtaPhi->Fill(detEta,detPhi);
      allmodsMap->Fill(xbin-1,ybin-1);
      bool anyerr=false;
      for (vector<string>::const_iterator it = meVec.begin(); it != meVec.end(); it++){
        if(anyerr) continue;
        string full_path = currDir + "/" + (*it);
        MonitorElement * me = bei->get(full_path);
        if (!me) continue;
        //use only presence of any FED error as error flag here:
        if(full_path.find("NErrors")!=string::npos && me->getEntries()>0){
	  full_path = full_path.replace(full_path.find("NErrors"),9,"errorType");
	  me = bei->get(full_path);
	  if(!me) anyerr=true;
	  else{
	    bool type30=false;
	    for(int jj=1; jj<16; jj++){
	      if(me->getBinContent(jj)>0.){
	        if(jj!=6) anyerr=true;
		else type30=true;
	      }
	    }
	    if(type30){
	      full_path = full_path.replace(full_path.find("errorType"),10,"TBMMessage");
	      me = bei->get(full_path);
	      if(!me) anyerr=true;
	      else{
	        for(int kk=1; kk<9; kk++){
		  if(me->getBinContent(kk)>0.){
		    if(kk!=6 && kk!=7) anyerr=true;
		    else anyerr=false;
		  }
		}
	      }
	    }
	  }
	}// if NErrors
	      
	      
        //use QTest results for error flag here:
        //if(me->hasError()||me->hasWarning()||me->hasOtherReport()) anyerr=true;
      }
      if(anyerr){
        errcount++;
        //errmodsEtaPhi->Fill(detEta,detPhi);
        errmodsMap->Fill(xbin-1,ybin-1);
      }
    }
  }

  vector<string> subDirVec = bei->getSubdirs();  
  for (vector<string>::const_iterator ic = subDirVec.begin();
       ic != subDirVec.end(); ic++) {
    bei->cd(*ic);
    init=false;
    fillGlobalQualityPlot(bei,init,eSetup);
    bei->goUp();
  }
  
  SummaryReportMap = bei->get("Pixel/EventInfo/reportSummaryMap");
  if(SummaryReportMap){ 
    float contents=0.;
    //for(int i=1; i!=61; i++)for(int j=1; j!=65; j++){
    for(int i=1; i!=29; i++)for(int j=1; j!=49; j++){
      //contents = (allmodsEtaPhi->GetBinContent(i,j))-(errmodsEtaPhi->GetBinContent(i,j));
      contents = (allmodsMap->GetBinContent(i,j))-(errmodsMap->GetBinContent(i,j));
      //cout<<"all: "<<allmodsMap->GetBinContent(i,j)<<" , error: "<<errmodsMap->GetBinContent(i,j)<<" , contents: "<<contents<<endl;
      //goodmodsEtaPhi->SetBinContent(i,j,contents);
      goodmodsMap->SetBinContent(i,j,contents);
      //if(allmodsEtaPhi->GetBinContent(i,j)>0){
      if(allmodsMap->GetBinContent(i,j)>0){
	//contents = (goodmodsEtaPhi->GetBinContent(i,j))/(allmodsEtaPhi->GetBinContent(i,j));
	contents = (goodmodsMap->GetBinContent(i,j))/(allmodsMap->GetBinContent(i,j));
      }else{
        contents = -1.;
      }
      //cout<<"MAP: "<<i<<","<<j<<","<<contents<<endl;
      SummaryReportMap->setBinContent(i,j,contents);
    }
  }
  if(allmodsMap) allmodsMap->Clear();
  if(goodmodsMap) goodmodsMap->Clear();
  if(errmodsMap) errmodsMap->Clear();
  //cout<<"counters: "<<count<<" , "<<errcount<<endl;
}

void SiPixelInformationExtractor::findNoisyPixels(DQMStore * bei)
{
  bei->cd();
  bei->setCurrentFolder("Pixel/EventInfo");
  MonitorElement * me0 = bei->get("Pixel/EventInfo/iEvent");
  int nevents=-1;
  if(me0) nevents = me0->getIntValue(); 
  if(nevents==0) nevents=-1;
  //cout<<"nevents: "<<nevents<<endl;
  //nevents=1000;
  bei->cd();
  bei->setCurrentFolder("Pixel/Endcap");
  EndcapNdigisFREQProjection = bei->book1D("endcapNdigisFREQProjection","Endcap: Digi event rate per module",1000,0.,1.);
  EndcapNdigisFREQProjection = bei->get("Pixel/Endcap/endcapNdigisFREQProjection");
  MonitorElement * me1 = bei->get("Pixel/Endcap/SUMDIG_ndigisFREQ_Endcap");
  if(me1){
    for(int i=1; i!=673; i++){
      //cout<<"entries: "<<me1->getBinContent(i)/float(nevents)<<" in bin "<<i<<endl;
      EndcapNdigisFREQProjection->Fill(me1->getBinContent(i)/float(nevents));
      //if(i>63&&i<90) cout<<"noisy edges: "<<me1->getBinContent(i)/float(nevents)<<endl;
    }
  }
  
  bei->cd();
  bei->setCurrentFolder("Pixel/Barrel");
  BarrelNdigisFREQProjection = bei->book1D("barrelNdigisFREQProjection","Barrel: Digi event rate per module",1000,0.,1.);
  BarrelNdigisFREQProjection = bei->get("Pixel/Barrel/barrelNdigisFREQProjection");
  MonitorElement * me2 = bei->get("Pixel/Barrel/SUMDIG_ndigisFREQ_Barrel");
  if(me2){
    for(int i=1; i!=769; i++){
      BarrelNdigisFREQProjection->Fill(me2->getBinContent(i)/float(nevents));
    }
  }
}


//
// -- Create Images 
//
void SiPixelInformationExtractor::createImages(DQMStore* bei){
  histoPlotter_->createPlots(bei);
}

//
// -- Set HTML Header in xgi output
//
void SiPixelInformationExtractor::setHTMLHeader(xgi::Output * out) {
  out->getHTTPResponseHeader().addHeader("Content-Type", "text/html");
  out->getHTTPResponseHeader().addHeader("Pragma", "no-cache");   
  out->getHTTPResponseHeader().addHeader("Cache-Control", "no-store, no-cache, must-revalidate,max-age=0");
  out->getHTTPResponseHeader().addHeader("Expires","Mon, 26 Jul 1997 05:00:00 GMT");
}
//
// -- Set XML Header in xgi output
//
void SiPixelInformationExtractor::setXMLHeader(xgi::Output * out) {
  out->getHTTPResponseHeader().addHeader("Content-Type", "text/xml");
  out->getHTTPResponseHeader().addHeader("Pragma", "no-cache");   
  out->getHTTPResponseHeader().addHeader("Cache-Control", "no-store, no-cache, must-revalidate,max-age=0");
  out->getHTTPResponseHeader().addHeader("Expires","Mon, 26 Jul 1997 05:00:00 GMT");
  *out << "<?xml version=\"1.0\" ?>" << std::endl;

}
//
// -- Set Plain Header in xgi output
//
void SiPixelInformationExtractor::setPlainHeader(xgi::Output * out) {
  out->getHTTPResponseHeader().addHeader("Content-Type", "text/plain");
  out->getHTTPResponseHeader().addHeader("Pragma", "no-cache");   
  out->getHTTPResponseHeader().addHeader("Cache-Control", "no-store, no-cache, must-revalidate,max-age=0");
  out->getHTTPResponseHeader().addHeader("Expires","Mon, 26 Jul 1997 05:00:00 GMT");

}
