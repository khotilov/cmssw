#include "DQM/SiPixelMonitorClient/interface/SiPixelInformationExtractor.h"
#include "DQM/SiPixelMonitorClient/interface/SiPixelUtility.h"
#include "DQMServices/Core/interface/DaqMonitorBEInterface.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "DQMServices/WebComponents/interface/CgiReader.h"
#include "DQM/SiStripCommon/interface/ExtractTObject.h"

#include "TText.h"
#include "TROOT.h"
#include "TPad.h"
#include "TStyle.h"
#include "TSystem.h"
#include "TString.h"
#include "TImage.h"
#include "TPaveText.h"
#include "TImageDump.h"

#include <iostream>
using namespace std;

//
// -- Constructor
// 
SiPixelInformationExtractor::SiPixelInformationExtractor() {
  edm::LogInfo("SiPixelInformationExtractor") << 
    " Creating SiPixelInformationExtractor " << "\n" ;
}
//
// --  Destructor
// 
SiPixelInformationExtractor::~SiPixelInformationExtractor() {
  edm::LogInfo("SiPixelInformationExtractor") << 
    " Deleting SiPixelInformationExtractor " << "\n" ;
  //  if (theCanvas) delete theCanvas;
}
//
// --  Fill Histo and Module List
// 
void SiPixelInformationExtractor::fillModuleAndHistoList(MonitorUserInterface * mui, vector<string>& modules, vector<string>& histos) {
//cout<<"entering SiPixelInformationExtractor::fillModuleAndHistoList"<<endl;
  string currDir = mui->pwd();
  if (currDir.find("Module_") != string::npos)  {
    if (histos.size() == 0) {
      //cout<<"currDir="<<currDir<<endl;
      vector<string> contents = mui->getMEs();    
      for (vector<string>::const_iterator it = contents.begin();
	   it != contents.end(); it++) {
	string hname = (*it).substr(0, (*it).find("_module_"));
	//cout<<"hname="<<hname<<endl;
        histos.push_back(hname);
        string mId=" ";
	if(hname.find("adc")!=string::npos) mId = (*it).substr((*it).find("adc_module_")+11, 9);
        if(mId!=" ") modules.push_back(mId);
        //cout<<"mId="<<mId<<endl;
      }    
    }
  } else {  
    vector<string> subdirs = mui->getSubdirs();
    for (vector<string>::const_iterator it = subdirs.begin();
	 it != subdirs.end(); it++) {
      mui->cd(*it);
      fillModuleAndHistoList(mui, modules, histos);
      mui->goUp();
    }
  }
//  fillBarrelList(mui, modules, histos);
//cout<<"leaving SiPixelInformationExtractor::fillModuleAndHistoList"<<endl;
}
/*void SiPixelInformationExtractor::fillModuleAndHistoList(MonitorUserInterface * mui, vector<string>& modules, vector<string>& histos) {
cout<<"entering SiPixelInformationExtractor::fillModuleAndHistoList"<<endl;
  mui->cd();
  fillBarrelList(mui, modules, histos); 
  mui->cd();

cout<<"leaving SiPixelInformationExtractor::fillModuleAndHistoList"<<endl;
}*/
void SiPixelInformationExtractor::createModuleTree(MonitorUserInterface* mui) {
//cout<<"entering SiPixelInformationExtractor::createModuleTree..."<<endl;
  string structure_name;
  vector<string> me_names;
  if (!configParser_->getMENamesForTree(structure_name, me_names)){
    cout << "SiPixelInformationExtractor::createModuleTree: Failed to read Tree configuration parameters!! ";
    return;
  }
  mui->cd();
  fillBarrelList(mui, structure_name, me_names);
  mui->cd();
  fillEndcapList(mui, structure_name, me_names);
  mui->cd();
  actionExecutor_->createLayout(mui);
  string fname = "test1.xml";
  configWriter_->write(fname);
  if (configWriter_) delete configWriter_;
  configWriter_ = 0;
//cout<<"leaving SiPixelInformationExtractor::createModuleTree..."<<endl;
}
void SiPixelInformationExtractor::fillBarrelList(MonitorUserInterface* mui,
                               string dir_name,vector<string>& me_names) {
  //cout<<"entering SiPixelInformationExtractor::fillBarrelList..."<<endl;
  string currDir = mui->pwd();
  if (currDir.find(dir_name) != string::npos)  {
    vector<MonitorElement*> mod_mes;
    vector<string> contents = mui->getMEs(); 
    for (vector<string>::const_iterator iv = me_names.begin();
	 iv != me_names.end(); iv++) {
      for (vector<string>::const_iterator im = contents.begin();
	   im != contents.end(); im++) {
        string sname = (*iv);
        string tname = sname.substr(8,(sname.find("_",8)-8));
	if (((*im)).find(tname) == 0) {
	  string fullpathname = mui->pwd() + "/" + (*im); 
          MonitorElement* me = getModuleME(mui, fullpathname);
	}
      }
    }
  } else {  
    vector<string> subdirs = mui->getSubdirs();
    for (vector<string>::const_iterator it = subdirs.begin();
       it != subdirs.end(); it++) {
      if((*it).find("PixelEndcap")!=string::npos) continue;
      mui->cd(*it);
      fillBarrelList(mui, dir_name, me_names);
      mui->goUp();
    }
  }
  //cout<<"...leaving SiPixelActionExecutor::fillBarrelSummary!"<<endl;
}
void SiPixelInformationExtractor::fillEndcapList(MonitorUserInterface* mui,
                               string dir_name,vector<string>& me_names) {
  //cout<<"entering SiPixelInformationExtractor::fillEndcapList..."<<endl;
  string currDir = mui->pwd();
  if (currDir.find(dir_name) != string::npos)  {
    vector<MonitorElement*> mod_mes;
    vector<string> contents = mui->getMEs(); 
    for (vector<string>::const_iterator iv = me_names.begin();
	 iv != me_names.end(); iv++) {
      for (vector<string>::const_iterator im = contents.begin();
	   im != contents.end(); im++) {
        string sname = (*iv);
        string tname = sname.substr(8,(sname.find("_",8)-8));
	if (((*im)).find(tname) == 0) {
	  string fullpathname = mui->pwd() + "/" + (*im); 
          MonitorElement* me = getModuleME(mui, fullpathname);
	}
      }
    }
  } else {  
    vector<string> subdirs = mui->getSubdirs();
    for (vector<string>::const_iterator it = subdirs.begin();
       it != subdirs.end(); it++) {
      if((mui->pwd()).find("PixelBarrel")!=string::npos) mui->goUp();
      mui->cd((*it));
      if((*it).find("PixelBarrel")!=string::npos) continue;
      fillEndcapList(mui, dir_name, me_names);
      mui->goUp();
    }
  }
  //cout<<"...leaving SiPixelActionExecutor::fillBarrelSummary!"<<endl;
}
//
// -- Get Summary ME
//
MonitorElement* SiPixelInformationExtractor::getModuleME(MonitorUserInterface* mui,string me_name) {
//cout<<"Entering SiPixelInformationExtractor::getModuleME..."<<endl;
  MonitorElement* me = 0;
  // If already booked
  vector<string> contents = mui->getMEs();    
  for (vector<string>::const_iterator it = contents.begin();
       it != contents.end(); it++) {
    if ((*it).find(me_name) == 0) {
      string fullpathname = mui->pwd() + "/" + (*it); 
      me = mui->get(fullpathname);
      if (me) {
	MonitorElementT<TNamed>* obh1 = dynamic_cast<MonitorElementT<TNamed>*> (me);
	if (obh1) {
	  TH1F * root_obh1 = dynamic_cast<TH1F *> (obh1->operator->());
	  if (root_obh1) root_obh1->Reset();        
	}
	return me;
      }
    }
  }
  //cout<<"...leaving SiPixelInformationExtractor::getModuleME!"<<endl;
}

/*void SiPixelInformationExtractor::fillBarrelList(MonitorUserInterface* mui,
                               vector<string>& modules, vector<string>& histos) {
  cout<<"entering SiPixelInformationExtractor::fillBarrelList..."<<endl;
  string currDir = mui->pwd();
  if (currDir.find("Ladder_") != string::npos)  {
    vector<string> subdirs = mui->getSubdirs();
    int ndet = 0;
    for (vector<string>::const_iterator it = subdirs.begin();
       it != subdirs.end(); it++) {
      if ( (*it).find("Module_") == string::npos) continue;
      mui->cd(*it);
      ndet++;
      vector<string> contents = mui->getMEs(); 
      for (vector<string>::const_iterator im = contents.begin();
	   im != contents.end(); im++) {
	string hname = (*im).substr(0, (*im).find("_module_"));
//	cout<<"hname="<<hname<<endl;
//        histos.push_back(hname);
	string fullpathname = mui->pwd() + "/" + (*im); 
	MonitorElement *  me = mui->get(fullpathname);
//        string mId=" ";
//	if(hname.find("adc")!=string::npos) mId = fullpathname;
//	if(hname.find("col")!=string::npos) mId = fullpathname;
//	if(hname.find("row")!=string::npos) mId = fullpathname;
//	if(hname.find("ndigis")!=string::npos) mId = fullpathname;
//        if(mId!=" ") modules.push_back(mId);
//	cout<<"mId="<<mId<<endl;
      }      
      mui->goUp();
    }
  } else {  
    vector<string> subdirs = mui->getSubdirs();
    for (vector<string>::const_iterator it = subdirs.begin();
       it != subdirs.end(); it++) {
      if((*it).find("PixelEndcap")!=string::npos) continue;
      mui->cd(*it);
      fillBarrelList(mui, modules, histos);
      mui->goUp();
    }
    //fillGrandBarrelSummaryHistos(mui, me_names);
  }
  cout<<"...leaving SiPixelInformationExtractor::fillBarrelList!"<<endl;
}
*/




//
// --  Fill Module Histo List
// 
void SiPixelInformationExtractor::printModuleHistoList(MonitorUserInterface * mui, ostringstream& str_val){
//cout<<"entering SiPixelInformationExtractor::printModuleHistoList"<<endl;
  static string indent_str = "";

  string currDir = mui->pwd();
  string dname = currDir.substr(currDir.find_last_of("/")+1);
  //if (dname.find("_module_") ==0) return;
  str_val << "<li><a href=\"#\" id=\"" 
          << currDir << "\">" << dname << "</a>" << endl;
  vector<string> meVec = mui->getMEs(); 
  vector<string> subDirVec = mui->getSubdirs();
  if ( meVec.size()== 0  && subDirVec.size() == 0 ) {
    str_val << "</li> "<< endl;    
    return;
  }
  str_val << "<ul>" << endl; 
  for (vector<string>::const_iterator it = meVec.begin();
       it != meVec.end(); it++) {
    if ((*it).find("_module_")!=string::npos) {
      str_val << "<li class=\"dhtmlgoodies_sheet.gif\"><a href=\"javascript:DrawSingleHisto('"
           << currDir << "/"<< (*it) << "')\">" << (*it) << "</a></li>" << endl;
    }
  }

  for (vector<string>::const_iterator ic = subDirVec.begin();
       ic != subDirVec.end(); ic++) {
    mui->cd(*ic);
    printModuleHistoList(mui, str_val);
    mui->goUp();
  }
  str_val << "</ul> "<< endl;  
  str_val << "</li> "<< endl;  
//cout<<"leaving SiPixelInformationExtractor::printModuleHistoList"<<endl;
}
//
// --  Fill Summary Histo List
// 
void SiPixelInformationExtractor::printSummaryHistoList(MonitorUserInterface * mui, ostringstream& str_val){
//cout<<"entering SiPixelInformationExtractor::printSummaryHistoList"<<endl;
  static string indent_str = "";

  string currDir = mui->pwd();
  string dname = currDir.substr(currDir.find_last_of("/")+1);
  if (dname.find("Module_") ==0) return;
  str_val << "<li><a href=\"#\" id=\"" 
          << currDir << "\">" << dname << "</a>" << endl;
  vector<string> meVec = mui->getMEs(); 
  vector<string> subDirVec = mui->getSubdirs();
  if ( meVec.size()== 0  && subDirVec.size() == 0 ) {
    str_val << "</li> "<< endl;    
    return;
  }
  str_val << "<ul>" << endl;      
  for (vector<string>::const_iterator it = meVec.begin();
       it != meVec.end(); it++) {
    if ((*it).find("Summary") == 0) {
      str_val << "<li class=\"dhtmlgoodies_sheet.gif\"><a href=\"javascript:DrawSingleHisto('"
           << currDir << "/"<< (*it) << "')\">" << (*it) << "</a></li>" << endl;
    }
  }

  for (vector<string>::const_iterator ic = subDirVec.begin();
       ic != subDirVec.end(); ic++) {
    mui->cd(*ic);
    printSummaryHistoList(mui, str_val);
    mui->goUp();
  }
  str_val << "</ul> "<< endl;  
  str_val << "</li> "<< endl;  
//cout<<"leaving SiPixelInformationExtractor::printSummaryHistoList"<<endl;
}
//
// --  Fill Alarm List
// 
void SiPixelInformationExtractor::printAlarmList(MonitorUserInterface * mui, ostringstream& str_val){
//cout<<"entering SiPixelInformationExtractor::printAlarmList"<<endl;
  static string indent_str = "";

  string currDir = mui->pwd();
  string dname = currDir.substr(currDir.find_last_of("/")+1);
  string image_name;
  selectImage(image_name,mui->getStatus(currDir));
  str_val << "<li><a href=\"#\" id=\"" 
          << currDir << "\">" << dname << "</a> <img src=\"" 
          << image_name << "\">" << endl;
  vector<string> subDirVec = mui->getSubdirs();
  vector<string> meVec = mui->getMEs(); 
  if (subDirVec.size() == 0 && meVec.size() == 0) {
    str_val << "</li> "<< endl;    
    return;
  }
  str_val << "<ul>" << endl;
  if (dname.find("Module_") != string::npos) {
    if (meVec.size() > 0) {
      for (vector<string>::const_iterator it = meVec.begin();
	   it != meVec.end(); it++) {
        string full_path = currDir + "/" + (*it);
	MonitorElement * me = mui->get(full_path);
	if (!me) continue;
        dqm::qtests::QR_map my_map = me->getQReports();
        if (my_map.size() > 0) {
	  string image_name1;
	  selectImage(image_name1,my_map);
	  str_val << "<li class=\"dhtmlgoodies_sheet.gif\"><a href=\"javascript:ReadStatus('"
		<< full_path<< "')\">" << (*it) << "</a><img src=\""
		<< image_name1 << "\""<< "</li>" << endl;
        }
      }
    }
  }
  for (vector<string>::const_iterator ic = subDirVec.begin();
       ic != subDirVec.end(); ic++) {
    mui->cd(*ic);
    printAlarmList(mui, str_val);
    mui->goUp();
  }
  str_val << "</ul> "<< endl;  
  str_val << "</li> "<< endl;  
//cout<<"leaving SiPixelInformationExtractor::printAlarmList"<<endl;
}
//
// --  Get Selected Monitor Elements
// 
void SiPixelInformationExtractor::selectSingleModuleHistos(MonitorUserInterface * mui, string mid, vector<string>& names, vector<MonitorElement*>& mes) {
//cout<<"entering SiPixelInformationExtractor::selectSingleModuleHistos"<<endl;
  string currDir = mui->pwd();
  //cout<<"currDir="<<currDir<<endl;
//  if (currDir.find("Module_") != string::npos &&
//      currDir.find(mid) != string::npos )  {
  if (currDir.find("Module_") != string::npos)  {
    vector<string> contents = mui->getMEs();    
    for (vector<string>::const_iterator it = contents.begin();
	 it != contents.end(); it++) {
      if((*it).find(mid) != string::npos){
        for (vector<string>::const_iterator ih = names.begin();
	   ih != names.end(); ih++) {
	  string temp_s = (*it).substr(0, (*it).find("_module_"));
	  //	if ((*it).find((*ih)) != string::npos) {
	  if (temp_s == (*ih)) {
	    string full_path = currDir + "/" + (*it);
	    //cout<<"full_path="<<full_path<<endl;
	    MonitorElement * me = mui->get(full_path.c_str());
	    if (me) mes.push_back(me);
	  }  
        }
      }
    }
    if (mes.size() >0) return;
  } else {  
    vector<string> subdirs = mui->getSubdirs();
    for (vector<string>::const_iterator it = subdirs.begin();
	 it != subdirs.end(); it++) {
      mui->cd(*it);
      selectSingleModuleHistos(mui, mid, names, mes);
      mui->goUp();
    }
  }
//cout<<"leaving SiPixelInformationExtractor::selectSingleModuleHistos"<<endl;
}
//
// --  Plot Selected Monitor Elements
// 
void SiPixelInformationExtractor::plotSingleModuleHistos(MonitorUserInterface* mui, multimap<string, string>& req_map) {
//cout<<"entering SiPixelInformationExtractor::plotSingleModuleHistos"<<endl;
  vector<string> item_list;  

  string mod_id = getItemValue(req_map,"ModId");
  //cout<<"mod_id in plotSingleModuleHistos:"<<mod_id<<endl;
  if (mod_id.size() < 9) return;
  item_list.clear();     
  getItemList(req_map,"histo", item_list); // item_list holds all histos to plot
  vector<MonitorElement*> me_list;

  mui->cd();
  selectSingleModuleHistos(mui, mod_id, item_list, me_list);
  mui->cd();

  plotHistos(req_map,me_list);
//cout<<"leaving SiPixelInformationExtractor::plotSingleModuleHistos"<<endl;
}
//
// -- plot a Histogram
//
void SiPixelInformationExtractor::plotSingleHistogram(MonitorUserInterface * mui,
		       std::multimap<std::string, std::string>& req_map){
//cout<<"entering SiPixelInformationExtractor::plotSingleHistogram"<<endl;
  vector<string> item_list;  

  string path_name = getItemValue(req_map,"Path");
  if (path_name.size() == 0) return;
  
  MonitorElement* me = mui->get(path_name);
  vector<MonitorElement*> me_list;
  if (me) {
    me_list.push_back(me);
    plotHistos(req_map,me_list);
  }
//cout<<"leaving SiPixelInformationExtractor::plotSingleHistogram"<<endl;
}
//
//  plot Histograms in a Canvas
//
void SiPixelInformationExtractor::plotHistos(multimap<string,string>& req_map, 
 
  			   vector<MonitorElement*> me_list){
//cout<<"entering SiPixelInformationExtractor::plotHistos"<<endl;
  int nhist = me_list.size();
  if (nhist == 0) return;
  int width = 600;
  int height = 600;
  TCanvas canvas("TestCanvas", "Test Canvas");
  canvas.Clear();
  gROOT->Reset(); gStyle->SetPalette(1);
  int ncol, nrow;
 
  float xlow = -1.0;
  float xhigh = -1.0;
  
  if (nhist == 1) {
    if (hasItem(req_map,"xmin")) xlow = atof(getItemValue(req_map,"xmin").c_str());
    if (hasItem(req_map,"xmax")) xhigh = atof(getItemValue(req_map,"xmax").c_str()); 
    ncol = 1;
    nrow = 1;
  } else {
    ncol = atoi(getItemValue(req_map, "cols").c_str());
    nrow = atoi(getItemValue(req_map, "rows").c_str());
    if (ncol*nrow < nhist) {
      if (nhist == 2) {
	ncol = 1;
	nrow = 2;
      } else if (nhist == 3) {
	ncol = 1;
	nrow = 3;
      } else if (nhist == 4) {
	ncol = 2;
	nrow = 3;
      } else if (nhist > 4 && nhist <= 10) {
        ncol = 2;
	nrow = nhist/ncol+1;
      } else if (nhist > 10 && nhist <= 20) {
        ncol = 3;
	nrow = nhist/ncol+1;
      } else if (nhist > 20 && nhist <= 40) {
         ncol = 4;
	 nrow = nhist/ncol+1;
      } 		

    }
  }

  if (hasItem(req_map,"width")) 
              width = atoi(getItemValue(req_map, "width").c_str());    
  if (hasItem(req_map,"height"))
              height = atoi(getItemValue(req_map, "height").c_str());

  canvas.SetWindowSize(width,height);
  canvas.Divide(ncol, nrow);
  int i=0;
  for (vector<MonitorElement*>::const_iterator it = me_list.begin();
       it != me_list.end(); it++) {
    i++;
    int istat =  SiPixelUtility::getStatus((*it));
    string tag;
    int icol;
    SiPixelUtility::getStatusColor(istat, icol, tag);
  
    MonitorElementT<TNamed>* ob = 
      dynamic_cast<MonitorElementT<TNamed>*>((*it));
    if (ob) {
//    cout<<"Has a ME to plot!"<<endl;
      canvas.cd(i);
      //      TAxis* xa = ob->operator->()->GetXaxis();
      //      xa->SetRangeUser(xlow, xhigh);
      if(hasItem(req_map,"colpal")){
        gROOT->Reset(); gStyle->SetPalette(1); gStyle->SetOptStat(0);
        ob->operator->()->Draw("colz");
      }else{
        gStyle->SetOptStat(1);
        ob->operator->()->Draw();
      }
      if (icol != 1) {
	TText tt;
	tt.SetTextSize(0.12);
	tt.SetTextColor(icol);
	tt.DrawTextNDC(0.5, 0.5, tag.c_str());
      }
      if (hasItem(req_map,"logy")) {
	  gPad->SetLogy(1);
      }
    }
  }
  gStyle->SetPalette(1);
  canvas.Update();
  fillImageBuffer(canvas);
  canvas.Clear();
//cout<<"leaving SiPixelInformationExtractor::plotHistos"<<endl;
}
//
// read the Module And HistoList
//
void SiPixelInformationExtractor::readModuleAndHistoList(MonitorUserInterface* mui, xgi::Output * out, bool coll_flag) {
//cout<<"entering SiPixelInformationExtractor::readModuleAndHistoList"<<endl;
   std::vector<std::string> hnames;
   std::vector<std::string> mod_names;
   if (coll_flag)  mui->cd("Collector/Collated");
   fillModuleAndHistoList(mui, mod_names, hnames);
   //for (std::vector<std::string>::iterator im = mod_names.begin();
   //     im != mod_names.end(); im++) cout<<"mod_names="<<*im<<endl;
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

   for (std::vector<std::string>::iterator ih = hnames.begin();
        ih != hnames.end(); ih++) {
     *out << "<Histo>" << *ih << "</Histo>" << endl;     

   }
   *out << "</HistoList>" << endl;
   *out << "</ModuleAndHistoList>" << endl;
   if (coll_flag)  mui->cd();
//cout<<"leaving SiPixelInformationExtractor::readModuleAndHistoList"<<endl;
}
//
// read the Structure And ModuleHisto Tree
//
void SiPixelInformationExtractor::readModuleHistoTree(MonitorUserInterface* mui, string& str_name, xgi::Output * out, bool coll_flag) {
//cout<<"entering  SiPixelInformationExtractor::readModuleHistoTree"<<endl;
  ostringstream modtree;
  if (goToDir(mui, str_name, coll_flag)) {
    modtree << "<ul id=\"dhtmlgoodies_tree\" class=\"dhtmlgoodies_tree\">" << endl;
    printModuleHistoList(mui,modtree);
    modtree <<"</ul>" << endl;   
  } else {
    modtree << "Desired Directory does not exist";
  }
  out->getHTTPResponseHeader().addHeader("Content-Type", "text/plain");
  *out << modtree.str();
   mui->cd();
//cout<<"leaving  SiPixelInformationExtractor::readModuleHistoTree"<<endl;
}
//
// read the Structure And SummaryHistogram List
//
void SiPixelInformationExtractor::readSummaryHistoTree(MonitorUserInterface* mui, string& str_name, xgi::Output * out, bool coll_flag) {
//cout<<"entering  SiPixelInformationExtractor::readSummaryHistoTree"<<endl;
  ostringstream sumtree;
  if (goToDir(mui, str_name, coll_flag)) {
    sumtree << "<ul id=\"dhtmlgoodies_tree\" class=\"dhtmlgoodies_tree\">" << endl;
    printSummaryHistoList(mui,sumtree);
    sumtree <<"</ul>" << endl;   
  } else {
    sumtree << "Desired Directory does not exist";
  }
  out->getHTTPResponseHeader().addHeader("Content-Type", "text/plain");
  *out << sumtree.str();
   mui->cd();
//cout<<"leaving  SiPixelInformationExtractor::readSummaryHistoTree"<<endl;
}
//
// read the Structure And Alarm Tree
//
void SiPixelInformationExtractor::readAlarmTree(MonitorUserInterface* mui, 
                  string& str_name, xgi::Output * out, bool coll_flag){
//cout<<"entering SiPixelInformationExtractor::readAlarmTree"<<endl;
  ostringstream alarmtree;
  if (goToDir(mui, str_name, coll_flag)) {
    alarmtree << "<ul id=\"dhtmlgoodies_tree\" class=\"dhtmlgoodies_tree\">" << endl;
    printAlarmList(mui,alarmtree);
    alarmtree <<"</ul>" << endl; 
  } else {
    alarmtree << "Desired Directory does not exist";
  }
  out->getHTTPResponseHeader().addHeader("Content-Type", "text/plain");
  *out << alarmtree.str();
   mui->cd();
//cout<<"leaving SiPixelInformationExtractor::readAlarmTree"<<endl;
}
//
// Get elements from multi map
//
void SiPixelInformationExtractor::getItemList(multimap<string, string>& req_map, 
                      string item_name,vector<string>& items) {
//cout<<"entering SiPixelInformationExtractor::getItemList"<<endl;
  items.clear();
  for (multimap<string, string>::const_iterator it = req_map.begin();
       it != req_map.end(); it++) {
    //cout<<"....item_name="<<item_name<<endl;
    //cout<<"....first="<<it->first<<" ....second="<<it->second<<endl;
    if (it->first == item_name) {
      items.push_back(it->second);
    }
  }
//cout<<"leaving SiPixelInformationExtractor::getItemList"<<endl;
}
//
//  check a specific item in the map
//
bool SiPixelInformationExtractor::hasItem(multimap<string,string>& req_map,
					  string item_name){
//cout<<"entering SiPixelInformationExtractor::hasItem"<<endl;
  multimap<string,string>::iterator pos = req_map.find(item_name);
  if (pos != req_map.end()) return true;
  return false;  
//cout<<"leaving SiPixelInformationExtractor::hasItem"<<endl;
}
//
// check the value of an item in the map
//  
string SiPixelInformationExtractor::getItemValue(multimap<string,string>& req_map,
						 std::string item_name){
//cout<<"entering SiPixelInformationExtractor::getItemValue"<<endl;
  multimap<string,string>::iterator pos = req_map.find(item_name);
  string value = " ";
  if (pos != req_map.end()) {
    value = pos->second;
  }
  return value;
//cout<<"leaving SiPixelInformationExtractor::getItemValue"<<endl;
}
//
// write the canvas in a string
//
void SiPixelInformationExtractor::fillImageBuffer(TCanvas& c1) {
//cout<<"entering SiPixelInformationExtractor::fillImageBuffer"<<endl;
  c1.SetFixedAspectRatio(kTRUE);
  c1.SetCanvasSize(530, 440);
  gStyle->SetPalette(1);
  // Now extract the image
  // 114 - stands for "no write on Close"
  TImageDump imgdump("tmp.png", 114);
  c1.Paint();

 // get an internal image which will be automatically deleted
 // in the imgdump destructor
  TImage *image = imgdump.GetImage();

  char *buf;
  int sz;
  image->GetImageBuffer(&buf, &sz);         /* raw buffer */
  pictureBuffer_.str("");
  for (int i = 0; i < sz; i++)
    pictureBuffer_ << buf[i];
  
  delete [] buf;
//cout<<"leaving SiPixelInformationExtractor::fillImageBuffer"<<endl;
}
//
// get the plot
//
const ostringstream&  SiPixelInformationExtractor::getImage() const {
//cout<<"entering SiPixelInformationExtractor::getImage"<<endl;
  return pictureBuffer_;
//cout<<"leaving SiPixelInformationExtractor::getImage"<<endl;
}
//
// go to a specific directory after scanning
//
bool SiPixelInformationExtractor::goToDir(MonitorUserInterface* mui, string& sname, bool flg){ 
//cout<<"entering SiPixelInformationExtractor::goToDir"<<endl;
  mui->cd();
  mui->cd("Collector");
//  cout << mui->pwd() << endl;
  vector<string> subdirs;
  subdirs = mui->getSubdirs();
  if (subdirs.size() == 0) return false;
  
  if (flg) mui->cd("Collated");
  else mui->cd(subdirs[0]);
//  cout << mui->pwd() << endl;
  subdirs.clear();
  subdirs = mui->getSubdirs();
  if (subdirs.size() == 0) return false;
  mui->cd(sname);
  string dirName = mui->pwd();
  if (dirName.find(sname) != string::npos) return true;
  else return false;  
//cout<<"leaving SiPixelInformationExtractor::goToDir"<<endl;
}
//
// -- Get Image name from status
//
void SiPixelInformationExtractor::selectImage(string& name, int status){
//cout<<"entering SiPixelInformationExtractor::selectImage"<<endl;
  if (status == dqm::qstatus::STATUS_OK) name="images/LI_green.gif";
  else if (status == dqm::qstatus::WARNING) name="images/LI_yellow.gif";
  else if (status == dqm::qstatus::ERROR) name="images/LI_red.gif";
  else if (status == dqm::qstatus::OTHER) name="images/LI_orange.gif";
  else  name="images/LI_blue.gif";
//cout<<"leaving SiPixelInformationExtractor::selectImage"<<endl;
}
//
// -- Get Image name from ME
//
void SiPixelInformationExtractor::selectImage(string& name, dqm::qtests::QR_map& test_map){
//cout<<"entering SiPixelInformationExtractor::selectImage"<<endl;
  int istat = 999;
  int status = 0;
  for (dqm::qtests::QR_map::const_iterator it = test_map.begin(); it != test_map.end();
       it++) {
    status = it->second->getStatus();
    if (status > istat) istat = status;
  }
  selectImage(name, status);
//cout<<"leaving SiPixelInformationExtractor::selectImage"<<endl;
}
//
// -- Get Warning/Error Messages
//
void SiPixelInformationExtractor::readStatusMessage(MonitorUserInterface* mui, string& path,xgi::Output * out) {
//cout<<"entering SiPixelInformationExtractor::readStatusMessage"<<endl;
  MonitorElement* me = mui->get(path);
  ostringstream test_status;
  if (!me) {
    test_status << " ME Does not exist ! ";
  } else {
    dqm::qtests::QR_map test_map = me->getQReports();
    for (dqm::qtests::QR_map::const_iterator it = test_map.begin(); it != test_map.end();
	 it++) {
      int status = it->second->getStatus();
      if (status == dqm::qstatus::WARNING) test_status << " Warning : ";
      else if (status == dqm::qstatus::ERROR) test_status << " Error : ";
      else if (status == dqm::qstatus::STATUS_OK) test_status << " Ok : ";
      else if (status == dqm::qstatus::OTHER) test_status << " Other(" << status << ") : ";
      test_status << it->second->getMessage();
    }      
  }
  out->getHTTPResponseHeader().addHeader("Content-Type", "text/plain");
  *out << test_status.str();
//cout<<"leaving SiPixelInformationExtractor::readStatusMessage"<<endl;
}

