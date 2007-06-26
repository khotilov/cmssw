// -*- C++ -*-
//
// Package:    SiPixelGainCalibration
// Class:      SiPixelIsAliveCalibration
// 
/**\class SiPixelIsAliveCalibration SiPixelIsAliveCalibration.cc CalibTracker/SiPixelGainCalibration/src/SiPixelIsAliveCalibration.cc

 Description: <one line class summary>

 Implementation:
     <Notes on implementation>
*/
//
// Original Author:  Freya Blekman
//         Created:  Thu Jun 14 18:06:29 CEST 2007
// $Id$
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "CalibTracker/SiPixelIsAliveCalibration/interface/SiPixelIsAliveCalibration.h"
#include "DataFormats/SiPixelDetId/interface/PXBDetId.h"
#include "DataFormats/SiPixelDetId/interface/PXFDetId.h"
#include "DataFormats/SiPixelDigi/interface/PixelDigi.h"
#include "DataFormats/DetId/interface/DetId.h"
#include "DataFormats/SiPixelDigi/interface/PixelDigiCollection.h"
#include "DataFormats/Common/interface/DetSetVector.h"
#include "Geometry/TrackerGeometryBuilder/interface/PixelGeomDetUnit.h"
#include "Geometry/TrackerGeometryBuilder/interface/PixelGeomDetType.h"
#include "Geometry/TrackerGeometryBuilder/interface/TrackerGeometry.h"
#include "Geometry/Records/interface/TrackerDigiGeometryRecord.h" 


#include "PhysicsTools/UtilAlgos/interface/TFileDirectory.h"
#include "PhysicsTools/UtilAlgos/interface/TFileService.h"

#include "TH2F.h"
#include "TStyle.h"
#include "TObjArray.h"
#include <string>
#include <iostream>
#include <sstream>
#include <fstream>
//
// class decleration
//

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
SiPixelIsAliveCalibration::SiPixelIsAliveCalibration(const edm::ParameterSet& iConfig):use_calib_(iConfig.getUntrackedParameter<bool>("useCalibFile",false)),
								 inputconfigfile_( iConfig.getUntrackedParameter<std::string>( "inputFileName","/afs/cern.ch/cms/Tracker/Pixel/forward/ryd/calib_070106d.dat" ) ),eventno_counter(0),
  src_( iConfig.getUntrackedParameter<std::string>("src","source"))
{
   //now do what ever initialization is needed
  theHistos_=new TObjArray();
  if(use_calib_)
    calib_ = new PixelCalib(inputconfigfile_);
}


SiPixelIsAliveCalibration::~SiPixelIsAliveCalibration()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called to for each event  ------------
void
SiPixelIsAliveCalibration::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;
   iSetup.get<TrackerDigiGeometryRecord>().get( geom_ );
   eventno_counter++;
   edm::Handle< edm::DetSetVector<PixelDigi> > pixelDigis;
   iEvent.getByLabel( src_, pixelDigis );
   edm::DetSetVector<PixelDigi>::const_iterator digiIter;
   edm::LogInfo("DEBUG") << "SiPixelGainCalibrationAnalysis:found " << pixelDigis->size() << " digi hits..." << std::endl;
   int ntriggers = 1;
   if(use_calib_)
     ntriggers = calib_->nTriggersPerPattern();
   //   std::cout << "ntriggeers = " << ntriggers << std::endl;
   for (digiIter = pixelDigis->begin() ; digiIter != pixelDigis->end(); digiIter++){// loop over detector units...
     unsigned int detid = digiIter->id; 
     edm::DetSet<PixelDigi>::const_iterator ipix;
     for(ipix = digiIter->data.begin(); ipix!=digiIter->end(); ipix++){
       //       unsigned int fed_channel =ipix->channel();
       edm::LogInfo("DEBUG") << "SiPixelIsAliveCalibration: looking at hit with detid " << detid <<", row " << ipix->row() << " column " << ipix->column() << std::endl;
       if(ipix->adc()>0)
	 fill(detid,ipix->adc(), ipix->column(), ipix->row(),ntriggers);
     }
   }
}
// ------------ method to write a .C file to plot the histograms
void SiPixelIsAliveCalibration::writeOutRootMacro(void){
  TString filename =  "readTheHistos.C";// therootfileservice_->GetName();
  std::cout << filename << std::endl;
  //  filename.Replace(".root",".C");
  ofstream outfile;
  outfile.open(filename.Data());
  outfile << "//macro (unnamed) to read histograms. Automatically generated by SiPixelIsAliveCalibration" << std::endl;
  outfile <<"// By Freya Blekman, fblekman@lepp.cornell.edu" << std::endl;
  outfile << "{" << std::endl;
  outfile << "gROOT->Reset();" << std::endl;
  outfile << "TFile *thefile = new TFile(\"histogramsPixelAlive_fromDigi.root\",\"read\"); // check that this is the correct file you're reading in..." << std::endl;
  outfile << "thefile->cd(); " << std::endl;
  for(int i=0; i<theHistos_->GetEntries(); i++){
    TH2F *hist = (TH2F*)theHistos_->At(i);
    if(hist){
    outfile << "TH2F *hist"<< i << "= (TH2F*) thefile->Get(\"SiPixelIsAliveCalibration/" << hist->GetName() << "\"); // modify the file if you've changed module label" << std::endl;
    outfile << "if(hist"<<i<<"){" << std::endl;
    outfile << "\tif(hist"<<i<<"->GetEntries()>0){" << std::endl;
    outfile << "\t \tTCanvas *canv"<< i<<"= new TCanvas(\"canv"<<i<<"\",\"canv"<<i << "\","<< 200*hist->GetNbinsX()/52.<< ","<<200*hist->GetNbinsY()/80.<<");" << std::endl;
    outfile << "\t \t hist"<<i << "->Scale(1/10.);" << std::endl;
    outfile << "\t \thist"<<i << "->Draw(\"colz\");" << std::endl;
    if(i==0)
      outfile << "\t\t canv"<<i<<"->Print(\"pixelalive.pdf(\");"<< std::endl;
    else if(i==theHistos_->GetEntries()-1)
      outfile << "\t\t canv"<<i<<"->Print(\"pixelalive.pdf)\");"<< std::endl;
    else
      outfile << "\t\t canv"<<i<<"->Print(\"pixelalive.pdf\");"<< std::endl;
    outfile << "\t }" << std::endl;
    outfile << "}" << std::endl;
}
  }
  
  outfile << "}" << std::endl;
  outfile.close();
}

// ------------ method called once each job just before starting event loop  ------------
void 
SiPixelIsAliveCalibration::beginJob(const edm::EventSetup&)
{
}
//------------------------------------------
// function that adds data to histograms.
void SiPixelIsAliveCalibration::fill(unsigned int detid, unsigned int adc, unsigned int icol, unsigned int irow, unsigned int ntimes){
  if(!theHistos_->FindObject(makeHistName(detid)))
    init(detid);
  
  TH2F *histo = (TH2F*) theHistos_->FindObject(makeHistName(detid));
  if(histo)
    histo->Fill(icol,irow,1./(float)ntimes);
  
}
TString SiPixelIsAliveCalibration::makeHistName(unsigned int detid){
  DetId detId(detid);
  TString result = "";
  std::ostringstream nameofhist;  
  if(detid==0)
    return result;

  if(detId.subdetId()==1){
    nameofhist << "BPIX_";
    nameofhist << detid;
  }
  else if(detId.subdetId()==2){
    nameofhist << "FPIX_";
    nameofhist << detid;
  }
  else{
    nameofhist << "Unknown PIX_";
    nameofhist<<detid;
  }
  result = nameofhist.str().c_str();
  return result;

}
// function to create new histograms...
void SiPixelIsAliveCalibration::init(unsigned int detid){
  TString nameofhist = makeHistName(detid); 
  if(theHistos_->FindObject(nameofhist))
    return;
  std::ostringstream labelofhist;
  edm::LogInfo("DEBUG") << "looking at det ID : " << detid << std::endl;
  DetId detId(detid);

  if(detid==0)
    return;
  if(detId.subdetId()==1){
    PXBDetId pdetId = PXBDetId(detid);
    labelofhist << "BPIX, ";
    labelofhist << "layer ";
    labelofhist << pdetId.layer();
    labelofhist << ", ladder";
    labelofhist << pdetId.ladder();
    labelofhist << "" ;
						
  }
  else if(detId.subdetId()==2){
    PXFDetId pdetId = PXFDetId(detid);
    labelofhist << "FPIX, ";
    labelofhist<< "disk ";
    labelofhist<<pdetId.disk();
    labelofhist<<", blade ";
    labelofhist<<pdetId.blade();
    labelofhist<<", module ";
    labelofhist<<pdetId.module();
    labelofhist<<", side ";
    labelofhist<<pdetId.side();
    labelofhist<<", panel ";
    labelofhist<< pdetId.panel();

  }
  else if(detId.det()==1){
    labelofhist<<"unknown det id ";
  }
  else return;
  
  std::cout << nameofhist << " " << labelofhist.str() << std::endl;

  const TrackerGeometry& theTracker(*geom_);

  const PixelGeomDetUnit *theGeomDet = dynamic_cast<const PixelGeomDetUnit*> ( theTracker.idToDet(detId) );   
  unsigned int ncols=theGeomDet->specificTopology().ncolumns();
  unsigned int nrows=theGeomDet->specificTopology().nrows();
  //  create histo 
  TH2F *histo = therootfileservice_->make<TH2F>(nameofhist.Data(),labelofhist.str().c_str(),ncols,0,ncols,nrows,0,nrows);
  gStyle->SetOptStat(0);
  histo->SetXTitle("columns");
  histo->SetYTitle("rows");
  histo->SetDrawOption("colz");
  histo->SetMinimum(0.);
  
  theHistos_->Add(histo);
  std::string msg="adding histogram named ";
  msg+=nameofhist;
  msg+=" with bins: ";
  msg+=ncols;
  msg+=",";
  msg+=nrows;
  edm::LogInfo("DEBUG") << msg << std::endl;
}

// ------------ method called once each job just after ending the event loop  ------------
 void 
SiPixelIsAliveCalibration::endJob() {
  //  std::cout << "in endjob..." << std::endl;
  writeOutRootMacro();
}


