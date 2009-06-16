#include "VisMonitoring/DQMServer/interface/DQMRenderPlugin.h"
#include "utils.h"

#include "TH1F.h"
#include "TH2F.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TGaxis.h"
#include "TColor.h"
#include "TROOT.h"
#include "TGraph.h"
#include "TLine.h"

#include <iostream>
#include <math.h>


class ESRenderPlugin : public DQMRenderPlugin {

 public:

  // virtual void initialise( int argc, char **argv );

  virtual bool applies( const DQMNet::CoreObject &o, const VisDQMImgInfo &i );

  virtual void preDraw( TCanvas *c, const DQMNet::CoreObject &o, const VisDQMImgInfo &i, VisDQMRenderInfo&  r);

  virtual void postDraw( TCanvas *c, const DQMNet::CoreObject &o, const VisDQMImgInfo &i );

 private:

  void preDrawTH2F( TCanvas *c, const DQMNet::CoreObject &o );
//  void preDrawTH1F( TCanvas *c, const DQMNet::CoreObject &o );

  void postDrawTH2F( TCanvas *c, const DQMNet::CoreObject &o );
//  void postDrawTH1F( TCanvas *c, const DQMNet::CoreObject &o );

  double NEntries;

};


bool ESRenderPlugin::applies( const DQMNet::CoreObject &o, const VisDQMImgInfo &i ) {

#ifdef DEBUG
	std::cout << "ESRenderPlugin:applies " << o.name << std::endl;
#endif

	if( o.name.find( "EcalPreshower" ) != std::string::npos ) {
		if( o.name.find( "ESOccupancyTask" ) != std::string::npos ){
			return true;
		}
		if( o.name.find( "ESIntegrityClient" ) != std::string::npos ){
			return true;
		}
	}

	return false;

}


void ESRenderPlugin::preDraw( TCanvas *c, const DQMNet::CoreObject &o, const VisDQMImgInfo &i, VisDQMRenderInfo &r ) {

#ifdef DEBUG
	std::cout << "ESRenderPlugin:preDraw " << o.name << std::endl;
#endif

	c->cd();

	gStyle->Reset("Default");

	gStyle->SetCanvasColor(10);
	gStyle->SetPadColor(10);
	gStyle->SetFillColor(10);
	gStyle->SetFrameFillColor(10);
	gStyle->SetStatColor(10);
	gStyle->SetTitleFillColor(10);

	TGaxis::SetMaxDigits(4);

	gStyle->SetOptTitle(kTRUE);
	gStyle->SetTitleBorderSize(0);

	gStyle->SetOptStat(kFALSE);
	gStyle->SetStatBorderSize(1);

	gStyle->SetOptFit(kFALSE);

	gROOT->ForceStyle();

	if( dynamic_cast<TH2F*>( o.object ) ) {
		preDrawTH2F( c, o );
	}

	r.drawOptions = "";

#ifdef DEBUG
	std::cout << "done" << std::endl;
#endif

}

void ESRenderPlugin::postDraw( TCanvas *c, const DQMNet::CoreObject &o, const VisDQMImgInfo &i ) {

#ifdef DEBUG
	std::cout << "ESRenderPlugin:postDraw " << o.name << std::endl;
#endif

	c->cd();

	//  if( dynamic_cast<TH2F*>( o.object ) ) {
	//    postDrawTH2F( c, o );
	//  }

#ifdef DEBUG
	std::cout << "done" << std::endl;
#endif

}



void ESRenderPlugin::preDrawTH2F( TCanvas *c, const DQMNet::CoreObject &o ) {

	TH2F* obj = dynamic_cast<TH2F*>( o.object );

	assert( obj );

	std::string name = o.name.substr(o.name.rfind("/")+1);

	//  int nbx = obj->GetNbinsX();
	//  int nby = obj->GetNbinsY();
	int colorbar[6] = {1,2,3,4,5,6};

	gStyle->SetPaintTextFormat();

	gStyle->SetOptStat(kFALSE);
	obj->SetStats(kFALSE);
	gPad->SetLogy(kFALSE);

	if( name.find( "Integrity Summary" ) != std::string::npos ) {
		gStyle->SetPalette(6,colorbar);
		obj->SetMinimum(0.5);
		obj->SetMaximum(6.5);
		obj->SetOption("colz");
		gPad->SetRightMargin(0.15);
		gStyle->SetPaintTextFormat("+g");
		return;
	}


	if( name.find( "RecHit 2D Occupancy" ) != std::string::npos ) {
		gStyle->SetPalette(1);
      		obj->GetXaxis()->SetRange(0,39);
      		obj->GetYaxis()->SetRange(0,39);
		NEntries = obj->GetBinContent(40,40);
		obj->Scale(1/NEntries);
		obj->SetMaximum(33);
		obj->SetOption("colz");
		gPad->SetRightMargin(0.15);
		gStyle->SetPaintTextFormat("+g");
		return;
	}

	if( name.find( "DigiHit 2D Occupancy" ) != std::string::npos ) {
		gStyle->SetPalette(1);
      		obj->GetXaxis()->SetRange(0,39);
      		obj->GetYaxis()->SetRange(0,39);
		NEntries = obj->GetBinContent(40,40);
		obj->Scale(1/NEntries);
//		obj->SetMaximum(33);
		obj->SetOption("colz");
		gPad->SetRightMargin(0.15);
		gStyle->SetPaintTextFormat("+g");
		return;
	}



}

void ESRenderPlugin::postDrawTH2F( TCanvas *c, const DQMNet::CoreObject &o ) {

	TH2F* obj = dynamic_cast<TH2F*>( o.object );

	assert( obj );

	std::string name = o.name.substr(o.name.rfind("/")+1);

	if( name.find( "2D Occupancy" ) != std::string::npos ) return;

}


static ESRenderPlugin instance;
