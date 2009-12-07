/*!
  \File BeamRenderPlugin
  \Display Plugin for BeamSpot DQM Histograms
  \author 
  \version $Revision: 1.2 $
  \date $Date: 2009/11/19 04:43:33 $
*/

#include "VisMonitoring/DQMServer/interface/DQMRenderPlugin.h"
#include "utils.h"

#include "TProfile2D.h"
#include "TProfile.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TColor.h"
#include "TText.h"
#include <cassert>

class BeamRenderPlugin : public DQMRenderPlugin {

public:
  virtual bool applies( const VisDQMObject &o, const VisDQMImgInfo & ) {
    if ((o.name.find( "BeamMonitor/" ) == std::string::npos) &&
	(o.name.find( "BeamMonitor_PixelLess/" ) == std::string::npos))
      return false;
      
    if (o.name.find( "/EventInfo/" ) != std::string::npos)
      return true;
      
    if (o.name.find( "/Fit/" ) != std::string::npos)
      return true;

    return false;

  }

  virtual void preDraw( TCanvas *c, const VisDQMObject &o, const VisDQMImgInfo &, VisDQMRenderInfo & ) {
    c->cd();

    if ( dynamic_cast<TH2F*>( o.object ) ) {
      preDrawTH2F( c, o );
    }

    if ( dynamic_cast<TH1F*>( o.object ) ) {
      preDrawTH1F( c, o );
    }

  }

  virtual void postDraw( TCanvas *c, const VisDQMObject &o, const VisDQMImgInfo & ) {
    c->cd();
    
    if ( dynamic_cast<TH2F*>( o.object ) ) {
      postDrawTH2F( c, o );
    }

  }
  
private:
  
  void preDrawTH2F( TCanvas *c, const VisDQMObject &o ) {
    
    TH2F* obj = dynamic_cast<TH2F*>( o.object );
    assert( obj );

    // This applies to all
    gStyle->SetCanvasBorderMode( 0 );
    gStyle->SetPadBorderMode( 0 );
    gStyle->SetPadBorderSize( 0 );
    
    TAxis* xa = obj->GetXaxis();
    TAxis* ya = obj->GetYaxis();

    xa->SetTitleOffset(0.7);
    xa->SetTitleSize(0.05);
    xa->SetLabelSize(0.04);

    ya->SetTitleOffset(0.7);
    ya->SetTitleSize(0.05);
    ya->SetLabelSize(0.04);

    if ( o.name.find( "reportSummaryMap" )  != std::string::npos) {
      obj->SetStats( kFALSE );
      dqm::utils::reportSummaryMapPalette(obj);
      obj->SetOption("colz");
      return;
    }

    if ( o.name.find( "trk_vx_vy" )  != std::string::npos) {
      gStyle->SetOptStat(11);
      obj->SetOption("colz");
      return;
    }

    if ( o.name.find( "fitResults" )  != std::string::npos) {
      c->SetGrid();
      obj->SetStats( kFALSE );
      return;
    }

  }

  void preDrawTH1F( TCanvas *, const VisDQMObject &o ) {
    TH1F* obj = dynamic_cast<TH1F*>( o.object );
    assert( obj );

    TAxis* ya = obj->GetYaxis();
    ya->SetTitleOffset(1.1);
    ya->SetTitleSize(0.04);
    ya->SetLabelSize(0.03);

    if( o.name.find( "nTrk_lumi" )  != std::string::npos) {
      gStyle->SetOptStat(11);
      return;
    }
  }

  void postDrawTH2F( TCanvas *c, const VisDQMObject &o ) {

    TH2F* obj = dynamic_cast<TH2F*>( o.object );
    assert( obj );

    std::string name = o.name.substr(o.name.rfind("/")+1);

    if ( name.find( "reportSummaryMap" ) != std::string::npos ) {
      c->SetGridx();
      c->SetGridy();
      return;
    }

  }
  
};

static BeamRenderPlugin instance;
