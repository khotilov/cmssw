/*!
  \file SiPixelRenderPlugin
  \brief Display Plugin for Pixel DQM Histograms
  \author P.Merkel
  \version $Revision: 1.8 $
  \date $Date: 2009/02/16 11:50:04 $
*/

#include "TProfile2D.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TColor.h"
#include "TText.h"
#include <cassert>

#include "DQM/RenderPlugins/src/SiPixelRenderPlugin.h"
#include "DQM/RenderPlugins/src/utils.h"

bool SiPixelRenderPlugin::applies( const DQMNet::CoreObject &o, const VisDQMImgInfo &i ) {
 
  if( o.name.find( "Pixel/" ) != std::string::npos ) {
    return true;
  } 

  return false;

}

void SiPixelRenderPlugin::preDraw( TCanvas *c, const DQMNet::CoreObject &o, const VisDQMImgInfo &i, VisDQMRenderInfo &r ) {

  c->cd();

  if( dynamic_cast<TH2*>( o.object ) ) {
    preDrawTH2( c, o );
  }
  else if( dynamic_cast<TH1*>( o.object ) ) {
    preDrawTH1( c, o );
  }

}


void SiPixelRenderPlugin::preDrawTH2( TCanvas *c, const DQMNet::CoreObject &o ) {

  TH2* obj = dynamic_cast<TH2*>( o.object );

  assert( obj );

  // This applies to all
  gStyle->SetCanvasBorderMode( 0 );
  gStyle->SetPadBorderMode( 0 );
  gStyle->SetPadBorderSize( 0 );
  //    (data->pad)->SetLogy( 0 );;
  gStyle->SetOptStat( 0 );
  obj->SetStats( kFALSE );

  TAxis* xa = obj->GetXaxis();
  TAxis* ya = obj->GetYaxis();

  xa->SetTitleOffset(0.7);
  xa->SetTitleSize(0.05);
  xa->SetLabelSize(0.04);

  ya->SetTitleOffset(0.7);
  ya->SetTitleSize(0.05);
  ya->SetLabelSize(0.04);

  if( o.name.find( "hitmap" ) != std::string::npos  ||
      o.name.find( "Occupancy" ) != std::string::npos) {
    gStyle->SetPalette(1);
    obj->SetOption("colz");
    return;
  }
    
  TH2F* obj2 = dynamic_cast<TH2F*>( o.object );

   if( o.name.find( "reportSummaryMap" ) != std::string::npos ) {
     dqm::utils::reportSummaryMapPalette(obj2);
     return;
   }  

  return;
}

void SiPixelRenderPlugin::preDrawTH1( TCanvas *c, const DQMNet::CoreObject &o ) {

  TH1* obj = dynamic_cast<TH1*>( o.object );

  assert( obj );

  // This applies to all
  gStyle->SetOptStat(0111);
  if ( obj->GetMaximum(1.e5) > 0. ) {
    gPad->SetLogy(1);
  } else {
   gPad->SetLogy(0);
  }
  
  if( o.name.find( "adcCOMB" ) != std::string::npos && obj->GetEntries() > 0. ) gPad->SetLogy(1);
  if( o.name.find( "size_siPixelClusters" ) != std::string::npos && obj->GetEntries() > 0. ) gPad->SetLogx(1);
  if( o.name.find( "OnTrack" ) != std::string::npos && o.name.find( "charge" ) != std::string::npos ) obj->SetTitle("ClusterCharge_OnTrack");
  if( o.name.find( "OnTrack" ) != std::string::npos && o.name.find( "size" ) != std::string::npos ) obj->SetTitle("ClusterSize_OnTrack");
  if( o.name.find( "OffTrack" ) != std::string::npos && o.name.find( "charge" ) != std::string::npos ) obj->SetTitle("ClusterCharge_OffTrack");
  if( o.name.find( "OffTrack" ) != std::string::npos && o.name.find( "size" ) != std::string::npos ) obj->SetTitle("ClusterSize_OffTrack");

}

void SiPixelRenderPlugin::postDraw( TCanvas *c, const DQMNet::CoreObject &o, const VisDQMImgInfo &i ) {

#ifdef DEBUG
  std::cout << "SiPixelRenderPlugin:postDraw " << o.name << std::endl;
#endif

  c->cd();

  if( dynamic_cast<TH1*>( o.object ) ) {
    postDrawTH1( c, o );
  }

#ifdef DEBUG
  std::cout << "done" << std::endl;
#endif

}

void SiPixelRenderPlugin::postDrawTH1( TCanvas *c, const DQMNet::CoreObject &o ) {

  TText tt;
  tt.SetTextSize(0.12);
  if (o.flags == 0) return;
  else {
/*    if (o.flags & DQMNet::DQM_FLAG_REPORT_ERROR) {
      tt.SetTextColor(2);
      tt.DrawTextNDC(0.5, 0.5, "Error");
    } else if (o.flags & DQMNet::DQM_FLAG_REPORT_WARNING) {
      tt.SetTextColor(5);
      tt.DrawTextNDC(0.5, 0.5, "Warning");
    } else if (o.flags & DQMNet::DQM_FLAG_REPORT_OTHER) { 
      tt.SetTextColor(1);
      tt.DrawTextNDC(0.5, 0.5, "Other ");      
    } else {
      tt.SetTextColor(3);
      tt.DrawTextNDC(0.5, 0.5, "Ok ");
    }     */ 
  }  

  return;

}
