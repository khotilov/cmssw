/*!
  \file SiStripRenderPlugin
  \brief Display Plugin for SiStrip DQM Histograms
  \author S. Dutta 
  \version $Revision: 1.7 $
  \date $Date: 2008/06/18 20:00:50 $
*/

#include "TProfile2D.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TColor.h"
#include "TText.h"
#include <cassert>

#include "DQM/RenderPlugins/src/SiStripRenderPlugin.h"
#include "DQM/RenderPlugins/src/utils.h"

bool SiStripRenderPlugin::applies( const DQMNet::CoreObject &o, const VisDQMImgInfo &i ) {
 
  if (o.name.find( "SiStrip/" )   == std::string::npos) return false;

  if( o.name.find( "/EventInfo/" ) != std::string::npos ) {
    return true;
  } 

  if( o.name.find( "/MechanicalView/" ) != std::string::npos ) {
    return true;
  } 

  if( o.name.find( "/ReadoutView/" ) != std::string::npos ) {
    return true;
  } 

  if( o.name.find( "/Tracks/" ) != std::string::npos ) {
    return true;
  } 

  return false;

}

void SiStripRenderPlugin::preDraw( TCanvas *c, const DQMNet::CoreObject &o, const VisDQMImgInfo &i, VisDQMRenderInfo &r ) {

  c->cd();

  if( dynamic_cast<TH2F*>( o.object ) ) {
    preDrawTH2F( c, o );
  }
  else if( dynamic_cast<TH1F*>( o.object ) ) {
    preDrawTH1F( c, o );
  }
}


void SiStripRenderPlugin::preDrawTH2F( TCanvas *c, const DQMNet::CoreObject &o ) {

  TH2F* obj = dynamic_cast<TH2F*>( o.object );

  assert( obj );

  // This applies to all
  gStyle->SetCanvasBorderMode( 0 );
  gStyle->SetPadBorderMode( 0 );
  gStyle->SetPadBorderSize( 0 );
  //    (data->pad)->SetLogy( 0 );;
  //  gStyle->SetOptStat( 0 );


  TAxis* xa = obj->GetXaxis();
  TAxis* ya = obj->GetYaxis();

  xa->SetTitleOffset(0.7);
  xa->SetTitleSize(0.05);
  xa->SetLabelSize(0.04);

  ya->SetTitleOffset(0.7);
  ya->SetTitleSize(0.05);
  ya->SetLabelSize(0.04);

  
  if( o.name.find( "PedsEvolution" ) != std::string::npos) {
    gStyle->SetOptStat( 1111 );
    obj->SetStats( kTRUE );
    obj->SetOption( "lego2" );
    return;
  } 
  if( o.name.find( "CMDistribution " )  != std::string::npos){
    obj->GetXaxis()->LabelsOption("d");
    obj->SetOption( "lego2" );
    return;
  }
  if( o.name.find( "CMSlopeDistribution " )  != std::string::npos) {
    obj->GetXaxis()->LabelsOption("d");
    obj->SetOption( "lego2" );
    return;
  }

  if( o.name.find( "PedestalDistribution " )  != std::string::npos) {
    obj->GetXaxis()->LabelsOption("d");
    obj->SetOption( "lego" );
    return;
  }

  if( o.name.find( "reportSummaryMap" )  != std::string::npos) {
    obj->SetStats( kFALSE );
    dqm::utils::reportSummaryMapPalette(obj);
    obj->SetOption("colztext");
    return;
  }

  if( o.name.find( "SummaryOfCabling" )  != std::string::npos) {
    obj->SetStats( kFALSE );
    obj->SetOption("text");
    return;
  }
  return;
}

void SiStripRenderPlugin::preDrawTH1F( TCanvas *c, const DQMNet::CoreObject &o ) {

  TH1F* obj = dynamic_cast<TH1F*>( o.object );

  assert( obj );

  // This applies to all
  gStyle->SetOptStat(1110);
//  if ( obj->GetMaximum(1.e5) > 0. ) {
//    gPad->SetLogy(1);
//  } else {
//   gPad->SetLogy(0);
//  }
  if( o.name.find( "Summary_MeanNumberOfDigis" )  != std::string::npos) {
    obj->SetStats( kFALSE );
    obj->SetMaximum(0.5);
    obj->SetMinimum(-0.01);
    return;
  }

  if( o.name.find( "Summary_MeanNumberOfClusters" )  != std::string::npos) {
    obj->SetStats( kFALSE );
    obj->SetMaximum(0.05);
    obj->SetMinimum(-0.001);
    return;
  }

  if( o.name.find( "Summary_MeanClusterWidth" )  != std::string::npos) {
    obj->SetStats( kFALSE );
    obj->SetMaximum(20.0);
    obj->SetMinimum(-0.1);
    return;
  }

}

void SiStripRenderPlugin::postDraw( TCanvas *c, const DQMNet::CoreObject &o, const VisDQMImgInfo &i ) {

#ifdef DEBUG
  std::cout << "SiStripRenderPlugin:postDraw " << o.name << std::endl;
#endif

  c->cd();

  if( dynamic_cast<TH1F*>( o.object ) ) {
    postDrawTH1F( c, o );
  }

  if( dynamic_cast<TH2F*>( o.object ) ) {
    postDrawTH2F( c, o );
  }


#ifdef DEBUG
  std::cout << "done" << std::endl;
#endif

}

void SiStripRenderPlugin::postDrawTH1F( TCanvas *c, const DQMNet::CoreObject &o ) {

  TText tt;
  tt.SetTextSize(0.12);
  if (o.flags == 0) return;
  else {
    if (o.flags & DQMNet::DQM_FLAG_REPORT_ERROR) {
      tt.SetTextColor(2);
      tt.DrawTextNDC(0.5, 0.5, "Error");
    } else if (o.flags & DQMNet::DQM_FLAG_REPORT_WARNING) {
      tt.SetTextColor(5);
      tt.DrawTextNDC(0.5, 0.5, "Warning");
    } else if (o.flags & DQMNet::DQM_FLAG_REPORT_OTHER) { 
      tt.SetTextColor(1);
      tt.DrawTextNDC(0.5, 0.5, "Other ");      
    } 
  }  

  return;

}
void SiStripRenderPlugin::postDrawTH2F( TCanvas *c, const DQMNet::CoreObject &o ) {

  TH2F* obj = dynamic_cast<TH2F*>( o.object );

  assert( obj );

  std::string name = o.name.substr(o.name.rfind("/")+1);

  if( name.find( "reportSummaryMap" ) != std::string::npos ) {
    c->SetGridx();
    c->SetGridy();
    return;
  }

  if( name.find( "SummaryOfCabling" ) != std::string::npos ) {
     c->SetGridx();
     c->SetGridy();
     return;
   }
}
