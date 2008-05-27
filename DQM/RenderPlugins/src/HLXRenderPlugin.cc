// -*- C++ -*-
//
// Package:     RenderPlugins
// Class  :     HLXRenderPlugin
// 
// Implementation:
//     <Notes on implementation>
//
// Original Author:  
//         Created:  Sat Apr 19 20:02:57 CEST 2008
// $Id: HLXRenderPlugin.cc,v 1.4 2008/05/13 07:13:53 neadam Exp $
//

// system include files
#include <iostream>
#include <cassert>

// user include files
#include "DQM/RenderPlugins/src/HLXRenderPlugin.h"

#include "TH1F.h"
#include "TProfile.h"

#include "TStyle.h"
#include "TCanvas.h"
#include "TColor.h"

#include "TGraph.h"
#include "TLine.h"
#include "TROOT.h"


//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// member functions
//
void HLXRenderPlugin::initialise( int argc, char **argv ) 
{

}

bool HLXRenderPlugin::applies( const DQMNet::CoreObject &o, const VisDQMImgInfo &i ) 
{
  if( o.name.find( "HLX/Luminosity" )  == 0 ){
    return true;
  }
  if( o.name.find( "HLX/HFMinus" ) == 0 ){
    return true;
  }
  if( o.name.find( "HLX/HFPlus" ) == 0 ){
    return true;
  }
  if( o.name.find( "HLX/HFCompare" ) == 0 ){
    return true;
  }
  if( o.name.find( "HLX/Average" ) == 0 ){ 
    return true;
  }
  if( o.name.find( "HLX/CheckSums" ) == 0 ){ 
    return true;
  }

  return false;
}

void HLXRenderPlugin::preDraw( TCanvas *c, const DQMNet::CoreObject &o, 
			       const VisDQMImgInfo &i, VisDQMRenderInfo &r ) 
{

   c->cd();

   gStyle->Reset("Plain");
   gStyle->SetErrorX(0);
   gStyle->SetEndErrorSize(10);
   gStyle->SetCanvasColor(10);
   gStyle->SetPadColor(10);
   gStyle->SetFillColor(10);
   gStyle->SetOptTitle(kTRUE);
   gStyle->SetTitleBorderSize(0);
   gStyle->SetOptStat(kTRUE);
   gStyle->SetStatBorderSize(1);
   gStyle->SetStatColor(kWhite);
   gStyle->SetTitleFillColor(kWhite);
   gStyle->SetPadTickX(1);
   gStyle->SetPadTickY(1);
   gStyle->SetPadLeftMargin(0.15);

   gROOT->ForceStyle();

   if( dynamic_cast<TProfile*>( o.object ) ) {
      preDrawTProfile( c, o );
   }
   else if( dynamic_cast<TH1F*>( o.object ) ) {
      preDrawTH1F( c, o );
   }

   r.drawOptions = "";

}

void HLXRenderPlugin::preDrawTProfile( TCanvas *c, const DQMNet::CoreObject &o ) 
{
   TProfile* obj = dynamic_cast<TProfile*>( o.object );

   assert( obj );

   // Average plots for Etsum/Tower Occupancy
   obj->SetStats(kTRUE);
   obj->GetYaxis()->SetTitleOffset(1.6);

   if( o.name.find("EtSum") != std::string::npos && o.name.find("Lumi") == std::string::npos )
   {
      if( obj->GetMaximum() > 0 ) gPad->SetLogy();
   }

   return;
}

void HLXRenderPlugin::preDrawTH1F( TCanvas *c, const DQMNet::CoreObject &o ) 
{

   TH1F* obj = dynamic_cast<TH1F*>( o.object );

   assert( obj );

   // Average plots for Etsum/Tower Occupancy
   obj->SetStats(kTRUE);
   
   obj->SetFillColor(kBlue);
   obj->SetLineColor(kBlue);
   obj->SetFillStyle(1001);
   obj->GetYaxis()->SetTitleOffset(1.6);

   if( o.name.find("EtSum") != std::string::npos || o.name.find("ETSum") != std::string::npos )
   {
      if( obj->GetMaximum() > 0 ) gPad->SetLogy();
   }
   
   return;
}

// Anything that needs to be done after the histogram is drawn
void HLXRenderPlugin::postDraw( TCanvas *c, const DQMNet::CoreObject &o, 
			       const VisDQMImgInfo &i ) 
{

   // Do nothing here for now ...

//    c->cd();
//    if( dynamic_cast<TProfile*>( o.object ) ) {
//       postDrawTProfile2D( c, o );
//    }
//    else if( dynamic_cast<TH1F*>( o.object ) ) {
//       postDrawTH1F( c, o );
//    }


}

//
// const member functions
//

//
// static member functions
//
