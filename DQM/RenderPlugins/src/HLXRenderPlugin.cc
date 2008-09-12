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
// $Id: HLXRenderPlugin.cc,v 1.9 2008/08/22 11:52:02 lat Exp $
//

// system include files
#include <iostream>
#include <cassert>

// user include files
#include "DQM/RenderPlugins/src/HLXRenderPlugin.h"
#include "DQM/RenderPlugins/src/utils.h"

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
  if( o.name.find( "HLX/Luminosity" )  != std::string::npos ){
    return true;
  }
  if( o.name.find( "HLX/HFMinus" ) != std::string::npos ){
    return true;
  }
  if( o.name.find( "HLX/HFPlus" ) != std::string::npos ){
    return true;
  }
  if( o.name.find( "HLX/HFCompare" ) != std::string::npos ){
    return true;
  }
  if( o.name.find( "HLX/Average" ) != std::string::npos ){ 
    return true;
  }
  if( o.name.find( "HLX/CheckSums" ) != std::string::npos ){ 
    return true;
  }
  if( o.name.find( "HLX/History" ) != std::string::npos ){ 
    return true;
  }
  if( o.name.find( "HLX/EventInfo" ) != std::string::npos ){ 
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

   r.drawOptions = "";
   if( dynamic_cast<TProfile*>( o.object ) ) {
      preDrawTProfile( c, o );
   }
   else if( dynamic_cast<TH1F*>( o.object ) ) {
      preDrawTH1F( c, o );
   }
   else if( dynamic_cast<TH2F*>( o.object ) ) {
      preDrawTH2F( c, o );
      //r.drawOptions = "colz";
   }


}

void HLXRenderPlugin::preDrawTProfile( TCanvas *c, const DQMNet::CoreObject &o ) 
{
   //std::cout << "NADIA: Here in predraw PROFILE!!!!" << std::endl;

   TProfile* obj = dynamic_cast<TProfile*>( o.object );

   assert( obj );

   // Average plots for Etsum/Tower Occupancy
   obj->SetStats(kTRUE);
   obj->GetYaxis()->SetTitleOffset(1.6);

   if( o.name.find("EtSum") != std::string::npos && 
       o.name.find("Lumi") == std::string::npos && 
       o.name.find("Hist") == std::string::npos )
   {
      int maxBin = obj->GetMaximumBin();
      int minBin = obj->GetMinimumBin();

      double objMax = obj->GetBinContent(maxBin);
      double objMin = obj->GetBinContent(minBin);
      if( objMax > 0 && objMin > 0 ) 
      {
	 if( (objMax/objMin) > 1e4 ) gPad->SetLogy();
      }

      obj->SetMaximum(objMax*1.02);
      obj->SetMinimum(objMin*0.98);
   }

   if( o.name.find("SumAllOcc") != std::string::npos )
   {
      obj->SetMaximum( 1.25*obj->GetMaximum() );
      obj->SetMinimum( 0.75*obj->GetMinimum() );
   }

   // History histograms ...
   if( o.name.find("Hist") != std::string::npos )
   {
      obj->SetMarkerStyle(kFullCircle);
      obj->SetMarkerSize(0.8);
      // Loop over the bins and find the first empty one
      int firstEmptyBin = obj->GetNbinsX();
      for( int iBin = 1; iBin<=obj->GetNbinsX(); ++iBin )
      {
	 if( obj->GetBinContent(iBin) == 0 )
	 {
	    firstEmptyBin = iBin;
	    break;
	 }
      }
      if( firstEmptyBin < 10 ) firstEmptyBin = 10;
      // Now set the new range
      obj->GetXaxis()->SetRange(0,firstEmptyBin+1);
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
      int maxBin = obj->GetMaximumBin();
      int minBin = obj->GetMinimumBin();

      double objMax = obj->GetBinContent(maxBin);
      double objMin = obj->GetBinContent(minBin);
      if( objMax > 0 && objMin > 0 ) 
      {
	 if( (objMax/objMin) > 1e4 ) gPad->SetLogy();
      }

      obj->SetMaximum(objMax*1.02);
      obj->SetMinimum(objMin*0.98);
   }
   
   return;
}

void HLXRenderPlugin::preDrawTH2F( TCanvas *c, const DQMNet::CoreObject &o )
{

   TH2F* obj = dynamic_cast<TH2F*>( o.object );

   assert( obj );
   obj->SetStats(kFALSE);
   obj->SetMinimum(-.001);
   obj->SetMaximum(1.02);

   dqm::utils::reportSummaryMapPalette(obj);
   obj->SetOption("colztext");

//    unsigned int Number = 3;
//    double       Red[]   = { 0.00, 0.00, 1.00};
//    double       Green[] = { 0.00, 1.00, 0.00};
//    double       Blue[]  = { 1.00, 0.00, 0.00};
//    double       Stops[] = { 0.00, 0.50, 1.00 };
//    int          nb = 100;
//    int          i = TColor::CreateGradientColorTable(Number,Stops,Blue,Green,Red,nb);
//    //gStyle->SetPalette(1);

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
