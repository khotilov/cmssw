// $Id: EBRenderPlugin.cc,v 1.31 2007/11/29 14:42:56 dellaric Exp $

/*!
  \file EBRenderPlugin
  \brief Display Plugin for Quality Histograms
  \author G. Della Ricca
  \author B. Gobbo 
  \version $Revision: 1.31 $
  \date $Date: 2007/11/29 14:42:56 $
*/

#include <TH1F.h>
#include <TH2F.h>
#include <TH3F.h>
#include <TProfile.h>
#include <TProfile2D.h>

#include <TStyle.h>
#include <TCanvas.h>
#include <TColor.h>

#include <iostream>

#include "DQMServices/Core/interface/MonitorElement.h"

#include "DQM/EcalCommon/interface/ColorPalette.h"
#include <DQM/EcalCommon/interface/Numbers.h>
#include "EBRenderPlugin.h"

static bool  first = true;
static TH2C* t1 = new TH2C( "eb_text1", "text1", 85, 0,  85, 20,   0, 20 );
static TH2C* t2 = new TH2C( "eb_text2", "text2", 17, 0,  17,  4,   0,  4 );
static TH2C* t3 = new TH2C( "eb_text3", "text3", 10, 0,  10,  5,   0,  5 );
static TH2C* t4 = new TH2C( "eb_text4", "text4",  2, 0,   2,  1,   0,  1 );
static TH2C* t6 = new TH2C( "eb_text6", "text6", 18, 0, 360,  2, -85, 85 );
static TH2C* t7 = new TH2C( "eb_text7", "text7", 18, -M_PI*(9+1.5)/9, M_PI*(9-1.5)/9, 2, -1.479, 1.479);
static TH2C* t8 = new TH2C( "eb_text8", "text8", 18, 0., 72., 2, -17., 17. );

static TGaxis* ax1 = new TGaxis(-M_PI*(9+1.5)/9, -1.479, M_PI*(9-1.5)/9, -1.479, -M_PI*(9+1.5)/9, M_PI*(9-1.5)/9, 40306, "N");

void EBRenderPlugin::initialise( int argc, char **argv, DaqMonitorBEInterface *bei ) {

  if( ! first ) return;

  first = false;

  for( short i=0; i<6; i++ ) {
    TColor* color;
    if( ! gROOT->GetColor( 301+i )) {
      color = new TColor( 301+i, ecdqm::rgb[i][0], ecdqm::rgb[i][1], ecdqm::rgb[i][2], "" );
    }
    else {
      color = gROOT->GetColor( 301+i );
      color->SetRGB( ecdqm::rgb[i][0], ecdqm::rgb[i][1], ecdqm::rgb[i][2] );
    }
  }
  for( short i=0; i<6; i++ ) pCol3[i]  = i+301;

  for( short i=0; i<10; i++ ) {
    TColor* color;
    if( ! gROOT->GetColor( 401+i )) {
      color = new TColor( 401+i, ecdqm::rgb2[i][0], ecdqm::rgb2[i][1], ecdqm::rgb2[i][2], "" );
    }
    else {
      color = gROOT->GetColor( 401+i );
      color->SetRGB( ecdqm::rgb2[i][0], ecdqm::rgb2[i][1], ecdqm::rgb2[i][2] );
    }
  }
  for( short i=0; i<10; i++ ) pCol4[i] = i+401;

  text1 = t1;
  text2 = t2;
  text3 = t3;
  text4 = t4;
  text6 = t6;
  text7 = t7;
  text8 = t8;

  text1->SetMinimum(   0.10 );
  text2->SetMinimum(   0.10 );
  text3->SetMinimum(   0.10 );
  text4->SetMinimum(   0.10 );
  text6->SetMinimum( -18.01 );
  text7->SetMinimum( -18.01 );
  text8->SetMinimum( -18.01 );

  for( short i=0; i<68; i++ ) {
    text1->Fill( 2+(i/4)*5, 2+(i%4)*5, i+1 );
    text2->Fill( i/4, i%4, i+1 );
  }
  for( short i=0; i<2; i++ ) {
    text3->Fill( 2+i*5, 2, i+1+68 );
    text4->Fill( i, 0., i+1+68 );
  }
  for ( short i=0; i<36; i++ ) {
    int x = 1 + i%18;
    int y = 1 + i/18;
    text6->SetBinContent(x, y, Numbers::iEB(i+1));
  }
  for ( short i=0; i<36; i++ ) {
    int x = 1 + i%18;
    int y = 2 - i/18;
    int z = x + 8;
    if ( z > 18 ) z = z - 18;
    if ( y == 1 ) {
      text7->SetBinContent(x, y, -z);
    } else {
      text7->SetBinContent(x, y, +z);
    }
  }
  for ( short i=0; i<36; i++ ) {
    int x = 1 + i%18;
    int y = 1 + i/18;
    text8->SetBinContent(x, y, Numbers::iEB(i+1));
  }

  text1->SetMarkerSize( 2 );
  text2->SetMarkerSize( 2 );
  text3->SetMarkerSize( 2 );
  text4->SetMarkerSize( 2 );
  text6->SetMarkerSize( 2 );
  text7->SetMarkerSize( 2 );
  text8->SetMarkerSize( 2 );

}

bool EBRenderPlugin::applies( const ObjInfo &o, const ImgInfo &i ) {
 
  if( o.name.find( "EcalBarrel/EB" ) == 0 ) {
    return true;
  } 

  if( o.name.find( "EcalBarrel/EcalInfo" ) == 0 ) {
    return true;
  }

  return false;

}

void EBRenderPlugin::preDraw( TCanvas *c, const ObjInfo &o, const ImgInfo &i, RenderInfo &r ) {

#ifdef DEBUG
  std::cout << "EBRenderPlugin:preDraw " << o.name << std::endl;
#endif

  c->cd();

  if( dynamic_cast<TProfile2D*>( o.object ) ) {
    preDrawTProfile2D( c, o );
  }
  else if( dynamic_cast<TProfile*>( o.object ) ) {
    preDrawTProfile( c, o );
  }
  else if( dynamic_cast<TH3F*>( o.object ) ) {
    preDrawTH3F( c, o );
  }
  else if( dynamic_cast<TH2F*>( o.object ) ) {
    preDrawTH2F( c, o );
  }
  else if( dynamic_cast<TH1F*>( o.object ) ) {
    preDrawTH1F( c, o );
  }

#ifdef DEBUG
  std::cout << "done" << std::endl;
#endif

}

void EBRenderPlugin::preDrawTProfile2D( TCanvas *c, const ObjInfo &o ) {

  TProfile2D* obj = dynamic_cast<TProfile2D*>( o.object );

  assert( obj );

  gStyle->SetPaintTextFormat();

  gStyle->SetOptStat(0);
  obj->SetStats(kFALSE);
  gPad->SetLogy(0);

  // Occupancy-like (10 x grays) plots
  if( o.name.find( "EBCLT" ) < o.name.size() ) {
    gPad->SetGridx();
    gPad->SetGridy();
    obj->GetXaxis()->SetNdivisions(40118, kFALSE);
    obj->GetYaxis()->SetNdivisions(170102, kFALSE);
    obj->SetMinimum(0.0);
    gStyle->SetPalette(10, pCol4);
    obj->SetOption("colz");
    gStyle->SetPaintTextFormat("+g");
    return;
  }

  int nbx = obj->GetNbinsX();
  int nby = obj->GetNbinsY();

  if( nbx == 72 && nby == 34 ) {
    gPad->SetGridx();
    gPad->SetGridy();
    obj->GetXaxis()->SetNdivisions(18, kFALSE);
    obj->GetYaxis()->SetNdivisions(2, kFALSE);
  }

  if( nbx == 85 && nby == 20 ) {
    gPad->SetGridx();
    gPad->SetGridy();
    obj->GetXaxis()->SetNdivisions(17);
    obj->GetYaxis()->SetNdivisions(4);
  }

  if( nbx == 17 && nby == 4 ) {
    gPad->SetGridx();
    gPad->SetGridy();
    obj->GetXaxis()->SetNdivisions(17);
    obj->GetYaxis()->SetNdivisions(4);
  }

  // Occupancy-like (10 x grays) plots
  obj->SetMinimum(0.0);
  gStyle->SetPalette(10, pCol4);
  obj->SetOption("colz");
  return;

}

void EBRenderPlugin::preDrawTProfile( TCanvas *c, const ObjInfo &o ) {

  TProfile* obj = dynamic_cast<TProfile*>( o.object );

  assert( obj );

  gStyle->SetOptStat("euomr");
  obj->SetStats(kTRUE);
  gPad->SetLogy(0);

  return;

}

void EBRenderPlugin::preDrawTH3F( TCanvas *c, const ObjInfo &o ) {
  
  TH3F* obj = dynamic_cast<TH3F*>( o.object );
  
  assert( obj );

  gStyle->SetOptStat(0);
  obj->SetStats( kFALSE );
  gPad->SetLogy(0);

  return;
  
} 

void EBRenderPlugin::preDrawTH2F( TCanvas *c, const ObjInfo &o ) {

  TH2F* obj = dynamic_cast<TH2F*>( o.object );

  assert( obj );

  gStyle->SetPaintTextFormat();

  gStyle->SetOptStat(0);
  obj->SetStats( kFALSE );
  gPad->SetLogy(0);

  // Occupancy-like (10 x grays) plots
  if( o.name.find( "EBCLT" ) < o.name.size() ) {
    gPad->SetGridx();
    gPad->SetGridy();
    obj->GetXaxis()->SetNdivisions(40118, kFALSE);
    obj->GetYaxis()->SetNdivisions(170102, kFALSE);
    obj->SetMinimum(0.0);
    gStyle->SetPalette(10, pCol4);
    obj->SetOption("colz");
    gStyle->SetPaintTextFormat("+g");
    return;
  }

  int nbx = obj->GetNbinsX();
  int nby = obj->GetNbinsY();

  if( nbx == 85 && nby == 20 ) {
    gPad->SetGridx();
    gPad->SetGridy();
    obj->GetXaxis()->SetNdivisions(17);
    obj->GetYaxis()->SetNdivisions(4);
  }

  if( nbx == 17 && nby == 4 ) {
    gPad->SetGridx();
    gPad->SetGridy();
    obj->GetXaxis()->SetNdivisions(17);
    obj->GetYaxis()->SetNdivisions(4);
  }

  if( nbx == 10 && nby == 5 ) {
    gPad->SetGridx();
    gPad->SetGridy();
    obj->GetXaxis()->SetNdivisions(10);
    obj->GetYaxis()->SetNdivisions(1);
  }

  if( nbx == 2 && nby == 1 ) {
    gPad->SetGridx();
    gPad->SetGridy();
    obj->GetXaxis()->SetNdivisions(2);
    obj->GetYaxis()->SetNdivisions(1);
  }

  if( nbx == 360 && nby == 170 ) { 
    gPad->SetGridx();
    gPad->SetGridy();
    obj->GetXaxis()->SetNdivisions(18, kFALSE);
    obj->GetYaxis()->SetNdivisions(2, kFALSE);
  }

  if( nbx == 90 && nby == 20 ) { 
    gPad->SetGridx();
    gPad->SetGridy();
    obj->GetXaxis()->SetNdivisions(18, kFALSE);
    obj->GetYaxis()->SetNdivisions(2, kFALSE);
  }

  if( nbx == 72 && nby == 34 ) {
    gPad->SetGridx();
    gPad->SetGridy();
    obj->GetXaxis()->SetNdivisions(18, kFALSE);
    obj->GetYaxis()->SetNdivisions(2, kFALSE);
  }

  // Occupancy-like (10 x grays) plots
  if( o.name.find( "EBIT" ) < o.name.size() &&
      o.name.find( "quality" ) >= o.name.size() ) {
    obj->SetMinimum(0.0);
    gStyle->SetPalette(10, pCol4);
    obj->SetOption("colz");
    gStyle->SetPaintTextFormat("+g");
    return;
  }

  // Occupancy-like (10 x grays) plots
  if( o.name.find( "EBTTT" ) < o.name.size() &&
      o.name.find( "quality" ) >= o.name.size() ) {
    obj->SetMinimum(0.0);
    gStyle->SetPalette(10, pCol4);
    obj->SetOption("colz");
    gStyle->SetPaintTextFormat("+g");
    return;
  }

  // Occupancy-like (10 x grays) plots
  if( o.name.find( "EBOT" ) < o.name.size() ) {
    obj->SetMinimum(0.0);
    gStyle->SetPalette(10, pCol4);
    obj->SetOption("colz");
    gStyle->SetPaintTextFormat("+g");
    return;
  }

  // Occupancy-like (10 x grays) plots
  if( o.name.find( "EBCT" ) < o.name.size() ) {
    obj->SetMinimum(0.0);
    gStyle->SetPalette(10, pCol4);
    obj->SetOption("colz");
    gStyle->SetPaintTextFormat("+g");
    return;
  }

  // Quality-like (green, yellow, red) plots
  if( o.name.find( "summary" ) < o.name.size() ) {
    gStyle->SetOptStat(" ");
    obj->SetMinimum(-0.00000001);
    obj->SetMaximum(6.0);
    gStyle->SetPalette(6, pCol3);
    obj->SetOption("col");
    gStyle->SetPaintTextFormat("+g");
    return;
  }

  // Quality-like (green, yellow, red) plots
  if( o.name.find( "quality" ) < o.name.size() ) {
    obj->SetMinimum(-0.00000001);
    obj->SetMaximum(6.0);
    gStyle->SetPalette(6, pCol3);
    obj->SetOption("col");
    gStyle->SetPaintTextFormat("+g");
    return;
  }

  // Occupancy-like (10 x grays) plots
  if( o.name.find( "EBMM event" ) < o.name.size() ) {
    obj->SetMinimum(0.0);
    gStyle->SetPalette(10, pCol4);
    obj->SetOption("colz");
    gStyle->SetPaintTextFormat("+g");
    return;
  }

}

void EBRenderPlugin::preDrawTH1F( TCanvas *c, const ObjInfo &o ) {

  TH1F* obj = dynamic_cast<TH1F*>( o.object );

  assert( obj );

  gStyle->SetOptStat("euomr");
  obj->SetStats(kTRUE);
  gPad->SetLogy(0);

  int nbx = obj->GetNbinsX();

  if ( obj->GetMaximum() > 0. ) gPad->SetLogy(1);

  if ( nbx == 10 ) gPad->SetLogy(0);
  if ( nbx == 1700 ) gPad->SetLogy(0);

  if( o.name.find( "EVTTYPE" ) < o.name.size() ) {
   gPad->SetBottomMargin(0.4);
   obj->GetXaxis()->LabelsOption("v");
  }

  if( o.name.find( "EBMM DCC" ) < o.name.size() ) {
   gPad->SetBottomMargin(0.2);
   obj->GetXaxis()->LabelsOption("v");
  }

  return;

}

void EBRenderPlugin::postDraw( TCanvas *c, const ObjInfo &o, const ImgInfo &i ) {

#ifdef DEBUG
  std::cout << "EBRenderPlugin:postDraw " << o.name << std::endl;
#endif

  c->cd();

  if( dynamic_cast<TProfile2D*>( o.object ) ) {
    postDrawTProfile2D( c, o );
  }
  else if( dynamic_cast<TH3F*>( o.object ) ) {
    postDrawTH3F( c, o );
  }
  else if( dynamic_cast<TH2F*>( o.object ) ) {
    postDrawTH2F( c, o );
  }

#ifdef DEBUG
  std::cout << "done" << std::endl;
#endif

}

void EBRenderPlugin::postDrawTProfile2D( TCanvas *c, const ObjInfo &o ) {

  TProfile2D* obj = dynamic_cast<TProfile2D*>( o.object );

  assert( obj );

  if( o.name.find( "EBCLT" ) < o.name.size() ) {
    int x1 = text7->GetXaxis()->FindBin(obj->GetXaxis()->GetXmin());
    int x2 = text7->GetXaxis()->FindBin(obj->GetXaxis()->GetXmax());
    int y1 = text7->GetYaxis()->FindBin(obj->GetYaxis()->GetXmin());
    int y2 = text7->GetYaxis()->FindBin(obj->GetYaxis()->GetXmax());
    text7->GetXaxis()->SetRange(x1, x2);
    text7->GetYaxis()->SetRange(y1, y2);
    text7->Draw("text,same");
    //ax1->Draw();
    return;
  }

  int nbx = obj->GetNbinsX();
  int nby = obj->GetNbinsY();

  if( nbx == 17 && nby == 4 ) {
    int x1 = text2->GetXaxis()->FindBin(obj->GetXaxis()->GetXmin());
    int x2 = text2->GetXaxis()->FindBin(obj->GetXaxis()->GetXmax());
    int y1 = text2->GetYaxis()->FindBin(obj->GetYaxis()->GetXmin());
    int y2 = text2->GetYaxis()->FindBin(obj->GetYaxis()->GetXmax());
    text2->GetXaxis()->SetRange(x1, x2);
    text2->GetYaxis()->SetRange(y1, y2);
    text2->Draw("text,same");
    return;
  }

  if( nbx == 72 && nby == 34 ) {
    int x1 = text8->GetXaxis()->FindBin(obj->GetXaxis()->GetXmin());
    int x2 = text8->GetXaxis()->FindBin(obj->GetXaxis()->GetXmax());
    int y1 = text8->GetYaxis()->FindBin(obj->GetYaxis()->GetXmin());
    int y2 = text8->GetYaxis()->FindBin(obj->GetYaxis()->GetXmax());
    text8->GetXaxis()->SetRange(x1, x2);
    text8->GetYaxis()->SetRange(y1, y2);
    text8->Draw("text,same");
    return;
  }

  int x1 = text1->GetXaxis()->FindBin(obj->GetXaxis()->GetXmin());
  int x2 = text1->GetXaxis()->FindBin(obj->GetXaxis()->GetXmax());
  int y1 = text1->GetYaxis()->FindBin(obj->GetYaxis()->GetXmin());
  int y2 = text1->GetYaxis()->FindBin(obj->GetYaxis()->GetXmax());
  text1->GetXaxis()->SetRange(x1, x2);
  text1->GetYaxis()->SetRange(y1, y2);
  text1->Draw("text,same");
  return;

}

void EBRenderPlugin::postDrawTH3F( TCanvas *c, const ObjInfo &o ) {

  TH3F* obj = dynamic_cast<TH3F*>( o.object );

  assert( obj );

  return;

}

void EBRenderPlugin::postDrawTH2F( TCanvas *c, const ObjInfo &o ) {

  TH2F* obj = dynamic_cast<TH2F*>( o.object );

  assert( obj );

  if( o.name.find( "EBCLT" ) < o.name.size() ) {
    int x1 = text7->GetXaxis()->FindBin(obj->GetXaxis()->GetXmin());
    int x2 = text7->GetXaxis()->FindBin(obj->GetXaxis()->GetXmax());
    int y1 = text7->GetYaxis()->FindBin(obj->GetYaxis()->GetXmin());
    int y2 = text7->GetYaxis()->FindBin(obj->GetYaxis()->GetXmax());
    text7->GetXaxis()->SetRange(x1, x2);
    text7->GetYaxis()->SetRange(y1, y2);
    text7->Draw("text,same");
    //ax1->Draw();
    return;
  }

  int nbx = obj->GetNbinsX();
  int nby = obj->GetNbinsY();

  if( nbx == 85 && nby == 20 ) {
    int x1 = text1->GetXaxis()->FindBin(obj->GetXaxis()->GetXmin());
    int x2 = text1->GetXaxis()->FindBin(obj->GetXaxis()->GetXmax());
    int y1 = text1->GetYaxis()->FindBin(obj->GetYaxis()->GetXmin());
    int y2 = text1->GetYaxis()->FindBin(obj->GetYaxis()->GetXmax());
    text1->GetXaxis()->SetRange(x1, x2);
    text1->GetYaxis()->SetRange(y1, y2);
    text1->Draw("text,same");
    return;
  }

  if( nbx == 17 && nby == 4 ) {
    int x1 = text2->GetXaxis()->FindBin(obj->GetXaxis()->GetXmin());
    int x2 = text2->GetXaxis()->FindBin(obj->GetXaxis()->GetXmax());
    int y1 = text2->GetYaxis()->FindBin(obj->GetYaxis()->GetXmin());
    int y2 = text2->GetYaxis()->FindBin(obj->GetYaxis()->GetXmax());
    text2->GetXaxis()->SetRange(x1, x2);
    text2->GetYaxis()->SetRange(y1, y2);
    text2->Draw("text,same");
    return;
  }

  if( nbx == 10 && nby == 5 ) {
    int x1 = text3->GetXaxis()->FindBin(obj->GetXaxis()->GetXmin());
    int x2 = text3->GetXaxis()->FindBin(obj->GetXaxis()->GetXmax());
    int y1 = text3->GetYaxis()->FindBin(obj->GetYaxis()->GetXmin());
    int y2 = text3->GetYaxis()->FindBin(obj->GetYaxis()->GetXmax());
    text3->GetXaxis()->SetRange(x1, x2);
    text3->GetYaxis()->SetRange(y1, y2);
    text3->Draw("text,same");
    return;
  }

  if( nbx == 2 && nby == 1 ) {
    int x1 = text4->GetXaxis()->FindBin(obj->GetXaxis()->GetXmin());
    int x2 = text4->GetXaxis()->FindBin(obj->GetXaxis()->GetXmax());
    int y1 = text4->GetYaxis()->FindBin(obj->GetYaxis()->GetXmin());
    int y2 = text4->GetYaxis()->FindBin(obj->GetYaxis()->GetXmax());
    text4->GetXaxis()->SetRange(x1, x2);
    text4->GetYaxis()->SetRange(y1, y2);
    text4->Draw("text,same");
    return;
  }

  if( nbx == 72 && nby == 34 ) {
    int x1 = text8->GetXaxis()->FindBin(obj->GetXaxis()->GetXmin());
    int x2 = text8->GetXaxis()->FindBin(obj->GetXaxis()->GetXmax());
    int y1 = text8->GetYaxis()->FindBin(obj->GetYaxis()->GetXmin());
    int y2 = text8->GetYaxis()->FindBin(obj->GetYaxis()->GetXmax());
    text8->GetXaxis()->SetRange(x1, x2);
    text8->GetYaxis()->SetRange(y1, y2);
    text8->Draw("text,same");
    return;
  }

  if( o.name.find( "summary" ) < o.name.size() ) {
    int x1 = text6->GetXaxis()->FindBin(obj->GetXaxis()->GetXmin());
    int x2 = text6->GetXaxis()->FindBin(obj->GetXaxis()->GetXmax());
    int y1 = text6->GetYaxis()->FindBin(obj->GetYaxis()->GetXmin());
    int y2 = text6->GetYaxis()->FindBin(obj->GetYaxis()->GetXmax());
    text6->GetXaxis()->SetRange(x1, x2);
    text6->GetYaxis()->SetRange(y1, y2);
    text6->Draw("text,same");
    return;
  }

}

