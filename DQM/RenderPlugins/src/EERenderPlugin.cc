// $Id: EERenderPlugin.cc,v 1.54 2008/01/20 17:11:56 dellaric Exp $

/*!
  \file EERenderPlugin
  \brief Display Plugin for Quality Histograms
  \author G. Della Ricca
  \author B. Gobbo 
  \version $Revision: 1.54 $
  \date $Date: 2008/01/20 17:11:56 $
*/

#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"
#include "TProfile.h"
#include "TProfile2D.h"

#include "TStyle.h"
#include "TCanvas.h"
#include "TColor.h"

#include "TGraph.h"
#include "TLine.h"
#include "TROOT.h"

#include <iostream>

#include "DQM/EcalCommon/interface/ColorPalette.h"
#include "DQM/EcalCommon/interface/Numbers.h"

#include "EERenderPlugin.h"

static bool  first = true;

void EERenderPlugin::initialise( int argc, char **argv ) {

  if( ! first ) return;

  first = false;

  for( int i=0; i<6; i++ ) {
    TColor* color = gROOT->GetColor( 301+i );
    if ( ! color ) color = new TColor( 301+i, 0, 0, 0, "");
    color->SetRGB( ecdqm::rgb[i][0], ecdqm::rgb[i][1], ecdqm::rgb[i][2] );
  }

  for( int i=0; i<10; i++ ) {
    TColor* color = gROOT->GetColor( 401+i );
    if ( ! color ) color = new TColor( 401+i, 0, 0, 0, "");
    color->SetRGB( ecdqm::rgb2[i][0], ecdqm::rgb2[i][1], ecdqm::rgb2[i][2] );
  }
  
  for( int i=0; i<10; i++ ) {
    TColor* color = gROOT->GetColor( 501+i );
    if ( ! color ) color = new TColor( 501+i, 0, 0, 0, "");
    color->SetRGB( ecdqm::rgb2[i][1], 0, 0 );
  }

  for( short i=0; i<6; i++ ) pCol3[i]  = i+301;
  for( short i=0; i<10; i++ ) pCol4[i] = i+401;
  for( short i=0; i<10; i++ ) pCol5[i] = i+501;

  text1 = new TH2S( "ee_text1", "text1", 100, -2., 98., 100, -2., 98.);
  text3 = new TH2C( "ee_text3", "text3", 10, 0,  10,  5,   0,  5 );
  text4 = new TH2C( "ee_text4", "text4",  2, 0,   2,  1,   0,  1 );
  text6 = new TH2C( "ee_text6", "text6", 10, 0., 100., 10, 0., 100. );
  text7 = new TH2C( "ee_text7", "text7", 10, 0., 100., 10, 0., 100. );
  text8 = new TH2C( "ee_text8", "text8", 10, -150., 150., 10, -150., 150. );
  text9 = new TH2C( "ee_text9", "text9", 10, -150., 150., 10, -150., 150. );

  text1->SetMinimum(  0.01 );
  text3->SetMinimum(  0.01 );
  text4->SetMinimum(  0.01 );
  text6->SetMinimum( -9.01 );
  text7->SetMinimum( +0.01 );
  text8->SetMinimum( -9.01 );
  text9->SetMinimum( +0.01 );

  text6->SetMaximum( -0.01 );
  text7->SetMaximum( +9.01 );
  text8->SetMaximum( -0.01 );
  text9->SetMaximum( +9.01 );

  for ( short j=0; j<400; j++ ) {
    int x = 5*(1 + j%20);
    int y = 5*(1 + j/20);
    text1->SetBinContent(x, y, Numbers::inTowersEE[j]);
  }

  for ( short i=1; i<=10; i++) {
    for ( int j=1; j<=10; j++) {
      text6->SetBinContent(i, j, -10);
      text7->SetBinContent(i, j, -10);
      text8->SetBinContent(i, j, -10);
      text9->SetBinContent(i, j, -10);
    }
  }
  for( short i=0; i<2; i++ ) {
    text3->Fill( 2+i*5, 2, i+1+68 );
    text4->Fill( i, 0., i+1+68 );
  }

  text6->SetBinContent(2, 5, -7);
  text6->SetBinContent(2, 7, -8);
  text6->SetBinContent(4, 9, -9);
  text6->SetBinContent(7, 9, -1);
  text6->SetBinContent(9, 7, -2);
  text6->SetBinContent(9, 5, -3);
  text6->SetBinContent(8, 3, -4);
  text6->SetBinContent(6, 2, -5);
  text6->SetBinContent(3, 3, -6);

  text7->SetBinContent(2, 5, +3);
  text7->SetBinContent(2, 7, +2);
  text7->SetBinContent(4, 9, +1);
  text7->SetBinContent(7, 9, +9);
  text7->SetBinContent(9, 7, +8);
  text7->SetBinContent(9, 5, +7);
  text7->SetBinContent(8, 3, +6);
  text7->SetBinContent(5, 2, +5);
  text7->SetBinContent(3, 3, +4);

  text8->SetBinContent(2, 5, -7);
  text8->SetBinContent(2, 7, -8);
  text8->SetBinContent(4, 9, -9);
  text8->SetBinContent(7, 9, -1);
  text8->SetBinContent(9, 7, -2);
  text8->SetBinContent(9, 5, -3);
  text8->SetBinContent(8, 3, -4);
  text8->SetBinContent(6, 2, -5);
  text8->SetBinContent(3, 3, -6);

  text9->SetBinContent(2, 5, +3);
  text9->SetBinContent(2, 7, +2);
  text9->SetBinContent(4, 9, +1);
  text9->SetBinContent(7, 9, +9);
  text9->SetBinContent(9, 7, +8);
  text9->SetBinContent(9, 5, +7);
  text9->SetBinContent(8, 3, +6);
  text9->SetBinContent(5, 2, +5);
  text9->SetBinContent(3, 3, +4);

  text1->SetMarkerSize( 2 );
  text3->SetMarkerSize( 2 );
  text4->SetMarkerSize( 2 );
  text6->SetMarkerSize( 2 );
  text7->SetMarkerSize( 2 );
  text8->SetMarkerSize( 2 );
  text9->SetMarkerSize( 2 );

}

bool EERenderPlugin::applies( const ObjInfo &o, const ImgInfo &i ) {
 
  if( o.name.find( "EcalEndcap/EE" ) == 0 ) {
    return true;
  } 

  if( o.name.find( "EcalEndcap/EcalInfo" ) == 0 ) {
    return true;
  }

  return false;

}

void EERenderPlugin::preDraw( TCanvas *c, const ObjInfo &o, const ImgInfo &i, RenderInfo &r ) {

#ifdef DEBUG
  std::cout << "EERenderPlugin:preDraw " << o.name << std::endl;
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

  r.drawOptions = "";

#ifdef DEBUG
  std::cout << "done" << std::endl;
#endif

}

void EERenderPlugin::preDrawTProfile2D( TCanvas *c, const ObjInfo &o ) {

  TProfile2D* obj = dynamic_cast<TProfile2D*>( o.object );

  assert( obj );

  gStyle->SetPaintTextFormat();

  gStyle->SetOptStat(0);
  obj->SetStats(kFALSE);

  // Occupancy-like (10 x grays) plots
  if( o.name.find( "EECLT" ) < o.name.size() ) {
    gPad->SetGridx();
    gPad->SetGridy();
    obj->GetXaxis()->SetNdivisions(10, kFALSE);
    obj->GetYaxis()->SetNdivisions(10, kFALSE);
    obj->SetMinimum(0.0);
    gStyle->SetPalette(10, pCol4);
    obj->SetOption("colz");
    gStyle->SetPaintTextFormat("+g");
    return;
  }

  int nbx = obj->GetNbinsX();
  int nby = obj->GetNbinsY();

  if( nbx == 50 && nby == 50 ) {
    gPad->SetGridx();
    gPad->SetGridy();
    obj->GetXaxis()->SetNdivisions(10);
    obj->GetYaxis()->SetNdivisions(10);
  }

  if( nbx == 100 && nby == 100 ) {
    gPad->SetGridx();
    gPad->SetGridy();
    obj->GetXaxis()->SetNdivisions(10);
    obj->GetYaxis()->SetNdivisions(10);
  }

  // Occupancy-like (10 x grays) plots
  obj->SetMinimum(0.0);
  gStyle->SetPalette(10, pCol4);
  obj->SetOption("colz");
  return;

}

void EERenderPlugin::preDrawTProfile( TCanvas *c, const ObjInfo &o ) {

  TProfile* obj = dynamic_cast<TProfile*>( o.object );

  assert( obj );

  gStyle->SetOptStat("euomr");
  obj->SetStats(kTRUE);
  gPad->SetLogy(0);

  return;

}

void EERenderPlugin::preDrawTH3F( TCanvas *c, const ObjInfo &o ) {

  TH3F* obj = dynamic_cast<TH3F*>( o.object );

  assert( obj );

  gStyle->SetOptStat(0);
  obj->SetStats( kFALSE );
  gPad->SetLogy(0);

  return;

}

void EERenderPlugin::preDrawTH2F( TCanvas *c, const ObjInfo &o ) {

  TH2F* obj = dynamic_cast<TH2F*>( o.object );

  assert( obj );

  gStyle->SetPaintTextFormat();

  gStyle->SetOptStat(0);
  obj->SetStats( kFALSE );
  gPad->SetLogy(0);

  // Occupancy-like (10 x grays) plots
  if( o.name.find( "EECLT" ) < o.name.size() ) {
    gPad->SetGridx();
    gPad->SetGridy();
    obj->GetXaxis()->SetNdivisions(10, kFALSE);
    obj->GetYaxis()->SetNdivisions(10, kFALSE);
    obj->SetMinimum(0.0);
    gStyle->SetPalette(10, pCol4);
    obj->SetOption("colz");
    gStyle->SetPaintTextFormat("+g");
    return;
  }

  int nbx = obj->GetNbinsX();
  int nby = obj->GetNbinsY();

  if( nbx == 50 && nby == 50 ) {
    gPad->SetGridx();
    gPad->SetGridy();
    obj->GetXaxis()->SetNdivisions(10);
    obj->GetYaxis()->SetNdivisions(10);
  }

  if( nbx == 20 && nby == 20 ) {
    gPad->SetGridx();
    gPad->SetGridy();
    obj->GetXaxis()->SetNdivisions(20);
    obj->GetYaxis()->SetNdivisions(20);
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

  if( nbx == 100 && nby == 100 ) {
    gPad->SetGridx();
    gPad->SetGridy();
    obj->GetXaxis()->SetNdivisions(10);
    obj->GetYaxis()->SetNdivisions(10);
  }

  if( nbx == 90 && nby == 20 ) {
    gPad->SetGridx();
    gPad->SetGridy();
    obj->GetXaxis()->SetNdivisions(18);
    obj->GetYaxis()->SetNdivisions(2);
  }

  // Occupancy-like (10 x grays) plots
  if( o.name.find( "EEIT" ) < o.name.size() &&
      o.name.find( "quality" ) >= o.name.size() ) {
    obj->SetMinimum(0.0);
    gStyle->SetPalette(10, pCol5);
    obj->SetOption("colz");
    gStyle->SetPaintTextFormat("+g");
    return;
  }

  // Occupancy-like (10 x grays) plots
  if( o.name.find( "EETTT" ) < o.name.size() &&
      o.name.find( "quality" ) >= o.name.size() ) {
    obj->SetMinimum(0.0);
    if( o.name.find( "Error" ) >= o.name.size() ) {
      gStyle->SetPalette(10, pCol4);
    } else {
      gStyle->SetPalette(10, pCol5);
    }
    obj->SetOption("colz");
    gStyle->SetPaintTextFormat("+g");
    return;
  }

  // Occupancy-like (10 x grays) plots
  if( o.name.find( "EEOT" ) < o.name.size() ) {
    obj->SetMinimum(0.0);
    gStyle->SetPalette(10, pCol4);
    obj->SetOption("colz");
    gStyle->SetPaintTextFormat("+g");
    return;
  }

  // Occupancy-like (10 x grays) plots
  if( o.name.find( "EESFT" ) < o.name.size() &&
      o.name.find( "summary" ) >= o.name.size() ) {
    obj->SetMinimum(0.0);
    gStyle->SetPalette(10, pCol5);
    obj->SetOption("colz");
    gStyle->SetPaintTextFormat("+g");
    return;
  }

  // Occupancy-like (10 x grays) plots
  if( o.name.find( "EECT" ) < o.name.size() ) {
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
  if( o.name.find( "EEMM event" ) < o.name.size() ) {
    obj->SetMinimum(0.0);
    gStyle->SetPalette(10, pCol4);
    obj->SetOption("colz");
    gStyle->SetPaintTextFormat("+g");
    return;
  }

}

void EERenderPlugin::preDrawTH1F( TCanvas *c, const ObjInfo &o ) {

  TH1F* obj = dynamic_cast<TH1F*>( o.object );

  assert( obj );

  gStyle->SetOptStat("euomr");
  obj->SetStats(kTRUE);
  gPad->SetLogy(0);

  int nbx = obj->GetNbinsX();

  if ( obj->GetMaximum() > 0. ) gPad->SetLogy(1);

  if ( nbx == 10 ) gPad->SetLogy(0);
  if ( nbx == 850 ) gPad->SetLogy(0);

  if( o.name.find( "EVTTYPE" ) < o.name.size() ) {
   gPad->SetBottomMargin(0.4);
   obj->GetXaxis()->LabelsOption("v");
  }

  if( o.name.find( "EEMM DCC" ) < o.name.size() ) {
   gPad->SetBottomMargin(0.2);
   obj->GetXaxis()->LabelsOption("v");
  }

  if( o.name.find( "front-end status bits" ) < o.name.size() ) {
   gPad->SetBottomMargin(0.25);
   obj->GetXaxis()->LabelsOption("v");
  }

  if( o.name.find( "front-end status errors summary" ) < o.name.size() ) {
   gPad->SetBottomMargin(0.2);
   obj->GetXaxis()->LabelsOption("v");
  }

  if( o.name.find( "quality errors summary" ) < o.name.size() ) {
   gPad->SetBottomMargin(0.2);
   obj->GetXaxis()->LabelsOption("v");
  }

  if( o.name.find( "EEOT digi occupancy summary 1D" ) < o.name.size() ) {
   gPad->SetBottomMargin(0.2);
   obj->GetXaxis()->LabelsOption("v");
  }

  return;

}

void EERenderPlugin::postDraw( TCanvas *c, const ObjInfo &o, const ImgInfo &i ) {

#ifdef DEBUG
  std::cout << "EERenderPlugin:postDraw " << o.name << std::endl;
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
  else if( dynamic_cast<TH1F*>( o.object ) ) {
    preDrawTH1F( c, o );
  }

#ifdef DEBUG
  std::cout << "done" << std::endl;
#endif

}

void EERenderPlugin::postDrawTProfile2D( TCanvas *c, const ObjInfo &o ) {

  TProfile2D* obj = dynamic_cast<TProfile2D*>( o.object );

  assert( obj );

  c->SetBit(TGraph::kClipFrame);
  TLine l;
  l.SetLineWidth(1);
  for ( int i=0; i<201; i=i+1){
    if ( (Numbers::ixSectorsEE[i]!=0 || Numbers::iySectorsEE[i]!=0) && (Numbers::ixSectorsEE[i+1]!=0 || Numbers::iySectorsEE[i+1]!=0) ) {
      if( o.name.find( "EECLT" ) < o.name.size() ) {
        l.DrawLine(3.0*(Numbers::ixSectorsEE[i]-50), 3.0*(Numbers::iySectorsEE[i]-50), 3.0*(Numbers::ixSectorsEE[i+1]-50), 3.0*(Numbers::iySectorsEE[i+1]-50));
      } else {
        l.DrawLine(Numbers::ixSectorsEE[i], Numbers::iySectorsEE[i], Numbers::ixSectorsEE[i+1], Numbers::iySectorsEE[i+1]);
      }
    }
  }

  if( o.name.find( "EECLT" ) < o.name.size() ) {
    if( o.name.find( "EE -" ) < o.name.size() ) {
      int x1 = text8->GetXaxis()->FindFixBin(obj->GetXaxis()->GetXmin());
      int x2 = text8->GetXaxis()->FindFixBin(obj->GetXaxis()->GetXmax());
      int y1 = text8->GetYaxis()->FindFixBin(obj->GetYaxis()->GetXmin());
      int y2 = text8->GetYaxis()->FindFixBin(obj->GetYaxis()->GetXmax());
      text8->GetXaxis()->SetRange(x1, x2);
      text8->GetYaxis()->SetRange(y1, y2);
      text8->Draw("text,same");
    }
    if( o.name.find( "EE +" ) < o.name.size() ) {
      int x1 = text9->GetXaxis()->FindFixBin(obj->GetXaxis()->GetXmin());
      int x2 = text9->GetXaxis()->FindFixBin(obj->GetXaxis()->GetXmax());
      int y1 = text9->GetYaxis()->FindFixBin(obj->GetYaxis()->GetXmin());
      int y2 = text9->GetYaxis()->FindFixBin(obj->GetYaxis()->GetXmax());
      text9->GetXaxis()->SetRange(x1, x2);
      text9->GetYaxis()->SetRange(y1, y2);
      text9->Draw("text,same");
    }
    return;
  }

  int x1 = text1->GetXaxis()->FindFixBin(obj->GetXaxis()->GetXmin());
  int x2 = text1->GetXaxis()->FindFixBin(obj->GetXaxis()->GetXmax());
  int y1 = text1->GetYaxis()->FindFixBin(obj->GetYaxis()->GetXmin());
  int y2 = text1->GetYaxis()->FindFixBin(obj->GetYaxis()->GetXmax());
  text1->GetXaxis()->SetRange(x1, x2);
  text1->GetYaxis()->SetRange(y1, y2);
  text1->Draw("text,same");
  return;

}

void EERenderPlugin::postDrawTH3F( TCanvas *c, const ObjInfo &o ) {

  TH3F* obj = dynamic_cast<TH3F*>( o.object );

  assert( obj );

  return;

}

void EERenderPlugin::postDrawTH2F( TCanvas *c, const ObjInfo &o ) {

  TH2F* obj = dynamic_cast<TH2F*>( o.object );

  assert( obj );

  c->SetBit(TGraph::kClipFrame);
  TLine l;
  l.SetLineWidth(1);
  for ( int i=0; i<201; i=i+1){
    if ( (Numbers::ixSectorsEE[i]!=0 || Numbers::iySectorsEE[i]!=0) && (Numbers::ixSectorsEE[i+1]!=0 || Numbers::iySectorsEE[i+1]!=0) ) {
      if( o.name.find( "EECLT" ) < o.name.size() ) {
        l.DrawLine(3.0*(Numbers::ixSectorsEE[i]-50), 3.0*(Numbers::iySectorsEE[i]-50), 3.0*(Numbers::ixSectorsEE[i+1]-50), 3.0*(Numbers::iySectorsEE[i+1]-50));
      } else {
        l.DrawLine(Numbers::ixSectorsEE[i], Numbers::iySectorsEE[i], Numbers::ixSectorsEE[i+1], Numbers::iySectorsEE[i+1]);
      }
    }
  }

  if( o.name.find( "EECLT" ) < o.name.size() ) {
    if( o.name.find( "EE -" ) < o.name.size() ) {
      int x1 = text8->GetXaxis()->FindFixBin(obj->GetXaxis()->GetXmin());
      int x2 = text8->GetXaxis()->FindFixBin(obj->GetXaxis()->GetXmax());
      int y1 = text8->GetYaxis()->FindFixBin(obj->GetYaxis()->GetXmin());
      int y2 = text8->GetYaxis()->FindFixBin(obj->GetYaxis()->GetXmax());
      text8->GetXaxis()->SetRange(x1, x2);
      text8->GetYaxis()->SetRange(y1, y2);    
      text8->Draw("text,same");
    }
    if( o.name.find( "EE +" ) < o.name.size() ) {
      int x1 = text9->GetXaxis()->FindFixBin(obj->GetXaxis()->GetXmin());
      int x2 = text9->GetXaxis()->FindFixBin(obj->GetXaxis()->GetXmax());
      int y1 = text9->GetYaxis()->FindFixBin(obj->GetYaxis()->GetXmin());
      int y2 = text9->GetYaxis()->FindFixBin(obj->GetYaxis()->GetXmax());
      text9->GetXaxis()->SetRange(x1, x2);
      text9->GetYaxis()->SetRange(y1, y2);    
      text9->Draw("text,same");
    }
    return;
  }

  if( o.name.find( "EECLT" ) < o.name.size() ) {
    if( o.name.find( "EE -" ) < o.name.size() ) {
      int x1 = text6->GetXaxis()->FindFixBin(obj->GetXaxis()->GetXmin());
      int x2 = text6->GetXaxis()->FindFixBin(obj->GetXaxis()->GetXmax());
      int y1 = text6->GetYaxis()->FindFixBin(obj->GetYaxis()->GetXmin());
      int y2 = text6->GetYaxis()->FindFixBin(obj->GetYaxis()->GetXmax());
      text6->GetXaxis()->SetRange(x1, x2);
      text6->GetYaxis()->SetRange(y1, y2);
      text6->Draw("text,same");
    }
    if( o.name.find( "EE +" ) < o.name.size() ) {
      int x1 = text7->GetXaxis()->FindFixBin(obj->GetXaxis()->GetXmin());
      int x2 = text7->GetXaxis()->FindFixBin(obj->GetXaxis()->GetXmax());
      int y1 = text7->GetYaxis()->FindFixBin(obj->GetYaxis()->GetXmin());
      int y2 = text7->GetYaxis()->FindFixBin(obj->GetYaxis()->GetXmax());
      text7->GetXaxis()->SetRange(x1, x2);
      text7->GetYaxis()->SetRange(y1, y2);
      text7->Draw("text,same");
    }
    return;
  }

  int nbx = obj->GetNbinsX();
  int nby = obj->GetNbinsY();

  if( nbx == 50 && nby == 50 ) {
    int x1 = text1->GetXaxis()->FindFixBin(obj->GetXaxis()->GetXmin());
    int x2 = text1->GetXaxis()->FindFixBin(obj->GetXaxis()->GetXmax());
    int y1 = text1->GetYaxis()->FindFixBin(obj->GetYaxis()->GetXmin());
    int y2 = text1->GetYaxis()->FindFixBin(obj->GetYaxis()->GetXmax());
    text1->GetXaxis()->SetRange(x1, x2);
    text1->GetYaxis()->SetRange(y1, y2);
    text1->Draw("text,same");
    return;
  }

  if( nbx == 10 && nby == 5 ) {
    int x1 = text3->GetXaxis()->FindFixBin(obj->GetXaxis()->GetXmin());
    int x2 = text3->GetXaxis()->FindFixBin(obj->GetXaxis()->GetXmax());
    int y1 = text3->GetYaxis()->FindFixBin(obj->GetYaxis()->GetXmin());
    int y2 = text3->GetYaxis()->FindFixBin(obj->GetYaxis()->GetXmax());
    text3->GetXaxis()->SetRange(x1, x2);
    text3->GetYaxis()->SetRange(y1, y2);
    text3->Draw("text,same");
    return;
  }

  if( nbx == 2 && nby == 1 ) {
    int x1 = text4->GetXaxis()->FindFixBin(obj->GetXaxis()->GetXmin());
    int x2 = text4->GetXaxis()->FindFixBin(obj->GetXaxis()->GetXmax());
    int y1 = text4->GetYaxis()->FindFixBin(obj->GetYaxis()->GetXmin());
    int y2 = text4->GetYaxis()->FindFixBin(obj->GetYaxis()->GetXmax());
    text4->GetXaxis()->SetRange(x1, x2);
    text4->GetYaxis()->SetRange(y1, y2);
    text4->Draw("text,same");
    return;
  }

  if( o.name.find( "summary" ) < o.name.size() ) {
    if( o.name.find( "EE -" ) < o.name.size() ) {
      int x1 = text6->GetXaxis()->FindFixBin(obj->GetXaxis()->GetXmin());
      int x2 = text6->GetXaxis()->FindFixBin(obj->GetXaxis()->GetXmax());
      int y1 = text6->GetYaxis()->FindFixBin(obj->GetYaxis()->GetXmin());
      int y2 = text6->GetYaxis()->FindFixBin(obj->GetYaxis()->GetXmax());
      text6->GetXaxis()->SetRange(x1, x2);
      text6->GetYaxis()->SetRange(y1, y2);    
      text6->Draw("text,same");
    }
    if( o.name.find( "EE +" ) < o.name.size() ) {
      int x1 = text7->GetXaxis()->FindFixBin(obj->GetXaxis()->GetXmin());
      int x2 = text7->GetXaxis()->FindFixBin(obj->GetXaxis()->GetXmax());
      int y1 = text7->GetYaxis()->FindFixBin(obj->GetYaxis()->GetXmin());
      int y2 = text7->GetYaxis()->FindFixBin(obj->GetYaxis()->GetXmax());
      text7->GetXaxis()->SetRange(x1, x2);
      text7->GetYaxis()->SetRange(y1, y2);    
      text7->Draw("text,same");
    }
    return;
  }

}

void EERenderPlugin::postDrawTH1F( TCanvas *c, const ObjInfo &o ) {

  TH1F* obj = dynamic_cast<TH1F*>( o.object );

  assert( obj );

  return;

}

