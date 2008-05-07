// $Id: DTRenderPlugin.cc,v 1.10 2008/03/11 11:43:37 battilan Exp $

/*!
  \file EBRenderPlugin
  \brief Display Plugin for Quality Histograms
  \author G. Masetti
  \version $Revision: 1.10 $
  \date $Date: 2008/03/11 11:43:37 $
*/

#include "TProfile2D.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TColor.h"
#include <cassert>

//#include "DQM/EcalCommon/interface/ColorPalette.h"
//#include <DQM/EcalCommon/interface/Numbers.h>
#include "DQM/RenderPlugins/src/DTRenderPlugin.h"


bool DTRenderPlugin::applies( const DQMNet::CoreObject &o, const VisDQMImgInfo &i ) {
 
  if( o.name.find( "/DT/" ) < o.name.size() ) {
    return true;
  } 

  return false;

}

void DTRenderPlugin::preDraw( TCanvas *c, const DQMNet::CoreObject &o, const VisDQMImgInfo &i, VisDQMRenderInfo &r ) {

#ifdef DEBUG
  std::cout << "DTRenderPlugin:preDraw " << o.name << std::endl;
#endif

  c->cd();

//  gPad->SetFrameFillColor(10);
//  if (o.error) gPad->SetFillColor(2);
//  if (o.warning) gPad->SetFillColor(5);
//  if (o.other) gPad->SetFillColor(16);

  if( dynamic_cast<TProfile2D*>( o.object ) ) {
    preDrawTProfile2D( c, o );
  }
  else if( dynamic_cast<TProfile*>( o.object ) ) {
    preDrawTProfile( c, o );
  }
  else if( dynamic_cast<TH2*>( o.object ) ) {
    preDrawTH2( c, o );
  }
  else if( dynamic_cast<TH1*>( o.object ) ) {
    preDrawTH1( c, o );
  }

#ifdef DEBUG
  std::cout << "done" << std::endl;
#endif

}

void DTRenderPlugin::preDrawTProfile2D( TCanvas *c, const DQMNet::CoreObject &o ) {

  return;

}

void DTRenderPlugin::preDrawTProfile( TCanvas *c, const DQMNet::CoreObject &o ) {

  return;

}

void DTRenderPlugin::preDrawTH2( TCanvas *c, const DQMNet::CoreObject &o ) {

  TH2* obj = dynamic_cast<TH2*>( o.object );

  assert( obj );

  // This applies to all
  gStyle->SetCanvasBorderMode( 0 );
  gStyle->SetPadBorderMode( 0 );
  gStyle->SetPadBorderSize( 0 );
  //    (data->pad)->SetLogy( 0 );
  gStyle->SetOptStat( 0 );
  gStyle->SetPalette( 1 );
  obj->SetStats( kFALSE );
  
  //obj->SetOption( "box" );
  gStyle->SetPalette(1);
  obj->SetOption( "colz" );
  gPad->SetLogz(0);

  if( o.name.find( "FED770_EventLenght" ) < o.name.size() ) {
    gStyle->SetOptStat( 1111111 );
    gStyle->SetPalette( 1 );
    obj->SetStats( kTRUE );
    return;
  }

  if( o.name.find( "FED770_ROS" ) < o.name.size() ) {
    obj->GetXaxis()->LabelsOption("d");
    obj->GetXaxis()->SetLabelSize(0.05);
    obj->GetYaxis()->SetLabelSize(0.05);
    return;
  }

  if( o.name.find( "SCTriggerBX" ) < o.name.size() ) {
    obj->GetYaxis()->SetLabelSize(0.1);
    obj->GetXaxis()->SetTitle("Trigger BX");
    obj->GetYaxis()->SetRangeUser(0.,40.);
    return;
  }
  
  if( o.name.find( "SCTriggerQuality" ) < o.name.size() ) {
    obj->GetXaxis()->LabelsOption("h");
    obj->GetXaxis()->SetLabelSize(0.1);
    obj->GetYaxis()->SetLabelSize(0.1);
    return;
  }

  if (o.name.find( "TrigEff" ) < o.name.size() ) {
    obj->SetMinimum(0.);
    obj->SetMaximum(1.);
    return;
  }

  if (o.name.find( "BestQual" ) < o.name.size() ) {
    return;
  }

  if (o.name.find( "QualvsPhi" ) <o.name.size() ) {
    obj->GetYaxis()->SetRangeUser(-10.,30.);
    return;
  }   

  if (o.name.find( "QualvsPhi" ) <o.name.size() ||
      o.name.find( "QualDDUvsQualDCC" ) <o.name.size() ||
      o.name.find( "PositionvsQual" ) <o.name.size() ||
      o.name.find( "PosvsAngle" ) <o.name.size() ||
      o.name.find( "PhitkvsPhitrig" ) <o.name.size() ||
      o.name.find( "PhibtkvsPhibtrig" ) <o.name.size() ||
      o.name.find( "HitstkvsQualtrig" ) <o.name.size() ||
      o.name.find( "Flag1stvsQual" ) <o.name.size() ) {
    obj->SetOption( "box" );
    return;
  }

  return;


}

void DTRenderPlugin::preDrawTH1( TCanvas *c, const DQMNet::CoreObject &o ) {

  TH1* obj = dynamic_cast<TH1*>( o.object );

  assert( obj );

  // This applies to all
  gStyle->SetCanvasBorderMode( 0 );
  gStyle->SetPadBorderMode( 0 );
  gStyle->SetPadBorderSize( 0 );
  //      (data->pad)->SetLogy( 1 );
  gStyle->SetOptStat( 1111111 );
  obj->SetStats( kFALSE );

//   if ( obj->GetMaximum(1.e15) > 0. ) {
//     gPad->SetLogy(1);
//   } else {
//    gPad->SetLogy(0);
//   }

  if( o.name.find( "FED770TTSValues_Percent" ) < o.name.size() ) {
    gPad->SetLogy( 1 );
    return;
  }

  if( o.name.find( "FED770_ROSList" ) < o.name.size() ) {
    obj->GetXaxis()->SetTitle("# of ROS");
    return;
  }

  if( o.name.find( "hResDist" ) < o.name.size() ||
      o.name.find( "MeanTest" ) < o.name.size() ||
      o.name.find( "SigmaTest" ) < o.name.size() ||
      o.name.find( "xEfficiency" ) < o.name.size() ||
      o.name.find( "yEfficiency" ) < o.name.size() ||
      o.name.find( "Efficiency_" ) < o.name.size() ||
      o.name.find( "OccupancyDiff_" ) < o.name.size() ||
      o.name.find( "tTrigTest" ) < o.name.size() ) {
    
    TAttLine *line = dynamic_cast<TAttLine *> (o.object);
    assert (line);
    
    if (line) {
      if (o.flags & DQMNet::DQM_FLAG_REPORT_ERROR) {
	line->SetLineColor(TColor::GetColor("#CC0000"));
	//	  std::cout << name << " has error" << std::endl;
      }
      else if (o.flags & DQMNet::DQM_FLAG_REPORT_WARNING) {
	line->SetLineColor(TColor::GetColor("#993300"));
	//	  std::cout << name << " has worning" << std::endl;
      }
      else if (o.flags & DQMNet::DQM_FLAG_REPORT_OTHER) { 
	line->SetLineColor(TColor::GetColor("#FFCC00"));
	//	  std::cout << name << " has other report" << std::endl;
      }
      else {
	line->SetLineColor(TColor::GetColor("#000000"));
	//	  std::cout << name << " has nothing" << std::endl; 
      }  
    }   

    if ( o.name.find( "tTrigTest" ) < o.name.size() ) {
      obj->GetXaxis()->SetBinLabel(1,"SL1");
      obj->GetXaxis()->SetBinLabel(2,"SL2");
      obj->GetXaxis()->SetBinLabel(3,"SL3");
    }

    return;
  }
  
  if( o.name.find( "CorrectBX" ) < o.name.size() ) {
   gPad->SetLogy(0);
   obj->GetYaxis()->SetRangeUser(-10.,30.);
   obj->GetXaxis()->LabelsOption("hd");
   return;
  }

  if( o.name.find( "TrigEff" ) < o.name.size()  ||
      o.name.find( "CorrFrac" ) < o.name.size() ||
      o.name.find( "2ndFrac" ) < o.name.size()  ||
      o.name.find( "HFrac" ) < o.name.size() ) {
    gPad->SetLogy(0);
    obj->GetYaxis()->SetRangeUser(0.,1.1);
    obj->GetXaxis()->LabelsOption("hd");
    return;
  } 

  return;

}

void DTRenderPlugin::postDraw( TCanvas *c, const DQMNet::CoreObject &o, const VisDQMImgInfo &i ) {

#ifdef DEBUG
  std::cout << "DTRenderPlugin:postDraw " << o.name << std::endl;
#endif

  c->cd();

  if( dynamic_cast<TProfile2D*>( o.object ) ) {
    postDrawTProfile2D( c, o );
  }
  else if( dynamic_cast<TProfile*>( o.object ) ) {
    postDrawTProfile( c, o );
  }
  else if( dynamic_cast<TH2*>( o.object ) ) {
    postDrawTH2( c, o );
  }
  else if( dynamic_cast<TH1*>( o.object ) ) {
    postDrawTH2( c, o );
  }


#ifdef DEBUG
  std::cout << "done" << std::endl;
#endif

}

void DTRenderPlugin::postDrawTProfile2D( TCanvas *c, const DQMNet::CoreObject &o ) {

  return;

}

void DTRenderPlugin::postDrawTProfile( TCanvas *c, const DQMNet::CoreObject &o ) {

  return;

}

void DTRenderPlugin::postDrawTH2( TCanvas *c, const DQMNet::CoreObject &o ) {

  return;

}

void DTRenderPlugin::postDrawTH1( TCanvas *c, const DQMNet::CoreObject &o ) {

  return;

}

