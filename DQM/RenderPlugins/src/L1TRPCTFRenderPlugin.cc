/*!
  \file L1TRPCTFRenderPlugin.cc
  \\
  \\ Code shamelessly taken from, which was shamelessly borrowed 
  \\ from J. Temple's HcalRenderPlugin.cc code, which was shamelessly 
  \\ borrowed from S. Dutta's SiStripRenderPlugin.cc  code, G. Della Ricca 
  \\ and B. Gobbo's EBRenderPlugin.cc, and other existing subdetector plugins
*/

#include "DQM/RenderPlugins/src/L1TRPCTFRenderPlugin.h" 

#include "TProfile2D.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TColor.h"
#include "TText.h"

#include <cassert>
#include "DQM/RenderPlugins/src/utils.h"

void L1TRPCTFRenderPlugin::initialise (int argc, char ** argv)
{
  // same as RenderPlugin default for now (no special action taken)
  return;
}



bool L1TRPCTFRenderPlugin::applies(const DQMNet::CoreObject &o,
				  const VisDQMImgInfo &i)
{
  // determine whether core object is an L1TEMU object
#ifdef DEBUG 
  std::cout << "L1TRPCTFRenderPlugin:applies " << o.name << std::endl; 
#endif 
  if (o.name.find( "L1TRPCTF/" ) != std::string::npos )
    {
      return true;
    }
  
  return false;
}



void L1TRPCTFRenderPlugin::preDraw (TCanvas * c,
			      const DQMNet::CoreObject &o,
			      const VisDQMImgInfo &i,
			      VisDQMRenderInfo &r)
{

#ifdef DEBUG 
  std::cout << "L1TRPCTFRenderPlugin:preDraw " << o.name << std::endl; 
#endif 
  c->cd(); 
 
  // object is TH2 histogram
  if( dynamic_cast<TH2F*>( o.object ) )
    { 
      preDrawTH2F( c, o ); 
    } 

  // object is TH1 histogram
  else if( dynamic_cast<TH1F*>( o.object ) ) 
    { 
      preDrawTH1F( c, o ); 
    } 

  else
#ifdef DEBUG 
    std::cout << "L1TRPCTFRenderPlugin:preDraw  -- Cannot identify object " << o.name << std::endl; 
#endif 

  return;
} // L1TRPCTFRenderPlugin::preDraw(...)


void L1TRPCTFRenderPlugin::postDraw (TCanvas * c,
			       const DQMNet::CoreObject & o,
			       const VisDQMImgInfo & i)
{

#ifdef DEBUG 
  std::cout << "L1TRPCTFRenderPlugin:postDraw " << o.name << std::endl; 
#endif 

  // object is TH2 histogram
  if( dynamic_cast<TH2F*>( o.object ) )
    { 
      postDrawTH2F( c, o ); 
    } 

  // object is TH1 histogram
  else if( dynamic_cast<TH1F*>( o.object ) ) 
    { 
      postDrawTH1F( c, o ); 
    } 

  else
#ifdef DEBUG 
    std::cout << "L1TRenderPlugin:postDraw  -- Cannot identify object " << o.name << std::endl; 
#endif 

  return;
} // L1TRenderPlugin::postDraw(...)



////////////////////////////////////////////////////////////////////////

//  private functions


void L1TRPCTFRenderPlugin::preDrawTH1F ( TCanvas *c, const DQMNet::CoreObject &o )
{
  // Do we want to do anything special yet with TH1F histograms?

  TH1F* obj = dynamic_cast<TH1F*>( o.object ); 
  assert (obj); // checks that object indeed exists

  // Code used in SiStripRenderPlugin -- do we want similar defaults?

  /*
    gStyle->SetOptStat(0111); 
    if ( obj->GetMaximum(1.e5) > 0. ) { 
    gPad->SetLogy(1); 
    } else { 
    gPad->SetLogy(0); 
    } 
  */

  return;

} // preDrawTH1(...)




void L1TRPCTFRenderPlugin::preDrawTH2F ( TCanvas *c, const DQMNet::CoreObject &o )
{
  

  TH2F* obj = dynamic_cast<TH2F*>( o.object ); 
  assert( obj ); 

  /*
  //put in preDrawTH2F  
  if( o.name.find( "reportSummaryMap" )  != std::string::npos) {
    obj->SetStats( kFALSE );
    dqm::utils::reportSummaryMapPalette(obj);
    obj->SetOption("colz");


    obj->SetTitle("L1TEMU Report Summary Map");
    
    obj->GetXaxis()->SetNdivisions(1,true);
    obj->GetYaxis()->SetNdivisions(12,true);
    obj->GetXaxis()->CenterLabels();
    obj->GetYaxis()->CenterLabels();
    
    gPad->SetGrid(1,1);
    return;
  }

  gStyle->SetCanvasBorderMode( 0 ); 
  gStyle->SetPadBorderMode( 0 ); 
  gStyle->SetPadBorderSize( 0 ); 
  
  // I don't think we want to set stats to 0 for Hcal
  //gStyle->SetOptStat( 0 ); 
  //obj->SetStats( kFALSE ); 

  // Use same labeling format as SiStripRenderPlugin.cc
  TAxis* xa = obj->GetXaxis(); 
  TAxis* ya = obj->GetYaxis(); 
  
  xa->SetTitleOffset(0.7); 
  xa->SetTitleSize(0.05); 
  xa->SetLabelSize(0.04); 
  
  ya->SetTitleOffset(0.7); 
  ya->SetTitleSize(0.05); 
  ya->SetLabelSize(0.04); 

  // Now the important stuff -- set 2D hist drawing option to "colz"
  gStyle->SetPalette(1);
  obj->SetOption("colz");
  */

  return;
} // preDrawTH2F(...)




void L1TRPCTFRenderPlugin::postDrawTH1F( TCanvas *c, const DQMNet::CoreObject &o )
{

  // Add error/warning text to 1-D histograms.  Do we want this at this time?
  /*
  TText tt; 
  tt.SetTextSize(0.12); 

  if (o.flags == 0) return; 

  else 
    { 
      if (o.flags & DQMNet::DQM_FLAG_REPORT_ERROR) 
	{ 
	  tt.SetTextColor(2); // error color = RED 
	  tt.DrawTextNDC(0.5, 0.5, "Error"); 
	}  // DQM_FLAG_REPORT_ERROR

      else if (o.flags & DQMNet::DQM_FLAG_REPORT_WARNING) 
	{ 
	  tt.SetTextColor(5); 
	  tt.DrawTextNDC(0.5, 0.5, "Warning"); // warning color = YELLOW 
	} // DQM_FLAG_REPORT_WARNING

      else if (o.flags & DQMNet::DQM_FLAG_REPORT_OTHER) 
	{  
	  tt.SetTextColor(1); // other color = BLACK
	  tt.DrawTextNDC(0.5, 0.5, "Other ");       
	} // DQM_FLAG_REPORT_OTHER

      else 
	{ 
	  tt.SetTextColor(3); 
	  tt.DrawTextNDC(0.5, 0.5, "Ok "); 
	} //else
    } // else (  o.flags != 0  )
  */
  return;

} // postDrawTH1F(...)




void L1TRPCTFRenderPlugin::postDrawTH2F( TCanvas *c, const DQMNet::CoreObject &o )
{
  // nothing to put here just yet
  // in the future, we can add text output based on error status,
  // or set bin range based on filled histograms, etc.  
  // Maybe add a big "OK" sign to histograms with no entries (i.e., no errors)?

  return;

} // postDrawTH2F(...)
