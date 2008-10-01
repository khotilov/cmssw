/*!
  \file HLTRenderPlugin.cc
  \\
  \\ Code shamelessly borrowed from J. Temple's HcalRenderPlugin.cc code, 
  \\ which was shamelessly borrowed from S. Dutta's SiStripRenderPlugin.cc 
  \\ code, G. Della Ricca and B. Gobbo's EBRenderPlugin.cc, and other existing
  \\ subdetector plugins
  \\ preDraw and postDraw methods now check whether histogram was a TH1
  \\ or TH2, and call a private method appropriate for the histogram type
  $Id: HLTRenderPlugin.cc,v 1.4 2008/08/28 21:50:39 wittich Exp $
  $Log: HLTRenderPlugin.cc,v $
  Revision 1.4  2008/08/28 21:50:39  wittich
  Rate histos: Also put in low range minimums in case we start in the
  middle of a run

*/

#include "DQM/RenderPlugins/src/HLTRenderPlugin.h" 

#include "TProfile2D.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TColor.h"
#include "TText.h"

#include <cassert>
#include "DQM/RenderPlugins/src/utils.h"

void
HLTRenderPlugin::initialise (int argc, char ** argv)
{
  // same as RenderPlugin default for now (no special action taken)
  return;
}



bool HLTRenderPlugin::applies(const DQMNet::CoreObject &o,
				  const VisDQMImgInfo &i)
{
  // determine whether core object is an HLT object
#ifdef DEBUG 
  std::cout << "HLTRenderPlugin:applies " << o.name << std::endl; 
#endif 
  if (o.name.find( "HLT/" ) != std::string::npos  )
    {
      return true;
    }
  
  return false;
}



void HLTRenderPlugin::preDraw (TCanvas * c,
			      const DQMNet::CoreObject &o,
			      const VisDQMImgInfo &i,
			      VisDQMRenderInfo &r)
{

#ifdef DEBUG 
  std::cout << "HLTRenderPlugin:preDraw " << o.name << std::endl; 
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
    std::cout << "HLTRenderPlugin:preDraw  -- Cannot identify object " << o.name << std::endl; 
#endif 

  return;
} // HLTRenderPlugin::preDraw(...)


void HLTRenderPlugin::postDraw (TCanvas * c,
			       const DQMNet::CoreObject & o,
			       const VisDQMImgInfo & i)
{

#ifdef DEBUG 
  std::cout << "HLTRenderPlugin:postDraw " << o.name << std::endl; 
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
    std::cout << "HLTRenderPlugin:postDraw  -- Cannot identify object " << o.name << std::endl; 
#endif 

  return;
} // HLTRenderPlugin::postDraw(...)



////////////////////////////////////////////////////////////////////////

//  private functions


void HLTRenderPlugin::preDrawTH1F ( TCanvas *c, const DQMNet::CoreObject &o )
{
  // Do we want to do anything special yet with TH1F histograms?

  TH1F* obj = dynamic_cast<TH1F*>( o.object ); 
  assert (obj); // checks that object indeed exists


  // rate histograms
  if ( o.name.find("rate_p") != std::string::npos) {
    gStyle->SetOptStat(11);
    obj->GetXaxis()->SetTitle("Luminosity Segment Number");
    obj->GetYaxis()->SetTitle("Rate (Hz)");
    int nbins = obj->GetNbinsX();
    int maxRange = nbins;
    for ( int i = nbins; i > 0; --i ) {
      if ( obj->GetBinContent(i) != 0 ) {
	maxRange = i;
	break;
      }
    }
    int minRange = 0;
    for ( int i = 0; i <= nbins; ++i ) {
      if ( obj->GetBinContent(i) != 0 ) {
	minRange = i;
	break;
      }
    }

    obj->GetXaxis()->SetRange(minRange, maxRange);
  }

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




void HLTRenderPlugin::preDrawTH2F ( TCanvas *c, const DQMNet::CoreObject &o )
{
  

  TH2F* obj = dynamic_cast<TH2F*>( o.object ); 
  assert( obj ); 
  
  //put in preDrawTH2F  
    if( o.name.find( "reportSummaryMap" )  != std::string::npos) {
   obj->SetStats( kFALSE );
    dqm::utils::reportSummaryMapPalette(obj);
    obj->SetOption("colz");
    obj->SetTitle("HLT Report Summary Map");
    obj->GetXaxis()->SetNdivisions(1,true);
    obj->GetYaxis()->SetNdivisions(5,true);
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

  return;
} // preDrawTH2F(...)




void HLTRenderPlugin::postDrawTH1F( TCanvas *c, const DQMNet::CoreObject &o )
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




void HLTRenderPlugin::postDrawTH2F( TCanvas *c, const DQMNet::CoreObject &o )
{
  // nothing to put here just yet
  // in the future, we can add text output based on error status,
  // or set bin range based on filled histograms, etc.  
  // Maybe add a big "OK" sign to histograms with no entries (i.e., no errors)?

  return;

} // postDrawTH2F(...)
	 
