/*!
  \file SiPixelRenderPlugin
  \brief Display Plugin for Pixel DQM Histograms
  \author P.Merkel
  \version $Revision: 1.24 $
  \date $Date: 2009/12/01 12:29:46 $
*/

#include "VisMonitoring/DQMServer/interface/DQMRenderPlugin.h"
#include "utils.h"

#include "TProfile2D.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TGraphPolar.h"
#include "TColor.h"
#include "TText.h"
#include "TLine.h"
#include <cassert>

#define PI_12 0.261799
#define PI    3.141592
#define PI_2  1.570796

class SiPixelRenderPlugin : public DQMRenderPlugin
{
public:
  virtual bool applies( const VisDQMObject &o, const VisDQMImgInfo & )
    {
      if( o.name.find( "Pixel/" ) != std::string::npos )
        return true;

      return false;
    }

  virtual void preDraw( TCanvas *c, const VisDQMObject &o, const VisDQMImgInfo &, VisDQMRenderInfo & )
    {
      c->cd();
      if( dynamic_cast<TH2*>( o.object ) )
      {
        preDrawTH2( c, o );
      }
      else if( dynamic_cast<TH1*>( o.object ) )
      {
        preDrawTH1( c, o );
      }
    }

  virtual void postDraw( TCanvas *c, const VisDQMObject &o, const VisDQMImgInfo & )
    {
      c->cd();

      if( dynamic_cast<TH2*>( o.object ) )
      {
        postDrawTH2( c, o );
      }
      else if( dynamic_cast<TH1*>( o.object ) )
      {
        postDrawTH1( c, o );
      }
    }

private:
/*
void paletteGraph(Double_t min, Double_t max){
	gStyle->SetPalette(1);
	gStyle->SetOptStat(0);
	gStyle->SetNumberContours(255);
	TCanvas* pal = new TCanvas("Palette", "Palette", 150, 600);
	TH2F* palet = new TH2F("Palette", "Palette", 1, 0, 1, 255, min, max);
	Double_t dc = (max - min) / 255;
	for(Int_t i = 0; i < 255; i++)
		palet->SetBinContent(1, i + 1, min + i * dc);
	palet->Draw("COL");
	gPad->SetLeftMargin(0.5);
	palet->SetLabelSize(0.2,"Y");
	palet->SetLabelOffset(0.05, "Y");
	palet->SetLabelOffset(99, "X");
	palet->SetTitle("");
	palet->SetTickLength(0,"X");
	palet->SetLabelSize(0,"X");
	gPad->SetGrid(false,false);
}
*/
Int_t getColor(TH2F* fH, Int_t i, Int_t j, Double_t zmin, Double_t zmax){
	Int_t ncolors  = gStyle->GetNumberOfColors();
	
	Double_t zc;
    	Double_t dz = zmax - zmin;
	if(dz <= 0)
		return -1;
  
	Int_t ndiv = fH->GetContour();
	if (ndiv == 0 ) {
		ndiv = gStyle->GetNumberContours();
		fH->SetContour(ndiv);
	}

	Int_t ndivz  = TMath::Abs(ndiv);

	if (fH->TestBit(TH1::kUserContour) == 0)
		fH->SetContour(ndiv);

	Double_t scale = ndivz/dz;
	Int_t color;
	const Double_t z = fH->GetBinContent(i,j);
	
	if (z < zmin || z==0.)
		return -1;	
	if (fH->TestBit(TH1::kUserContour)) {
		zc = fH->GetContourLevelPad(0);
		if (z < zc)
			return -1;
		color = -1;
		for (Int_t k=0; k<ndiv; k++) {
			zc = fH->GetContourLevelPad(k);
			if (z < zc)
				continue;
			else
				color++;
		}
	} else
		color = Int_t(0.01+(z-zmin)*scale);

	Int_t theColor = Int_t((color+0.99)*Double_t(ncolors)/Double_t(ndivz));
	if (z >= zmax)
		theColor = ncolors-1;
	return theColor;
}

void polarGraph(TCanvas *c, TH2F* map){
	TGraphPolar *graphs[24][7];
	TColor *colors = new TColor();

	colors->SetPalette(1, 0);
	gStyle->SetPalette(1);
	gStyle->SetNumberContours(255);

		for(Int_t i = 0 ; i < 24; i++)
			for(Int_t j = 0 ; j < 7; j++)
			{
				graphs[i][j] = makeSlice(j + 1., j + 1.9999, PI_2-(i+1)*PI_12, PI_2-(i)*PI_12,
						colors->GetColorPalette(getColor(map,
i + 1, j + 1, map->GetMinimum(), map->GetMaximum())), 5,
						map->GetName());
				graphs[i][j]->Draw("F same");
			}
		c->Update();
		graphs[0][0]->GetPolargram()->SetToRadian();
		graphs[0][0]->SetMaximum(8);
		graphs[0][0]->GetPolargram()->SetNdivPolar(24);
		graphs[0][0]->GetPolargram()->SetNdivRadial(8);
	
		Int_t cnt(7);

		for(Int_t i = 0 ; i < 6 ; i++)
			graphs[0][0]->GetPolargram()->SetPolarLabel(i,Form("%d", cnt--));
		for(Int_t i = 6 ; i < 18 ; i++)
			graphs[0][0]->GetPolargram()->SetPolarLabel(i,Form("%d", i-5));
		cnt = 12;
		for(Int_t i = 19 ; i < 24 ; i++)
			graphs[0][0]->GetPolargram()->SetPolarLabel(i,Form("%d", cnt--));

}


TGraphPolar* makeSlice(Double_t rA, Double_t rB, Double_t phiA, Double_t phiB, Int_t
color, Int_t npts, std::string title){
	Double_t *r = new Double_t[2*npts+3];
	Double_t *th = new Double_t[2*npts+3];
	
	r[0] = rA;
	th[0] = phiA;
	for(Int_t i = 0 ; i <= npts ; i++){
		r[i + 1] = rB;
		th[i + 1] = i * (phiB - phiA) / npts + phiA;
		r[npts+2+i] = rA;
		th[npts+2+i] = phiB - i * (phiB - phiA) / npts;
	}	
	
	
	TGraphPolar *grP = new TGraphPolar(2*npts+3,th,r);
	
	grP->SetTitle(title.c_str());
	grP->SetFillColor(color);
	grP->SetLineColor(color);
	grP->SetLineWidth(2);
	return grP;
}

  void preDrawTH2( TCanvas *, const VisDQMObject &o )
    {
      TH2* obj = dynamic_cast<TH2*>( o.object );
      assert( obj );

      // This applies to all
      gStyle->SetCanvasBorderMode( 0 );
      gStyle->SetPadBorderMode( 0 );
      gStyle->SetPadBorderSize( 0 );
      //    (data->pad)->SetLogy( 0 );;
      gStyle->SetOptStat( 0 );
      obj->SetStats( kFALSE );
      if(o.name.find( "sizeYvsEta" ) != std::string::npos){
        gStyle->SetOptStat( 1111111 );
	obj->SetStats( kTRUE ); gPad->SetLogz(1);
      }

      if( o.name.find( "reportSummaryMap" ) == std::string::npos){
        TAxis* xa = obj->GetXaxis();
        TAxis* ya = obj->GetYaxis();
        xa->SetTitleOffset(0.7);
        xa->SetTitleSize(0.065);
        xa->SetLabelSize(0.065);
        ya->SetTitleOffset(0.75);
        ya->SetTitleSize(0.065);
        ya->SetLabelSize(0.065);
      }

      if( o.name.find( "endcapOccupancyMap" ) != std::string::npos ) obj->SetTitle("Endcap Digi Occupancy Map");
      if( o.name.find( "hitmap" ) != std::string::npos  ||
          o.name.find( "Occupancy" ) != std::string::npos ||
	  o.name.find( "position_siPixelClusters" ) != std::string::npos ||
	  (o.name.find( "TRKMAP" ) != std::string::npos && o.name.find( "Layer" ) != std::string::npos) ||
	  o.name.find( "FedChLErrArray" ) != std::string::npos ||
	  o.name.find( "FedChNErrArray" ) != std::string::npos ||
	  o.name.find( "FedETypeNErrArray" ) != std::string::npos || 
	  o.name.find( "sizeYvsEta" ) != std::string::npos) 
      {
        gStyle->SetPalette(1);
        obj->SetOption("colz");
        return;
      }
      
      
      
      TH2F* obj2 = dynamic_cast<TH2F*>( o.object );

      if( o.name.find( "reportSummaryMap" ) != std::string::npos )
      {
        gPad->SetGrid();
        if(obj->GetNbinsX()==7) gPad->SetLeftMargin(0.3);
        dqm::utils::reportSummaryMapPalette(obj2);
        return;
      }
    }

  void preDrawTH1( TCanvas *, const VisDQMObject &o )
    {
      TH1* obj = dynamic_cast<TH1*>( o.object );
      assert( obj );

      // This applies to all
      gStyle->SetOptStat(111);
//       if ( obj->GetMaximum(1.e5) > 0. )
//       {
//         gPad->SetLogy(1);
//       }
//       else
//       {
//         gPad->SetLogy(0);
//       }
      
      TAxis* xa = obj->GetXaxis();
      TAxis* ya = obj->GetYaxis();
      xa->SetTitleOffset(0.7);
      xa->SetTitleSize(0.065);
      xa->SetLabelSize(0.065);
      ya->SetTitleOffset(0.75);
      ya->SetTitleSize(0.065);
      ya->SetLabelSize(0.065);

      if( o.name.find( "adcCOMB" ) != std::string::npos && obj->GetEntries() > 0. ) gPad->SetLogy(1);
      if( o.name.find( "chargeCOMB" ) != std::string::npos && obj->GetEntries() > 0. ){ obj->GetXaxis()->SetRange(1,51); gPad->SetLogy(1);}
      if( o.name.find( "OnTrack/charge_siPixelClusters" ) != std::string::npos ){ obj->GetXaxis()->SetRange(1,51); }
      if( o.name.find( "OffTrack/charge_siPixelClusters" ) != std::string::npos ){ obj->GetXaxis()->SetRange(1,51); }
      if( o.name.find( "barrelEventRate" ) != std::string::npos && obj->GetEntries() > 0. ) {gPad->SetLogx(1); gPad->SetLogy(1);}
      if( o.name.find( "endcapEventRate" ) != std::string::npos && obj->GetEntries() > 0. ) {gPad->SetLogx(1); gPad->SetLogy(1);}
      if( o.name.find( "ALLMODS_chargeCOMB" ) != std::string::npos ) obj->GetXaxis()->SetRange(1,51);
      if( o.name.find( "FEDEntries" ) != std::string::npos ) gStyle->SetOptStat(0);
//       if( o.name.find( "size_siPixelClusters" ) != std::string::npos && obj->GetEntries() > 0. ) gPad->SetLogx(1);
      if( o.name.find( "OnTrack" ) != std::string::npos && o.name.find( "charge" ) != std::string::npos ) obj->SetTitle("ClusterCharge_OnTrack");
      if( o.name.find( "OnTrack" ) != std::string::npos && o.name.find( "size" ) != std::string::npos ) obj->SetTitle("ClusterSize_OnTrack");
      if( o.name.find( "OffTrack" ) != std::string::npos && o.name.find( "charge" ) != std::string::npos ) obj->SetTitle("ClusterCharge_OffTrack");
      if( o.name.find( "OffTrack" ) != std::string::npos && o.name.find( "size" ) != std::string::npos ) obj->SetTitle("ClusterSize_OffTrack");
   
     // prettify for shifters:
//       if( o.name.find( "SUMDIG_ndigis_" ) != std::string::npos || 
//           o.name.find( "SUMCLU_nclusters_" ) != std::string::npos || 
// 	  o.name.find( "SUMCLU_size_" ) != std::string::npos ) gPad->SetLogy(0);
//       if( o.name.find( "SUMOFF_ndigis_" ) != std::string::npos || 
//           o.name.find( "SUMOFF_nclusters_" ) != std::string::npos || 
// 	  o.name.find( "SUMOFF_size_" ) != std::string::npos ) gPad->SetLogy(0);
    }

  void postDrawTH1( TCanvas *, const VisDQMObject &o )
    {
      TText tt;
      tt.SetTextSize(0.12);
      if (o.flags == 0) return;
      else
      {
        /*    if (o.flags & DQMNet::DQM_PROP_REPORT_ERROR)
              {
              tt.SetTextColor(2);
              tt.DrawTextNDC(0.5, 0.5, "Error");
              }
              else if (o.flags & DQMNet::DQM_PROP_REPORT_WARN)
              {
              tt.SetTextColor(5);
              tt.DrawTextNDC(0.5, 0.5, "Warning");
              }
              else if (o.flags & DQMNet::DQM_PROP_REPORT_OTHER)
              {
              tt.SetTextColor(1);
              tt.DrawTextNDC(0.5, 0.5, "Other ");
              }
              else
              {
              tt.SetTextColor(3);
              tt.DrawTextNDC(0.5, 0.5, "Ok ");
              }
        */
      }
      
      TH1* obj = dynamic_cast<TH1*>( o.object );
      assert( obj );
      
      if( o.name.find( "SUMOFF_adc_Barrel" ) != std::string::npos ){ 
        TLine tl1; tl1.SetLineColor(4); tl1.DrawLine(0.,75.,192.,75.); 
        TLine tl2; tl2.SetLineColor(4); tl2.DrawLine(0.,100.,192.,100.); 
      }
      else if( o.name.find( "SUMDIG_adc_Barrel" ) != std::string::npos ){ 
        TLine tl1; tl1.SetLineColor(4); tl1.DrawLine(0.,60.,768.,60.); 
        TLine tl2; tl2.SetLineColor(4); tl2.DrawLine(0.,140.,768.,140.); 
      }
      else if( o.name.find( "SUMOFF_adc_Endcap" ) != std::string::npos ){ 
        TLine tl1; tl1.SetLineColor(4); tl1.DrawLine(0.,100.,96.,100.); 
        TLine tl2; tl2.SetLineColor(4); tl2.DrawLine(0.,130.,96.,130.); 
      }
      else if( o.name.find( "SUMDIG_adc_Endcap" ) != std::string::npos ){ 
        TLine tl1; tl1.SetLineColor(4); tl1.DrawLine(0.,90.,672.,90.); 
        TLine tl2; tl2.SetLineColor(4); tl2.DrawLine(0.,160.,672.,160.); 
      }
      else if( o.name.find( "SUMOFF_ndigis_Barrel" ) != std::string::npos ){ 
        TLine tl; tl.SetLineColor(4); tl.DrawLine(0.,7.,192.,7.); 
        TLine t2; t2.SetLineColor(4); t2.DrawLine(0.,50.,192.,50.); 
      }
      else if( o.name.find( "SUMDIG_ndigis_Barrel" ) != std::string::npos ){ 
        TLine tl; tl.SetLineColor(4); tl.DrawLine(0.,1.,768.,1.); 
        TLine t2; t2.SetLineColor(4); t2.DrawLine(0.,18.,768.,18.); 
      }
      else if( o.name.find( "SUMOFF_ndigis_Endcap" ) != std::string::npos ){ 
        obj->SetMaximum(8.);
        TLine tl; tl.SetLineColor(4); tl.DrawLine(0.,2.,96.,2.); 
        TLine t2; t2.SetLineColor(4); t2.DrawLine(0.,6.,96.,6.); 
      }
      else if( o.name.find( "SUMDIG_ndigis_Endcap" ) != std::string::npos ){ 
        TLine tl; tl.SetLineColor(4); tl.DrawLine(0.,1.,672.,1.); 
        TLine t2; t2.SetLineColor(4); t2.DrawLine(0.,7.,672.,7.); 
      }
      else if( o.name.find( "SUMOFF_charge_Barrel" ) != std::string::npos ){ 
        TLine tl1; tl1.SetLineColor(4); tl1.DrawLine(0.,40.,192.,40.); 
        TLine tl2; tl2.SetLineColor(4); tl2.DrawLine(0.,150.,192.,150.); 
      }
      else if( o.name.find( "SUMOFF_size_Barrel" ) != std::string::npos ){ 
        TLine tl1; tl1.SetLineColor(4); tl1.DrawLine(0.,4.,192.,4.); 
        TLine tl2; tl2.SetLineColor(4); tl2.DrawLine(0.,11.,192.,11.); 
      }
      else if( o.name.find( "SUMCLU_charge_Barrel" ) != std::string::npos ){ 
        TLine tl1; tl1.SetLineColor(4); tl1.DrawLine(0., 30.,768., 30.); 
        TLine tl2; tl2.SetLineColor(4); tl2.DrawLine(0.,110.,768.,110.); 
      }
      else if( o.name.find( "SUMCLU_nclusters_Barrel" ) != std::string::npos ){ 
        TLine tl1; tl1.SetLineColor(4); tl1.DrawLine(0., 0.5,768., 0.5); 
        TLine tl2; tl2.SetLineColor(4); tl2.DrawLine(0.,4.,768.,4.); 
      }
      else if( o.name.find( "SUMCLU_size_Barrel" ) != std::string::npos ){ 
        TLine tl1; tl1.SetLineColor(4); tl1.DrawLine(0.,1.5,768.,1.5); 
        TLine tl2; tl2.SetLineColor(4); tl2.DrawLine(0.,11.,768.,11.); 
      }
      else if( o.name.find( "SUMOFF_charge_Endcap" ) != std::string::npos ){ 
        TLine tl1; tl1.SetLineColor(4); tl1.DrawLine(0., 7.,96., 7.); 
        TLine tl2; tl2.SetLineColor(4); tl2.DrawLine(0.,50.,96.,50.); 
      }
      else if( o.name.find( "SUMOFF_size_Endcap" ) != std::string::npos ){ 
        TLine tl1; tl1.SetLineColor(4); tl1.DrawLine(0.,1.7,96.,1.7); 
        TLine tl2; tl2.SetLineColor(4); tl2.DrawLine(0.,2.5,96.,2.5); 
      }
      else if( o.name.find( "SUMCLU_charge_Endcap" ) != std::string::npos ){ 
        TLine tl1; tl1.SetLineColor(4); tl1.DrawLine(0., 20.,672., 20.); 
        TLine tl2; tl2.SetLineColor(4); tl2.DrawLine(0.,50.,672.,50.); 
      }
      else if( o.name.find( "SUMCLU_nclusters_Endcap" ) != std::string::npos ){ 
        TLine tl1; tl1.SetLineColor(4); tl1.DrawLine(0., 0.5,672., 0.5); 
        TLine tl2; tl2.SetLineColor(4); tl2.DrawLine(0.,3.,672.,3.); 
      }
      else if( o.name.find( "SUMCLU_size_Endcap" ) != std::string::npos ){ 
        TLine tl1; tl1.SetLineColor(4); tl1.DrawLine(0., 1.,672., 1.); 
        TLine tl2; tl2.SetLineColor(4); tl2.DrawLine(0.,4.,672.,4.); 
      }
      else if( o.name.find( "OnTrack/charge_siPixelClusters" ) != std::string::npos ){ 
        Int_t ibin = obj->GetMaximumBin(); 
        Double_t val = obj->GetBinContent(ibin); 
        TLine tl; tl.SetLineColor(4); tl.DrawLine(21.,0.,21.,val); 
      }
//       else if( o.name.find( "OnTrack/size_siPixelClusters" ) != std::string::npos ||
//                o.name.find( "OffTrack/size_siPixelClusters" ) != std::string::npos ){ 
//         Int_t ibin = obj->GetMaximumBin(); 
//         Double_t val = obj->GetBinContent(ibin); 
//         TLine tl; tl.SetLineColor(4); tl.DrawLine(10.,0.,10.,val); 
//       }
      
    }

  void postDrawTH2( TCanvas *c, const VisDQMObject &o )
{

      TH2* obj = dynamic_cast<TH2*>( o.object );
      assert( obj );			  	
      if( o.name.find( "TRKMAP" ) != std::string::npos && o.name.find( "Disc" ) != std::string::npos ){
        c->Clear();
	polarGraph(c, dynamic_cast<TH2F*>(obj));
      }

      if( o.name.find( "reportSummaryMap" ) != std::string::npos )
      {
	if(obj->GetNbinsX()==40){
	  TLine tl1; tl1.SetLineWidth(1); tl1.SetLineColor(1); tl1.DrawLine(32.,27.,32.,37.); 
	  TLine tl2; tl2.SetLineWidth(1); tl2.SetLineColor(1); tl2.DrawLine(32.,27.,40.,27.); 
	}  
      }

}

};

static SiPixelRenderPlugin instance;
