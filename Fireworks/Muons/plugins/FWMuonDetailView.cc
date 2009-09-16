// -*- C++ -*-
//
// Package:     Calo
// Class  :     FWMuonDetailView
// $Id: FWMuonDetailView.cc,v 1.11 2009/09/06 23:14:35 dmytro Exp $
//

#include "TEveLegoEventHandler.h"

// ROOT includes
#include "TLatex.h"
#include "TEveCalo.h"
#include "TEveStraightLineSet.h"
#include "TEvePointSet.h"
#include "TEveScene.h"
#include "TEveViewer.h"
#include "TGLViewer.h"
#include "TEveManager.h"
#include "TCanvas.h" 
#include "TEveCaloLegoOverlay.h"
#include "TRootEmbeddedCanvas.h"

// Fireworks includes
#include "Fireworks/Muons/plugins/FWMuonDetailView.h"
#include "Fireworks/Calo/interface/FWECALDetailViewBuilder.h"
#include "Fireworks/Core/interface/FWModelId.h"
#include "Fireworks/Core/interface/FWEventItem.h"

//
// constructors and destructor
//
FWMuonDetailView::FWMuonDetailView()
{
}

FWMuonDetailView::~FWMuonDetailView()
{
}

//
// member functions
//
void FWMuonDetailView::build(const FWModelId &id, const reco::Muon* iMuon, TEveWindowSlot* base)
{
   if(0==iMuon) return;
   TEveWindowPack* eveWindow = base->MakePack();
   eveWindow->SetShowTitleBar(kFALSE);

   TEveWindow* ew(0);
   TEveWindowSlot* slot(0);
   TEveScene*      scene(0);
   TEveViewer*     viewer(0);
   TGVerticalFrame* ediFrame(0);

   ////////////////////////////////////////////////////////////////////////
   //                              Sub-view 1
   ///////////////////////////////////////////////////////////////////////

   // prepare window
   slot = eveWindow->NewSlot();
   ew = FWDetailViewBase::makePackViewer(slot, ediFrame, viewer, scene);
   ew->SetElementName("Muon based view");
   FWDetailViewBase::setEveWindow(ew);

   // build ECAL objects
   double eta = iMuon->eta();
   double phi = iMuon->phi();
   
   if ( iMuon->isEnergyValid() )
     {
	eta = iMuon->calEnergy().ecal_position.eta();
	phi = iMuon->calEnergy().ecal_position.phi();
     }

   FWECALDetailViewBuilder builder(id.item()->getEvent(), id.item()->getGeom(),
				   eta, phi, 10);
   builder.showSuperClusters(kGreen+2, kGreen+4);
   if ( iMuon->isEnergyValid() ) {
      std::vector<DetId> ids;
      ids.push_back(iMuon->calEnergy().ecal_id);
      builder.setColor(kYellow,ids);
   }
   TEveCaloLego* lego = builder.build();
   scene->AddElement(lego);
   
   addInfo(iMuon, scene);

   // draw legend in latex
   TRootEmbeddedCanvas* ec = new TRootEmbeddedCanvas("Embeddedcanvas", ediFrame, 100, 100, 0);
   ediFrame->AddFrame(ec, new TGLayoutHints(kLHintsExpandX|kLHintsExpandY));
   ediFrame->MapSubwindows();
   ediFrame->Layout();
   ec->GetCanvas()->SetBorderMode(0);
   makeLegend(iMuon, id, ec->GetCanvas());
   
   // draw axis at the window corners
   // std::cout << "TEveViewer: " << viewer << std::endl;
   TGLViewer* glv =  viewer->GetGLViewer();
   TEveCaloLegoOverlay* overlay = new TEveCaloLegoOverlay();
   overlay->SetShowPlane(kFALSE);
   overlay->SetShowPerspective(kFALSE);
   overlay->SetCaloLego(lego);
   overlay->SetShowScales(0); // temporary
   glv->AddOverlayElement(overlay);

   // set event handler and flip camera to top view at beginning
   glv->SetCurrentCamera(TGLViewer::kCameraOrthoXOY);
   TEveLegoEventHandler* eh = 
     new TEveLegoEventHandler(lego,(TGWindow*)glv->GetGLWidget(), (TObject*)glv);
   glv->SetEventHandler(eh);
   glv->UpdateScene();
   glv->CurrentCamera().Reset();

   scene->Repaint(true);
   viewer->GetGLViewer()->RequestDraw(TGLRnrCtx::kLODHigh);
   gEve->Redraw3D();
   
   ////////////////////////////////////////////////////////////////////////
   //                              Sub-view 2
   ///////////////////////////////////////////////////////////////////////
   // slot = eveWindow->NewSlot();
   // ew = FWDetailViewBase::makePackViewer(slot, ediFrame, viewer, scene);
   // ew->SetElementName("View C");

   // eveWindow->GetTab()->SetTab(1);
}

void
FWMuonDetailView::makeLegend(const reco::Muon *muon,
				 const FWModelId& id,
				 TCanvas* textCanvas)
{
   textCanvas->cd();

   TLatex* latex = new TLatex(0.02, 0.970, "");
   const double textsize(0.05);
   latex->SetTextSize(2*textsize);

   float_t x = 0.02;
   float   y = 0.95;
   float fontsize = latex->GetTextSize()*0.6;
   
   latex->DrawLatex(x, y, id.item()->modelName(id.index()).c_str() );
   y -= fontsize;
   latex->SetTextSize(textsize);
   fontsize = latex->GetTextSize()*0.6;

   latex->DrawLatex(x, y, Form(" p_{T} = %.1f GeV, #eta = %0.2f, #varphi = %0.2f", 
			       muon->pt(), muon->eta(), muon->phi()) );
   y -= fontsize;
   // summary
   if (muon->charge() > 0)
     latex->DrawLatex(x, y, " charge = +1");
   else 
     latex->DrawLatex(x, y, " charge = -1");
   y -= fontsize;
   
   if (! muon->isEnergyValid() ) return; 
   // delta phi/eta in
   latex->DrawLatex(x, y, "ECAL energy in:");
   y -= fontsize;
   latex->DrawLatex(x, y,  Form(" crossed crystalls = %.3f",
				muon->calEnergy().em) );
   y -= fontsize;
   latex->DrawLatex(x, y,  Form(" 3x3 crystall shape = %.3f",
				muon->calEnergy().emS9) );
   y -= fontsize;
   latex->DrawLatex(x, y,  Form(" 5x5 crystall shape = %.3f",
				muon->calEnergy().emS25) );
}

void
FWMuonDetailView::addInfo(const reco::Muon *i, TEveElementList* tList)
{
   // muon direction at vertex
   
   bool barrel = fabs(i->eta())<1.5;
   if ( i->isEnergyValid() ) barrel = fabs(i->calEnergy().ecal_position.eta()) < 1.5;
   
   TEvePointSet *direction = new TEvePointSet("muon direction");
   direction->SetTitle("Muon direction at vertex");
   direction->SetPickable(kTRUE);
   direction->SetMarkerStyle(2);
   if (barrel) {
      direction->SetNextPoint(i->eta(), i->phi(), 0);
      direction->SetMarkerSize(0.01);
   }else{
      direction->SetNextPoint(310*fabs(tan(i->theta()))*cos(i->phi()), 
			      310*fabs(tan(i->theta()))*sin(i->phi()),
			      0);
      direction->SetMarkerSize(1);
   }
   direction->SetMarkerColor(kYellow);
   tList->AddElement(direction);
   
   if (! i->isEnergyValid() ) return;
   
   // ecal position
   Double_t x(0), y(0);
   TEvePointSet *ecalposition = new TEvePointSet("ecal position");
   ecalposition->SetPickable(kTRUE);
   ecalposition->SetTitle("Track position at ECAL surface");
   if ( barrel ) {
      x = i->calEnergy().ecal_position.eta();
      y = i->calEnergy().ecal_position.phi();
   } else {
      x = i->calEnergy().ecal_position.x();
      y = i->calEnergy().ecal_position.y();
   }
   ecalposition->SetNextPoint(x,y,0);
   ecalposition->SetMarkerSize(1);
   ecalposition->SetMarkerStyle(4);
   ecalposition->SetMarkerColor(kRed);
   tList->AddElement(ecalposition);

   // hcal position
   TEvePointSet *hcalposition = new TEvePointSet("hcal position");
   hcalposition->SetTitle("Track position at HCAL surface");
   hcalposition->SetPickable(kTRUE);
   if ( barrel ) {
      x = i->calEnergy().hcal_position.eta();
      y = i->calEnergy().hcal_position.phi();
      hcalposition->SetMarkerSize(0.01);
   } else {
      x = i->calEnergy().hcal_position.x();
      y = i->calEnergy().hcal_position.y();
      hcalposition->SetMarkerSize(1);
   }
   hcalposition->SetNextPoint(x, y, 0);
   hcalposition->SetMarkerStyle(2);
   hcalposition->SetMarkerColor(kRed);
   tList->AddElement(hcalposition);
   
   // draw a line connecting the two positions
   TEveStraightLineSet *lines = new TEveStraightLineSet("Muon trajectory in ECAL","Muon trajectory in ECAL");
   lines->SetPickable(kTRUE);
   if ( barrel ) {
      lines->AddLine(i->calEnergy().ecal_position.eta(),
		     i->calEnergy().ecal_position.phi(),
		     0,
		     i->calEnergy().hcal_position.eta(),
		     i->calEnergy().hcal_position.phi(),
		     0);
   } else {
      lines->AddLine(i->calEnergy().ecal_position.x(),
		     i->calEnergy().ecal_position.y(),
		     0,
		     i->calEnergy().hcal_position.x(),
		     i->calEnergy().hcal_position.y(),
		     0);
   }
   lines->SetLineColor(kRed);
   tList->AddElement(lines);
}

REGISTER_FWDETAILVIEW(FWMuonDetailView);
