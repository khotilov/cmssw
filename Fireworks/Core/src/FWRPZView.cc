// -*- C++ -*-
//
// Package:     Core
// Class  :     FWRPZView
//
// Implementation:
//     <Notes on implementation>
//
// Original Author:  Chris Jones
//         Created:  Tue Feb 19 10:33:25 EST 2008
// $Id: FWRPZView.cc,v 1.37 2011/02/24 19:15:31 matevz Exp $
//

// system include files
#include <stdexcept>
#include <boost/bind.hpp>
#include <boost/shared_ptr.hpp>

#include "TGLViewer.h"
#include "TGLScenePad.h"
#include "TEveManager.h"
#include "TEveElement.h"
#include "TEveScene.h"
#include "TEveProjections.h"
#include "TEveProjectionAxes.h"
#include "TGLabel.h"

//!!! FIXME add get/sets for TEveCalo2D for CellIDs
//          add option in FWEveProjectionManager::SetCenter()
#define protected public  
#include "TEveCalo.h"
#include "TEveProjectionManager.h"
#undef protected

// user include files
#include "Fireworks/Core/interface/FWRPZView.h"
#include "Fireworks/Core/interface/FWRPZViewGeometry.h"
#include "Fireworks/Core/interface/FWBeamSpot.h"
#include "Fireworks/Core/interface/Context.h"
#include "Fireworks/Core/interface/fwLog.h"
#include "Fireworks/Core/interface/FWViewContext.h"
#include "Fireworks/Core/interface/FWViewContext.h"
#include "Fireworks/Core/interface/FWViewEnergyScale.h"
#include "Fireworks/Core/interface/CmsShowViewPopup.h"

FWRPZViewGeometry* FWRPZView::s_geometryList = 0;
const float FWRPZView::s_distortF = 0.001;
const float FWRPZView::s_distortFInv = 1000;
//
// constructors and destructor
//
FWRPZView::FWRPZView(TEveWindowSlot* iParent, FWViewType::EType id) :
   FWEveView(iParent, id),
   m_calo(0),
   m_fishEyeDistortion(this,"Distortion",0., 0., 100.),
   m_fishEyeR(this,"FixedRadius",(double)fireworks::Context::caloR1(), 0.0, 150.0),

   m_caloDistortion(this,"Calo compression",1.0,0.01,10.),
   m_muonDistortion(this,"Muon compression",0.2,0.01,10.),
   m_showProjectionAxes(this,"Show projection axes", false),
   m_compressMuon(this,"Compress detectors",false),
   m_showHF(0),
   m_showEndcaps(0)
{

   TEveProjection::EPType_e projType = (id == FWViewType::kRhoZ) ? TEveProjection::kPT_RhoZ : TEveProjection::kPT_RPhi;

   m_projMgr = new TEveProjectionManager(projType);
   m_projMgr->IncDenyDestroy();
   m_projMgr->SetImportEmpty(kTRUE);

   m_projMgr->GetProjection()->SetDistortion(m_fishEyeDistortion.value()*s_distortF);
   m_projMgr->GetProjection()->SetFixR(m_fishEyeR.value());

   if ( id == FWViewType::kRhoPhi || id == FWViewType::kRhoPhiPF) {
      m_projMgr->GetProjection()->AddPreScaleEntry(0, fireworks::Context::caloR1(), 1.0);
      m_projMgr->GetProjection()->AddPreScaleEntry(0, 300, 0.2);
   } else {
      m_projMgr->GetProjection()->AddPreScaleEntry(0, fireworks::Context::caloR1(), 1.0);
      m_projMgr->GetProjection()->AddPreScaleEntry(1, 310, 1.0);
      m_projMgr->GetProjection()->AddPreScaleEntry(0, 370, 0.2);
      m_projMgr->GetProjection()->AddPreScaleEntry(1, 580, 0.2);
   }
 

   // camera  
   viewerGL()->SetCurrentCamera(TGLViewer::kCameraOrthoXOY);
   if ( TGLOrthoCamera* camera = dynamic_cast<TGLOrthoCamera*>( &(viewerGL()->CurrentCamera()) ) ) {
      camera->SetZoomMax(1e6);
   }
   geoScene()->GetGLScene()->SetSelectable(kFALSE);

   m_axes = new TEveProjectionAxes(m_projMgr);
   m_axes->SetRnrState(m_showProjectionAxes.value());
   m_showProjectionAxes.changed_.connect(boost::bind(&FWRPZView::showProjectionAxes,this));
   eventScene()->AddElement(m_axes);

   if ( id != FWViewType::kRhoZ ) {
      m_showEndcaps = new FWBoolParameter(this,"Include EndCaps", true);
      m_showEndcaps->changed_.connect(  boost::bind(&FWRPZView::setEtaRng, this) );
      m_showHF = new FWBoolParameter(this,"Include HF", true);
      m_showHF->changed_.connect(  boost::bind(&FWRPZView::setEtaRng, this) );
   }

   m_fishEyeDistortion.changed_.connect(boost::bind(&FWRPZView::doFishEyeDistortion,this));
   m_fishEyeR.changed_.connect(boost::bind(&FWRPZView::doFishEyeDistortion,this));

   m_caloDistortion.changed_.connect(boost::bind(&FWRPZView::doPreScaleDistortion,this));
   m_muonDistortion.changed_.connect(boost::bind(&FWRPZView::doPreScaleDistortion,this));
   m_compressMuon.changed_.connect(boost::bind(&FWRPZView::doCompression,this,_1));
}

FWRPZView::~FWRPZView()
{
   m_calo->Destroy();
   m_projMgr->DecDenyDestroy();
}

//
// member functions
//

TEveCaloViz*
FWRPZView::getEveCalo() const
{
   return static_cast<TEveCaloViz*>(m_calo);
}

void
FWRPZView::setContext(const fireworks::Context& ctx)
{
   FWEveView::setContext(ctx);

   if (!s_geometryList)
   {
      s_geometryList = new  FWRPZViewGeometry(ctx);
      s_geometryList->IncDenyDestroy();
   }
   m_projMgr->ImportElements(s_geometryList->getGeoElements(typeId()), geoScene());

   TEveCaloData* data = context().getCaloData();

   TEveCalo3D* calo3d = new TEveCalo3D(data);

   m_calo = static_cast<TEveCalo2D*> (m_projMgr->ImportElements(calo3d, eventScene()));

   if (typeId() == FWViewType::kRhoPhiPF)
      m_calo->SetBarrelRadius(177);
   else
      m_calo->SetBarrelRadius(context().caloR1(false));

   m_calo->SetEndCapPos(context().caloZ1(false));
   m_calo->SetAutoRange(false);
   m_calo->SetScaleAbs(true);
}

void
FWRPZView::eventBegin()
{  
   if (context().getBeamSpot())
   {
      FWBeamSpot& b = *(context().getBeamSpot());

      // AMT: should call directly TEveProjectionManager::SetCenter()
      // and add an option to skip EveProjectionManager::ProjectChildren() 
    
      fwLog(fwlog::kDebug) << Form("%s::eventBegin Set projection center (%f, %f, %f) \n", typeName().c_str(), b.x0(), b.y0(), b.z0());
      m_projMgr->fCenter.Set(b.x0(), b.y0(), b.z0());
      m_projMgr->fProjection->SetCenter(m_projMgr->fCenter);

      // camera move
      TGLCamera& cam = viewerGL()->CurrentCamera();
      cam.SetExternalCenter(true);
      if (typeId() != FWViewType::kRhoZ)
      {
         double r = TMath::Sqrt( b.x0()*b.x0() +  b.y0()*b.y0());
         cam.SetCenterVec(b.z0(), TMath::Sign(r, b.y0()), 0);
      }
      else
      {
         cam.SetCenterVec(b.x0(), b.y0(), b.z0());
      }
   }
}

void
FWRPZView::doFishEyeDistortion()
{
   TEveProjection* p = m_projMgr->GetProjection();
    if (p->GetDistortion() != m_fishEyeDistortion.value()*s_distortFInv)
   p->SetDistortion(m_fishEyeDistortion.value()*s_distortF);
   if (p->GetFixR() != m_fishEyeR.value())
      p->SetFixR(m_fishEyeR.value());

   m_projMgr->ProjectChildren();
   gEve->Redraw3D();
}

void
FWRPZView::doPreScaleDistortion()
{
   if ( typeId() == FWViewType::kRhoPhi || typeId() == FWViewType::kRhoPhiPF ) {
      m_projMgr->GetProjection()->ChangePreScaleEntry(0,1,m_caloDistortion.value());
      m_projMgr->GetProjection()->ChangePreScaleEntry(0,2,m_muonDistortion.value());
   } else {
      m_projMgr->GetProjection()->ChangePreScaleEntry(0,1,m_caloDistortion.value());
      m_projMgr->GetProjection()->ChangePreScaleEntry(0,2,m_muonDistortion.value());
      m_projMgr->GetProjection()->ChangePreScaleEntry(1,1,m_caloDistortion.value());
      m_projMgr->GetProjection()->ChangePreScaleEntry(1,2,m_muonDistortion.value());
   }
   m_projMgr->UpdateName();
   m_projMgr->ProjectChildren();
   gEve->Redraw3D();
}

void
FWRPZView::doCompression(bool flag)
{
   m_projMgr->GetProjection()->SetUsePreScale(flag);
   m_projMgr->UpdateName();
   m_projMgr->ProjectChildren();
   gEve->Redraw3D();
}

void
FWRPZView::importElements(TEveElement* iChildren, float iLayer, TEveElement* iProjectedParent)
{
   float oldLayer = m_projMgr->GetCurrentDepth();
   m_projMgr->SetCurrentDepth(iLayer);
   //make sure current depth is reset even if an exception is thrown
   boost::shared_ptr<TEveProjectionManager> sentry(m_projMgr,
                                                   boost::bind(&TEveProjectionManager::SetCurrentDepth,
                                                               _1,oldLayer));
   m_projMgr->ImportElements(iChildren,iProjectedParent);
}


void
FWRPZView::addTo(FWConfiguration& iTo) const
{
   FWEveView::addTo(iTo);
   TGLOrthoCamera* camera = dynamic_cast<TGLOrthoCamera*>( &(viewerGL()->CurrentCamera()) );
   if (camera) addToOrthoCamera(camera, iTo);
}

void
FWRPZView::setFrom(const FWConfiguration& iFrom)
{
   FWEveView::setFrom(iFrom);
   
   TGLOrthoCamera* camera = dynamic_cast<TGLOrthoCamera*>( &(viewerGL()->CurrentCamera()) );
   if (camera) setFromOrthoCamera(camera, iFrom);
}

void
FWRPZView::setEtaRng()
{
   if (typeId() != FWViewType::kRhoZ)
   {
      // rng controllers only in RhoPhi
      double eta_range = context().caloMaxEta();
      if (!m_showHF->value() ) eta_range = 3.0;
      if (!m_showEndcaps->value() ) eta_range = context().caloTransEta();
      m_calo->SetEta(-eta_range,eta_range);
   }

   FWEveView::setupEnergyScale();
}

void
FWRPZView::voteCaloMaxVal()
{
   if (! m_calo->GetData()->Empty())
   {
      m_calo->AssertCellIdCache();
      Float_t sumEt, sumE;
      TEveCaloData::CellData_t cellData;
      typedef std::vector<TEveCaloData::vCellId_t*>           vBinCells_t;
      typedef std::vector<TEveCaloData::vCellId_t*>::iterator vBinCells_i;

      vBinCells_t   cellLists = m_calo->fCellLists;
      for (vBinCells_i it = cellLists.begin(); it != cellLists.end(); it++)
      {
         TEveCaloData::vCellId_t* binCells = *it;
         if (binCells) {
            sumEt =  0; sumE = 0;
            TEveCaloData::vCellId_i a = binCells->end();
      
            for (TEveCaloData::vCellId_i k = binCells->begin(); k != a; ++k)
            { 
               m_calo->GetData()->GetCellData((*k), cellData);
               sumEt += cellData.Value(true);
               sumE  += cellData.Value(false);
            } 
            // printf("vote sum %f %f \n", sumEt , sumE);
            context().voteMaxEtAndEnergy(sumEt, sumE);
         }
      }
   }
}

void FWRPZView::showProjectionAxes( )
{   
   if ( m_showProjectionAxes.value() )
      m_axes->SetRnrState(kTRUE);
   else
      m_axes->SetRnrState(kFALSE);
   gEve->Redraw3D();
}

void 
FWRPZView::populateController(ViewerParameterGUI& gui) const
{
   FWEveView::populateController(gui);

   gui.requestTab("Projection");
   TGCompositeFrame* f = gui.getTabContainer();

   f->AddFrame(new TGLabel(f, "FishEye:"));
   gui.addParam(&m_fishEyeDistortion).
      addParam(&m_fishEyeR).
      separator();

   f->AddFrame(new TGLabel(f, "PreScales:"));

   gui.requestTab("Projection").
      addParam(&m_compressMuon).
      addParam(&m_muonDistortion).
      addParam(&m_caloDistortion).separator();

   gui.requestTab("Projection").
      addParam(&m_showProjectionAxes);

   if (typeId() == FWViewType::kRhoPhi || typeId() == FWViewType::kRhoPhiPF) 
   {
      gui.requestTab("Calo").
         addParam(m_showHF).
         addParam(m_showEndcaps);
   }
}

