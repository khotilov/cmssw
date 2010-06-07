// -*- C++ -*-
//
// Package:     Core
// Class  :     Context
//
// Implementation:
//     <Notes on implementation>
//
// Original Author:  Chris Jones
//         Created:  Tue Sep 30 14:57:12 EDT 2008
// $Id: Context.cc,v 1.19 2010/05/31 13:01:25 amraktad Exp $
//

// system include files

// user include files
#include "TH2.h"
#include "TMath.h"
#include "TEveTrackPropagator.h"
#include "TEveCaloData.h"
#include "Fireworks/Core/interface/fw3dlego_xbins.h"

#include "Fireworks/Core/interface/Context.h"
#include "Fireworks/Core/interface/FWMagField.h"

using namespace fireworks;
//
// constants, enums and typedefs
//

//
// static data member definitions
//

const float Context::s_ecalR = 126;
const float Context::s_ecalZ = 306;
const float Context::s_transitionAngle = atan( s_ecalR/s_ecalZ );

//
// constructors and destructor
//
Context::Context(FWModelChangeManager* iCM,
                 FWSelectionManager* iSM,
                 FWEventItemsManager* iEM,
                 FWColorManager* iColorM
                 ) :
   m_changeManager(iCM),
   m_selectionManager(iSM),
   m_eventItemsManager(iEM),
   m_colorManager(iColorM),
   m_geom(0),
   m_propagator(0),
   m_trackerPropagator(0),
   m_muonPropagator(0),
   m_magField(0),
   m_caloData(0),
   m_caloDataHF(0)
{
}

Context::~Context()
{
   delete m_magField;
}

void
Context::initEveElements()
{
   m_magField = new FWMagField();

   // common propagator, helix stepper
   m_propagator = new TEveTrackPropagator();
   m_propagator->SetMagFieldObj(m_magField, false);
   m_propagator->SetMaxR(123.0);
   m_propagator->SetMaxZ(300.0);
   m_propagator->SetProjTrackBreaking(TEveTrackPropagator::kPTB_UseLastPointPos);
   m_propagator->SetRnrPTBMarkers(kTRUE);
   m_propagator->IncDenyDestroy();
   // tracker propagator
   m_trackerPropagator = new TEveTrackPropagator();
   m_trackerPropagator->SetStepper( TEveTrackPropagator::kRungeKutta );
   m_trackerPropagator->SetMagFieldObj(m_magField, false);
   m_trackerPropagator->SetDelta(0.02);
   m_trackerPropagator->SetMaxR(123.0);
   m_trackerPropagator->SetMaxZ(300.0);
   m_trackerPropagator->SetProjTrackBreaking(TEveTrackPropagator::kPTB_UseLastPointPos);
   m_trackerPropagator->SetRnrPTBMarkers(kTRUE);
   m_trackerPropagator->IncDenyDestroy();
   // muon propagator
   m_muonPropagator = new TEveTrackPropagator();
   m_muonPropagator->SetStepper( TEveTrackPropagator::kRungeKutta );
   m_muonPropagator->SetMagFieldObj(m_magField, false);
   m_muonPropagator->SetDelta(0.05);
   m_muonPropagator->SetMaxR(850.f);
   m_muonPropagator->SetMaxZ(1100.f);
   m_muonPropagator->SetProjTrackBreaking(TEveTrackPropagator::kPTB_UseLastPointPos);
   m_muonPropagator->SetRnrPTBMarkers(kTRUE);
   m_muonPropagator->IncDenyDestroy();

   // general calo data
   {
      m_caloData = new TEveCaloDataHist();
      m_caloData->IncDenyDestroy();
      
      Bool_t status = TH1::AddDirectoryStatus();
      TH1::AddDirectory(kFALSE); //Keeps histogram from going into memory
      TH2F* dummy = new TH2F("background",
                             "background",
                             82, fw3dlego::xbins,
                             72, -1*TMath::Pi(), TMath::Pi());
      
      TH1::AddDirectory(status);
      Int_t sliceIndex = m_caloData->AddHistogram(dummy);
      (m_caloData)->RefSliceInfo(sliceIndex).Setup("background", 0., 0);
   }
   // HF calo data
   {
      m_caloDataHF = new TEveCaloDataVec(1);
      m_caloDataHF->IncDenyDestroy();
      
      m_caloDataHF->RefSliceInfo(0).Setup("bg", 0.3, kRed);
      m_caloDataHF->SetEtaBins(new TAxis(82, fw3dlego::xbins));
      m_caloDataHF->SetPhiBins(new TAxis( 72, -1*TMath::Pi(), TMath::Pi()));
   }
}

void
Context::deleteEveElements()
{
   m_propagator->DecDenyDestroy();
   m_trackerPropagator->DecDenyDestroy();
   m_muonPropagator->DecDenyDestroy();
   m_caloData->DecDenyDestroy();
   m_caloDataHF->DecDenyDestroy();
}

//
// static member functions
//


//
// static data member definitions
//
//
// const member functions
//
//
// member functions
