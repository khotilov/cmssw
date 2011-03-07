// -*- C++ -*-
//
// Package:     Core
// Class  :     CmsShowCommon
// 
// Implementation:
//     [Notes on implementation]
//
// Original Author:  Alja Mrak-Tadel
//         Created:  Fri Sep 10 14:50:32 CEST 2010
// $Id: CmsShowCommon.cc,v 1.12 2010/11/27 22:08:24 amraktad Exp $
//

// system include files
#include <boost/bind.hpp>

// user include files
#include "TEveManager.h"
#include "TEveTrackPropagator.h"

#include "Fireworks/Core/interface/CmsShowCommon.h"
#include "Fireworks/Core/interface/FWEveView.h"

#include "Fireworks/Core/interface/Context.h"


//
// constructors and destructor
//
CmsShowCommon::CmsShowCommon(fireworks::Context* c):
   FWConfigurableParameterizable(2),
   m_trackBreak(this, "     ", 1l, 0l, 2l), // do not want to render text at setter
   m_drawBreakPoints(this, "Show y=0 points as markers", true),
   m_context(c),
   m_backgroundColor(this, "backgroundColIdx", 1l, 0l, 1000l),
   m_gamma(this, "Brightness", 0l, -15l, 15l),
   m_geomTransparency2D(this, "Transparency 2D", long(colorManager()->geomTransparency(true)), 0l, 100l),
   m_geomTransparency3D(this, "Transparency 3D", long(colorManager()->geomTransparency(false)), 0l, 100l),
   m_energyScale(new FWViewEnergyScale("global", 2))
{
   // projections 
   m_trackBreak.addEntry(0, "Jump to proper hemisphere");
   m_trackBreak.addEntry(1, "Stay on first point side");
   m_trackBreak.addEntry(2, "Stay on last point side");

   // colors
   char name[32];
   for (int i = 0; i < kFWGeomColorSize; ++i)
   {
      snprintf(name, 31, "GeometryColor %d ", i);
      m_geomColors[i] = new FWLongParameter(this, name   , long(colorManager()->geomColor(FWGeomColorIndex(i))), 1000l, 1100l);
   }  

   m_trackBreak.changed_.connect(boost::bind(&CmsShowCommon::setTrackBreakMode, this));
   m_drawBreakPoints.changed_.connect(boost::bind(&CmsShowCommon::setDrawBreakMarkers, this));
   m_gamma.changed_.connect(boost::bind(&CmsShowCommon::setGamma, this));
}

CmsShowCommon::~CmsShowCommon()
{
}


const FWColorManager* CmsShowCommon::colorManager() const
{
   return m_context->colorManager();
}
//
// member functions
//

void
CmsShowCommon::setDrawBreakMarkers()
{
   m_context->getTrackPropagator()->SetRnrPTBMarkers(m_drawBreakPoints.value());
   m_context->getTrackerTrackPropagator()->SetRnrPTBMarkers(m_drawBreakPoints.value());
   m_context->getMuonTrackPropagator()->SetRnrPTBMarkers(m_drawBreakPoints.value());
   gEve->Redraw3D();
}

void
CmsShowCommon::setTrackBreakMode()
{
   m_context->getTrackPropagator()->SetProjTrackBreaking(m_trackBreak.value());
   m_context->getTrackerTrackPropagator()->SetProjTrackBreaking(m_trackBreak.value());
   m_context->getMuonTrackPropagator()->SetProjTrackBreaking(m_trackBreak.value());
   gEve->Redraw3D();
}

void
CmsShowCommon::setGamma()
{
   m_context->colorManager()->setBrightness(m_gamma.value());
}

void
CmsShowCommon::switchBackground()
{
   m_context->colorManager()->switchBackground();
   m_backgroundColor.set(colorManager()->background());
}

void
CmsShowCommon::setGeomColor(FWGeomColorIndex cidx, Color_t iColor)
{
   m_geomColors[cidx]->set(iColor);
   m_context->colorManager()->setGeomColor(cidx, iColor);
}

void
CmsShowCommon::setGeomTransparency(int iTransp, bool projected)
{
   if (projected)
      m_geomTransparency2D.set(iTransp);
   else
      m_geomTransparency3D.set(iTransp);

   m_context->colorManager()->setGeomTransparency(iTransp, projected);
}

//______________________________________________________________________________

void
CmsShowCommon::addTo(FWConfiguration& oTo) const
{
   m_backgroundColor.set(int(colorManager()->background()));

   FWConfigurableParameterizable::addTo(oTo);
   m_energyScale->addTo(oTo);
}

void
CmsShowCommon::setFrom(const FWConfiguration& iFrom)
{  
   for(const_iterator it =begin(), itEnd = end();
       it != itEnd;
       ++it) {
      (*it)->setFrom(iFrom);      
   }  
 
   // handle old and new energy scale configuration if existing
   if (iFrom.valueForKey("ScaleMode"))
   {
      long mode  = atol(iFrom.valueForKey("ScaleMode")->value().c_str());      

      float convert;
      if (iFrom.valueForKey("EnergyToLength [GeV/m]"))
         convert = atof(iFrom.valueForKey("EnergyToLength [GeV/m]")->value().c_str());
      else
         convert = atof(iFrom.valueForKey("ValueToHeight [GeV/m]")->value().c_str());

      float maxH;
      if (iFrom.valueForKey("MaximumLength [m]"))
         maxH = atof(iFrom.valueForKey("MaximumLength [m]")->value().c_str());
      else
         maxH = atof(iFrom.valueForKey("MaxTowerH [m]")->value().c_str());
         
      int et = atoi(iFrom.valueForKey("PlotEt")->value().c_str());
      m_energyScale->SetFromCmsShowCommonConfig(mode, convert, maxH, et);
   }
      
   // background
   FWColorManager* cm =  m_context->colorManager();
   cm->setBackgroundAndBrightness( FWColorManager::BackgroundColorIndex(m_backgroundColor.value()), m_gamma.value());
   
   // geom colors
   cm->setGeomTransparency( m_geomTransparency2D.value(), true);
   cm->setGeomTransparency( m_geomTransparency3D.value(), false);

   for (int i = 0; i < kFWGeomColorSize; ++i)
      cm->setGeomColor(FWGeomColorIndex(i), m_geomColors[i]->value());
}
