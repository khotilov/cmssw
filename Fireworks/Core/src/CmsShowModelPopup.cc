// -*- C++ -*-
//
// Package:     Core
// Class  :     CmsShowModelPopup
// 
// Implementation:
//     <Notes on implementation>
//
// Original Author:  
//         Created:  Fri Jun 27 11:23:08 EDT 2008
// $Id: CmsShowModelPopup.cc,v 1.4 2008/07/07 02:14:20 chrjones Exp $
//

// system include file
#include <iostream>
#include <set>
#include <sigc++/sigc++.h>
#include <boost/bind.hpp>
#include "TClass.h"
#include "TGFrame.h"
#include "TGButton.h"
#include "TGLabel.h"
#include "TGString.h"
#include "TColor.h"
#include "TG3DLine.h"
#include "TGFont.h"

// user include files
#include "Fireworks/Core/interface/CmsShowModelPopup.h"
#include "Fireworks/Core/interface/FWDisplayProperties.h"
#include "Fireworks/Core/src/FWListModel.h"
#include "Fireworks/Core/src/FWColorSelect.h"
#include "Fireworks/Core/interface/FWModelChangeSignal.h"
#include "Fireworks/Core/interface/FWModelChangeManager.h"
#include "Fireworks/Core/interface/FWEventItem.h"
//#include "Fireworks/Core/src/FWListModel.h"
#include "Fireworks/Core/interface/FWModelId.h"
#include "Fireworks/Core/interface/FWSelectionManager.h"
#include "Fireworks/Core/interface/FWDetailViewManager.h"

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
CmsShowModelPopup::CmsShowModelPopup(FWDetailViewManager* iManager,
                                     FWSelectionManager* iSelMgr,
                                     const TGWindow* p, UInt_t w, UInt_t h):
TGTransientFrame(gClient->GetDefaultRoot(),p,w,h),
m_detailViewManager(iManager)
{
  m_changes = iSelMgr->selectionChanged_.connect(boost::bind(&CmsShowModelPopup::fillModelPopup, this, _1));
   
  SetCleanup(kDeepCleanup);
  TGHorizontalFrame* objectFrame = new TGHorizontalFrame(this);
  m_modelLabel = new TGLabel(objectFrame, " ");
  TGFont* defaultFont = gClient->GetFontPool()->GetFont(m_modelLabel->GetDefaultFontStruct());
  m_modelLabel->SetTextFont(gClient->GetFontPool()->GetFont(defaultFont->GetFontAttributes().fFamily, 14, defaultFont->GetFontAttributes().fWeight + 2, defaultFont->GetFontAttributes().fSlant));
  m_modelLabel->SetTextJustify(kTextLeft);
  objectFrame->AddFrame(m_modelLabel, new TGLayoutHints(kLHintsExpandX));
  AddFrame(objectFrame, new TGLayoutHints(kLHintsExpandX, 2, 2, 0, 0));
  TGHorizontal3DLine* nameObjectSeperator = new TGHorizontal3DLine(this, 200, 5);
  AddFrame(nameObjectSeperator, new TGLayoutHints(kLHintsNormal, 0, 0, 5, 5));
  TGHorizontalFrame* colorSelectFrame = new TGHorizontalFrame(this, 200, 100);
  TGLabel* colorSelectLabel = new TGLabel(colorSelectFrame, "Color:");
  colorSelectFrame->AddFrame(colorSelectLabel, new TGLayoutHints(kLHintsNormal, 0, 50, 0, 0));
  TGString* graphicsLabel = new TGString(" ");
  Pixel_t selection = gVirtualX->GetPixel(kRed);
  std::vector<Pixel_t> colors;
  colors.push_back((Pixel_t)gVirtualX->GetPixel(kRed));
  colors.push_back((Pixel_t)gVirtualX->GetPixel(kBlue));
  colors.push_back((Pixel_t)gVirtualX->GetPixel(kYellow));
  colors.push_back((Pixel_t)gVirtualX->GetPixel(kGreen));
  colors.push_back((Pixel_t)gVirtualX->GetPixel(kCyan));
  colors.push_back((Pixel_t)gVirtualX->GetPixel(kMagenta));
  colors.push_back((Pixel_t)gVirtualX->GetPixel(kOrange));
  colors.push_back((Pixel_t)gVirtualX->GetPixel(TColor::GetColorDark(kRed)));
  colors.push_back((Pixel_t)gVirtualX->GetPixel(TColor::GetColorDark(kBlue)));
  colors.push_back((Pixel_t)gVirtualX->GetPixel(TColor::GetColorDark(kYellow)));
  colors.push_back((Pixel_t)gVirtualX->GetPixel(TColor::GetColorDark(kGreen)));
  colors.push_back((Pixel_t)gVirtualX->GetPixel(TColor::GetColorDark(kCyan)));
  colors.push_back((Pixel_t)gVirtualX->GetPixel(TColor::GetColorDark(kMagenta)));
  colors.push_back((Pixel_t)gVirtualX->GetPixel(TColor::GetColorDark(kOrange)));
  bool haveColor = false;
  for (std::vector<Pixel_t>::const_iterator iCol = colors.begin(); iCol != colors.end(); ++iCol) {
    if (*iCol == selection) haveColor = true;
  }
  if(!haveColor) {
    printf("Error: Color is not present in palette!\n");
    colors.push_back(selection);
  }
  m_colorSelectWidget = new FWColorSelect(colorSelectFrame, graphicsLabel, selection, colors, -1);
  m_colorSelectWidget->SetEnabled(kFALSE);
  colorSelectFrame->AddFrame(m_colorSelectWidget);
  AddFrame(colorSelectFrame);
  TGHorizontal3DLine* colorVisSeperator = new TGHorizontal3DLine(this, 200, 5);
  AddFrame(colorVisSeperator, new TGLayoutHints(kLHintsNormal, 0, 0, 5, 5));
  m_isVisibleButton = new TGCheckButton(this, "Visible");
  m_isVisibleButton->SetState(kButtonDown, kFALSE);
  m_isVisibleButton->SetEnabled(kFALSE);
  AddFrame(m_isVisibleButton);
  AddFrame(new TGHorizontal3DLine(this, 200, 5), new TGLayoutHints(kLHintsNormal, 0, 0, 5, 5));
  m_openDetailedViewButton = new TGTextButton(this,"Open Detailed View");
  m_openDetailedViewButton->SetEnabled(kFALSE);
  AddFrame(m_openDetailedViewButton);
  m_openDetailedViewButton->Connect("Clicked()","CmsShowModelPopup", this, "openDetailedView()");
  this->Connect("CloseWindow()","CmsShowModelPopup", this, "windowClosing()");
  SetWindowName("Model Inspector");
  Resize(GetDefaultSize());
  MapSubwindows();
  Layout();
  MapWindow();
}

// CmsShowModelPopup::CmsShowModelPopup(const CmsShowModelPopup& rhs)
// {
//    // do actual copying here;
// }

CmsShowModelPopup::~CmsShowModelPopup()
{
   m_changes.disconnect();
   disconnectAll();
}

//
// assignment operators
//
// const CmsShowModelPopup& CmsShowModelPopup::operator=(const CmsShowModelPopup& rhs)
// {
//   //An exception safe implementation is
//   CmsShowModelPopup temp(rhs);
//   swap(rhs);
//
//   return *this;
// }

//
// member functions
//
void
CmsShowModelPopup::fillModelPopup(const FWSelectionManager& iSelMgr) {
  disconnectAll();
  if (iSelMgr.selected().size() > 0) {
    bool multipleNames(false);
    bool multipleColors(false);
    bool multipleVis(false);
    m_models = iSelMgr.selected();
    FWModelId id;
    FWModelId prevId;
    const FWEventItem* item = 0;
    const FWEventItem* prevItem = 0;
    for (std::set<FWModelId>::iterator it_mod = m_models.begin(); it_mod != m_models.end(); ++it_mod) {
      if (it_mod != m_models.begin()) {
	item = (*it_mod).item();
	if (item->name() != prevItem->name()) multipleNames = true;
	if (item->modelInfo((*it_mod).index()).displayProperties().color() != prevItem->modelInfo(prevId.index()).displayProperties().color()) 
	  multipleColors = true;
	if (item->modelInfo((*it_mod).index()).displayProperties().isVisible() != prevItem->modelInfo(prevId.index()).displayProperties().isVisible())
	  multipleVis = true;
      }
      prevId = *it_mod;
      prevItem = (*it_mod).item();
    }
    id = *m_models.begin();
    item = (*(m_models.begin())).item();
    if (multipleNames) 
      m_modelLabel->SetText("Multiple Items");
    else 
      m_modelLabel->SetText(item->name().c_str());
    if(m_models.size()==1) {
       m_openDetailedViewButton->SetEnabled(m_detailViewManager->haveDetailViewFor(id));
    }
    m_colorSelectWidget->SetColor(gVirtualX->GetPixel(item->modelInfo(id.index()).displayProperties().color()));
    m_isVisibleButton->SetDisabledAndSelected(item->modelInfo(id.index()).displayProperties().isVisible());
    m_colorSelectWidget->SetEnabled(kTRUE);
    m_isVisibleButton->SetEnabled(kTRUE);
    if (!(m_colorSelectWidget->HasConnection("ColorSelected(Pixel_t)")))
      m_colorSelectWidget->Connect("ColorSelected(Pixel_t)", "CmsShowModelPopup", this, "changeModelColor(Pixel_t)");
    if (!(m_isVisibleButton->HasConnection("Toggled(Bool_T)")))
      m_isVisibleButton->Connect("Toggled(Bool_t)", "CmsShowModelPopup", this, "toggleModelVisible(Bool_t)");
    //    m_displayChangedConn = m_item->defaultDisplayPropertiesChanged_.connect(boost::bind(&CmsShowEDI::updateDisplay, this));
    m_modelChangedConn = item->changed_.connect(boost::bind(&CmsShowModelPopup::updateDisplay, this));
    //    m_selectionChangedConn = m_selectionManager->selectionChanged_.connect(boost::bind(&CmsShowEDI::updateSelection, this));
    m_destroyedConn = item->goingToBeDestroyed_.connect(boost::bind(&CmsShowModelPopup::disconnectAll, this));
    Layout();
  }    
}

void
CmsShowModelPopup::updateDisplay() {
  const FWEventItem* item;
  for (std::set<FWModelId>::iterator it_mod = m_models.begin(); it_mod != m_models.end(); ++it_mod) {
    item = (*it_mod).item();
    m_colorSelectWidget->SetColor(gVirtualX->GetPixel(item->modelInfo((*it_mod).index()).displayProperties().color()));
    m_isVisibleButton->SetState(item->modelInfo((*it_mod).index()).displayProperties().isVisible() ? kButtonDown : kButtonUp, kFALSE);
  }
}

void
CmsShowModelPopup::disconnectAll() {
  m_modelChangedConn.disconnect();
  m_destroyedConn.disconnect();
  m_colorSelectWidget->Disconnect("ColorSelected(Pixel_t)", this, "changeModelColor(Pixel_t)");
  m_isVisibleButton->Disconnect("Toggled(Bool_t)", this, "toggleModelVisible(Bool_t)");
  //  m_item = 0;
  //  m_model = 0;
  m_modelLabel->SetText(" ");
  m_colorSelectWidget->SetColor(gVirtualX->GetPixel(kRed));
  m_isVisibleButton->SetDisabledAndSelected(kTRUE);
  m_colorSelectWidget->SetEnabled(kFALSE);
  m_isVisibleButton->SetEnabled(kFALSE);
  m_openDetailedViewButton->SetEnabled(kFALSE);
}

void
CmsShowModelPopup::changeModelColor(Pixel_t pixel) {
  Color_t color(TColor::GetColor(pixel));
  const FWEventItem* item;
  for (std::set<FWModelId>::iterator it_mod = m_models.begin(); it_mod != m_models.end(); ++it_mod) {
    item = (*it_mod).item();
    const FWDisplayProperties changeProperties(color, item->modelInfo((*it_mod).index()).displayProperties().isVisible());
    item->setDisplayProperties((*it_mod).index(), changeProperties);
  }
}

void
CmsShowModelPopup::toggleModelVisible(Bool_t on) {
  const FWEventItem* item;
  for (std::set<FWModelId>::iterator it_mod = m_models.begin(); it_mod != m_models.end(); ++it_mod) {
    item = (*it_mod).item();
    const FWDisplayProperties changeProperties(item->modelInfo((*it_mod).index()).displayProperties().color(), on);
    item->setDisplayProperties((*it_mod).index(), changeProperties);
  }
  //  const FWDisplayProperties changeProperties(m_item->modelInfo(m_model->index()).displayProperties().color(), on);
  //  m_item->setDisplayProperties(m_model->index(), changeProperties);
}

void 
CmsShowModelPopup::openDetailedView()
{
   m_detailViewManager->openDetailViewFor( *(m_models.begin()) );
}

void
CmsShowModelPopup::windowClosing()
{
   delete this;
}

//
// const member functions
//

//
// static member functions
//
