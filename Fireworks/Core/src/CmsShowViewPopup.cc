// -*- C++ -*-
//
// Package:     Core
// Class  :     CmsShowViewPopup
//
// Implementation:
//     <Notes on implementation>
//
// Original Author:
//         Created:  Wed Jun 25 15:15:04 EDT 2008
// $Id: CmsShowViewPopup.cc,v 1.25 2010/08/09 08:04:48 eulisse Exp $
//

// system include files
#include <iostream>
#include <boost/bind.hpp>
#include "TGLabel.h"
#include "TGButton.h"
#include "TG3DLine.h"
#include "TGFrame.h"
#include "TGTab.h"
#include "TG3DLine.h"
#include "TEveWindow.h"

// user include files
#include "Fireworks/Core/interface/CmsShowViewPopup.h"
#include "Fireworks/Core/interface/FWViewBase.h"
#include "Fireworks/Core/interface/FWParameterSetterBase.h"
#include "Fireworks/Core/src/FWDialogBuilder.h"
#include "Fireworks/Core/interface/FWColorManager.h"


//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
CmsShowViewPopup::CmsShowViewPopup(const TGWindow* p, UInt_t w, UInt_t h, FWColorManager* iCMgr, FWViewBase* vb, TEveWindow* ew) :
   TGTransientFrame(gClient->GetDefaultRoot(),p, w, h),
   m_mapped(kFALSE),
   m_viewLabel(0),
   m_paramGUI(0),
   m_saveImageButton(0),
   m_changeBackground(0),
   m_colorManager(iCMgr),
   m_viewBase(0),
   m_eveWindow(0)
{
   m_colorManager->colorsHaveChanged_.connect(boost::bind(&CmsShowViewPopup::backgroundColorWasChanged,this));

   SetCleanup(kDeepCleanup);

   // label
   TGHorizontalFrame* viewFrame = new TGHorizontalFrame(this);
   m_viewLabel = new TGLabel(viewFrame, "No view selected");
   try
   {
      TGFont* defaultFont = gClient->GetFontPool()->GetFont(m_viewLabel->GetDefaultFontStruct());
      m_viewLabel->SetTextFont(gClient->GetFontPool()->GetFont(defaultFont->GetFontAttributes().fFamily, 14, defaultFont->GetFontAttributes().fWeight + 2, defaultFont->GetFontAttributes().fSlant));
   }
   catch(...)
   {
      // FIXME: looks like under certain conditions (e.g. in full framework)
      // GetFontPool() throws when the default font is not found. This is a
      // quick workaround, but we should probably investigate more.
   }

   m_viewLabel->SetTextJustify(kTextLeft);
   viewFrame->AddFrame(m_viewLabel, new TGLayoutHints(kLHintsExpandX));
   AddFrame(viewFrame, new TGLayoutHints(kLHintsExpandX, 2, 2, 0, 0));
   // background
   m_changeBackground = new TGTextButton(this,"Change Background Color");
   backgroundColorWasChanged();
   AddFrame(m_changeBackground, new TGLayoutHints(kLHintsNormal, 0, 0, 5, 5));
   m_changeBackground->Connect("Clicked()","CmsShowViewPopup",this,"changeBackground()");
   // save image
   m_saveImageButton= new TGTextButton(this,"Save Image ...");
   AddFrame(m_saveImageButton);
   m_saveImageButton->Connect("Clicked()","CmsShowViewPopup",this,"saveImage()");

  // content frame
   AddFrame(new TGHorizontal3DLine(this), new TGLayoutHints(kLHintsExpandX, 0, 0, 5, 5));
   m_paramGUI = new ViewerParameterGUI(this);
   AddFrame(m_paramGUI,new TGLayoutHints(kLHintsExpandX|kLHintsExpandY));

   SetWindowName("View Controller");
}

// CmsShowViewPopup::CmsShowViewPopup(const CmsShowViewPopup& rhs)
// {
//    // do actual copying here;
// }

CmsShowViewPopup::~CmsShowViewPopup()
{
}

void
CmsShowViewPopup::reset(FWViewBase* vb, TEveWindow* ew)
{
   m_viewBase = vb;
   m_eveWindow = ew;

   m_paramGUI->reset();

   // fill content
   if(m_viewBase) {
      m_saveImageButton->SetEnabled(kTRUE);
      m_viewLabel->SetText(m_viewBase->typeName().c_str());
      m_viewBase->populateController(*m_paramGUI);
      m_paramGUI->populateComplete();
      fMain = m_eveWindow->GetEveFrame();
   }
   else {
      fMain = 0;
      m_viewLabel->SetText("No view selected");
      m_saveImageButton->SetEnabled(kFALSE);
   }

   MapSubwindows();
   Resize(GetDefaultSize());
   Layout();
   if (fMain)
   {
      CenterOnParent(kTRUE, TGTransientFrame::kTopRight);
   }
}

void
CmsShowViewPopup::CloseWindow()
{
   UnmapWindow();
   closed_.emit();
}

void
CmsShowViewPopup::MapWindow()
{
   TGWindow::MapWindow();
   m_mapped = true;
}

void
CmsShowViewPopup::UnmapWindow()
{
   TGWindow::UnmapWindow();
   m_mapped = false;
}

void
CmsShowViewPopup::saveImage()
{
   if(m_viewBase)
      m_viewBase->promptForSaveImageTo(this);
}

void
CmsShowViewPopup::changeBackground()
{
   m_colorManager->setBackgroundColorIndex( FWColorManager::kBlackIndex == m_colorManager->backgroundColorIndex()?
                                            FWColorManager::kWhiteIndex:
                                            FWColorManager::kBlackIndex);
}

void
CmsShowViewPopup::backgroundColorWasChanged()
{
   if(FWColorManager::kBlackIndex == m_colorManager->backgroundColorIndex()) {
      m_changeBackground->SetText("Change Background Color to White");
   } else {
      m_changeBackground->SetText("Change Background Color to Black");
   }
}

//==============================================================================

ViewerParameterGUI::ViewerParameterGUI(const TGFrame* p):
   TGCompositeFrame(p),
   m_tab(0)
{
   SetCleanup(kDeepCleanup);
}

void
ViewerParameterGUI::reset()
{ 
   m_setters.clear();
   if (m_tab) RemoveFrame(m_tab);
   m_tab = 0;
}


ViewerParameterGUI&
ViewerParameterGUI::requestTab(const char* name)
{
   if (!m_tab)
   {
      m_tab = new TGTab(this);
      AddFrame(m_tab,  new TGLayoutHints(kLHintsExpandX));
   }

   TGCompositeFrame* cont = 0;
   cont = m_tab->GetTabContainer(name);
   if (!cont)
      cont = m_tab->AddTab(name);
   m_tab->SetTab(name);
   return *this;
}

/* Add parameter setter in the current tab.*/
ViewerParameterGUI&
ViewerParameterGUI::addParam( const FWParameterBase* param)
{
   boost::shared_ptr<FWParameterSetterBase> ptr( FWParameterSetterBase::makeSetterFor((FWParameterBase*)param) );
   ptr->attach((FWParameterBase*)param, this);
   TGCompositeFrame* parent = m_tab->GetCurrentContainer();

   TGFrame* pframe = ptr->build(parent);
   parent->AddFrame(pframe, new TGLayoutHints(kLHintsExpandX));
   m_setters.push_back(ptr);

   pframe->MapWindow();
   pframe->MapSubwindows();
   pframe->Layout();
   parent->MapSubwindows();
   parent->Layout();
   m_tab->Layout();
   parent->Resize(parent->GetDefaultSize());
   return *this;
}
   


/* Add separator in current tab. */
ViewerParameterGUI&
ViewerParameterGUI::separator()
{
   assert(m_tab);
   TGHorizontal3DLine* s = new TGHorizontal3DLine(m_tab->GetCurrentContainer());
   m_tab->GetCurrentContainer()->AddFrame(s, new TGLayoutHints(kLHintsExpandX,4 ,4, 2, 2));

   return *this;
}

/* Setup after finish gui build. */
void
ViewerParameterGUI::populateComplete()
{
   if (m_tab) m_tab->SetTab(0);
}
