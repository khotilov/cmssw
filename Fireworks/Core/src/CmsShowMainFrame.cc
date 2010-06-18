
// -*- C++ -*-
//
// Package:     Core
// Class  :     CmsShowMainFrame
//
// Implementation:
//     <Notes on implementation>
//
// Original Author:  Chris Jones
//         Created:  Thu May 29 20:58:23 CDT 2008
// $Id: CmsShowMainFrame.cc,v 1.94 2010/06/08 18:30:46 matevz Exp $
//
// hacks
#include "DataFormats/FWLite/interface/Event.h"

// system include files
#include <TCollection.h>
#include <TApplication.h>
#include <TGClient.h>
#include <TGSplitFrame.h>
#include <TGLayout.h>
#include <TGButton.h>
#include <TGMenu.h>
#include <TGLabel.h>
#include <TGTab.h>
#include <TGStatusBar.h>
#include <TGNumberEntry.h>
#include <KeySymbols.h>
#include <TGTextEntry.h>
#include <TGSlider.h>

#include <TSystem.h>
#include <TImage.h>
// user include files
#include "DataFormats/Provenance/interface/EventID.h"
#include "Fireworks/Core/interface/CSGAction.h"
#include "Fireworks/Core/interface/CSGContinuousAction.h"
#include "Fireworks/Core/interface/CmsShowMainFrame.h"
#include "Fireworks/Core/interface/ActionsList.h"
#include "Fireworks/Core/interface/BuilderUtils.h"

#include "Fireworks/Core/interface/FWGUIManager.h"
#include "Fireworks/Core/interface/FWCustomIconsButton.h"

#include "Fireworks/Core/interface/FWIntValueListener.h"
#include "Fireworks/Core/interface/fwLog.h"
#include "Fireworks/Core/src/FWCheckBoxIcon.h"

#include <fstream>

//
// constants, enums and typedefs
//

//
// static data member definitions
//


//
// constructors and destructor
//
CmsShowMainFrame::CmsShowMainFrame(const TGWindow *p,UInt_t w,UInt_t h,FWGUIManager *m) :
   TGMainFrame(p, w, h),
   m_filterEnableBtn(),
   m_filterShowGUIBtn(),
   m_runEntry(0),
   m_eventEntry(0),
   m_delaySliderListener(0),
   m_fworksAbout(0)
{
   const unsigned int backgroundColor=0x2f2f2f;
   const unsigned int textColor= 0xb3b3b3;

   m_manager = m;
   CSGAction *openData   = new CSGAction(this, cmsshow::sOpenData.c_str());
   CSGAction *appendData = new CSGAction(this, cmsshow::sAppendData.c_str());
   CSGAction *loadConfig = new CSGAction(this, cmsshow::sLoadConfig.c_str());
   CSGAction *saveConfig   = new CSGAction(this, cmsshow::sSaveConfig.c_str());
   CSGAction *saveConfigAs = new CSGAction(this, cmsshow::sSaveConfigAs.c_str());
   CSGAction *exportImage  = new CSGAction(this, cmsshow::sExportImage.c_str());
   CSGAction *exportImages = new CSGAction(this, cmsshow::sExportAllImages.c_str());
   CSGAction *quit = new CSGAction(this, cmsshow::sQuit.c_str());
   CSGAction *undo = new CSGAction(this, cmsshow::sUndo.c_str());
   undo->disable(); //NOTE: All disables happen again later in this routine
   CSGAction *redo  = new CSGAction(this, cmsshow::sRedo.c_str());
   redo->disable(); //NOTE: All disables happen again later in this routine
   CSGAction *cut   = new CSGAction(this, cmsshow::sCut.c_str());
   cut->disable();  //NOTE: All disables happen again later in this routine
   CSGAction *copy  = new CSGAction(this, cmsshow::sCopy.c_str());
   copy->disable(); //NOTE: All disables happen again later in this routine
   CSGAction *paste = new CSGAction(this, cmsshow::sPaste.c_str());
   paste->disable();//NOTE: All disables happen again later in this routine
   CSGAction *goToFirst = new CSGAction(this, cmsshow::sGotoFirstEvent.c_str());
   CSGAction *goToLast = new CSGAction(this, cmsshow::sGotoLastEvent.c_str());

   CSGAction *showBrightnessInsp = new CSGAction(this, cmsshow::sShowBrightnessInsp.c_str());
   CSGAction *nextEvent          = new CSGAction(this, cmsshow::sNextEvent.c_str());
   CSGAction *previousEvent      = new CSGAction(this, cmsshow::sPreviousEvent.c_str());

   CSGContinuousAction *playEvents     = new CSGContinuousAction(this, cmsshow::sPlayEvents.c_str());
   CSGContinuousAction *playEventsBack = new CSGContinuousAction(this, cmsshow::sPlayEventsBack.c_str());
   CSGContinuousAction *loop           = new CSGContinuousAction(this, cmsshow::sAutoRewind.c_str());
   CSGAction *showObjInsp          = new CSGAction(this, cmsshow::sShowObjInsp.c_str());
   CSGAction *showEventDisplayInsp = new CSGAction(this, cmsshow::sShowEventDisplayInsp.c_str());
   CSGAction *showMainViewCtl      = new CSGAction(this, cmsshow::sShowMainViewCtl.c_str());
   CSGAction *showAddCollection    = new CSGAction(this, cmsshow::sShowAddCollection.c_str());
   CSGAction *help                 = new CSGAction(this, cmsshow::sHelp.c_str());
   CSGAction *searchFiles          = new CSGAction(this, cmsshow::sSearchFiles.c_str());

   CSGAction *keyboardShort      = new CSGAction(this, cmsshow::sKeyboardShort.c_str());
   CSGAction *helpGL             = new CSGAction(this, cmsshow::sHelpGL.c_str());
   CSGAction *colorset           = new CSGAction(this, cmsshow::sBackgroundColor.c_str());

   m_nextEvent = nextEvent;
   m_previousEvent = previousEvent;
   m_goToFirst = goToFirst;
   m_goToLast = goToLast;
   m_playEvents = playEvents;
   m_playEventsBack = playEventsBack;
   m_loopAction = loop;

   goToFirst->setToolTip("Goto first event");
   goToLast->setToolTip("Goto last event");
   previousEvent->setToolTip("Goto previous event");
   nextEvent->setToolTip("Goto next event");
   playEvents->setToolTip("Play events");
   playEventsBack->setToolTip("Play events backwards");

   TGCompositeFrame *menuTopFrame = new TGCompositeFrame(this, 1, 1, kHorizontalFrame, backgroundColor);

   TGMenuBar *menuBar = new TGMenuBar(menuTopFrame, this->GetWidth(), 28, kHorizontalFrame);

   TGPopupMenu *fileMenu = new TGPopupMenu(gClient->GetRoot());
   menuBar->AddPopup("File", fileMenu, new TGLayoutHints(kLHintsTop | kLHintsLeft, 0, 4, 2, 0));
   
   openData->createMenuEntry(fileMenu);
   appendData->createMenuEntry(fileMenu);
   searchFiles->createMenuEntry(fileMenu);
   //searchFiles->disable();
   loadConfig->createMenuEntry(fileMenu);
   saveConfig->createMenuEntry(fileMenu);
   saveConfigAs->createMenuEntry(fileMenu);
   fileMenu->AddSeparator();
    
   exportImage->createMenuEntry(fileMenu);
   exportImages->createMenuEntry(fileMenu);
   fileMenu->AddSeparator();

   quit->createMenuEntry(fileMenu);

   openData->createShortcut(kKey_O, "CTRL", GetId());
   loadConfig->createShortcut(kKey_L, "CTRL", GetId());
   saveConfig->createShortcut(kKey_S, "CTRL", GetId());
   saveConfigAs->createShortcut(kKey_S, "CTRL+SHIFT", GetId());
   exportImage->createShortcut(kKey_P, "CTRL", GetId());
   exportImages->createShortcut(kKey_P, "CTRL+SHIFT", GetId());
   quit->createShortcut(kKey_Q, "CTRL", GetId());

   TGPopupMenu *editMenu = new TGPopupMenu(gClient->GetRoot());
   menuBar->AddPopup("Edit", editMenu, new TGLayoutHints(kLHintsTop | kLHintsLeft, 0, 4, 2, 0));
   undo->createMenuEntry(editMenu);
   undo->createShortcut(kKey_Z, "CTRL", GetId());
   redo->createMenuEntry(editMenu);
   redo->createShortcut(kKey_Z, "CTRL+SHIFT", GetId());
   editMenu->AddSeparator();

   cut->createMenuEntry(editMenu);
   cut->createShortcut(kKey_X, "CTRL", GetId());
   copy->createMenuEntry(editMenu);
   copy->createShortcut(kKey_C, "CTRL", GetId());
   paste->createMenuEntry(editMenu);
   paste->createShortcut(kKey_V, "CTRL", GetId());

   TGPopupMenu *viewMenu = new TGPopupMenu(gClient->GetRoot());
   menuBar->AddPopup("View", viewMenu, new TGLayoutHints(kLHintsTop | kLHintsLeft, 0, 4, 2, 0));  

   m_newViewerMenu = new TGPopupMenu(gClient->GetRoot());
   viewMenu->AddPopup("New Viewer", m_newViewerMenu);
   viewMenu->AddSeparator();

   colorset->createMenuEntry(viewMenu);
   colorset->createShortcut(kKey_B, "CTRL", GetId());
 
   showBrightnessInsp->createMenuEntry(viewMenu);

   viewMenu->AddSeparator();

   nextEvent->createMenuEntry(viewMenu);
   nextEvent->createShortcut(kKey_Right, "CTRL", GetId());
   previousEvent->createMenuEntry(viewMenu);
   previousEvent->createShortcut(kKey_Left, "CTRL", GetId());
   goToFirst->createMenuEntry(viewMenu);
   goToLast->createMenuEntry(viewMenu);
   playEvents->createMenuEntry(viewMenu);
   playEvents->createShortcut(kKey_Space, "CTRL", GetId());
   playEventsBack->createMenuEntry(viewMenu);
   playEventsBack->createShortcut(kKey_Space, "CTRL+SHIFT", GetId());
   loop->createMenuEntry(viewMenu);

   TGPopupMenu* windowMenu = new TGPopupMenu(gClient->GetRoot());
   menuBar->AddPopup("Window", windowMenu, new TGLayoutHints(kLHintsTop | kLHintsLeft, 0, 4, 2, 0));

   showObjInsp->createMenuEntry(windowMenu);
   showObjInsp->createShortcut(kKey_I, "CTRL", GetId());
   showEventDisplayInsp->createMenuEntry(windowMenu);
   showAddCollection->createMenuEntry(windowMenu);
   showMainViewCtl->createMenuEntry(windowMenu);

   TGPopupMenu *helpMenu = new TGPopupMenu(gClient->GetRoot());
   menuBar->AddPopup("Help", helpMenu, new TGLayoutHints(kLHintsTop | kLHintsLeft, 0, 4, 2, 0));
   help->createMenuEntry(helpMenu);
   keyboardShort->createMenuEntry(helpMenu);
   helpGL->createMenuEntry(helpMenu);

   // colors
   menuBar->SetBackgroundColor(backgroundColor);
   TIter next(menuBar->GetTitles());
   TGMenuTitle *title;
   while ((title = (TGMenuTitle *)next()))
      title->SetTextColor(textColor);

   menuTopFrame->AddFrame(menuBar, new TGLayoutHints(kLHintsExpandX, 0, 0, 0, 0));
   AddFrame(menuTopFrame, new TGLayoutHints(kLHintsExpandX, 0, 0, 0, 0));

   // !!!! MT Line separating menu from other window components.
   // I would even remove it and squeeze the navigation buttons up.
   AddFrame(new TGFrame(this, 1, 1, kChildFrame, 0x503020),
            new TGLayoutHints(kLHintsExpandX, 0, 0, 0, 0));

   TGHorizontalFrame *fullbar = new TGHorizontalFrame(this, this->GetWidth(), 30,0, backgroundColor);
   m_statBar = new TGStatusBar(this, this->GetWidth(), 12);
   AddFrame(m_statBar, new TGLayoutHints(kLHintsBottom | kLHintsExpandX));

   /**************************************************************************/
   // controls

   TGCompositeFrame* controlFrame = new TGVerticalFrame(fullbar, 10, 20, 0, backgroundColor);

   TGCompositeFrame* buttonFrame = new TGHorizontalFrame(controlFrame, 10, 10, 0, backgroundColor);
   TImage *imgBtn  = TImage::Open(FWCheckBoxIcon::coreIcondir()+"slider-bg-up.png");
   buttonFrame->SetBackgroundPixmap(imgBtn->GetPixmap());


   goToFirst->createCustomIconsButton(buttonFrame,
                                      fClient->GetPicture(FWCheckBoxIcon::coreIcondir()+"button-gotofirst.png"),
                                      fClient->GetPicture(FWCheckBoxIcon::coreIcondir()+"button-gotofirst-over.png"),
                                      fClient->GetPicture(FWCheckBoxIcon::coreIcondir()+"button-gotofirst-disabled.png"),
                                      new TGLayoutHints(kLHintsCenterY| kLHintsLeft, 4, 3, 10, 0));

   playEventsBack->createCustomIconsButton(buttonFrame,
                                           fClient->GetPicture(FWCheckBoxIcon::coreIcondir()+"button-backward.png"),
                                           fClient->GetPicture(FWCheckBoxIcon::coreIcondir()+"button-backward-over.png"),
                                           fClient->GetPicture(FWCheckBoxIcon::coreIcondir()+"button-backward-disabled.png"),
                                           fClient->GetPicture(FWCheckBoxIcon::coreIcondir()+"button-pause.png"),
                                           fClient->GetPicture(FWCheckBoxIcon::coreIcondir()+"button-pause-over.png"),
                                           new TGLayoutHints(kLHintsCenterY| kLHintsLeft, 2, 3, 10, 0));

   previousEvent->createCustomIconsButton(buttonFrame,
                                          fClient->GetPicture(FWCheckBoxIcon::coreIcondir()+"button-stepback.png"),
                                          fClient->GetPicture(FWCheckBoxIcon::coreIcondir()+"button-stepback-over.png"),
                                          fClient->GetPicture(FWCheckBoxIcon::coreIcondir()+"button-stepback-disabled.png"),
                                          new TGLayoutHints(kLHintsCenterY| kLHintsLeft, 2, 3, 10, 0));

   nextEvent->createCustomIconsButton(buttonFrame,
                                      fClient->GetPicture(FWCheckBoxIcon::coreIcondir()+"button-stepforward.png"),
                                      fClient->GetPicture(FWCheckBoxIcon::coreIcondir()+"button-stepforward-over.png"),
                                      fClient->GetPicture(FWCheckBoxIcon::coreIcondir()+"button-stepforward-disabled.png"),
                                      new TGLayoutHints(kLHintsCenterY| kLHintsLeft, 2, 3, 10, 0));


   playEvents->createCustomIconsButton(buttonFrame,
                                       fClient->GetPicture(FWCheckBoxIcon::coreIcondir()+"button-forward.png"),
                                       fClient->GetPicture(FWCheckBoxIcon::coreIcondir()+"button-forward-over.png"),
                                       fClient->GetPicture(FWCheckBoxIcon::coreIcondir()+"button-forward-disabled.png"),
                                       fClient->GetPicture(FWCheckBoxIcon::coreIcondir()+"button-pause.png"),
                                       fClient->GetPicture(FWCheckBoxIcon::coreIcondir()+"button-pause-over.png"),
                                       new TGLayoutHints(kLHintsCenterY| kLHintsLeft, 2, 3, 10, 0));

   goToLast->createCustomIconsButton(buttonFrame,
                                     fClient->GetPicture(FWCheckBoxIcon::coreIcondir()+"button-gotolast.png"),
                                     fClient->GetPicture(FWCheckBoxIcon::coreIcondir()+"button-gotolast-over.png"),
                                     fClient->GetPicture(FWCheckBoxIcon::coreIcondir()+"button-gotolast-disabled.png"),
                                     new TGLayoutHints(kLHintsCenterY| kLHintsLeft, 2, 3, 10, 0));



   controlFrame->AddFrame(buttonFrame, new TGLayoutHints(kLHintsTop | kLHintsLeft, 10, 0, 0, 0));

   /**************************************************************************/

   TGHorizontalFrame* sliderFrame = new TGHorizontalFrame(controlFrame, 10, 10, 0, backgroundColor);
   TImage *imgSld  = TImage::Open(FWCheckBoxIcon::coreIcondir()+"slider-bg-down.png");
   sliderFrame->SetBackgroundPixmap(imgSld->GetPixmap());
   TString sldBtn = FWCheckBoxIcon::coreIcondir() +"slider-button.png";

   m_delaySlider = new TGHSlider(sliderFrame, 109, kSlider1 | kScaleNo);
   sliderFrame->AddFrame(m_delaySlider, new TGLayoutHints(kLHintsTop | kLHintsLeft, 39, 8, 1, 3));
   m_delaySlider->SetRange(0, 10000);
   m_delaySlider->SetPosition(0);
   m_delaySlider->SetBackgroundColor(0x1a1a1a);
   m_delaySlider->ChangeSliderPic(sldBtn);

   controlFrame->AddFrame(sliderFrame, new TGLayoutHints(kLHintsTop | kLHintsLeft, 10, 0, 0, 0));

   fullbar->AddFrame(controlFrame, new TGLayoutHints(kLHintsLeft, 2, 2, 5, 5));

   m_delaySliderListener =  new FWIntValueListener();
   TQObject::Connect(m_delaySlider, "PositionChanged(Int_t)", "FWIntValueListenerBase",  m_delaySliderListener, "setValue(Int_t)");

   /**************************************************************************/
   // delay label
   TGVerticalFrame* delayFrame = new TGVerticalFrame(fullbar, 60, 10, 0, backgroundColor);
   TGLabel *label = new TGLabel(delayFrame, "Delay");
   label->SetTextJustify(kTextCenterX);
   label->SetTextColor(0xb3b3b3);
   label->SetBackgroundColor(backgroundColor);
   delayFrame->AddFrame(label, new TGLayoutHints(kLHintsTop | kLHintsCenterX, 0, 0, 22, 0));

   TGHorizontalFrame *labFixed = new TGHorizontalFrame(delayFrame, 70, 20, kFixedSize, backgroundColor);
   m_delayLabel = new TGLabel(labFixed, "0.0s");
   m_delayLabel->SetBackgroundColor(backgroundColor);
   m_delayLabel->SetTextJustify(kTextCenterX);
   m_delayLabel->SetTextColor(0xffffff);
   labFixed->AddFrame(m_delayLabel, new TGLayoutHints(kLHintsTop | kLHintsCenterX |kLHintsExpandX, 0, 0, 0, 0));
   delayFrame->AddFrame(labFixed, new TGLayoutHints(kLHintsLeft, 0, 4, 0, 0));

   fullbar->AddFrame(delayFrame, new TGLayoutHints(kLHintsTop | kFixedSize, 0, 0, 0, 0));

   /**************************************************************************/
   // text/num entries

   Int_t maxW =  fullbar->GetWidth() - controlFrame->GetWidth();
   TGVerticalFrame *texts = new TGVerticalFrame(fullbar, 400, 48, kFixedSize, backgroundColor);
   Int_t entryHeight = 20;

   // upper row
   TGHorizontalFrame *runInfo = new TGHorizontalFrame(texts, maxW, entryHeight, 0);
   runInfo->SetBackgroundColor(backgroundColor);
   TGHorizontalFrame *rLeft = new TGHorizontalFrame(runInfo, 200, 20);
   makeFixedSizeLabel(rLeft, "Run", backgroundColor, 0xffffff);
   m_runEntry = new TGNumberEntryField(rLeft, -1, 0, TGNumberFormat::kNESInteger);
   rLeft->AddFrame(m_runEntry, new TGLayoutHints(kLHintsNormal | kLHintsExpandX, 0,0,0,0));
   runInfo->AddFrame(rLeft, new TGLayoutHints(kLHintsLeft));

   TGHorizontalFrame *rRight = new TGHorizontalFrame(runInfo, 200, 20);
   makeFixedSizeLabel(rRight, "Event", backgroundColor, 0xffffff);
   m_eventEntry = new TGNumberEntryField(rRight, -1, 0, TGNumberFormat::kNESInteger);
   rRight->AddFrame(m_eventEntry, new TGLayoutHints(kLHintsNormal | kLHintsExpandX, 0,2,0,0));

   runInfo->AddFrame(rRight, new TGLayoutHints(kLHintsRight));

   texts->AddFrame(runInfo, new TGLayoutHints(kLHintsNormal | kLHintsExpandX, 0,0,0,1));

   // lower row
   TGHorizontalFrame *filterFrame = new TGHorizontalFrame(texts, maxW, entryHeight, 0, backgroundColor);
   
   TGCompositeFrame *lframe = new TGHorizontalFrame(filterFrame, 50, entryHeight, kFixedSize, backgroundColor);

   // filter state Off
   m_filterIcons[0] = fClient->GetPicture("unchecked_t.xpm");
   m_filterIcons[1] = fClient->GetPicture("unchecked_t.xpm");
   m_filterIcons[2] = fClient->GetPicture("unchecked_dis_t.xpm");
   
   // filter state On
   m_filterIcons[3] = fClient->GetPicture("checked_t.xpm");
   m_filterIcons[4] = fClient->GetPicture("checked_t.xpm");
   m_filterIcons[5] = fClient->GetPicture("checked_dis_t.xpm");
   
   // filter withdrawn
   m_filterIcons[6] = fClient->GetPicture(FWCheckBoxIcon::coreIcondir() + "icon-alert-ltgraybg.png");
   m_filterIcons[7] = fClient->GetPicture(FWCheckBoxIcon::coreIcondir() + "icon-alert-ltgraybg-over.png");
   m_filterIcons[8] = fClient->GetPicture(FWCheckBoxIcon::coreIcondir() + "icon-alert-ltgraybg.png");
   
   m_filterEnableBtn = new FWCustomIconsButton(lframe, m_filterIcons[0], m_filterIcons[1], m_filterIcons[2]);
   m_filterEnableBtn->SetBackgroundColor(backgroundColor);
   m_filterEnableBtn->SetToolTipText("Enable/disable event filtering");
   lframe->AddFrame(m_filterEnableBtn, new TGLayoutHints(kLHintsRight | kLHintsCenterY,0,0,2,2));
   filterFrame->AddFrame(lframe, new TGLayoutHints(kLHintsLeft|kLHintsCenterY,0,0,0,0));

   m_filterShowGUIBtn = new TGTextButton(filterFrame,"Event filtering is OFF");
   m_filterShowGUIBtn->SetBackgroundColor(backgroundColor);
   m_filterShowGUIBtn->SetTextColor(0xFFFFFF);
   m_filterShowGUIBtn->SetToolTipText("Edit filters");
   filterFrame->AddFrame(m_filterShowGUIBtn,new TGLayoutHints(kLHintsExpandX|kLHintsLeft|kLHintsTop,4,1,2,2));

   texts->AddFrame(filterFrame, new TGLayoutHints(kLHintsNormal | kLHintsExpandX, 0,0,1,0));
   fullbar->AddFrame(texts, new TGLayoutHints(kLHintsNormal| kLHintsCenterY, 20, 5, 5, 5));

   /**************************************************************************/
   TGVerticalFrame *texts2 = new TGVerticalFrame(fullbar, 200, 44, kFixedSize, backgroundColor);

   // Lumi
   m_lumiBlock = new TGLabel(texts2, "Lumi block id: ");
   m_lumiBlock->SetTextJustify(kTextLeft);
   m_lumiBlock->SetTextColor(0xffffff);
   m_lumiBlock->SetBackgroundColor(backgroundColor);
   texts2->AddFrame(m_lumiBlock, new TGLayoutHints(kLHintsNormal | kLHintsExpandX| kLHintsBottom, 0,0,3,1));

   // time
   m_timeText = new TGLabel(texts2, "...");
   m_timeText->SetTextJustify(kTextLeft);
   m_timeText->SetTextColor(0xffffff);
   m_timeText->SetBackgroundColor(backgroundColor);
   texts2->AddFrame(m_timeText, new TGLayoutHints(kLHintsNormal | kLHintsExpandX| kLHintsBottom, 0,0,0,1));

   fullbar->AddFrame(texts2, new TGLayoutHints(kLHintsNormal| kLHintsCenterY, 6, 5, 5, 5));

   /**************************************************************************/
   //  logo
   {
      TGVerticalFrame* parentLogoFrame = new TGVerticalFrame(fullbar, 70, 53, kFixedSize); 
      parentLogoFrame->SetBackgroundColor(backgroundColor);
      fullbar->AddFrame(parentLogoFrame, new TGLayoutHints(kLHintsRight| kLHintsCenterY, 0, 0, 0 ));

      TGVerticalFrame* logoFrame = new TGVerticalFrame(parentLogoFrame, 53, 53, kFixedSize);
      TImage *logoImg  = TImage::Open(FWCheckBoxIcon::coreIcondir() + "CMSRedOnBlackThick.png");
      logoFrame->SetBackgroundPixmap(logoImg->GetPixmap());
      parentLogoFrame->AddFrame(logoFrame, new TGLayoutHints(kLHintsRight | kLHintsCenterY, 0, 17, 0, 0));
   }
   {
      TGCompositeFrame *logoFrame = new TGCompositeFrame(this, 61, 23, kFixedSize | kHorizontalFrame, backgroundColor);
      FWCustomIconsButton *infoBut =
         new FWCustomIconsButton(logoFrame, fClient->GetPicture(FWCheckBoxIcon::coreIcondir()+"fireworksSmallGray.png"),
                                            fClient->GetPicture(FWCheckBoxIcon::coreIcondir()+"fireworksSmallGray-green.png"),
                                            fClient->GetPicture(FWCheckBoxIcon::coreIcondir()+"fireworksSmallGray-red.png"),
                                            fClient->GetPicture(FWCheckBoxIcon::coreIcondir()+"fireworksSmallGray-red.png"));
      logoFrame->AddFrame(infoBut);
      infoBut->Connect("Clicked()", "CmsShowMainFrame", this, "showFWorksInfo()");
      //TImage *logoImg  = TImage::Open( FWCheckBoxIcon::coreIcondir() + "fireworksSmallGray.png");
      //logoFrame->SetBackgroundPixmap(logoImg->GetPixmap());
      menuTopFrame->AddFrame(logoFrame, new TGLayoutHints(kLHintsRight | kLHintsBottom, 0, 13, 3, 1));
   }
  
   /**************************************************************************/
   AddFrame(fullbar, new TGLayoutHints(kLHintsExpandX, 0, 0, 0, 0));

   //Start disabled
   goToFirst->disable();
   goToLast->disable();
   previousEvent->disable();
   nextEvent->disable();
   playEvents->disable();
   playEventsBack->disable();
   loop->disable();
   
   //NOTE: There appears to be a bug in ROOT such that creating a menu item and setting it as
   // disabled immediately is ignored.  Therefore we have to wait till here to actually get ROOT
   // to disable these menu items
   undo->disable();
   redo->disable();
   cut->disable();
   copy->disable();
   paste->disable();
   
   
   TGSplitFrame *csArea = new TGSplitFrame(this, this->GetWidth(), this->GetHeight()-42);
   csArea->VSplit(200);
   csArea->GetFirst()->AddFrame(m_manager->createList(csArea->GetFirst()), new TGLayoutHints(kLHintsLeft | kLHintsExpandX | kLHintsExpandY));
   TGTab *tabFrame = new TGTab(csArea->GetSecond(), csArea->GetSecond()->GetWidth(), csArea->GetSecond()->GetHeight());

   m_manager->createViews(tabFrame);

   csArea->GetSecond()->AddFrame(tabFrame, new TGLayoutHints(kLHintsLeft | kLHintsExpandX | kLHintsExpandY));
   AddFrame(csArea,new TGLayoutHints(kLHintsTop | kLHintsExpandX | kLHintsExpandY,2,2,0,2));
   SetWindowName("cmsShow");
}

// CmsShowMainFrame::CmsShowMainFrame(const CmsShowMainFrame& rhs)
// {
//    // do actual copying here;
// }

CmsShowMainFrame::~CmsShowMainFrame() {
   Cleanup();
}

//
// assignment operators
//
// const CmsShowMainFrame& CmsShowMainFrame::operator=(const CmsShowMainFrame& rhs)
// {
//   //An exception safe implementation is
//   CmsShowMainFrame temp(rhs);
//   swap(rhs);
//
//   return *this;
// }

//
// member functions
//

CSGAction*
CmsShowMainFrame::createNewViewerAction(const std::string& iActionName)
{
   CSGAction* action(new CSGAction(this, iActionName.c_str()));
   action->createMenuEntry(m_newViewerMenu);
   return action;
}



void CmsShowMainFrame::loadEvent(const fwlite::Event& event) {

   if (event.id().run() != static_cast<unsigned int>(m_runEntry->GetIntNumber()))
      m_runEntry->SetIntNumber(event.id().run());

   if (event.id().event() != static_cast<unsigned int>(m_eventEntry->GetIntNumber()))
      m_eventEntry->SetIntNumber(event.id().event());

   m_timeText->SetText( fw::getLocalTime( event ).c_str() );
   char title[128];
   snprintf(title,128,"Lumi block id: %d", event.eventAuxiliary().luminosityBlock());
   m_lumiBlock->SetText( title );
}

void CmsShowMainFrame::enableNavigatorControls()
{
   m_nextEvent->enable();
   m_previousEvent->enable();
   m_goToFirst->enable();
   m_goToLast->enable();
   m_playEvents->enable();
   m_playEventsBack->enable();
   m_loopAction->enable();
}

void  CmsShowMainFrame::CloseWindow()
{
   getAction(cmsshow::sQuit)->activated();
}

void CmsShowMainFrame::quit() {
   getAction(cmsshow::sQuit)->activated();
}

void
CmsShowMainFrame::enableActions(bool enable)
{
   CSGActionSupervisor::enableActions(enable);

   m_runEntry->SetEnabled(enable);
   m_eventEntry->SetEnabled(enable);
   m_filterEnableBtn->SetEnabled(enable);
   m_filterShowGUIBtn->SetEnabled(enable);
}

void
CmsShowMainFrame::enablePrevious(bool enable)
{
   if (m_previousEvent != 0) {
      if (enable) {
         m_previousEvent->enable();
         m_goToFirst->enable();
         m_playEventsBack->enable();
      } else {
         m_previousEvent->disable();
         m_goToFirst->disable();
         m_playEventsBack->disable();
         m_playEventsBack->stop();
      }
   }
}

void
CmsShowMainFrame::enableNext(bool enable)
{
   if (m_nextEvent != 0) {
      if (enable) {
         m_nextEvent->enable();
         m_goToLast->enable();
         m_playEvents->enable();
      } else {
         m_nextEvent->disable();
         m_goToLast->disable();
         m_playEvents->disable();
         m_playEvents->stop();
      }
   }
}

bool
CmsShowMainFrame::nextIsEnabled()
{
   return m_nextEvent->isEnabled();
}

bool
CmsShowMainFrame::previousIsEnabled()
{
   return m_previousEvent->isEnabled();
}

void CmsShowMainFrame::updateStatusBar(const char* status) {
   m_statBar->SetText(status, 0);
   //force the status bar to update its image
   gClient->ProcessEventsFor(m_statBar);
}

void CmsShowMainFrame::clearStatusBar()
{
   m_statBar->SetText("", 0);
   //don't process immediately since we want this on the event queue
   // since results of the last action may still be happening
}

void CmsShowMainFrame::HandleMenu(Int_t id) {
   switch(id) {
      case 1:
      {
         gApplication->Terminate(0);
      }
      break;
      default:
         fwLog(fwlog::kInfo) << "Invalid menu id\n";
         break;
   }
}

Bool_t CmsShowMainFrame::HandleKey(Event_t *event) {
   if (event->fType == kGKeyPress) {
      const std::vector<CSGAction*>& alist = getListOfActions();
      std::vector<CSGAction*>::const_iterator it_act;
      Int_t keycode;
      Int_t modcode;
      for (it_act = alist.begin(); it_act != alist.end(); ++it_act) {
         keycode = (*it_act)->getKeycode();
         modcode = (*it_act)->getModcode();
         if ((event->fCode == (UInt_t)keycode) &&
             ((event->fState == (UInt_t)modcode) ||
              (event->fState == (UInt_t)(modcode | kKeyMod2Mask)) ||
              (event->fState == (UInt_t)(modcode | kKeyLockMask)) ||
              (event->fState == (UInt_t)(modcode | kKeyMod2Mask | kKeyLockMask)))) {
            (*it_act)->activated.emit();
            return kTRUE;
         }
      }
   }
   return kFALSE;
}

void
CmsShowMainFrame::setPlayDelayGUI(Float_t val, Bool_t sliderChanged)
{
   m_delayLabel->SetText(Form("%.1fs", val));
   if (sliderChanged)
      m_delaySlider->SetPosition(Int_t(val*1000));
}

void
CmsShowMainFrame::makeFixedSizeLabel(TGHorizontalFrame* p, const char* txt, UInt_t bgCol,  UInt_t txtCol)
{
   // Utility function.

   Int_t labW = 50;
   Int_t labH = 20;

   p->SetBackgroundColor(bgCol);
   TGCompositeFrame *lframe = new TGHorizontalFrame(p, labW, labH, kFixedSize, bgCol);
   TGLabel* label = new TGLabel(lframe, txt);
   label->SetBackgroundColor(bgCol);
   label->SetTextColor(txtCol);
   lframe->AddFrame(label,     new TGLayoutHints(kLHintsRight | kLHintsBottom));
   p->AddFrame(lframe, new TGLayoutHints(kLHintsLeft  | kLHintsBottom, 0, 4, 0, 0));
}

class InfoFrame : public TGMainFrame {
public:
   InfoFrame(const TGWindow* p, UInt_t w, UInt_t h, UInt_t opts) : TGMainFrame(p, w, h, opts) {}
   virtual ~InfoFrame() {}
   
   virtual void CloseWindow()
   {
      UnmapWindow();  
   }
};

void
CmsShowMainFrame::showFWorksInfo()
{
   if (m_fworksAbout == 0)
   {
      TString infoFileName("$(CMSSW_BASE)/src/Fireworks/Core/data/version.txt");
      gSystem->ExpandPathName(infoFileName);
      
      ifstream infoFile(infoFileName);
      TString infoText;
      infoText.ReadFile(infoFile);
      infoFile.close();
 
      const UInt_t ww = 280, hh = 180;
      
      m_fworksAbout = new InfoFrame(gClient->GetRoot(), ww, hh, kVerticalFrame | kFixedSize);
      m_fworksAbout->SetWMSizeHints(ww, hh, ww, hh, 0, 0);
      m_fworksAbout->SetBackgroundColor(0x2f2f2f);
      
      TGFrame* logoFrame = new TGFrame(m_fworksAbout, 140, 48, kFixedSize);
      TImage *logoImg  = TImage::Open(FWCheckBoxIcon::coreIcondir()+"logo-fireworks.png");
      logoFrame->SetBackgroundPixmap(logoImg->GetPixmap());
      m_fworksAbout->AddFrame(logoFrame, new TGLayoutHints(kLHintsTop | kLHintsCenterX, 0, 0, 16, 0));
      
      TGLabel* label = new TGLabel(m_fworksAbout, infoText);
      label->SetBackgroundColor(0x2f2f2f);
      label->SetForegroundColor(0xffffff);
      m_fworksAbout->AddFrame(label, new TGLayoutHints(kLHintsCenterX | kLHintsCenterY));      
            
      TGTextButton* btn = new TGTextButton(m_fworksAbout, "  OK  ");
      btn->SetBackgroundColor(0x2f2f2f);
      btn->SetForegroundColor(0xffffff);
      m_fworksAbout->AddFrame(btn, new TGLayoutHints(kLHintsBottom | kLHintsCenterX, 0, 0, 0, 12));
      btn->Connect("Clicked()", "TGMainFrame", m_fworksAbout, "CloseWindow()");

      m_fworksAbout->MapSubwindows();
      m_fworksAbout->Layout();
   }
   
   m_fworksAbout->MapRaised();
}


