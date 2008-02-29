// -*- C++ -*-
//
// Package:     Core
// Class  :     FWGUIManager
// 
// Implementation:
//     <Notes on implementation>
//
// Original Author:  Chris Jones
//         Created:  Mon Feb 11 11:06:40 EST 2008
// $Id: FWGUIManager.cc,v 1.6 2008/02/24 20:39:05 dmytro Exp $
//

// system include files
#include <boost/bind.hpp>
#include <stdexcept>
#include <iostream>

#include "TGButton.h"
#include "TGComboBox.h"
#include "TGTextEntry.h"
#include "TApplication.h"
#include "TROOT.h"
#include "TSystem.h"
#include "TGSplitFrame.h"
#include "TGTab.h"
#include "TGListTree.h"
//EVIL, no accessor for the editor yet
//#define protected public
#include "TEveBrowser.h"
//#undef protected
#include "TEveManager.h"
#include "TEveGedEditor.h"

// user include files
#include "Fireworks/Core/interface/FWGUIManager.h"
#include "Fireworks/Core/interface/FWGUISubviewArea.h"

#include "Fireworks/Core/interface/FWSelectionManager.h"
#include "Fireworks/Core/interface/FWModelExpressionSelector.h"
#include "Fireworks/Core/interface/FWEventItemsManager.h"
#include "Fireworks/Core/interface/FWViewBase.h"

#include "Fireworks/Core/src/FWListEventItem.h"

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
FWGUIManager::FWGUIManager(FWSelectionManager* iSelMgr,
                           FWEventItemsManager* iEIMgr,
                           bool iDebugInterface
):
m_selectionManager(iSelMgr),
m_eiManager(iEIMgr),
m_continueProcessingEvents(false),
m_waitForUserAction(true),
m_code(0)
{
   m_selectionManager->selectionChanged_.connect(boost::bind(&FWGUIManager::selectionChanged,this,_1));
   m_eiManager->newItem_.connect(boost::bind(&FWGUIManager::newItem,
                                             this, _1) );

   // These are only needed temporarilty to work around a problem which 
   // Matevz has patched in a later version of the code
   TApplication::NeedGraphicsLibs();
   gApplication->InitializeGraphics();
   
   TEveManager::Create();
   TEveBrowser* browser = gEve->GetBrowser();
   // TGFrame* f = (TGFrame*) gClient->GetDefaultRoot();
   // browser->MoveResize(f->GetX(), f->GetY(), f->GetWidth(), f->GetHeight());
   // browser->Resize( gClient->GetDisplayWidth(), gClient->GetDisplayHeight() );
   
   //should check to see if already has our tab
   {
      browser->StartEmbedding(TRootBrowser::kLeft);
      {
         TGMainFrame* frmMain=new TGMainFrame(gClient->GetRoot(),
                                              1000,
                                              600);
         frmMain->SetWindowName("GUI");
         frmMain->SetCleanup(kDeepCleanup);
         
         TGHorizontalFrame* hf = new TGHorizontalFrame(frmMain);
         //We need an error handling system which can properly report
         // errors and decide what to do
         // given that we are an interactive system we need to leave
         // the code in a good state so that users can decided to 
         // continue or not
         {
            if(0==gSystem->Getenv("ROOTSYS")) {
               std::cerr<<"environment variable ROOTSYS is not set" <<
               std::endl;
               throw std::runtime_error("ROOTSYS environment variable not set");
            }
            TString icondir(Form("%s/icons/",gSystem->Getenv("ROOTSYS")));
            
            //m_homeButton= new TGPictureButton(hf, gClient->GetPicture(icondir+"GoHome.gif"));
            m_homeButton= new TGPictureButton(hf, gClient->GetPicture(icondir+"first_t.xpm"));
            const unsigned int kButtonSize = 30;
            m_homeButton->SetToolTipText("Go back to first event");
            m_homeButton->SetMinHeight(kButtonSize);
            m_homeButton->SetMinWidth(kButtonSize);
            m_homeButton->SetHeight(kButtonSize);
            m_homeButton->SetWidth(kButtonSize);
            hf->AddFrame(m_homeButton);
            m_homeButton->Connect("Clicked()", "FWGUIManager", this, "goHome()");
            
            
            //m_backwardButton= new TGPictureButton(hf, gClient->GetPicture(icondir+"GoBack.gif"));
            m_backwardButton= new TGPictureButton(hf, gClient->GetPicture(icondir+"previous_t.xpm"));
            m_backwardButton->SetToolTipText("Go back one event");
            m_backwardButton->SetMinHeight(kButtonSize);
            m_backwardButton->SetMinWidth(kButtonSize);
            m_backwardButton->SetHeight(kButtonSize);
            m_backwardButton->SetWidth(kButtonSize);
            hf->AddFrame(m_backwardButton);
            m_backwardButton->Connect("Clicked()", "FWGUIManager", this, "goBack()");
            
            //m_advanceButton= new TGPictureButton(hf, gClient->GetPicture(icondir+"GoForward.gif"));
            m_advanceButton= new TGPictureButton(hf, gClient->GetPicture(icondir+"next_t.xpm"));
            m_advanceButton->SetToolTipText("Go to next event");
            const unsigned int kExpand = 10;
            m_advanceButton->SetMinHeight(kButtonSize+kExpand);
            m_advanceButton->SetMinWidth(kButtonSize+kExpand);
            m_advanceButton->SetHeight(kButtonSize+kExpand);
            m_advanceButton->SetWidth(kButtonSize+kExpand);
            hf->AddFrame(m_advanceButton);
            m_advanceButton->Connect("Clicked()", "FWGUIManager", this, "goForward()");
            
            //m_stopButton= new TGPictureButton(hf, gClient->GetPicture(icondir+"StopLoading.gif"));
            m_stopButton= new TGPictureButton(hf, gClient->GetPicture(icondir+"stop_t.xpm"));
            m_stopButton->SetToolTipText("Stop looping over event");
            m_stopButton->SetMinHeight(kButtonSize);
            m_stopButton->SetMinWidth(kButtonSize);
            m_stopButton->SetHeight(kButtonSize);
            m_stopButton->SetWidth(kButtonSize);
            hf->AddFrame(m_stopButton);
            m_stopButton->Connect("Clicked()", "FWGUIManager", this, "stop()");
            
         }
         frmMain->AddFrame(hf);
         //frmMain->SetEditable();
         TEveGListTreeEditorFrame* ltf = new TEveGListTreeEditorFrame(frmMain);
         //frmMain->SetEditable(kFALSE);
         frmMain->AddFrame(ltf);
         m_listTree = ltf->GetListTree();
         m_eventObjects =  new TEveElementList("Physics Objects");
         m_listTree->OpenItem(m_eventObjects->AddIntoListTree(m_listTree,
                                                              reinterpret_cast<TGListTreeItem*>(0))
                              );
         m_views =  new TEveElementList("Views");
         m_views->AddIntoListTree(m_listTree,reinterpret_cast<TGListTreeItem*>(0));
         m_editor = ltf->GetEditor();
         m_editor->DisplayElement(0);

         {
            //m_listTree->Connect("mouseOver(TGListTreeItem*, UInt_t)", "FWGUIManager",
              //                 this, "itemBelowMouse(TGListTreeItem*, UInt_t)");
            m_listTree->Connect("Clicked(TGListTreeItem*, Int_t, UInt_t, Int_t, Int_t)", "FWGUIManager",
                               this, "itemClicked(TGListTreeItem*, Int_t, UInt_t, Int_t, Int_t)");
            m_listTree->Connect("DoubleClicked(TGListTreeItem*, Int_t)", "FWGUIManager",
                               this, "itemDblClicked(TGListTreeItem*, Int_t)");
            m_listTree->Connect("KeyPressed(TGListTreeItem*, ULong_t, ULong_t)", "FWGUIManager",
                               this, "itemKeyPress(TGListTreeItem*, UInt_t, UInt_t)");
         }
         
         TGGroupFrame* vf = new TGGroupFrame(frmMain,"Selection",kVerticalFrame);
         {
            
            TGGroupFrame* vf2 = new TGGroupFrame(vf,"Expression");
            m_selectionItemsComboBox = new TGComboBox(vf2,200);
            m_selectionItemsComboBox->Resize(200,20);
            vf2->AddFrame(m_selectionItemsComboBox, new TGLayoutHints(kLHintsTop | kLHintsLeft,0,5,5,5));
            m_selectionExpressionEntry = new TGTextEntry(vf2,"$.pt() > 10");
            vf2->AddFrame(m_selectionExpressionEntry, new TGLayoutHints(kLHintsExpandX,0,5,5,5));
            m_selectionRunExpressionButton = new TGTextButton(vf2,"Select by Expression");
            vf2->AddFrame(m_selectionRunExpressionButton);
            m_selectionRunExpressionButton->Connect("Clicked()","FWGUIManager",this,"selectByExpression()");
            vf->AddFrame(vf2);
            
            m_unselectAllButton = new TGTextButton(vf,"Unselect All");
            m_unselectAllButton->Connect("Clicked()", "FWGUIManager",this,"unselectAll()");
            vf->AddFrame(m_unselectAllButton);
            m_unselectAllButton->SetEnabled(kFALSE);
          
         }
         frmMain->AddFrame(vf);
         frmMain->MapSubwindows();
         frmMain->Resize();
         frmMain->MapWindow();
      }
      browser->StopEmbedding();
      browser->SetTabTitle("Fireworks",TRootBrowser::kLeft);
   }
   {
      //pickup our other icons
      const char* cmspath = gSystem->Getenv("CMSSW_BASE");
      if(0 == cmspath) {
         throw std::runtime_error("CMSSW_BASE environment variable not set");
      }
      TString coreIcondir(Form("%s/src/Fireworks/Core/icons/",gSystem->Getenv("CMSSW_BASE")));
      
      
      
      browser->StartEmbedding(TRootBrowser::kRight);
      {
         m_mainFrame = new TGMainFrame(gClient->GetRoot(),600,450);
         m_splitFrame = new TGSplitFrame(m_mainFrame, 800, 600);
         m_mainFrame->AddFrame(m_splitFrame, new TGLayoutHints(kLHintsExpandX | kLHintsExpandY));
         // split it once
         m_splitFrame->HSplit(434);
         // then split each part again (this will make four parts)
         m_splitFrame->GetSecond()->VSplit(400);

         TGSplitFrame* sf = m_splitFrame->GetFirst();
         m_viewFrames.push_back(sf);

         sf = m_splitFrame->GetSecond()->GetFirst();
         TGCompositeFrame* hf = new FWGUISubviewArea(sf,m_splitFrame);
         m_viewFrames.push_back(hf);
         (sf)->AddFrame(hf,new TGLayoutHints(kLHintsExpandX | 
                                             kLHintsExpandY) );

         
         sf=m_splitFrame->GetSecond()->GetSecond();
         hf = new FWGUISubviewArea(sf,m_splitFrame);
         m_viewFrames.push_back(hf);
         (sf)->AddFrame(hf,new TGLayoutHints(kLHintsExpandX | 
                                             kLHintsExpandY) );
         m_nextFrame = m_viewFrames.begin();

         m_mainFrame->MapSubwindows();
         m_mainFrame->Resize();
         m_mainFrame->MapWindow();
         
      }
      browser->StopEmbedding();
      browser->SetTabTitle("Views",TRootBrowser::kRight);
   }
   if(not iDebugInterface) {
      browser->GetTabLeft()->RemoveTab(0);
      browser->GetTabLeft()->RemoveTab(0);
      browser->GetTabRight()->RemoveTab(0);
   }
   
}

// FWGUIManager::FWGUIManager(const FWGUIManager& rhs)
// {
//    // do actual copying here;
// }

FWGUIManager::~FWGUIManager()
{
}

//
// assignment operators
//
// const FWGUIManager& FWGUIManager::operator=(const FWGUIManager& rhs)
// {
//   //An exception safe implementation is
//   FWGUIManager temp(rhs);
//   swap(rhs);
//
//   return *this;
// }

//
// member functions
//
void 
FWGUIManager::addFrameHoldingAView(TGFrame* iChild)
{
   (*m_nextFrame)->AddFrame(iChild,new TGLayoutHints(kLHintsExpandX | 
                                                     kLHintsExpandY) );
   
   m_mainFrame->MapSubwindows();
   m_mainFrame->Resize();
   iChild->Resize();
   m_mainFrame->MapWindow();
   
   ++m_nextFrame;
}

TGFrame* 
FWGUIManager::parentForNextView()
{
   assert(m_nextFrame != m_viewFrames.end());
   return *m_nextFrame;
}


void 
FWGUIManager::registerViewBuilder(const std::string& iName, 
                                  ViewBuildFunctor& iBuilder)
{
   m_nameToViewBuilder[iName]=iBuilder;
}

void 
FWGUIManager::createView(const std::string& iName)
{
   NameToViewBuilder::iterator itFind = m_nameToViewBuilder.find(iName);
   if(itFind == m_nameToViewBuilder.end()) {
      throw std::runtime_error(std::string("Unable to create view named ")+iName+" because it is unknown");
   }
   addFrameHoldingAView((itFind->second(parentForNextView()))->frame());
   
   TEveElementList* lst = new TEveElementList(iName.c_str(),"");
   lst->AddIntoListTree(m_listTree,m_views);
   
}



void
FWGUIManager::goForward()
{
   m_continueProcessingEvents = true;
   m_code = 1;
}

void
FWGUIManager::goBack()
{
   m_continueProcessingEvents = true;
   m_code = -1;
}

void
FWGUIManager::goHome()
{
   m_continueProcessingEvents = true;
   m_code = -2;
}

void
FWGUIManager::stop()
{
   m_continueProcessingEvents = true;
   m_code = -3;
}

void
FWGUIManager::waitForUserAction()
{
   m_waitForUserAction = true;
}

void
FWGUIManager::doNotWaitForUserAction()
{
   m_waitForUserAction = false;
}

void 
FWGUIManager::selectByExpression()
{
   FWModelExpressionSelector selector;
   selector.select(*(m_eiManager->begin()+m_selectionItemsComboBox->GetSelected()),
                   m_selectionExpressionEntry->GetText());
}

void 
FWGUIManager::unselectAll()
{
   m_selectionManager->clearSelection();
}

void 
FWGUIManager::selectionChanged(const FWSelectionManager& iSM)
{
   m_unselectAllButton->SetEnabled( 0 !=iSM.selected().size() );
}

void 
FWGUIManager::processGUIEvents()
{
   gSystem->ProcessEvents();
}

void
FWGUIManager::newItem(const FWEventItem* iItem)
{
   m_selectionItemsComboBox->AddEntry(iItem->name().c_str(),iItem->id());
   if(iItem->id()==0) {
      m_selectionItemsComboBox->Select(0);
   }
   //TEveElementList* lst = new TEveElementList(iItem->name().c_str(),"",kTRUE);
   //lst->SetMainColor(iItem->defaultDisplayProperties().color());
   //NEED TO CHANGE THE SIGNATURE OF THE SIGNAL
   TEveElementList* lst = new FWListEventItem( const_cast<FWEventItem*>(iItem) );
   lst->AddIntoListTree(m_listTree,m_eventObjects);
}


bool
FWGUIManager::waitingForUserAction() const
{
   return m_waitForUserAction;
}

//
// const member functions
//

namespace {
   //guarantee that no matter how we go back to Cint that
   // we have disabled these buttons
   struct EnableButton {
      EnableButton( TGButton* iButton):
      m_button(iButton)
      {
         if(0!=m_button) {
            m_button->SetEnabled();
         }
      }
      ~EnableButton()
      {
         m_button->SetEnabled(kFALSE);
         gSystem->DispatchOneEvent(kFALSE);
      }
      
   private:
      TGButton* m_button;
   };
   
}

int
FWGUIManager::allowInteraction()
{
   //need to reset
   m_continueProcessingEvents = false;
   EnableButton homeB(m_homeButton);
   EnableButton advancedB(m_advanceButton);
   EnableButton backwardB(m_backwardButton);
   EnableButton stopB(m_stopButton);
   //Unselect all doesn't need this since the selection manager will 
   // properly update this button
   //EnableButton stopUnselect(m_unselectAllButton);
   EnableButton stopSelect(m_selectionRunExpressionButton);
   
   //m_viewManager->newEventAvailable();
   
   //check for input at least once
   gSystem->ProcessEvents();
   while(not gROOT->IsInterrupted() and
         m_waitForUserAction and 
         not m_continueProcessingEvents) {
      // gSystem->ProcessEvents();
      gSystem->DispatchOneEvent(kFALSE);
   }
   return m_code;
}


void 
FWGUIManager::itemChecked(TObject* obj, Bool_t state)
{
}
void 
FWGUIManager::itemClicked(TGListTreeItem *item, Int_t btn,  UInt_t mask, Int_t x, Int_t y)
{
   TEveElement* el = static_cast<TEveElement*>(item->GetUserData());
   m_editor->DisplayElement(el);
}
void 
FWGUIManager::itemDblClicked(TGListTreeItem* item, Int_t btn)
{
}
void 
FWGUIManager::itemKeyPress(TGListTreeItem *entry, UInt_t keysym, UInt_t mask)
{
}

void 
FWGUIManager::itemBelowMouse(TGListTreeItem* item, UInt_t)
{
}

//
// static member functions
//
