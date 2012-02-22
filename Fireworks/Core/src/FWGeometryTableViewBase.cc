#include <iostream>

#include <boost/bind.hpp>
#include <boost/regex.hpp>

#include "Fireworks/Core/interface/FWGeometryTableViewBase.h"
#include "Fireworks/Core/interface/FWGeoTopNode.h"
#include "Fireworks/Core/interface/fwLog.h"
#include "Fireworks/Core/interface/FWGeometryTableManagerBase.h"
#include "Fireworks/TableWidget/interface/FWTableWidget.h"
#include "Fireworks/Core/interface/FWGUIManager.h"
#include "Fireworks/Core/interface/FWColorManager.h"
#include "Fireworks/Core/interface/FWParameterSetterBase.h"
#include "Fireworks/Core/src/FWColorSelect.h"
#include "Fireworks/Core/src/FWPopupMenu.cc"
#include "Fireworks/Core/src/FWGeoTopNodeScene.h"

#include "TGFileDialog.h"
#include "TGeoNode.h"
#include "TGeoMatrix.h"
#include "TGStatusBar.h"
#include "TGButton.h"
#include "TGLabel.h"
#include "TGLPhysicalShape.h"
#include "TGMenu.h"
#include "TGComboBox.h"
// #define PERFTOOL_BROWSER
#include "TEvePointSet.h"
#include "TGeoShape.h"
#include "TGeoBBox.h"
#include "TEveManager.h"
#include "TEveGeoNode.h"
#include "TEveScene.h"
#include "TEveSceneInfo.h"
#include "TEveViewer.h"
#include "TGLViewer.h"
#include "TGLCamera.h"
#include "TEveSelection.h"
#ifdef PERFTOOL_BROWSER 
#include <google/profiler.h>
#endif

//______________________________________________________________________________
//==============================================================================
//==============================================================================
//==============================================================================
//==============================================================================


Bool_t FWGeometryTableViewBase::FWViewCombo::HandleButton(Event_t* event) 
{
   if (event->fType == kButtonPress)
   {
      bool map = false;

      FWPopupMenu* m_viewPopup = new FWPopupMenu(0);

      TEveElementList* views = gEve->GetViewers();
      int idx = 0;

      for (TEveElement::List_i it = views->BeginChildren(); it != views->EndChildren(); ++it)
      { 
         TEveViewer* v = ((TEveViewer*)(*it));
         if (strstr( v->GetElementName(), "3D") )
         {     
            bool added = false;          
            m_viewPopup->AddEntry(v->GetElementName(), idx);

            for (TEveElement::List_i eit = v->BeginChildren(); eit != v->EndChildren(); ++eit )
            {
               TEveScene* s = ((TEveSceneInfo*)*eit)->GetScene();
               if (m_el && s->HasChildren() && s->FirstChild() == m_el) {
                  added = true;
                  break;
               }
            }
            map = true;
            if (added)
               m_viewPopup->CheckEntry(idx);
         }
         ++idx;
      }

      if (map) {

         Window_t wdummy;
         Int_t ax,ay;
         gVirtualX->TranslateCoordinates(GetId(),
                                         gClient->GetDefaultRoot()->GetId(),
                                         event->fX, event->fY, //0,0 in local coordinates
                                         ax,ay, //coordinates of screen
                                         wdummy);


         m_viewPopup->PlaceMenu(ax, ay, true,true);
         m_viewPopup->Connect("Activated(Int_t)",
                              "FWGeometryTableViewBase",
                              const_cast<FWGeometryTableViewBase*>(m_tableView),
                              "selectView(Int_t)");
      }
      else
      {
         fwLog(fwlog::kInfo) << "No 3D View added. \n";
      }
   }
   return true;
}
//==============================================================================
//==============================================================================
//==============================================================================
//==============================================================================
//==============================================================================


FWGeometryTableViewBase::FWGeometryTableViewBase(TEveWindowSlot* iParent,FWViewType::EType type, FWColorManager* colMng )
   : FWViewBase(type),
     m_topNodeIdx(this, "TopNodeIndex", -1l, 0, 1e7),
     m_autoExpand(this,"ExpandList:", 1l, 0l, 100l),
     m_enableHighlight(this,"EnableHighlight", true),
     m_colorManager(colMng),
     m_colorPopup(0),
     m_eveWindow(0),
     m_frame(0),
     m_viewBox(0),
     m_viewersConfig(0),
     m_enableRedraw(true),
     m_marker(0),
     m_eveTopNode(0),
     m_eveScene(0),
     m_tableRowIndexForColorPopup(-1)
{
   m_eveWindow = iParent->MakeFrame(0);
   TGCompositeFrame* xf = m_eveWindow->GetGUICompositeFrame();

   m_frame = new TGVerticalFrame(xf);


   xf->AddFrame(m_frame, new TGLayoutHints(kLHintsExpandX | kLHintsExpandY));   
}

void FWGeometryTableViewBase::postConst()
{
   m_tableWidget = new FWTableWidget(getTableManager(), m_frame); 
   m_frame->AddFrame(m_tableWidget,new TGLayoutHints(kLHintsExpandX|kLHintsExpandY,2,2,0,0));
   m_tableWidget->SetBackgroundColor(0xffffff);
   m_tableWidget->SetLineSeparatorColor(0x000000);
   m_tableWidget->SetHeaderBackgroundColor(0xececec);
   m_tableWidget->Connect("cellClicked(Int_t,Int_t,Int_t,Int_t,Int_t,Int_t)",
                          "FWGeometryTableViewBase",this,
                          "cellClicked(Int_t,Int_t,Int_t,Int_t,Int_t,Int_t)");
   m_tableWidget->disableGrowInWidth();
   //   resetSetters();


   m_frame->MapSubwindows();
   m_frame->Layout();
  m_eveWindow->GetGUICompositeFrame()->Layout();
   m_frame->MapWindow();
}
//______________________________________________________________________________

FWGeometryTableViewBase::~FWGeometryTableViewBase()
{
   // take out composite frame and delete it directly (zwithout the timeout)
   TGCompositeFrame *frame = m_eveWindow->GetGUICompositeFrame();
   frame->RemoveFrame( m_frame );
   delete m_frame;



   m_eveWindow->DestroyWindowAndSlot();
   delete getTableManager();
}


namespace {
TEveScene* getMarkerScene(TEveViewer* v)
{
  TEveElement* si = v->FindChild(Form("SI - GeoScene %s", v->GetElementName()));
  if(si) 
    return ((TEveSceneInfo*)(si))->GetScene();
  else
    return 0;
}
}
//==============================================================================


void FWGeometryTableViewBase::cdNode(int idx)
{
   std::string p;
   getTableManager()->getNodePath(idx, p);
   setPath(idx, p);
}

void FWGeometryTableViewBase::cdTop()
{
   std::string path = "/" ;
   path += getTableManager()->refEntries().at(0).name();
   setPath(-1, path ); 
}

void FWGeometryTableViewBase::cdUp()
{   
   if (getTopNodeIdx() != -1)
   {
      int pIdx = getTableManager()->refEntries()[getTopNodeIdx()].m_parent;
      std::string p;
      getTableManager()->getNodePath(pIdx, p);
      setPath(pIdx, p);
   }
}

void FWGeometryTableViewBase::setPath(int parentIdx, std::string&)
{
   m_eveTopNode->clearSelection();

   printf("set path %d \n", parentIdx);
   m_topNodeIdx.set(parentIdx);
   // getTableManager()->refEntries().at(getTopNodeIdx()).setBitVal(FWGeometryTableManagerBase::kVisNodeSelf,!m_disableTopNode.value() );
   getTableManager()->setLevelOffset(getTableManager()->refEntries().at(getTopNodeIdx()).m_level);
 

   checkExpandLevel();
   refreshTable3D(); 
}

//------------------------------------------------------------------------------

void  FWGeometryTableViewBase::checkExpandLevel()
{
   // check expand state
   int ae = m_autoExpand.value();
   if ( m_topNodeIdx.value() > 0) 
      ae += getTableManager()->refEntries().at(m_topNodeIdx.value()).m_level;

   for (FWGeometryTableManagerBase::Entries_i i = getTableManager()->refEntries().begin(); i !=  getTableManager()->refEntries().end(); ++i)
   {
      if (i->m_level  < ae)
         i->setBit(FWGeometryTableManagerBase::kExpanded);
      else
         i->resetBit(FWGeometryTableManagerBase::kExpanded);
   } 
}

//==============================================================================

void
FWGeometryTableViewBase::populate3DViewsFromConfig()
{
   // post-config 
   if (m_viewersConfig) {
      TEveElementList* viewers = gEve->GetViewers();
      const FWConfiguration::KeyValues* keyVals = m_viewersConfig->keyValues();

      if(0!=keyVals)  
      {
         for(FWConfiguration::KeyValuesIt it = keyVals->begin(); it!= keyVals->end(); ++it) {
    
            TString sname = it->first;
            TEveViewer* v = dynamic_cast<TEveViewer*>(viewers->FindChild(sname.Data()));
            if (!v)
            {
               fwLog(fwlog::kError)  << "FWGeometryTableViewBase::populate3DViewsFromConfig no viewer found " << it->first << std::endl;
               return;
            }
            v->AddScene(m_eveScene);  
            m_viewBox->setElement(m_eveTopNode);
            if (m_marker) getMarkerScene(v)->AddElement(m_marker);

            gEve->FullRedraw3D(false, true);
         }   
      }
   }
}

//==============================================================================

void 
FWGeometryTableViewBase::selectView(int idx)
{
   // callback from sleclect view popup menu

   m_viewBox->setElement(m_eveTopNode);

   TEveElement::List_i it = gEve->GetViewers()->BeginChildren();
   std::advance(it, idx);
   TEveViewer* v = (TEveViewer*)(*it);

   for (TEveElement::List_i eit = v->BeginChildren(); eit != v->EndChildren(); ++eit )
   {
      if ((((TEveSceneInfo*)(*eit))->GetScene()) == m_eveScene)
      {
        v->RemoveElement(*eit);
        if (m_marker) getMarkerScene(v)->RemoveElement(m_marker);
        gEve->Redraw3D();
        return;
      }
   }

   if (m_marker) getMarkerScene(v)->AddElement(m_marker); 
   v->AddScene(m_eveScene);
   gEve->Redraw3D();
}

//==============================================================================

void 
FWGeometryTableViewBase::cellClicked(Int_t iRow, Int_t iColumn, Int_t iButton, Int_t iKeyMod, Int_t x, Int_t y)
{
   // getTableManager()->setSelection(iRow, iColumn, iButton);
   int idx = getTableManager()->rowToIndex()[iRow];
   FWGeometryTableManagerBase::NodeInfo& ni = getTableManager()->refEntries()[idx];

   if (iButton == kButton1) 
   {
      if (iColumn == 0)
      {
         Window_t wdummy;
         Int_t xLoc,yLoc;
         gVirtualX->TranslateCoordinates(gClient->GetDefaultRoot()->GetId(), m_tableWidget->GetId(),  x, y, xLoc, yLoc, wdummy);

         bool sel = getTableManager()->firstColumnClicked(iRow, xLoc);

         if (sel) {
            int idx =getTableManager()->rowToIndex()[iRow];
	          // printf("cell clicled top node %p\n", (void*)m_eveTopNode);
            if (gEve->GetSelection()->HasChild( m_eveTopNode))
               gEve->GetSelection()->RemoveElement( m_eveTopNode);

            if (gEve->GetHighlight()->HasChild( m_eveTopNode))
               gEve->GetHighlight()->RemoveElement( m_eveTopNode);

            // reset bits and sets for old selected table entry
            m_eveTopNode->UnSelected();
            m_eveTopNode->UnHighlighted();


            if (m_eveTopNode->selectPhysicalFromTable(idx))
               gEve->GetSelection()->AddElement(m_eveTopNode);

            ni.setBit(FWGeometryTableManagerBase::kSelected);
            getTableManager()->redrawTable();
            gEve->Redraw3D();
	 }
      }
      else if (iColumn == 1)
      { 
         std::vector<Color_t> colors;
         m_colorManager->fillLimitedColors(colors);
      
         if (!m_colorPopup) {
            m_colorPopup = new FWColorPopup(gClient->GetDefaultRoot(), colors.front());
            m_colorPopup->InitContent("", colors);
            m_colorPopup->Connect("ColorSelected(Color_t)","FWGeometryTableViewBase", const_cast<FWGeometryTableViewBase*>(this), "nodeColorChangeRequested(Color_t");
         }
         m_tableRowIndexForColorPopup = idx;
         m_colorPopup->SetName("Selected");
         m_colorPopup->ResetColors(colors, m_colorManager->backgroundColorIndex()==FWColorManager::kBlackIndex);
         m_colorPopup->PlacePopup(x, y, m_colorPopup->GetDefaultWidth(), m_colorPopup->GetDefaultHeight());
         return;
      }
      else
      {
         bool elementChanged = false;
         if (iColumn == 2)
         {
            ni.switchBit(FWGeometryTableManagerBase::kVisNodeSelf);
            elementChanged = true;
         }
         else if (iColumn == 3)
         { 
            ni.switchBit(FWGeometryTableManagerBase::kVisNodeChld); 
            elementChanged = true;
         }
         else if (iColumn == 5)
         {
            // used in overlaps for RnrMarker column
            ni.switchBit(BIT(5));
            elementChanged = true;
         }

         if (elementChanged)
         {
            refreshTable3D();
         }
      }

      getTableManager()->dataChanged();

   }
   else if (iColumn == 0)
   {
      m_eveTopNode->popupMenu(x, y);
   }
}


void FWGeometryTableViewBase::setBackgroundColor()
{
   bool backgroundIsWhite = m_colorManager->backgroundColorIndex()==FWColorManager::kWhiteIndex;
   if(backgroundIsWhite) {
      m_tableWidget->SetBackgroundColor(0xffffff);
      m_tableWidget->SetLineSeparatorColor(0x000000);
   } else {
      m_tableWidget->SetBackgroundColor(0x000000);
      m_tableWidget->SetLineSeparatorColor(0xffffff);
   }
   getTableManager()->setBackgroundToWhite(backgroundIsWhite);
   gClient->NeedRedraw(m_tableWidget);
}

//______________________________________________________________________________

void FWGeometryTableViewBase::nodeColorChangeRequested(Color_t col)
{
   //   printf("color change %d \n", m_tableRowIndexForColorPopup);
   if (m_tableRowIndexForColorPopup >= 0) {
      FWGeometryTableManagerBase::NodeInfo& ni = getTableManager()->refEntries()[m_tableRowIndexForColorPopup];
      ni.m_color = col;
      ni.m_node->GetVolume()->SetLineColor(col);
      refreshTable3D();
      m_tableRowIndexForColorPopup = -1;
   }
}

//______________________________________________________________________________

void FWGeometryTableViewBase::refreshTable3D()
{
   if (gEve->GetSelection()->HasChild(m_eveTopNode))
      gEve->GetSelection()->RemoveElement(m_eveTopNode);

   if (gEve->GetHighlight()->HasChild(m_eveTopNode))
      gEve->GetHighlight()->RemoveElement(m_eveTopNode);

   m_eveTopNode->fSceneJebo->PadPaint(m_eveTopNode->fSceneJebo->GetPad());
   gEve->Redraw3D(); 

   getTableManager()->redrawTable();
}

//______________________________________________________________________________

void FWGeometryTableViewBase::addTo(FWConfiguration& iTo) const
{
   FWConfigurableParameterizable::addTo(iTo);

   FWConfiguration viewers(1);
   FWConfiguration tempArea;

   for(TEveElement::List_i k = gEve->GetViewers()->BeginChildren(); k!= gEve->GetViewers()->EndChildren(); ++k)
   {
      for (TEveElement::List_i eit = (*k)->BeginChildren(); eit != (*k)->EndChildren(); ++eit )
      {
         TEveScene* s = ((TEveSceneInfo*)*eit)->GetScene();
         if (s->GetGLScene() == m_eveTopNode->fSceneJebo)
         {
            viewers.addKeyValue( (*k)->GetElementName(), tempArea);
            break;
         }
      }
   }

   iTo.addKeyValue("Viewers", viewers, true);
}
  
//______________________________________________________________________________

void FWGeometryTableViewBase::setFrom(const FWConfiguration& iFrom)
{ 
   m_enableRedraw = false;
   for(const_iterator it =begin(), itEnd = end();
       it != itEnd;
       ++it) {

      //      printf("set from %s \n",(*it)->name().c_str() );
      (*it)->setFrom(iFrom);

   }  
   m_viewersConfig = iFrom.valueForKey("Viewers");

   m_enableRedraw = true;
   refreshTable3D();
}
