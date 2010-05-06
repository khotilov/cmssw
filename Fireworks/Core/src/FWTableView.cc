// -*- C++ -*-
//
// Package:     Core
// Class  :     FWTableView
//
// Implementation:
//     <Notes on implementation>
//
// Original Author:  Chris Jones
//         Created:  Thu Feb 21 11:22:41 EST 2008
// $Id: FWTableView.cc,v 1.20 2010/04/22 17:29:52 amraktad Exp $
//

// system include files
#include <stdlib.h>
#include <algorithm>
#include <boost/bind.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/numeric/conversion/converter.hpp>
#include <iostream>
#include <sstream>
#include <stdexcept>

#include "TRootEmbeddedCanvas.h"
#include "THStack.h"
#include "TCanvas.h"
#include "TClass.h"
#include "TH2F.h"
#include "TView.h"
#include "TColor.h"
#include "TEveScene.h"
#include "TGLViewer.h"
#include "TSystem.h"
#include "TGComboBox.h"
#include "TGLabel.h"
#include "TGTextView.h"
#include "TGTextEntry.h"
#include "TEveViewer.h"
#include "TEveManager.h"
#include "TEveWindow.h"
#include "TEveElement.h"
#include "TEveCalo.h"
#include "TEveElement.h"
#include "TEveLegoEventHandler.h"
#include "TGLWidget.h"
#include "TGLScenePad.h"
#include "TGLFontManager.h"
#include "TEveTrans.h"
#include "TGeoTube.h"
#include "TEveGeoNode.h"
#include "TEveStraightLineSet.h"
#include "TEveText.h"
#include "TGeoArb8.h"

// user include files
#include "Fireworks/Core/interface/FWColorManager.h"
#include "Fireworks/Core/interface/FWCustomIconsButton.h"
#include "Fireworks/Core/interface/FWModelChangeManager.h"
#include "Fireworks/Core/interface/FWSelectionManager.h"
#include "Fireworks/Core/interface/FWTableView.h"
#include "Fireworks/Core/interface/FWTableViewManager.h"
#include "Fireworks/Core/interface/FWEventItem.h"
#include "Fireworks/Core/interface/FWEveValueScaler.h"
#include "Fireworks/Core/interface/FWConfiguration.h"
#include "Fireworks/Core/interface/BuilderUtils.h"
#include "Fireworks/Core/interface/FWExpressionEvaluator.h"
#include "Fireworks/Core/interface/FWTableViewTableManager.h"
#include "Fireworks/Core/interface/fwLog.h"
#include "Fireworks/Core/src/FWGUIValidatingTextEntry.h"
#include "Fireworks/Core/src/FWExpressionValidator.h"
#include "Fireworks/TableWidget/interface/FWTableWidget.h"

static const TString& coreIcondir() 
{
   static TString path = Form("%s/src/Fireworks/Core/icons/",gSystem->Getenv("CMSSW_BASE"));
   if ( gSystem->AccessPathName(path.Data()) ){ // cannot find directory
	assert(gSystem->Getenv("CMSSW_RELEASE_BASE"));
	path = Form("%s/src/Fireworks/Core/icons/",gSystem->Getenv("CMSSW_RELEASE_BASE"));
   }
   return path;
}

/*
static 
const TGPicture* filtered(bool iBackgroundIsBlack)
{
   if(iBackgroundIsBlack) {
      static const TGPicture* s = gClient->GetPicture(coreIcondir()+"filtered-blackbg.png");
      return s;
   }
   static const TGPicture* s = gClient->GetPicture(coreIcondir()+"filtered-whitebg.png");
   return s;
   
}

static 
const TGPicture* filtered_over(bool iBackgroundIsBlack)
{
   if(iBackgroundIsBlack) {
      static const TGPicture* s = gClient->GetPicture(coreIcondir()+"filtered-whitebg-over.png");
      return s;
   }
   static const TGPicture* s = gClient->GetPicture(coreIcondir()+"filtered-whitebg-over.png");
   return s;
}
*/
/*
static 
const TGPicture* alert_over()
{
   static const TGPicture* s = gClient->GetPicture(coreIcondir()+"alert-blackbg-over.png");
   return s;
}

static 
const TGPicture* alert()
{
   static const TGPicture* s = gClient->GetPicture(coreIcondir()+"alert-blackbg.png");
   return s;
}
*/

/*
static 
const TGPicture* unfiltered(bool iBackgroundIsBlack)
{
   if(iBackgroundIsBlack) {
      static const TGPicture* s = gClient->GetPicture(coreIcondir()+"unfiltered-blackbg.png");
      return s;
   }
   static const TGPicture* s = gClient->GetPicture(coreIcondir()+"unfiltered-whitebg.png");
   return s;
}
static 
const TGPicture* unfiltered_over(bool iBackgroundIsBlack)
{
   if(iBackgroundIsBlack) {
      static const TGPicture* s = gClient->GetPicture(coreIcondir()+"unfiltered-blackbg-over.png");
      return s;
   }
   static const TGPicture* s = gClient->GetPicture(coreIcondir()+"unfiltered-whitebg-over.png");
   return s;   
}

static
const TGPicture* info(bool iBackgroundIsBlack)
{
   if(iBackgroundIsBlack) {
      static const TGPicture* s = gClient->GetPicture(coreIcondir()+"info2-blackbg.png");
      return s;
   }
   static const TGPicture* s = gClient->GetPicture(coreIcondir()+"info2-whitebg.png");
   return s;   
}

static
const TGPicture* info_over(bool iBackgroundIsBlack)
{
   if(iBackgroundIsBlack) {
      static const TGPicture* s = gClient->GetPicture(coreIcondir()+"info2-blackbg-over.png");
      return s;
   }
   static const TGPicture* s = gClient->GetPicture(coreIcondir()+"info2-whitebg-over.png");
   return s;
}

static
const TGPicture* info_disabled(bool iBackgroundIsBlack)
{
   if(iBackgroundIsBlack) {
      static const TGPicture* s = gClient->GetPicture(coreIcondir()+"info2-blackbg-disabled.png");
      return s;
   }
   static const TGPicture* s = gClient->GetPicture(coreIcondir()+"info2-whitebg-disabled.png");
   return s;
}
*/
static
const TGPicture* arrow_right(bool iBackgroundIsBlack)
{
   if(iBackgroundIsBlack) {
      static const TGPicture* s = gClient->GetPicture(coreIcondir()+"arrow-white-right-blackbg.png");
      return s;
   }
   static const TGPicture* s = gClient->GetPicture(coreIcondir()+"arrow-black-right-whitebg.png");
   return s;
}

static
const TGPicture* arrow_right_disabled(bool iBackgroundIsBlack)
{
   if(iBackgroundIsBlack) {
      static const TGPicture* s = gClient->GetPicture(coreIcondir()+"arrow-white-right-disabled-blackbg.png");
      return s;
   }
   static const TGPicture* s = gClient->GetPicture(coreIcondir()+"arrow-black-right-disabled-whitebg.png");
   return s;
}

static
const TGPicture* arrow_down(bool iBackgroundIsBlack)
{
   if(iBackgroundIsBlack) {
      static const TGPicture* s = gClient->GetPicture(coreIcondir()+"arrow-white-down-blackbg.png");
      return s;
   }
   static const TGPicture* s = gClient->GetPicture(coreIcondir()+"arrow-black-down-whitebg.png");
   return s;
}

static
const TGPicture* arrow_down_disabled(bool iBackgroundIsBlack)
{
   if(iBackgroundIsBlack) {
      static const TGPicture* s = gClient->GetPicture(coreIcondir()+"arrow-white-down-disabled-blackbg.png");
      return s;
   }
   static const TGPicture* s = gClient->GetPicture(coreIcondir()+"arrow-black-down-disabled-whitebg.png");
   return s;
}

//
// constants, enums and typedefs
//
static const std::string kTableView = "TableView";
static const std::string kCollection = "collection";
static const std::string kColumns = "columns";
static const std::string kSortColumn = "sortColumn";
static const std::string kDescendingSort = "descendingSort";

//
// constructors and destructor
//
FWTableView::FWTableView (TEveWindowSlot* iParent, FWTableViewManager *manager)
     : m_iColl(-1),
       m_manager(manager),
       m_tableManager(new FWTableViewTableManager(this)),
       m_tableWidget(0),
       m_showColumnUI(false),
       m_validator(new FWExpressionValidator),
       m_currentColumn(-1),
       m_useColumnsFromConfig(false)

{
     m_eveWindow = iParent->MakeFrame(0);
     TGCompositeFrame *frame = m_eveWindow->GetGUICompositeFrame();
//      TGHorizontalFrame *buttons = new TGHorizontalFrame(frame);
//      frame->AddFrame(buttons, new TGLayoutHints(kLHintsTop | kLHintsExpandX));

//      m_collection = new TGComboBox(buttons);
     m_vert = new TGVerticalFrame(frame);
     frame->AddFrame(m_vert, new TGLayoutHints(kLHintsExpandX | kLHintsExpandY));
     TGHorizontalFrame *header = new TGHorizontalFrame(m_vert);
     m_vert->AddFrame(header, new TGLayoutHints(kLHintsTop | kLHintsExpandX));
     const bool bgIsBlack = m_manager->colorManager().background() == kBlack;
     m_columnUIButton = new FWCustomIconsButton(header, 
						arrow_right(bgIsBlack),
						arrow_right_disabled(bgIsBlack),
						arrow_right_disabled(bgIsBlack));
     m_columnUIButton->Connect("Clicked()", "FWTableView", this, "toggleShowHide()");
     header->AddFrame(m_columnUIButton, new TGLayoutHints(kLHintsCenterY | kLHintsLeft,6,10));

     TGCompositeFrame *labfr = new TGHorizontalFrame(header, 60, 25, kFixedSize);
     TGLabel *label = new TGLabel(labfr, "Collection");
     labfr->AddFrame(label,  new TGLayoutHints(kLHintsLeft | kLHintsCenterY, 1,3,0,0));
     header->AddFrame(labfr, new TGLayoutHints(kLHintsLeft));

     m_collection = new TGComboBox(header);
     updateItems();
     header->AddFrame(m_collection, new TGLayoutHints(kLHintsLeft | kLHintsExpandX | kLHintsExpandY));
     m_collection->Connect("Selected(Int_t)", "FWTableView", this, "selectCollection(Int_t)");
     m_collection->Select(2, true);
     m_column_control = new TGVerticalFrame(m_vert);
     m_vert->AddFrame(m_column_control, new TGLayoutHints(kLHintsExpandX));
     TGLabel *column_control_label = new TGLabel(m_column_control, "Column editor");
//      column_control_label->SetBackgroundColor(bgIsBlack ? kBlack : kWhite);
//      column_control_label->SetForegroundColor(bgIsBlack ? kWhite : kBlack);
//      column_control_label->SetTextColor(bgIsBlack ? kWhite : kBlack);
     m_column_control->AddFrame(column_control_label, new TGLayoutHints(kLHintsExpandX));
     TGHorizontalFrame *column_control_fields = new TGHorizontalFrame(m_column_control);
     m_column_control->AddFrame(column_control_fields, new TGLayoutHints(kLHintsExpandX));
     m_column_name_field = new TGTextEntry(column_control_fields);
     m_column_name_field->SetMaxWidth(10);
     m_column_expr_field = new FWGUIValidatingTextEntry(column_control_fields);
//      m_column_expr_field->SetEnabled(kFALSE);
     m_column_expr_field->setValidator(m_validator);
     m_column_prec_field = new TGTextEntry(column_control_fields);
     m_column_prec_field->SetMaxWidth(10);
     TGLabel *name_label = new TGLabel(column_control_fields, "Title");
     TGLabel *expr_label = new TGLabel(column_control_fields, "Expression");
     TGLabel *prec_label = new TGLabel(column_control_fields, "Precision");
     column_control_fields->AddFrame(name_label, new TGLayoutHints(kLHintsBottom, 1, 1, 2, 2));
     column_control_fields->AddFrame(m_column_name_field, new TGLayoutHints(kLHintsExpandX));
     column_control_fields->AddFrame(expr_label, new TGLayoutHints(kLHintsBottom, 1, 1, 2, 2));
     column_control_fields->AddFrame(m_column_expr_field, new TGLayoutHints(kLHintsExpandX));
     column_control_fields->AddFrame(prec_label, new TGLayoutHints( kLHintsBottom, 1, 1, 2, 2)); 
     column_control_fields->AddFrame(m_column_prec_field, new TGLayoutHints(kLHintsExpandX));
     TGTextButton *add_button = new TGTextButton(column_control_fields, "Add");
     TGTextButton *del_button = new TGTextButton(column_control_fields, "Delete");
     TGTextButton *mod_button = new TGTextButton(column_control_fields, "Modify");
     add_button->Connect("Clicked()", "FWTableView", this, "addColumn()");
     del_button->Connect("Clicked()", "FWTableView", this, "deleteColumn()");
     mod_button->Connect("Clicked()", "FWTableView", this, "modifyColumn()");
     column_control_fields->AddFrame(add_button, new TGLayoutHints);
     column_control_fields->AddFrame(del_button, new TGLayoutHints);
     column_control_fields->AddFrame(mod_button, new TGLayoutHints);
     m_tableWidget = new FWTableWidget(m_tableManager, m_vert);
     resetColors(m_manager->colorManager());
     m_tableWidget->SetHeaderBackgroundColor(gVirtualX->GetPixel(kWhite));
     m_tableWidget->Connect("rowClicked(Int_t,Int_t,Int_t,Int_t,Int_t)", "FWTableView",
			    this, "modelSelected(Int_t,Int_t,Int_t,Int_t,Int_t)");
     m_tableWidget->Connect("columnClicked(Int_t,Int_t,Int_t)", "FWTableView",
			    this, "columnSelected(Int_t,Int_t,Int_t)");
     m_vert->AddFrame(m_tableWidget, new TGLayoutHints(kLHintsExpandX | kLHintsExpandY));
     frame->MapSubwindows();
     m_vert->HideFrame(m_column_control);
     frame->Layout();
     frame->MapWindow();
}

FWTableView::~FWTableView()
{
     // take out composite frame and delete it directly ( without the timeout) 
     TGCompositeFrame *frame = m_eveWindow->GetGUICompositeFrame();
     frame->RemoveFrame(m_vert);
     delete m_vert;

     m_eveWindow->DestroyWindowAndSlot();
     delete m_tableManager;
     delete m_validator;
}

void
FWTableView::setBackgroundColor(Color_t iColor) 
{
     m_tableWidget->SetBackgroundColor(gVirtualX->GetPixel(iColor));
//      m_tableWidget->SetBackgroundColor(TColor::Number2Pixel(iColor));
//    m_viewer->GetGLViewer()->SetClearColor(iColor);
}

void FWTableView::resetColors (const FWColorManager &manager)
{
     m_tableWidget->SetBackgroundColor(gVirtualX->GetPixel(manager.background()));
//      m_tableWidget->SetHeaderBackgroundColor(gVirtualX->GetPixel(manager.background()));
//      switch (manager.foreground()) {
//      case FWColorManager::kBlackIndex:
// 	  m_tableWidget->SetHeaderForegroundColor(gVirtualX->GetPixel(kBlack));
// 	  break;
//      default:
// 	  m_tableWidget->SetHeaderForegroundColor(0xffffff);
// 	  break;
//      }
     m_tableWidget->SetLineSeparatorColor(gVirtualX->GetPixel(manager.foreground()));
//      m_tableWidget->dataChanged();
}

//
// const member functions
//
TGFrame*
FWTableView::frame() const
{
     return 0;
//    return m_embeddedViewer->GetFrame();
}

const std::string&
FWTableView::typeName() const
{
   return staticTypeName();
}

void
FWTableView::addTo(FWConfiguration& iTo) const
{
     // are we the first FWTableView to go into the configuration?  If
     // we are, then we are responsible for writing out the list of
     // types (which we do by letting FWTableViewManager::addToImpl
     // write into our configuration)
     if (this == m_manager->m_views.front().get())
	  m_manager->addToImpl(iTo);
     // then there is the stuff we have to do anyway: remember what
     // collection we display
     FWConfiguration main(1);
     const std::string &collectionName = m_manager->items()[m_iColl]->name();
     FWConfiguration collection(collectionName);
     main.addKeyValue(kCollection, collection);
     FWConfiguration sortColumn(m_tableWidget->sortedColumn());
     main.addKeyValue(kSortColumn, sortColumn);
     FWConfiguration descendingSort(m_tableWidget->descendingSort());
     main.addKeyValue(kDescendingSort, descendingSort);
//      FWConfiguration columns(1);
//      for (std::vector<FWTableViewManager::TableEntry>::const_iterator 
// 	       i = m_tableManager->m_tableFormats->begin(),
// 	       iEnd = m_tableManager->m_tableFormats->end();
// 	  i != iEnd; ++i) {
// 	  columns.addValue(i->name);
// 	  columns.addValue(i->expression);
// 	  char prec[100];
// 	  snprintf(prec, 100, "%d", i->precision);
// 	  columns.addValue(prec);
//      }
//      main.addKeyValue(kColumns, columns);
     iTo.addKeyValue(kTableView, main);
     // take care of parameters
     FWConfigurableParameterizable::addTo(iTo);
}

void
FWTableView::setFrom(const FWConfiguration& iFrom)
{
     if (this == m_manager->m_views.front().get())
	  m_manager->setFrom(iFrom);
     try {
	  const FWConfiguration *main = iFrom.valueForKey(kTableView);
	  assert(main != 0);
	  // use the columns from the config, not the default columns for
	  // the collection type
//  	  m_useColumnsFromConfig = true;
// 	  m_tableManager->m_tableFormats->clear();
// 	  const FWConfiguration *columns = main->valueForKey(kColumns);
// 	  for (FWConfiguration::StringValuesIt it = columns->stringValues()->begin(),
// 		    itEnd = columns->stringValues()->end(); it != itEnd; ++it) {
// 	       const std::string &name = *it++;
// 	       const std::string &expr = *it++;
// 	       int prec = atoi(it->c_str());
// 	       FWTableViewManager::TableEntry e = { expr, name, prec };
// 	       m_tableManager->m_tableFormats->push_back(e);
// 	  }
	  const FWConfiguration *collection = main->valueForKey(kCollection);
	  const std::string &collectionName = collection->value();
	  // find item 
	  for (std::vector<const FWEventItem *>::const_iterator 
		    it = m_manager->items().begin(), 
		    itEnd = m_manager->items().end();
	       it != itEnd; ++it) {
	       if (*it && (*it)->name() == collectionName) {
		    m_collection->Select(it - m_manager->items().begin(), true);
		    break;
	       }
	  }
	  const FWConfiguration *sortColumn = main->valueForKey(kSortColumn);
	  const FWConfiguration *descendingSort = main->valueForKey(kDescendingSort);
	  if (sortColumn != 0 && descendingSort != 0) {
	       unsigned int sort = sortColumn->version();
	       bool descending = descendingSort->version();
	       if (sort < (( unsigned int) m_tableManager->numberOfColumns()))
		    m_tableWidget->sort(sort, descending);
	  }
     } catch (...) {
	  // configuration doesn't contain info for the table.  Be forgiving.
	  std::cerr << "This configuration file contains tables, but no column information.  "
	       "(It is probably old.)  Using defaults." << std::endl;
     }

//      main.addKeyValue(kCollection, collection);
//      FWConfiguration columns(1);
//      for (std::vector<FWTableViewManager::TableEntry>::const_iterator 
// 	       i = m_tableManager->m_tableFormats->begin(),
// 	       iEnd = m_tableManager->m_tableFormats->end();
// 	  i != iEnd; ++i) {
// 	  columns.addValue(i->name);
// 	  columns.addValue(i->expression);
// 	  columns.addValue(Form("%d", i->precision));
//      }
//      main.addKeyValue(kColumns, columns);
//      iTo.addKeyValue(kTableView, main);
//      // take care of parameters
//      FWConfigurableParameterizable::addTo(iTo);

     // take care of parameters
     FWConfigurableParameterizable::setFrom(iFrom);
}

void
FWTableView::saveImageTo(const std::string& iName) const
{
//    bool succeeded = m_viewer->GetGLViewer()->SavePicture(iName.c_str());
//    if(!succeeded) {
//       throw std::runtime_error("Unable to save picture");
//    }
}

void
FWTableView::toggleShowHide () 
{
     m_showColumnUI = not m_showColumnUI;
     const TGPicture* picture = 0;
     const TGPicture* down = 0;
     const TGPicture* disabled = 0;
     const bool bgIsBlack = m_manager->colorManager().background() == kBlack;
     if (m_showColumnUI) {
	  picture = arrow_down(bgIsBlack);
	  down = arrow_down_disabled(bgIsBlack);
	  disabled = arrow_down_disabled(bgIsBlack);
	  m_vert->ShowFrame(m_column_control);
     } else {
	  picture = arrow_right(bgIsBlack);
	  down = arrow_right_disabled(bgIsBlack);
	  disabled = arrow_right_disabled(bgIsBlack);
	  m_vert->HideFrame(m_column_control);
     }
     m_vert->Layout();
     m_columnUIButton->swapIcons(picture,down,disabled);
}

void FWTableView::updateItems ()
{
     int selected = m_collection->GetSelected();
     m_collection->RemoveAll();
     int index =0;
     for (std::vector<const FWEventItem *>::const_iterator it = m_manager->items().begin(), 
	       itEnd = m_manager->items().end();
	  it != itEnd; ++it,++index) {
        if(*it) {
           m_collection->AddEntry((*it)->name().c_str(), it - m_manager->items().begin());
        }
        if(m_iColl == index && 0 == *it) {
           //the collection we were showing is now gone
           m_iColl = -1;
           selected = -1;
        }
     }
     if (selected != -1 && selected < m_collection->GetNumberOfEntries())
	  m_collection->Select(selected, false);
}

void FWTableView::updateEvaluators ()
{
     m_tableManager->updateEvaluators();
}

const FWEventItem *FWTableView::item () const
{
     if (m_iColl == -1)
	  return 0;
     return m_manager->items()[m_iColl];
}

void FWTableView::dataChanged ()
{
//      const FWEventItem *item = m_manager->items()[m_iColl];
     updateEvaluators();
     m_tableManager->dataChanged();
//      std::vector<FWExpressionEvaluator> &ev = m_evaluators;
//      for (unsigned int i = 0; i < item->size(); ++i) {
// 	  for (unsigned int j = 0; j < ev.size(); ++j) {
// 	       printf("%s = %f\t", (*m_manager->tableFormats(item->modelType()->GetName())).second[j].name.c_str(),
// 		      ev[j].evalExpression(item->modelData(i)));
// 	  }
// 	  printf("\n");
//      }
//      fflush(stdout);
}


void FWTableView::selectCollection (Int_t i_coll)
{
//      printf("selected collection %d, ", i_coll);
     const FWEventItem *item = m_manager->items()[i_coll];
     assert(0!=item);
//      printf("%s\n", item->modelType()->GetName());
     m_iColl = i_coll;
//      m_validator = new FWExpressionValidator;
//      m_column_expr_field->setValidator(m_validator);
     if (m_validator != 0) {
// 	  std::cout << "setting validator to " << item->modelType()->GetName() << std::endl;
	  m_validator->setType(ROOT::Reflex::Type::ByTypeInfo(*(item->modelType()->GetTypeInfo())));
     } else {
// 	  std::cout << "can't set null validator\n";
     }
     if (not m_useColumnsFromConfig) {
	  if (m_manager->tableFormats(*item->modelType()) == m_manager->m_tableFormats.end()) {
               fwLog(fwlog::kInfo) << "No table format for objects of this type " << item->modelType()->GetName() << std::endl;
	       m_tableManager->m_tableFormats->clear();
	  } else {
	       m_tableManager->m_tableFormats = &m_manager->tableFormats(*item->modelType())->second;
	  }
     }
//      columnSelected(-1, 1, 0);
     dataChanged();
}

void FWTableView::modelSelected(Int_t iRow,Int_t iButton,Int_t iKeyMod,Int_t iGlobalX,Int_t iGlobalY)
{
     if(iKeyMod & kKeyControlMask) {      
	  item()->toggleSelect(iRow);
     } else {
	  FWChangeSentry sentry(*(item()->changeManager()));
	  item()->selectionManager()->clearSelection();
	  item()->select(iRow);
     }
   if(iButton == kButton3) {
      openSelectedModelContextMenu_(iGlobalX,iGlobalY);
   }
}

void FWTableView::columnSelected (Int_t iCol, Int_t iButton, Int_t iKeyMod)
{
     if (iButton == 1 || iButton == 3)
	  m_currentColumn = iCol;
     // update contents of the column editor
     if (m_currentColumn >= 0 && 
	 m_currentColumn < (int)m_tableManager->m_tableFormats->size()) {
	  const FWTableViewManager::TableEntry &entry = 
	       m_tableManager->m_tableFormats->at(m_currentColumn);
	  m_column_name_field->SetText(entry.name.c_str());
	  m_column_expr_field->SetText(entry.expression.c_str());
	  m_column_prec_field->SetText(Form("%d", entry.precision));
     } else {
	  m_column_name_field->SetText("");
	  m_column_expr_field->SetText("");
	  m_column_prec_field->SetText("");
     }
}

void FWTableView::addColumn ()
{
     std::string name = m_column_name_field->GetText();
     std::string expr = m_column_expr_field->GetText();
     // convert the precision to a long int
     char *endptr = 0;
     long int prec = strtol(m_column_prec_field->GetText(), &endptr, 0);
     if (name == "" || expr == "" || 
	 m_column_prec_field->GetText() == 0 || *endptr != 0) {
        fwLog(fwlog::kInfo) << "bad input\n";
	  fflush(stdout);
	  return;
     }
     fwLog(fwlog::kInfo) << "adding column "<<  name << ": " << expr << ", precision " << prec << std::endl;
     fflush(stdout);
//      m_manager->tableFormats(*item->modelType())
     FWTableViewManager::TableEntry e = { expr, name, prec };
     m_tableManager->m_tableFormats->push_back(e);
     m_currentColumn = (int)m_tableManager->m_tableFormats->size() + 1;
     // change needs to be propagated to all tables, because all
     // tables displaying objects of this type are affected
     m_manager->dataChanged();
}

void FWTableView::deleteColumn ()
{
     if (m_currentColumn >= 0 && 
	 m_currentColumn < (int)m_tableManager->m_tableFormats->size()) {
	  m_tableManager->m_tableFormats->erase(m_tableManager->m_tableFormats->begin() + 
					       m_currentColumn);
	  m_column_name_field->SetText("");
	  m_column_expr_field->SetText("");
	  m_column_prec_field->SetText("");
	  m_currentColumn = -1;
     }
     // change needs to be propagated to all tables, because all
     // tables displaying objects of this type are affected
     m_manager->dataChanged();
}

void FWTableView::modifyColumn ()
{
     std::string name = m_column_name_field->GetText();
     std::string expr = m_column_expr_field->GetText();
     // convert the precision to a long int
     char *endptr = 0;
     long int prec = strtol(m_column_prec_field->GetText(), &endptr, 0);
     if (name == "" || expr == "" || 
	 m_column_prec_field->GetText() == 0 || *endptr != 0) {
        fwLog(fwlog::kInfo) << "bad input\n";
	  fflush(stdout);
	  return;
     }
     fwLog(fwlog::kInfo) << "modify column "<<  name << ": " << expr << ", precision " << prec << std::endl;
     fflush(stdout);
//      m_manager->tableFormats(*item->modelType())
     FWTableViewManager::TableEntry e = { expr, name, prec };
     m_tableManager->m_tableFormats->at(m_currentColumn) = e;
     // change needs to be propagated to all tables, because all
     // tables displaying objects of this type are affected
     m_manager->dataChanged();
}

//
// static member functions
//
const std::string&
FWTableView::staticTypeName()
{
   static std::string s_name("Table");
   return s_name;
}
