#include "TGLabel.h"
#include "TGButtonGroup.h"

#include "Fireworks/Core/interface/FWGUIEventFilter.h"
#include "Fireworks/Core/interface/FWEventSelector.h"
#include "Fireworks/Core/interface/FWGUIEventSelector.h"
#include "Fireworks/Core/interface/FWCustomIconsButton.h"
#include "Fireworks/Core/interface/CmsShowNavigator.h"
#include "Fireworks/Core/interface/CSGAction.h"
#include "Fireworks/Core/src/FWCheckBoxIcon.h"

FWGUIEventFilter::FWGUIEventFilter(const TGWindow* parent):
   TGTransientFrame(gClient->GetRoot(), parent, m_width+4, m_height),
   m_applyAction(0),
   m_filterDisableAction(0),
   m_finishEditAction(0),

   m_origFilterMode(CmsShowNavigator::kOr),
   m_isOpen(false),
   m_filtersRemoved(false),

   m_validator(0),
   m_selectionFrameParent(0),
   m_selectionFrame(0),
   m_btnGroup(0),
   m_addBtn(0)
{
   SetWindowName("Event Filters");

   TGVerticalFrame* v1 = new TGVerticalFrame(this);
   AddFrame(v1, new TGLayoutHints(kLHintsExpandX |kLHintsExpandY));

   //-------------------- logical operations

   TGHorizontalFrame* headerFrame = new TGHorizontalFrame(v1, m_width, 2*m_entryHeight, 0);
   v1->AddFrame(headerFrame, new TGLayoutHints(kLHintsExpandX|kLHintsTop, 1, 1, 1, 1));

   m_btnGroup = new TGButtonGroup(headerFrame, "Outputs of enabled selectors are combined as the logical:",
                                  kHorizontalFrame | kFixedHeight);
   m_btnGroup->Resize(1, 38);
   new TGRadioButton(m_btnGroup, "OR");
   new TGRadioButton(m_btnGroup, "AND");
   m_btnGroup->SetLayoutHints(new TGLayoutHints(kLHintsLeft, 0, 16, 1, 0), 0);
   
   TQObject::Connect(m_btnGroup, "Clicked(Int_t)", "FWGUIEventFilter", this, "changeFilterMode(Int_t)");

   headerFrame->AddFrame(m_btnGroup, new TGLayoutHints(kLHintsTop | kLHintsCenterY | kLHintsExpandX, 2, 3, 0, 0));

   //-------------------- selection header

   m_selectionFrameParent =  new TGVerticalFrame(v1, m_width, m_entryHeight, 0);
   v1->AddFrame(m_selectionFrameParent, new TGLayoutHints(kLHintsExpandX|kLHintsTop, 2, 2, 0,0));

   // headers
   TGHorizontalFrame* selH = new TGHorizontalFrame(m_selectionFrameParent);
   m_selectionFrameParent->AddFrame(selH, new TGLayoutHints(kLHintsExpandX));

   {
      TGCompositeFrame *cfr = new TGHorizontalFrame(selH, 292, 22, kFixedSize);
      selH->AddFrame(cfr);
      cfr->AddFrame(new TGLabel(cfr, "Expression:"), new TGLayoutHints(kLHintsLeft|kLHintsBottom, 2, 2, 0, 0));
   }
   {
      TGCompositeFrame *cfr = new TGHorizontalFrame(selH, 122, 22, kFixedSize);
      selH->AddFrame(cfr);
      cfr->AddFrame(new TGLabel(cfr, "Comment:"), new TGLayoutHints(kLHintsLeft|kLHintsBottom, 2, 2, 0, 0));
   }
   {
      TGCompositeFrame *cfr = new TGHorizontalFrame(selH, 39, 22, kFixedSize);
      selH->AddFrame(cfr);
      cfr->AddFrame(new TGLabel(cfr, "Pass:"), new TGLayoutHints(kLHintsLeft|kLHintsBottom, 2, 2, 0, 0));
   }

   //-------------------- adding new selection

   TGHorizontalFrame* addBtnFrame = new TGHorizontalFrame(v1);
   v1->AddFrame(addBtnFrame, new TGLayoutHints(kLHintsRight));

   m_addBtn = new FWCustomIconsButton(addBtnFrame, fClient->GetPicture(FWCheckBoxIcon::coreIcondir() + "plus-sign.png"),
                                      fClient->GetPicture(FWCheckBoxIcon::coreIcondir() + "plus-sign-over.png"),
                                      fClient->GetPicture(FWCheckBoxIcon::coreIcondir() + "plus-sign-disabled.png"));

   addBtnFrame->AddFrame(m_addBtn, new TGLayoutHints(kLHintsRight|kLHintsExpandX|kLHintsExpandY, 0, 6, 4, 1));
   TQObject::Connect(m_addBtn, "Clicked()", "FWGUIEventFilter",  this, "newEntry()");

   //-------------------- external actions
   m_finishEditAction = new CSGAction(this, "Finish");

   TGHorizontalFrame* btnFrame = new TGHorizontalFrame(v1, 280, 30);
   v1->AddFrame(btnFrame, new TGLayoutHints(kLHintsCenterX | kLHintsExpandX | kLHintsBottom , 0, 0, 2, 4));

   m_applyAction = new CSGAction(this, "Apply Filters");

   TGTextButton* cancel = new TGTextButton(btnFrame," Close ");
   btnFrame->AddFrame(cancel, new TGLayoutHints(kLHintsLeft | kLHintsCenterY , 24, 4, 2, 4));
   cancel->Connect("Clicked()","FWGUIEventFilter", this, "CloseWindow()");


   m_applyBtn = new TGTextButton(btnFrame," Apply ");
   btnFrame->AddFrame(m_applyBtn, new TGLayoutHints(kLHintsRight | kLHintsCenterY, 4, 27, 2, 4));
   m_applyBtn->Connect("Clicked()","FWGUIEventFilter", this, "apply()");

   m_filterDisableAction = new CSGAction(this, " Disable Filtering ");
   m_filterDisableAction->createTextButton(btnFrame,new TGLayoutHints(kLHintsRight | kLHintsCenterY, 4, 18, 2, 4) );

}

void
FWGUIEventFilter::addSelector(FWEventSelector* sel)
{
   FWGUIEventSelector* es = new FWGUIEventSelector(m_selectionFrame, m_validator, sel);
   m_selectionFrame->AddFrame(es, new TGLayoutHints(kLHintsExpandX));
   TQObject::Connect(es, "removeSelector(FWGUIEventSelector*)", "FWGUIEventFilter",  this, "deleteEntry(FWGUIEventSelector*)");
   TQObject::Connect(es, "selectorChanged()", "FWGUIEventFilter",  this, "checkApplyButton()");

   m_guiSelectors.push_back(es);
}

void
FWGUIEventFilter::show( std::list<FWEventSelector*>* sels,  fwlite::Event* event, int filterMode)
{
   m_applyBtn->SetForegroundColor(0x000000);
   m_filtersRemoved = false;

   m_isOpen = true;
   m_validator = new FWHLTValidator(*event);

   m_origFilterMode = filterMode;

   // Button ids run from 1
   if (TGButton *btn = m_btnGroup->GetButton(filterMode)) {
      btn->SetDown();
      m_btnGroup->SetButton(filterMode);
   }

   assert(m_selectionFrame == 0);
   m_selectionFrame = new TGVerticalFrame(m_selectionFrameParent);
   m_selectionFrameParent->AddFrame(m_selectionFrame,  new TGLayoutHints(kLHintsExpandX));

   for(std::list<FWEventSelector*>::iterator i = sels->begin(); i != sels->end(); ++i)
      addSelector(*i);

   MapSubwindows();
   Layout();
   MapWindow();
}


///////////////////////////////////////////
//   Callbacks
///////////////////////////////////////////

void
FWGUIEventFilter::deleteEntry(FWGUIEventSelector* sel)
{
   m_filtersRemoved = true;

   m_guiSelectors.remove(sel);

   m_selectionFrame->RemoveFrame(sel);
   Layout();
   gClient->NeedRedraw(this);
}

void
FWGUIEventFilter::newEntry()
{
   addSelector(0);
   MapSubwindows();
   Layout();
}

void
FWGUIEventFilter::apply()
{
   m_applyAction->activate();

   m_origFilterMode = getFilterMode();
   m_filtersRemoved = false;
   m_applyBtn->SetForegroundColor(0x000000);
   fClient->NeedRedraw( this);

}

int
FWGUIEventFilter::getFilterMode()
{
   if (m_btnGroup->GetButton(1)->IsOn())
      return CmsShowNavigator::kOr;
   else
      return CmsShowNavigator::kAnd;
}

void
FWGUIEventFilter::CloseWindow()
{
   m_isOpen = false;
   m_selectionFrameParent->RemoveFrame(m_selectionFrame);
   m_selectionFrame = 0;

   FWGUIEventSelector* gs;
   for (std::list<FWGUIEventSelector*>::iterator i = m_guiSelectors.begin(); i != m_guiSelectors.end(); ++i)
   {
      gs = *i;
      delete gs;
   }

   m_guiSelectors.clear();

   delete m_validator;
   m_validator = 0;

   UnmapWindow();
   m_finishEditAction->activated();
}

void
FWGUIEventFilter::changeFilterMode(Int_t)
{
   checkApplyButton();
}

void
FWGUIEventFilter::checkApplyButton()
{
   bool changed = m_filtersRemoved;

   if (!changed)
   {
      changed = (getFilterMode() != m_origFilterMode);

      if (!changed)
      {
         std::list<FWGUIEventSelector*>::iterator i = m_guiSelectors.begin();
         while (i != m_guiSelectors.end())
         {
            if ((*i)->origSelector() == 0)
               break;

            if ( (*i)->guiSelector()->m_enabled    != (*i)->origSelector()->m_enabled  ||
                 (*i)->guiSelector()->m_expression != (*i)->origSelector()->m_expression )
               break;

            ++i;
         }
         changed = (i != m_guiSelectors.end());
      }
   }

   m_applyBtn->SetForegroundColor(changed ? 0x40FF80 : 0x000000);
   gClient->NeedRedraw( m_applyBtn);
}

