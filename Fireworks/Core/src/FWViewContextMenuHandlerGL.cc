#include "Fireworks/Core/interface/FWViewContextMenuHandlerGL.h"

#include "TClass.h"
#include "TEveViewer.h"
#include "TGLViewer.h"
#include "TGLAnnotation.h"
#include "TGLWidget.h"
#include "TEveVector.h"

#include "Fireworks/Core/interface/FWModelId.h"
#include "Fireworks/Core/interface/FWEventItem.h"
#include "Fireworks/Core/interface/FWEveView.h"


#include "Fireworks/Core/interface/FWItemValueGetter.h"
#include "Fireworks/Core/interface/FWSelectionManager.h"
#include "Fireworks/Core/interface/FWEventItem.h"
#include "Fireworks/Core/interface/Context.h"
#include "Fireworks/Core/interface/FWRPZView.h"

FWViewContextMenuHandlerGL::FWViewContextMenuHandlerGL( FWEveView* v):
m_view(v)
{
}

void 
FWViewContextMenuHandlerGL::init(FWViewContextMenuHandlerBase::MenuEntryAdder& adder, const FWModelId &id)
{
   adder.addEntry("Add Annotation", kAnnotate);
   if (FWViewType::isProjected(m_view->typeId()))
   { 
      const char* p  = id.item()->purpose().c_str();
      bool enabled = (strstr(p, "Beam Spot") || strstr(p, "Vertices") || strstr(p, "Vertex") );
      adder.addEntry("Use As Projection Origin", kCameraCenter, enabled);
      adder.addEntry("Reset Projection Origin", kResetCameraCenter, enabled);
   }
   else
   { 
      adder.addEntry("Set Camera Center", kCameraCenter);
      adder.addEntry("Reset Camera Center", kResetCameraCenter);
   }   
}

void 
FWViewContextMenuHandlerGL::select(int iEntryIndex, const FWModelId &id, int iX, int iY)
{
   TGLViewer* v = m_view->viewerGL();

   Window_t wdummy;
   Int_t x,y;
   gVirtualX->TranslateCoordinates(gClient->GetDefaultRoot()->GetId(), v->GetGLWidget()->GetId(), iX, iY, x, y, wdummy);   
   TGLVector3 pnt(x, y, 0.5*v->GetSelRec().GetMinZ());
   v->CurrentCamera().WindowToViewport(pnt);
   pnt = v->CurrentCamera().ViewportToWorld(pnt);

   switch (iEntryIndex)
   {
      case kAnnotate:
      {
         TGFrame* f = v->GetGLWidget();
         gVirtualX->TranslateCoordinates(gClient->GetDefaultRoot()->GetId(), f->GetId(), iX, iY, x, y, wdummy);

         std::string name = id.item()->modelName(id.index());
         if (id.item()->haveInterestingValue())
            name += ", " + id.item()->modelInterestingValueAsString(id.index());

         TGLAnnotation* an = new TGLAnnotation(v, name.c_str(),  x*1.f/f->GetWidth(), 1 - y*1.f/f->GetHeight(), pnt);
         an->SetUseColorSet(true);
         an->SetTextSize(0.03);
         break;
      }
      case kCameraCenter:
      {
         TEveVector center;
         if (FWViewType::isProjected(m_view->typeId()))
         {

            // AMT:: find way to get values without using FWItemValueGetter

            FWModelId mId = *(m_view->context().selectionManager()->selected().begin());
            const FWEventItem* item = mId.item();

            bool bs = strstr(item->purpose().c_str(), "Beam Spot");
            std::vector<std::pair<std::string,std::string> > func;
            func.push_back(std::pair<std::string,std::string>("", "cm"));
            {  
               func.back().first = bs ? "x0" : "x";
               FWItemValueGetter valueGetter(ROOT::Reflex::Type::ByTypeInfo(*(item->modelType()->GetTypeInfo())), func);
               center.fX = valueGetter.valueFor(item->modelData(mId.index())); 
            }
            {
               func.back().first = bs ? "y0" : "y";
               FWItemValueGetter valueGetter(ROOT::Reflex::Type::ByTypeInfo(*(item->modelType()->GetTypeInfo())), func);
               center.fY = valueGetter.valueFor(item->modelData(mId.index()));
            }
            {  
               func.back().first = bs ? "z0" : "z";
               FWItemValueGetter valueGetter(ROOT::Reflex::Type::ByTypeInfo(*(item->modelType()->GetTypeInfo())), func);
               center.fZ = valueGetter.valueFor(item->modelData(mId.index()));
            }

            FWRPZView* pv = static_cast<FWRPZView*>(m_view);
            pv->shiftOrigin(center);
         }

         else
         {
            v->CurrentCamera().SetCenterVec(pnt.X(), pnt.Y(), pnt.Z());
            v->CurrentCamera().SetExternalCenter(true);
            v->SetDrawCameraCenter(true);
         }

         break;
      }
      case kResetCameraCenter:
      { 
         if (FWViewType::isProjected(m_view->typeId()))
         {
            FWRPZView* pv = static_cast<FWRPZView*>(m_view);
            pv->resetOrigin();
         }

         v->CurrentCamera().SetExternalCenter(false);
         v->SetDrawCameraCenter(false);
         break;
      }
   }
}
