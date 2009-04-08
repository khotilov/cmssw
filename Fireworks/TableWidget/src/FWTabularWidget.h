#ifndef Fireworks_TableWidget_FWTabularWidget_h
#define Fireworks_TableWidget_FWTabularWidget_h
// -*- C++ -*-
//
// Package:     TableWidget
// Class  :     FWTabularWidget
// 
/**\class FWTabularWidget FWTabularWidget.h Fireworks/TableWidget/interface/FWTabularWidget.h

 Description: Widget that draws part of a table [Implementation detail of FWTableWidget]

 Usage:
    This class is used internally by FWTableWidget.

*/
//
// Original Author:  Chris Jones
//         Created:  Mon Feb  2 16:45:19 EST 2009
// $Id: FWTabularWidget.h,v 1.2 2009/03/04 16:13:05 chrjones Exp $
//

// system include files
#include <vector>
#include "TGFrame.h"

// user include files

// forward declarations
class FWTableManagerBase;

class FWTabularWidget : public TGFrame
{

   public:
      static const TGGC&  getDefaultGC();
      
      FWTabularWidget(FWTableManagerBase* iManager,const TGWindow* p=0, GContext_t context = getDefaultGC()());
      virtual ~FWTabularWidget();

      // ---------- const member functions ---------------------
      const std::vector<unsigned int>& widthOfTextInColumns() const { return m_widthOfTextInColumns;}
      UInt_t verticalOffset() const {return m_vOffset;}
      UInt_t horizontalOffset() const { return m_hOffset;}

      // ---------- static member functions --------------------

      // ---------- member functions ---------------------------
      void setWidthOfTextInColumns(const std::vector<unsigned int>& );
      void DoRedraw();
      TGDimension GetDefaultSize() const;

      void setVerticalOffset(UInt_t);
      void setHorizontalOffset(UInt_t);

      virtual Bool_t HandleButton(Event_t *event);

      void buttonPressed(Int_t row, Int_t column, Event_t* event, Int_t relX, Int_t relY); //*SIGNAL*
      void buttonReleased(Int_t row, Int_t column, Event_t* event, Int_t relX, Int_t relY); //*SIGNAL*

      void dataChanged();
      void needToRedraw();

      ClassDef(FWTabularWidget,0);
   
      void setLineContext(GContext_t iContext);
      void setBackgroundAreaContext(GContext_t iContext);
   
   private:
      //FWTabularWidget(const FWTabularWidget&); // stop default

      //const FWTabularWidget& operator=(const FWTabularWidget&); // stop default

      // ---------- member data --------------------------------
      static const int kTextBuffer=2;
      static const int kSeperatorWidth=1;

      void translateToRowColumn(Int_t iX, Int_t iY, Int_t& oRow, Int_t& oCol, Int_t&oRelX, Int_t& oRelY) const;

      FWTableManagerBase* m_table;
      std::vector<unsigned int> m_widthOfTextInColumns;
      int m_textHeight;
      int m_tableWidth;

      unsigned int m_vOffset;
      unsigned int m_hOffset;

      GContext_t m_normGC;
      GContext_t m_backgroundGC;
};


#endif
