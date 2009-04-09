#ifndef Fireworks_TableWidget_FWTableWidget_h
#define Fireworks_TableWidget_FWTableWidget_h
// -*- C++ -*-
//
// Package:     TableWidget
// Class  :     FWTableWidget
// 
/**\class FWTableWidget FWTableWidget.h Fireworks/TableWidget/interface/FWTableWidget.h

 Description: ROOT widget for representing data in a tabular form where the data is accessed via a FWTableManagerBase

 Usage:
    This widget creates a table made up of cells where the representation of the cells is controlled by the FWTableManagerBase instance
    passed to the FWTableWidget's constructor. If the data for the FWTableManagerBase changes, the view is automatically updated. See
    the documentation for FWTableManagerBase for further details.

*/
//
// Original Author:  Chris Jones
//         Created:  Mon Feb  2 16:45:47 EST 2009
// $Id: FWTableWidget.h,v 1.5 2009/04/08 19:04:15 chrjones Exp $
//

// system include files
#include "TGFrame.h"

// user include files

// forward declarations
class FWTableManagerBase;
class FWTabularWidget;
class TGVScrollBar;
class TGHScrollBar;
class TGGC;

class FWTableWidget : public TGCompositeFrame
{

   public:
      FWTableWidget(FWTableManagerBase* iManager,const TGWindow* p=0);
      virtual ~FWTableWidget();

      // ---------- const member functions ---------------------
      //the standard method gives the wrong values when using a TGTabularLayout
      TGDimension GetDefaultSize() const;

      // ---------- static member functions --------------------

      // ---------- member functions ---------------------------
      using  TGFrame::Resize;
      void Resize(UInt_t w, UInt_t h);
      void SetSize(const TGDimension &s);
      virtual void    MoveResize(Int_t x, Int_t y, UInt_t w = 0, UInt_t h = 0);
      void SetBackgroundColor(Pixel_t);
      void SetHeaderBackgroundColor(Pixel_t);
      void SetHeaderForegroundColor(Pixel_t);
      void SetLineSeparatorColor(Pixel_t);
   
      virtual Bool_t ProcessMessage(Long_t msg, Long_t parm1, Long_t parm2);

      /**Sorts the rows of data in the table based on values in column iColumn.
      If iDescendingSort is 'true' then the rows are sorted in descending order of the values.
      */
      void sort(UInt_t iColumn, bool iDescendingSort);

      void buttonReleasedInHeader(Int_t row, Int_t column, Event_t* event,Int_t, Int_t);
      void buttonReleasedInBody(Int_t row, Int_t column, Event_t* event, Int_t,Int_t);

      /**This signal is emitted if the mouse button is 'clicked' while the cursor
      was over a row. Arguments:
      iRow: the unsorted row number (natural ordering) of the row clicked
      iButton: the ROOT button value for the click event (e.g. says which button used)
      iKeyMod: the ROOT key modifier value for the click event (e.g. says if a keyboard key was being held)
      */
      void rowClicked(Int_t iRow, Int_t iButton, Int_t iKeyMod); //*SIGNAL*

      ClassDef(FWTableWidget,0);
   
      void dataChanged();

     void buttonPressedInRowHeader(Int_t row, Int_t column, Event_t* event, Int_t relX, Int_t relY);
     void buttonReleasedInRowHeader(Int_t row, Int_t column, Event_t* event, Int_t relX, Int_t relY);

   private:
      //FWTableWidget(const FWTableWidget&); // stop default

      //const FWTableWidget& operator=(const FWTableWidget&); // stop default

      // ---------- member data --------------------------------
      void handleResize(UInt_t w, UInt_t h);
      FWTableManagerBase* m_bodyTable;
      FWTableManagerBase* m_headerTable;
      FWTableManagerBase* m_rowHeaderTable;
      FWTabularWidget* m_header;
      FWTabularWidget* m_body;
      FWTabularWidget* m_rowHeader;
      TGVScrollBar* m_vSlider;
      TGHScrollBar* m_hSlider;
      bool m_showingVSlider;
      bool m_showingHSlider;

      int m_sortedColumn;
      bool m_descendingSort;
   
      TGGC* m_headerBackground;
      TGGC* m_headerForeground;
      TGGC* m_lineSeparator;

};


#endif
