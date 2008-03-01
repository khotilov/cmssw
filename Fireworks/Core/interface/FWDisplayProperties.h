#ifndef Fireworks_Core_FWDisplayProperties_h
#define Fireworks_Core_FWDisplayProperties_h
// -*- C++ -*-
//
// Package:     Core
// Class  :     FWDisplayProperties
// 
/**\class FWDisplayProperties FWDisplayProperties.h Fireworks/Core/interface/FWDisplayProperties.h

 Description: <one line class summary>

 Usage:
    <usage>

*/
//
// Original Author:  
//         Created:  Thu Jan  3 14:22:36 EST 2008
// $Id: FWDisplayProperties.h,v 1.2 2008/01/21 01:17:09 chrjones Exp $
//

// system include files
#include "Rtypes.h"

// user include files

// forward declarations

class FWDisplayProperties
{

   public:
      FWDisplayProperties(const Color_t& iColor = kWhite,
			  bool isVisible = true ):
	m_color(iColor),
	m_isVisible(isVisible) {}
      //virtual ~FWDisplayProperties();

      // ---------- const member functions ---------------------
      const Color_t& color() const {
	return m_color;
      }

      bool isVisible() const {
	return m_isVisible;
      }
   
      bool operator==(const FWDisplayProperties& iRHS) const {
         return m_color == iRHS.m_color && m_isVisible == iRHS.m_isVisible;
      }
      bool operator!=(const FWDisplayProperties& iRHS) const {
         return not (*this == iRHS);
      }
   
      // ---------- static member functions --------------------

      // ---------- member functions ---------------------------
      void setIsVisible(bool iSet) {
         m_isVisible = iSet;
      }
   private:
      //FWDisplayProperties(const FWDisplayProperties&); // stop default

      //const FWDisplayProperties& operator=(const FWDisplayProperties&); // stop default

      // ---------- member data --------------------------------
      Color_t m_color;
      bool m_isVisible;
};


#endif
