#ifndef Fireworks_Core_FWColorManager_h
#define Fireworks_Core_FWColorManager_h
// -*- C++ -*-
//
// Package:     Core
// Class  :     FWColorManager
// 
/**\class FWColorManager FWColorManager.h Fireworks/Core/interface/FWColorManager.h

 Description: <one line class summary>

 Usage:
    <usage>

*/
//
// Original Author:  Chris Jones
//         Created:  Tue Mar 24 10:07:58 CET 2009
// $Id: FWColorManager.h,v 1.4 2009/04/13 21:22:28 chrjones Exp $
//

// system include files
#include <vector>
#include "sigc++/signal.h"
#include "Rtypes.h"

// user include files

// forward declarations
class FWModelChangeManager;

enum FWGeomColorIndex {
   kFWMuonBarrelMainColorIndex,
   kFWMuonBarrelLineColorIndex,
   kFWMuonEndCapMainColorIndex,
   kFWMuonEndCapLineColorIndex,
   kFWTrackerColorIndex
};

class FWColorManager {

public:
   FWColorManager(FWModelChangeManager*);
   virtual ~FWColorManager();
   
   // ---------- const member functions ---------------------
   Color_t background() const {return m_background;}
   Color_t foreground() const {return m_foreground;}
   
   Color_t indexToColor(unsigned int) const;
   unsigned int numberOfIndicies() const;
   
   unsigned int colorToIndex(Color_t) const;

   //help with backward compatibility with old config files
   unsigned int oldColorToIndex(Color_t) const;
   
   bool colorHasIndex(Color_t) const;
   
   Color_t geomColor(FWGeomColorIndex) const;
   
   enum BackgroundColorIndex {kBlackIndex, kWhiteIndex};
   BackgroundColorIndex backgroundColorIndex() const;
   // ---------- static member functions --------------------
   
   // ---------- member functions ---------------------------
   void increaseBrightness();
   void decreaseBrightness();
   void setBackgroundColorIndex(BackgroundColorIndex);

   mutable sigc::signal<void> colorsHaveChanged_;
   //called after all the slots attached to colorsHaveChanged_ are done
   mutable sigc::signal<void> colorsHaveChangedFinished_;
private:
   FWColorManager(const FWColorManager&); // stop default
   
   const FWColorManager& operator=(const FWColorManager&); // stop default
   void updateBrightness();

   // ---------- member data --------------------------------

   Color_t m_background;
   Color_t m_foreground;
   FWModelChangeManager* m_changeManager;
   
   unsigned int m_startColorIndex;
   unsigned int m_startGeomColorIndex;
};


#endif
