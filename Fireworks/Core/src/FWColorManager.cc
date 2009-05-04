// -*- C++ -*-
//
// Package:     Core
// Class  :     FWColorManager
// 
// Implementation:
//     <Notes on implementation>
//
// Original Author:  Chris Jones
//         Created:  Tue Mar 24 10:10:01 CET 2009
// $Id: FWColorManager.cc,v 1.5 2009/04/14 15:45:25 chrjones Exp $
//

// system include files
#include <iostream>
#include <map>
#include <boost/shared_ptr.hpp>
#include "TColor.h"
#include "TROOT.h"
#include "TObjArray.h"
#include "TMath.h"

// user include files
#include "Fireworks/Core/interface/FWColorManager.h"
#include "Fireworks/Core/interface/FWModelChangeManager.h"

//
// constants, enums and typedefs
//

//
// static data member definitions
//
//static std::vector<Color_t>* s_forWhite=0;
//static std::vector<Color_t>* s_forBlack=0;

static Float_t m_gammaOff = 0;

enum {
   kFWRed =8,
   kFWBlue = 5,
   kFWCyan = 7,
   kFWGreen = 9,
   kFWMagenta = 1,
   kFWOrange = 4,
   kFWYellow = 0
/*
   kFWRed =0,
   kFWBlue = 1,
   kFWCyan = 2,
   kFWGreen = 3,
   kFWMagenta = 4,
   kFWOrange = 5,
   kFWYellow = 6,
   kFWGray = 7
 */
};

static const float s_forWhite[][3] ={
{0.82, 0.82, 0.17}, //yellow (made it a bit darker)
{0.53, 0.00, 0.69}, //purple
{0.98, 0.74, 0.01}, //yellowish-orange
{0.24, 0.00, 0.64}, //purplish-blue
{0.98, 0.60, 0.01}, //orange
{0.01, 0.28, 1.00}, //blue
{0.99, 0.33, 0.03}, //dark orange
{0.01, 0.83, 0.81}, //cyan
{1.00, 0.09, 0.00}, //red
{0.40, 0.69, 0.20}, //green
{0.65, 0.10, 0.29}, //burgundy
{0.65, 0.92, 0.17}, //lime
/*
{1.,0.,0.}, //red
{0.,0.,1.}, //blue
{0.,1.,1.},  //cyan
{0.,1.,0.}, //green
{1.,0.,1.}, //magenta
{1.,0.5,0.0}, //orange
{1.,1.,0.}, //yellow
{0.5,0.5,0.5}, //gray
*/
{0.99, 1.00, 0.46},
{0.89, 0.76, 0.93},
{0.99, 0.90, 0.64},
{0.82, 0.76, 0.92},
{1.00, 0.85, 0.64},
{0.75, 0.82, 0.99},
{1.00, 0.83, 0.76},
{0.75, 0.98, 0.96},
{0.99, 0.78, 0.74},
{0.80, 0.88, 0.70},
{0.92, 0.78, 0.82},
{0.72, 0.96, 0.58}
/*
{1.,0.3,0.3},
{0.3,0.3,1.},
{0.3,1.,1.},
{0.3,1.,0.3},
{1.,0.3,1.},
{1.,0.7,0.0},
{1.,1.,0.3},
{0.7,0.7,0.7}
 */
};

static const float s_forBlack[][3] ={
{1.00, 1.00, 0.20}, //yellow
{0.53, 0.00, 0.69}, //purple
{0.98, 0.74, 0.01}, //yellowish-orange
{0.24, 0.00, 0.64}, //purplish-blue
{0.98, 0.60, 0.01}, //orange
{0.01, 0.28, 1.00}, //blue
{0.99, 0.33, 0.03}, //dark orange
{0.01, 0.83, 0.81}, //cyan
{1.00, 0.09, 0.00}, //red
{0.40, 0.69, 0.20}, //green
{0.65, 0.10, 0.29}, //burgundy
{0.65, 0.92, 0.17}, //lime
/*
{1.,0.,0.}, //red
{0.,0.,1.}, //blue
{0.,1.,1.}, //cyan
{0.,1.,0.}, //green
{1.,0.,1.}, //magenta
{1.,0.5,0.0},  //orange
{1.,1.,0.}, //yellow
{0.5,0.5,0.5}, //gray
*/
{0.27, 0.27, 0.04},
{0.19, 0.00, 0.24},
{0.19, 0.15, 0.00},
{0.14, 0.00, 0.38},
{0.19, 0.11, 0.00},
{0.01, 0.10, 0.33},
{0.17, 0.05, 0.02},
{0.00, 0.33, 0.29},
{0.34, 0.03, 0.01},
{0.15, 0.24, 0.06},
{0.24, 0.02, 0.11},
{0.22, 0.30, 0.07}
/*
{0.7,0.0,0.0},
{0.0,0.0,0.7},
{0.0,.7,0.7},
{0.0,.7,0.},
{.7,0.,.7},
{.7,0.4,0.0},
{.7,.7,0.0},
{0.3,0.3,0.3}
 */
};

static const unsigned int s_size = sizeof(s_forBlack)/sizeof(s_forBlack[0]);

static const float s_geomForWhite[][3] ={
{1.,0.5,0.5},
{1.,0.5,0.5},
{0.45,0.68,1.},
{0.2,0.4,1.},
{0.57,1.,0.26}
};

static const float s_geomForBlack[][3] ={
{0x3f/256.,0.,0.},
{0x7f/256.,0.,0.},
{0.,0.,0x3f/256.},
{0.,0.,0x7f/256.},
{0.,0x7f/256.,0.}
};

static const unsigned int s_geomSize = sizeof(s_geomForBlack)/sizeof(s_geomForBlack[0]);

//==============================================================================

static
void resetColors(const float(* iColors)[3], unsigned int iSize, unsigned int iStart)
{
   TSeqCollection* colorTable = gROOT->GetListOfColors();
   
   TColor* c = static_cast<TColor*>(colorTable->At(iStart));
   unsigned int index = iStart;
   if(0==c || c->GetNumber() != static_cast<int>(iStart)) {
      TIter   next(colorTable);
      while( (c=static_cast<TColor*>( next() )) ) {
         if(c->GetNumber()==static_cast<int>(iStart)) {
            index = iStart;
            break;
         }
      }
   }
   assert(0!=c);
   
   for(unsigned int i = index; i< index+iSize; ++i,++iColors) {
      TColor* c = dynamic_cast<TColor*> (colorTable->At(i));
      float red = (*iColors)[0];
      float green = (*iColors)[1];
      float blue = (*iColors)[2];
     
      // apply brightness
      red     =  TMath::Power( red, (2.5 + m_gammaOff)/2.5);
      green  = TMath::Power(green, (2.5 + m_gammaOff)/2.5);
     blue    = TMath::Power(blue, (2.5 + m_gammaOff)/2.5);

      c->SetRGB(red,green,blue);
   }
}
//
// constructors and destructor
//
FWColorManager::FWColorManager(FWModelChangeManager* iManager):
m_background(kBlack),
m_foreground(kWhite),
m_changeManager(iManager)
{
   TObjArray* colorTable = dynamic_cast<TObjArray*>(gROOT->GetListOfColors());
   m_startColorIndex = static_cast<TColor*>(colorTable->Last())->GetNumber()+1;
   
   unsigned int index = m_startColorIndex;
   //std::cout <<"start color index "<<m_startColorIndex<<std::endl;
   
   const float(* itEnd)[3] = s_forBlack+s_size;
   for(const float(* it)[3] = s_forBlack;
       it != itEnd;
       ++it) {
      //NOTE: this constructor automatically places this color into the gROOT color list
      //std::cout <<" color "<< index <<" "<<(*it)[0]<<" "<<(*it)[1]<<" "<<(*it)[2]<<std::endl;
      new TColor(index++,(*it)[0],(*it)[1],(*it)[2]);
   }

   m_startGeomColorIndex = index;
   itEnd = s_geomForBlack+s_geomSize;
   for(const float(* it)[3] = s_geomForBlack;
       it != itEnd;
       ++it) {
      //NOTE: this constructor automatically places this color into the gROOT color list
      new TColor(index++,(*it)[0],(*it)[1],(*it)[2]);
   }   
}

// FWColorManager::FWColorManager(const FWColorManager& rhs)
// {
//    // do actual copying here;
// }

FWColorManager::~FWColorManager()
{
}

//
// assignment operators
//
// const FWColorManager& FWColorManager::operator=(const FWColorManager& rhs)
// {
//   //An exception safe implementation is
//   FWColorManager temp(rhs);
//   swap(rhs);
//
//   return *this;
// }

//
// member functions

void FWColorManager::updateBrightness()
{

   if(backgroundColorIndex() == kBlackIndex) {
      resetColors(s_forBlack,s_size,m_startColorIndex);
      resetColors(s_geomForBlack, s_geomSize, m_startGeomColorIndex);
   } else {
      resetColors(s_forWhite,s_size,m_startColorIndex);
      resetColors(s_geomForWhite, s_geomSize, m_startGeomColorIndex);
   }
   FWChangeSentry sentry(*m_changeManager);
   colorsHaveChanged_();
   colorsHaveChangedFinished_();
}

void
FWColorManager::decreaseBrightness()
{
   Float_t value =  m_gammaOff + 0.1;
   if (value < -0.5 || value > 0.5)
   {
      printf("Warning::Set brightness out of range.  value '%f' out of range [-0.5, 0.5]. \n", value);
   }
   
   m_gammaOff = value;
   updateBrightness();
}

void
FWColorManager::increaseBrightness()
{
   Float_t value =  m_gammaOff - 0.1;
   if (value < -0.5 || value > 0.5)
   {
      printf("Warning::Set brightness out of range.  value '%f' out of range [-0.5, 0.5].\n", value);
   }
   m_gammaOff = value;
   updateBrightness();
}

void 
FWColorManager::setBackgroundColorIndex(BackgroundColorIndex iIndex)
{
   if(backgroundColorIndex()!=iIndex) {
      if(backgroundColorIndex()==kBlackIndex) {
         m_background=kWhite;
         m_foreground=kBlack;
         resetColors(s_forWhite,s_size,m_startColorIndex);
         resetColors(s_geomForWhite, s_geomSize, m_startGeomColorIndex);
      } else {
         m_background=kBlack;
         m_foreground=kWhite;
         resetColors(s_forBlack,s_size,m_startColorIndex);
         resetColors(s_geomForBlack, s_geomSize, m_startGeomColorIndex);
      }
      FWChangeSentry sentry(*m_changeManager);
      colorsHaveChanged_();
      colorsHaveChangedFinished_();
   }
}

//
// const member functions
//
Color_t 
FWColorManager::indexToColor(unsigned int iIndex) const
{
   return m_startColorIndex+iIndex;
}
unsigned int 
FWColorManager::numberOfIndicies() const
{
   return s_size;
}

FWColorManager::BackgroundColorIndex 
FWColorManager::backgroundColorIndex() const
{
   if(m_background==kBlack) {
      return kBlackIndex;
   }
   return kWhiteIndex;
}

unsigned int 
FWColorManager::colorToIndex(Color_t iColor) const
{
   if(iColor < static_cast<int>(m_startColorIndex) ) {
      std::cerr <<"asked to convert a non-standard color "<<iColor<<". Will attempt to use old color scheme"<<std::endl;
      return oldColorToIndex(iColor);
   }
   return iColor - m_startColorIndex;
}

bool 
FWColorManager::colorHasIndex(Color_t iColor) const
{
   return iColor >= static_cast<int>(m_startColorIndex);
}


Color_t 
FWColorManager::geomColor(FWGeomColorIndex iIndex) const
{
   return m_startGeomColorIndex+iIndex;
}


static boost::shared_ptr<std::map<Color_t,unsigned int> > m_oldColorToIndexMap;

unsigned int 
FWColorManager::oldColorToIndex(Color_t iColor) const
{
   if(0==m_oldColorToIndexMap.get()) {
      m_oldColorToIndexMap = boost::shared_ptr<std::map<Color_t,unsigned int> >(new std::map<Color_t,unsigned int>());
      (*m_oldColorToIndexMap)[kRed]=kFWRed;
      (*m_oldColorToIndexMap)[kBlue]=kFWBlue;
      (*m_oldColorToIndexMap)[kYellow]=kFWYellow;
      (*m_oldColorToIndexMap)[kGreen]=kFWGreen;
      (*m_oldColorToIndexMap)[kCyan]=kFWCyan;
      (*m_oldColorToIndexMap)[kTeal]=kFWCyan;
      (*m_oldColorToIndexMap)[kMagenta]=kFWMagenta;
      (*m_oldColorToIndexMap)[kViolet]=kFWMagenta;
      (*m_oldColorToIndexMap)[kOrange]=kFWOrange;
      //(*m_oldColorToIndexMap)[kGray]=kFWGray;
      (*m_oldColorToIndexMap)[3]=kFWGreen;
      
   }
   return (*m_oldColorToIndexMap)[iColor];
}

//
// static member functions
//
