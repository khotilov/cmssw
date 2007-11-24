// $Id: EBMDisplayPlugins.h,v 1.7 2007/09/26 11:42:10 dellaric Exp $

#ifndef  EBMDisplayPlugins_H
# define EBMDisplayPlugins_H

/*!
  \file EBMDisplayPlugins
  \brief Display Plugin for Quality Histograms (2D)
  \author B. Gobbo 
  \version $Revision: 1.7 $
  \date $Date: 2007/09/26 11:42:10 $
*/

#include "VisMonitoring/VisDQMBase/interface/VisDQMDisplayPlugin.h"
#include <string>
#include <TH2C.h>

class EBMDisplayPlugins : public VisDQMDisplayPlugin {

 public:

  static const char * catalogLabel( void ) {
    return "Ecal Barrel Monitor Plots";	    
  }	
    
  EBMDisplayPlugins( IgState *state );

  virtual bool applies( DisplayData *data );

  virtual std::string preDraw( DisplayData *data );

  virtual void postDraw( DisplayData *data );

 private:

  int nbx;
  int nby;

  std::string name;

  // Temporary workaround due to erroneous multiple instantiation of this class...
  //static bool first;
  //static TH2C* t1;
  //static TH2C* t2;
  //static TH2C* t3;
  //static TH2C* t4;
  //static TH2C* t5;

  TH2C* text1;
  TH2C* text2;
  TH2C* text3;
  TH2C* text4;
  TH2C* text5;
  TH2C* text6;
  TH2C* text7;

  int pCol3[6]; 
  int pCol4[10];
    
  // private functions...
  std::string preDrawTProfile2D( DisplayData *data );
  std::string preDrawTProfile( DisplayData *data );
  std::string preDrawTH2F( DisplayData *data );
  std::string preDrawTH1F( DisplayData *data );
  void postDrawTProfile2D( DisplayData *data );
  void postDrawTH2F( DisplayData *data );
  template<class T> void adjustRange( T obj );

};

#endif //  EBMDisplayPlugins_H
