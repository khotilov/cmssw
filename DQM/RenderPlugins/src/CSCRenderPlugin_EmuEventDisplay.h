#ifndef CSCRenderPlugin_EmuEventDisplay_H
#define CSCRenderPlugin_EmuEventDisplay_H

/*
 * =====================================================================================
 *
 *       Filename:  CSCRenderPlugin_EmuEventDisplay.h
 *
 *    Description:  CSC Histogram Rendering Plugin
 *
 *        Version:  1.0
 *        Created:  12/02/2010 07:50:48 PM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Valdas Rapsevicius (VR), Valdas.Rapsevicius@cern.ch
 *        Company:  CERN, CH
 *
 * =====================================================================================
 */

#define CSC_RENDER_PLUGIN

#ifdef DQMLOCAL
#include "DQM/CSCMonitorModule/interface/CSCDQM_Detector.h"
#else
#include "CSCDQM_Detector.h"
#endif 

#include <string>
#include <iostream>
#include <vector>
#include <TH1.h>
#include <TH2.h>
#include <TBox.h>
#include <TText.h>
#include <TPRegexp.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TPolyLine.h>

/**
 * @class EmuEventDisplay
 * @brief Class that draws CSC Map diagram
 */
class EmuEventDisplay {
  
  private:

    cscdqm::Detector detector;

    TH2F* histo_zr;
    TH2F* histo_zphi;
    TH2F* histo_xy;

    std::vector<TPolyLine*> zrChambersUp;
    std::vector<TPolyLine*> zrChambersDown;
    std::vector<TPolyLine*> zpChambers;
    std::vector<TPolyLine*> xyChambers;

    std::vector<TBox*> zrHits;
    std::vector<TBox*> zrHitsUp;
    std::vector<TBox*> zrHitsDown;
    std::vector<TBox*> zpHits;
    std::vector<TBox*> xyHits;

    template<class T>
    void deleteItems(std::vector<T*> &v) {
      for (unsigned int i = 0; i < v.size(); i++) {
        T* t = v.at(i);
        delete t;
      }
      v.clear();
    }

  public:
    
    EmuEventDisplay();
    ~EmuEventDisplay();

    void drawEventDisplay_ZR(TH2* data, int histogramType);
    void drawEventDisplay_ZPhi(TH2* data);
    void drawEventDisplay_XY(TH2* data);

};

#endif
