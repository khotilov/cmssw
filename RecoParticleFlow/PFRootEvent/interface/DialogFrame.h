#ifndef PF_DialogFrame_h
#define PF_DialogFrame_h

#include <TGClient.h>
#include <TGButton.h>
#include <TGButtonGroup.h>
#include <TGSlider.h>
#include <TGNumberEntry.h>
#include <TGFrame.h>
#include <TGLabel.h>
#include <TButton.h>
#include <TCanvas.h>
#include <TMarker.h>
#include <TGraph.h>
#include <RQ_OBJECT.h>

#include "RecoParticleFlow/PFRootEvent/interface/PFRootEventManager.h"
#include "RecoParticleFlow/PFRootEvent/interface/DisplayManager.h"

#include <string>


class DialogFrame  : public TGMainFrame {
  RQ_OBJECT("DialogFrame")
  
    private:
  
  static const int EN=1;       // id threshold energy field
  static const int ENER=10;    // id threshold energy slider
  static const int MAXL=20;
  static const int PARTTYPE=30;
  static const int PFBLOCK=40;
   
    
  PFRootEventManager  *evMan_;
    
  DisplayManager      *display_;
  TGCompositeFrame    *mainFrame_;
  TGCompositeFrame    *cmdFrame_;
    

  TGCheckButton       *selectObject_[5];
  TGCheckButton       *printButton_[7];
  TGHSlider           *thresholdS_[5];
  TGNumberEntryField  *threshEntry_[5];
  TGNumberEntryField  *maxLineEntry_; 
  TGNumberEntryField  *particleTypeEntry_; 
  TGTextButton        *exitButton,*nextButton,*previousButton;
  TGTextButton        *reProcessButton;
  //int                  eventNr_;
  //int                  maxEvents_;
  //TButton             *Modify_;
  //TButton             *Cancel_;
  TCanvas             *attrView_;
  TMarker             *thisClusPattern_;
  TGraph              *trackPattern_;
  TGraph              *simplePartPattern_;   
    
  
 public:
  DialogFrame(PFRootEventManager *evman, DisplayManager *dm,const TGWindow *p,UInt_t w,UInt_t h);
  virtual ~DialogFrame(); 
  
  void closeCanvasAttr();   
  void createCmdFrame();
  void createCanvasAttr();
  void doDisplayGenParticle();
  void doLookForGenParticle();
  void doNextEvent();
  void doPreviousEvent();
  void doModifyOptions(unsigned obj);
  void doModifyPtThreshold(unsigned obj,long val);
  void isPFBlockVisible();
  void doPrint();
  void doPrintGenParticles();
  void doReProcessEvent();
  void selectPrintOption(int opt);
  //void modifyGraphicAttributes();
  void modifyAttr();
  void updateDisplay();
  void unZoom();

  virtual bool ProcessMessage(Long_t msg, Long_t parm1, Long_t parm2);
  virtual void CloseWindow();
     
};
#endif 
