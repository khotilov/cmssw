// -*- C++ -*-
#ifndef Fireworks_Core_CmsShowMain_h
#define Fireworks_Core_CmsShowMain_h
//
// Package:     Core
// Class  :     CmsShowMain
// 
/**\class CmsShowMain CmsShowMain.h Fireworks/Core/interface/CmsShowMain.h

 Description: Displays an fwlite::Event in ROOT

 Usage:
    <usage>

*/
//
// Original Author:  
//         Created:  Mon Dec  3 08:34:30 PST 2007
// $Id: CmsShowMain.h,v 1.6 2008/08/20 21:02:05 chrjones Exp $
//

// system include files
#include <vector>
#include <string>
#include <memory>
#include <boost/shared_ptr.hpp>
#include "Rtypes.h"

// user include files
#include "Fireworks/Core/interface/DetIdToMatrix.h"

// forward declarations
class TGPictureButton;
class TGComboBox;
class TGTextButton;
class TGTextEntry;
class FWEventItemsManager;
class FWViewManagerManager;
class FWModelChangeManager;
class FWSelectionManager;
class FWGUIManager;
class FWEventItem;
class FWPhysicsObjectDesc;
class FWConfigurationManager;
class FWTextView;
class TTimer;
class TMonitor;
class TSocket;

class CmsShowNavigator;
class CmsShowTaskExecutor;

namespace fwlite {
  class Event;
}

class FWDetailView;

class CmsShowMain
{
public:
  CmsShowMain(int argc, char *argv[]);
  virtual ~CmsShowMain();
  void resetInitialization();
  void draw(const fwlite::Event& event);
  void openData();
  void quit();

  // ---------- const member functions ---------------------
  const DetIdToMatrix& getIdToGeo() const { return m_detIdToGeo; }
   
  //  void writeConfigurationFile(const std::string& iFileName) const;
  // ---------- static member functions --------------------
  static void   setAutoFieldMode(bool state) { m_autoField = state; }
  static bool   isAutoField() { return m_autoField; }
  static double getMagneticField();
  static void   setMagneticField(double var);    
  static int    getFieldEstimates() { return m_numberOfFieldEstimates; }
  static void   guessFieldIsOn( bool guess );
  static double getCaloScale() { return m_caloScale; }
  static void   setCaloScale(double var) { m_caloScale = var; }
  static void setBrightness(unsigned int value);

  // ---------- member functions ---------------------------
  //  int draw(const fwlite::Event& );

  void registerPhysicsObject(const FWPhysicsObjectDesc&);
  void registerDetailView (const std::string &item_name, FWDetailView *view);
   
  void notified(TSocket*);
private:
  CmsShowMain(const CmsShowMain&); // stop default
  
  const CmsShowMain& operator=(const CmsShowMain&); // stop default
   
   void loadGeometry();
   void setupViewManagers();
   void setupConfiguration();
   void setupDetailedViewManagers();
   void setupDataHandling();
   void setupDebugSupport();
   void setupSocket(unsigned int);
  
   void playForward();
   void playBackward();
   void stopPlaying();
   void reachedEnd();
  // ---------- member data --------------------------------
  std::auto_ptr<FWConfigurationManager> m_configurationManager;
  std::auto_ptr<FWModelChangeManager> m_changeManager;
  std::auto_ptr<FWSelectionManager> m_selectionManager;
  std::auto_ptr<FWEventItemsManager> m_eiManager;
  std::auto_ptr<FWViewManagerManager> m_viewManager;
  std::auto_ptr<FWGUIManager> m_guiManager;
  std::auto_ptr<FWTextView> m_textView;

  CmsShowNavigator* m_navigator;

  DetIdToMatrix    m_detIdToGeo;
  std::string m_inputFileName;
  std::string m_configFileName;
  std::string m_geomFileName;
  static bool   m_autoField;                   // data derived magnetif field state
  static double m_magneticField;               
  static int    m_numberOfFieldEstimates;      
  static int    m_numberOfFieldIsOnEstimates;  
  static double m_caloScale;
   
   std::auto_ptr<CmsShowTaskExecutor> m_startupTasks;
   
   TTimer* m_playTimer;
   TTimer* m_playBackTimer;
   bool m_isPlaying;
   bool m_forward;
   Long_t m_playDelay;
   
   std::auto_ptr<TMonitor> m_monitor;
};


#endif
