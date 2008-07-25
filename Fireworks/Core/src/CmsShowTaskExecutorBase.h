#ifndef Fireworks_Core_CmsShowTaskExecutorBase_h
#define Fireworks_Core_CmsShowTaskExecutorBase_h
// -*- C++ -*-
//
// Package:     Core
// Class  :     CmsShowTaskExecutorBase
// 
/**\class CmsShowTaskExecutorBase CmsShowTaskExecutorBase.h Fireworks/Core/interface/CmsShowTaskExecutorBase.h

 Description: <one line class summary>

 Usage:
    <usage>

*/
//
// Original Author:  Chris Jones
//         Created:  Fri Jul 11 12:09:38 EDT 2008
// $Id: CmsShowTaskExecutorBase.h,v 1.2 2008/07/13 21:56:08 chrjones Exp $
//

// system include files
#include <RQ_OBJECT.h>
#include <sigc++/signal.h>

// user include files

// forward declarations
class TTimer;
class CmsShowTaskTimer;

class CmsShowTaskExecutorBase : public TQObject {
   RQ_OBJECT("CmsShowTaskExecutorBase");

public:
   CmsShowTaskExecutorBase();
   virtual ~CmsShowTaskExecutorBase();
   
   // ---------- const member functions ---------------------
   
   // ---------- static member functions --------------------
   
   // ---------- member functions ---------------------------
   void requestNextTask();
   void doNextTask();
   
   virtual void startDoingTasks()=0;

   sigc::signal<void> tasksCompleted_;

protected:
   virtual void doNextTaskImp() = 0;
   virtual bool moreTasksAvailable() = 0;
private:
   CmsShowTaskExecutorBase(const CmsShowTaskExecutorBase&); // stop default
   
   const CmsShowTaskExecutorBase& operator=(const CmsShowTaskExecutorBase&); // stop default
   
   // ---------- member data --------------------------------
   //TTimer* m_timer;
   CmsShowTaskTimer* m_taskTimer;
};


#endif
