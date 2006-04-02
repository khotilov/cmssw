#ifndef _MaterialBudgetAction_h
#define _MaterialBudgetAction_h
#include <string>
#include <vector>
#include <map>
 
// user include files
#include "Validation/Geometry/interface/MaterialBudgetTree.h"
#include "Validation/Geometry/interface/MaterialBudgetHistos.h"
#include "Validation/Geometry/interface/MaterialBudgetTxt.h"

#include "SimG4Core/Watcher/interface/SimProducer.h"
#include "SimG4Core/Notification/interface/Observer.h"

#include <CLHEP/Vector/LorentzVector.h>

using namespace std;

class BeginOfTrack;
class G4Step;
class EndOfTrack;
class G4StepPoint;

class MaterialBudgetAction : public SimProducer, 
			     public Observer<const BeginOfTrack*>,
			     public Observer<const G4Step*>,
			     public Observer<const EndOfTrack*>
{
   public:
      MaterialBudgetAction(const edm::ParameterSet&);
      virtual ~MaterialBudgetAction();

      void produce(edm::Event&, const edm::EventSetup&);


   private:
      MaterialBudgetAction(const MaterialBudgetAction&); // stop default

      const MaterialBudgetAction& operator=(const MaterialBudgetAction&); // stop default

      void update(const BeginOfTrack*);
      void update(const G4Step*);
      void update(const EndOfTrack*);

      void initRun();
      void processEvent( uint nEv );
      void endRun();
   
   private:
  void save( const G4Step* aStep );
  std::string getSubDetectorName( G4StepPoint* aStepPoint );
  std::string getPartName( G4StepPoint* aStepPoint );
  MaterialBudgetData* theData;
  MaterialBudgetTree* theTree;
  MaterialBudgetHistos* theHistos;
  MaterialBudgetTxt* theTxt;
  bool saveToTxt, saveToTree, saveToHistos;

};

#endif
