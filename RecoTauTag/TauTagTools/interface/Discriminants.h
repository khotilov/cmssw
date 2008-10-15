#ifndef RecoTauTag_TauTagTools_PFTauDiscriminants
#define RecoTauTag_TauTagTools_PFTauDiscriminants

#include "RecoTauTag/TauTagTools/interface/PFTauDiscriminantBase.h"

namespace PFTauDiscriminants {
using namespace std;

typedef reco::Particle::LorentzVector LorentzVector;

//forward declarations

class DecayMode : public DiscriminantBase<int> {
   public:
      DecayMode():DiscriminantBase<int>("DecayMode", "I", true, false, -1){};
      ~DecayMode(){};
   protected:
      void doComputation(PFTauDiscriminantManager* input, vector<int>& result);
};

class OutlierNCharged : public DiscriminantBase<int> {
   public:
      OutlierNCharged():DiscriminantBase<int>("OutlierNCharged", "I", true, false, -1){};
      ~OutlierNCharged(){};
   protected:
      void doComputation(PFTauDiscriminantManager* input, vector<int>& result);
};


class Pt : public DiscriminantBase<double>  {
   public:
      Pt():DiscriminantBase<double>("Pt", "D", true, false, 0.0){};
      ~Pt(){};
   protected:
      void doComputation(PFTauDiscriminantManager* input, vector<double>& result);
};

class Eta : public DiscriminantBase<double>  {
   public:
      Eta():DiscriminantBase<double>("Eta", "D", true, false, 0.0){};
      ~Eta(){};
   protected:
      void doComputation(PFTauDiscriminantManager* input, vector<double>& result);
};

class MainTrackPt : public DiscriminantBase<double>  {
   public:
      MainTrackPt():DiscriminantBase<double>("MainTrackPt", "D", true, false, -1){};
      ~MainTrackPt(){};
   protected:
      void doComputation(PFTauDiscriminantManager* input, vector<double>& result);
};

class MainTrackAngle : public DiscriminantBase<double>  {
   public:
      MainTrackAngle():DiscriminantBase<double>("MainTrackAngle", "D", true, false, -1){};
      ~MainTrackAngle(){};
   protected:
      void doComputation(PFTauDiscriminantManager* input, vector<double>& result);
};

class TrackPt : public DiscriminantBase<double> {
   public:
      TrackPt():DiscriminantBase<double>("TrackPt", "vector<double>", false, true, 0.0){};
      ~TrackPt(){};
   protected:
      void doComputation(PFTauDiscriminantManager* input, vector<double>& result);
};

class PiZeroPt : public DiscriminantBase<double> {
   public:
      PiZeroPt():DiscriminantBase<double>("PiZeroPt", "vector<double>", false, true, 0.0){};
      ~PiZeroPt(){};
   protected:
      void doComputation(PFTauDiscriminantManager* input, vector<double>& result);
};

class TrackAngle : public DiscriminantBase<double> {
   public:
      TrackAngle():DiscriminantBase<double>("TrackAngle", "vector<double>", false, true, 0.0){};
      ~TrackAngle(){};
   protected:
      void doComputation(PFTauDiscriminantManager* input, vector<double>& result);
};

class PiZeroAngle : public DiscriminantBase<double> {
   public:
      PiZeroAngle():DiscriminantBase<double>("PiZeroAngle", "vector<double>", false, true, 0.0){};
      ~PiZeroAngle(){};
   protected:
      void doComputation(PFTauDiscriminantManager* input, vector<double>& result);
};

class Dalitz : public DiscriminantBase<double> {
   public:
      Dalitz():DiscriminantBase<double>("Dalitz", "vector<double>", false, true, 0.0){};
      ~Dalitz(){};
   protected:
      void doComputation(PFTauDiscriminantManager* input, vector<double>& result);
};

// takes invariant mass of all objects in signal cone
class InvariantMassOfSignal : public DiscriminantBase<double> {
   public:
      InvariantMassOfSignal():DiscriminantBase<double>("InvariantMassOfSignal", "D", true, false, 0.0){};
      ~InvariantMassOfSignal(){};
   protected:
      void doComputation(PFTauDiscriminantManager* input, vector<double>& result);
};

// returns vector of invariant masses of larger and larger subsets of all signal objects e.g. result[2] is
// the invariant mass of the lead track with the next highest Pt object

class InvariantMass : public DiscriminantBase<double> {
   public:
      InvariantMass():DiscriminantBase<double>("InvariantMass", "vector<double>", false, true, 0.0){};
      ~InvariantMass(){};
   protected:
      void doComputation(PFTauDiscriminantManager* input, vector<double>& result);
};

class OutlierPt : public DiscriminantBase<double> {
   public:
      OutlierPt():DiscriminantBase<double>("OutlierPt", "vector<double>", false, true, 0.0){};
      ~OutlierPt(){};
   protected:
      void doComputation(PFTauDiscriminantManager* input, vector<double>& result);
};

class OutlierAngle : public DiscriminantBase<double> {
   public:
      OutlierAngle():DiscriminantBase<double>("OutlierAngle", "vector<double>", false, true, 0.0){};
      ~OutlierAngle(){};
   protected:
      void doComputation(PFTauDiscriminantManager* input, vector<double>& result);
};

class ChargedOutlierPt : public DiscriminantBase<double> {
   public:
      ChargedOutlierPt():DiscriminantBase<double>("ChargedOutlierPt", "vector<double>", false, true, 0.0){};
      ~ChargedOutlierPt(){};
   protected:
      void doComputation(PFTauDiscriminantManager* input, vector<double>& result);
};

class ChargedOutlierAngle : public DiscriminantBase<double> {
   public:
      ChargedOutlierAngle():DiscriminantBase<double>("ChargedOutlierAngle", "vector<double>", false, true, 0.0){};
      ~ChargedOutlierAngle(){};
   protected:
      void doComputation(PFTauDiscriminantManager* input, vector<double>& result);
};

class NeutralOutlierPt : public DiscriminantBase<double> {
   public:
      NeutralOutlierPt():DiscriminantBase<double>("NeutralOutlierPt", "vector<double>", false, true, 0.0){};
      ~NeutralOutlierPt(){};
   protected:
      void doComputation(PFTauDiscriminantManager* input, vector<double>& result);
};

class NeutralOutlierAngle : public DiscriminantBase<double> {
   public:
      NeutralOutlierAngle():DiscriminantBase<double>("NeutralOutlierAngle", "vector<double>", false, true, 0.0){};
      ~NeutralOutlierAngle(){};
   protected:
      void doComputation(PFTauDiscriminantManager* input, vector<double>& result);
};


}
#endif





