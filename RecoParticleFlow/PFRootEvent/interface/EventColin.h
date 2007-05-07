#ifndef Demo_PFRootEvent_EventColin_h
#define Demo_PFRootEvent_EventColin_h

#include <vector>

class EventColin {
  
 public:
  EventColin() : number_(-1), nTracks_(0) {}
  virtual ~EventColin() {}
  void reset() {
    number_ = -1;
    particles_.clear();
    clusters_.clear();
    clustersIsland_.clear();
    jetsMC_.clear();
    jetsEHT_.clear();
    jetsPF_.clear();
    nTracks_=0;
  }

  class Particle {
  public:
    Particle() : eta(0),phi(0),e(0) {}
    virtual ~Particle() {}
    double eta; 
    double phi;
    double e;
  };   

  class Cluster {
  public:
    Cluster() : eta(0),phi(0),e(0),layer(0),type(0) {}
    virtual ~Cluster() {}
    double eta; 
    double phi;
    double e;
    int layer;
    int type;
  };   

  class Jet {
  public:
    Jet() : eta(0),phi(0),et(0),e(0) {}
    virtual ~Jet() {}
    double eta;
    double phi;
    double et;
    double e;
  };

  

  void setNumber(int number) {number_ = number;}
  void setNTracks(int nTracks) {nTracks_ = nTracks;}

  void addParticle( const Particle& ptc ) {
    particles_.push_back(ptc);
  }

  void addCluster( const Cluster& ptc ) {
    clusters_.push_back(ptc);
  }

  void addClusterIsland( const Cluster& ptc ) {
    clustersIsland_.push_back(ptc);
  }

  void addJetMC( const Jet & jets ) {
    jetsMC_.push_back(jets);
  }

  void addJetEHT( const Jet & jets ) {
    jetsEHT_.push_back(jets);
  }

  void addJetPF( const Jet & jets ) {
    jetsPF_.push_back(jets);
  }


  const std::vector<EventColin::Particle>& particles() 
    {return particles_;}
  const std::vector<EventColin::Cluster>&  clusters() 
    {return clusters_;}
  const std::vector<EventColin::Cluster>&  clustersIsland() 
    {return clustersIsland_;}
  const std::vector<EventColin::Jet>&      jetsMC()
    {return jetsMC_;}
  const std::vector<EventColin::Jet>&      jetsEHT()
    {return jetsEHT_;}
  const std::vector<EventColin::Jet>&      jetsPF()
    {return jetsPF_;}

 private:

  int                               number_;
  int                               nTracks_;
  std::vector<EventColin::Particle> particles_;
  std::vector<EventColin::Cluster>  clusters_;
  std::vector<EventColin::Cluster>  clustersIsland_;
  std::vector<EventColin::Jet>      jetsMC_;
  std::vector<EventColin::Jet>      jetsEHT_;
  std::vector<EventColin::Jet>      jetsPF_;
};

#endif
