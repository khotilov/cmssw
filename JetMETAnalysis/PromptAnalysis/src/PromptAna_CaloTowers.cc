#include "JetMETAnalysis/PromptAnalysis/interface/PromptAna_CaloTowers.h"
#include "FWCore/Framework/interface/Event.h"

PromptAna_CaloTowers::PromptAna_CaloTowers(const edm::ParameterSet& iConfig) 
  : inputTag(iConfig.getParameter<edm::InputTag>("InputTag"))
  , prefix  (iConfig.getParameter<std::string>  ("Prefix"  ))
  , suffix  (iConfig.getParameter<std::string>  ("Suffix"  ))
{
  produces <std::vector<double> > ( prefix + "EmEt"  + suffix );
  produces <std::vector<double> > ( prefix + "HadEt"  + suffix );
  produces <std::vector<double> > ( prefix + "OuterEt"  + suffix );
  produces <std::vector<double> > ( prefix + "Eta"  + suffix );
  produces <std::vector<double> > ( prefix + "Phi"  + suffix );
  produces <std::vector<double> > ( prefix + "EcalTime"  + suffix );
  produces <std::vector<double> > ( prefix + "HcalTime"  + suffix );
  produces <std::vector<int> >    ( prefix + "Ieta"  + suffix );
  produces <std::vector<int> >    ( prefix + "Iphi"  + suffix );
  produces <std::vector<int> >    ( prefix + "TowerStatusWord"  + suffix );
}

void PromptAna_CaloTowers::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) 
{
  std::auto_ptr<std::vector<double> >  emEt               ( new std::vector<double>()  ) ;
  std::auto_ptr<std::vector<double> >  hadEt              ( new std::vector<double>()  ) ;
  std::auto_ptr<std::vector<double> >  outerEt            ( new std::vector<double>()  ) ;
  std::auto_ptr<std::vector<double> >  eta                ( new std::vector<double>()  ) ;
  std::auto_ptr<std::vector<double> >  phi                ( new std::vector<double>()  ) ;
  std::auto_ptr<std::vector<double> >  ecalTime           ( new std::vector<double>()  ) ;
  std::auto_ptr<std::vector<double> >  hcalTime           ( new std::vector<double>()  ) ;
  std::auto_ptr<std::vector<int> >     ieta               ( new std::vector<int>()  ) ;
  std::auto_ptr<std::vector<int> >     iphi               ( new std::vector<int>()  ) ;
  std::auto_ptr<std::vector<int> >     towerStatusWord    ( new std::vector<int>()  ) ;

  //Get the CaloTower Collection
  edm::Handle<CaloTowerCollection> calotowercollection;
  iEvent.getByLabel(inputTag, calotowercollection);
  
  //Fill the variables
  for(CaloTowerCollection::const_iterator it = calotowercollection->begin(); it != calotowercollection->end() ; ++it )
    {
      emEt                -> push_back(it->emEt());
      hadEt               -> push_back(it->hadEt());
      outerEt             -> push_back(it->outerEt());
      eta                 -> push_back(it->p4().eta());
      phi                 -> push_back(it->p4().phi());
      ieta                -> push_back(it->ieta());
      iphi                -> push_back(it->iphi());
      ecalTime            -> push_back(it->ecalTime());
      hcalTime            -> push_back(it->hcalTime());
      towerStatusWord     -> push_back(it->towerStatusWord());
    }

  iEvent.put( emEt                   ,  prefix + "EmEt" + suffix );
  iEvent.put( hadEt                  ,  prefix + "HadEt" + suffix );
  iEvent.put( outerEt                ,  prefix + "OuterEt" + suffix );
  iEvent.put( eta                    ,  prefix + "Eta" + suffix );
  iEvent.put( phi                    ,  prefix + "Phi" + suffix );
  iEvent.put( ieta                   ,  prefix + "Ieta"  + suffix );
  iEvent.put( iphi                   ,  prefix + "Iphi"  + suffix );
  iEvent.put( ecalTime               ,  prefix + "EcalTime"  + suffix );
  iEvent.put( hcalTime               ,  prefix + "HcalTime"  + suffix );
  iEvent.put( towerStatusWord        ,  prefix + "TowerStatusWord"  + suffix );
}
