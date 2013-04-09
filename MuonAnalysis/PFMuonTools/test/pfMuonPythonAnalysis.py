import ROOT
import itertools
from DataFormats.FWLite import Events, Handle

events = Events ('metOutput.root')

muonHandle  = Handle ('std::vector<reco::Muon>')
muonCleanedHandle  = Handle ('std::vector<reco::Muon>')
pfMuonHandle  = Handle ('std::vector<reco::PFCandidate>')
pfMuonCleanedHandle  = Handle ('std::vector<reco::PFCandidate>')
metHandle  = Handle ('std::vector<reco::PFMET>')
metCleanedHandle  = Handle ('std::vector<reco::PFMET>')

metOld = ROOT.TH1F('metOld','Met', 100,0,1000)
metNew = ROOT.TH1F('metNew','Met', 100,0,1000)
muOld = ROOT.TH1F('muOld','MuOld', 100,0,1000)
muNew = ROOT.TH1F('muNew','MuNew', 100,0,1000)


for event in events:
    event.getByLabel('muonsCleaned',muonCleanedHandle)
    event.getByLabel('muons',muonHandle)
    event.getByLabel('pfMuons',pfMuonCleanedHandle)
    event.getByLabel('oldPFMuons',pfMuonHandle)
    event.getByLabel('pfMet',metHandle)
    event.getByLabel('pfMetCleaned',metCleanedHandle)

    muons =muonHandle.product()
    met = metHandle.product()[0]
    muonsCleaned = muonCleanedHandle.product()
    metCleaned = metCleanedHandle.product()[0]
    pfMuons =pfMuonCleanedHandle.product()
    oldPFMuons =pfMuonHandle.product()



    metOld.Fill(met.pt())
    metNew.Fill(metCleaned.pt())

    for mu in muons:
        muOld.Fill(mu.pt())
    for mu in muonsCleaned:
        muNew.Fill(mu.pt())
        

    if abs(met.pt()-metCleaned.pt())/metCleaned.pt()>0.01:
        print '------DIFFERENCE-------'
        print 'Old MET',met.pt()
        print 'New MET',metCleaned.pt()
        print 'All Muons Before------'
        for muon in muons:
            print muon.pt(),muon.eta(),muon.phi(),muon.charge(),muon.isGlobalMuon(),muon.isTrackerMuon()
        print 'All Muons After reassignment------'
        for muon in muonsCleaned:
            print muon.pt(),muon.eta(),muon.phi(),muon.charge(),muon.isGlobalMuon(),muon.isTrackerMuon()

        print 'PF Muons Before------'
        for muon in oldPFMuons:
            print muon.pt(),muon.eta(),muon.phi(),muon.charge()
        print 'PF Muons After------'
        for muon in pfMuons:
            print muon.pt(),muon.eta(),muon.phi(),muon.charge()


        print 'Searching for split tracks  = SS low pass pairs'    
        sortedMuonsCleaned = sorted(muonsCleaned,key=lambda x: x.pt(),reverse=True)    
        if len(sortedMuonsCleaned)>1:
            print 'DiMuMassCleaned',(sortedMuonsCleaned[0].p4()+sortedMuonsCleaned[1].p4()).M()

        for muon1,muon2 in itertools.combinations(sortedMuonsCleaned,2):
            if muon1.charge()==muon2.charge():
                if abs(muon1.eta()-muon2.eta())<0.02:
                    if abs(muon1.phi()-muon2.phi())<0.02:
                        print 'FOUND SPLIT TRACK CANDIDATE',(muon1.p4()+muon2.p4()).M()

#        if metCleaned.pt()-met.pt()>200:
#            import pdb;pdb.set_trace()
    
f=ROOT.TFile("plots.root","RECREATE")
f.cd()
metOld.Write()
metNew.Write()
f.Close()
