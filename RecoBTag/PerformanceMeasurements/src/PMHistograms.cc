#include "RecoBTag/PerformanceMeasurements/interface/PMHistograms.h"

#include "TH2F.h"
#include "TString.h"

//______________
void PMHistograms::Add()
{

    std::map<std::string, int>  quark_color;
    quark_color[""] = 1;
    quark_color["b"] = 2;
    quark_color["c"] = 3;
    quark_color["uds"] = 4;
    quark_color["g"] = 6;

    TString           ftagger;
    TString           flevel;
    TString           fAwaytagger;
    TString           fAwaylevel;
    std::map< TString, float > fTrackCountingMap;
    TAxis             fJetPtAxis;
    TAxis             fJetEtaAxis;
    TAxis             fCorrPtAxis;
    TAxis             fCorrEtaAxis;

    ftagger = "TrackCounting";
    flevel  = "Loose";
    fAwaytagger = "TrackCounting";
    fAwaylevel = "Loose";

    fTrackCountingMap["Loose"]  = 2.0; // use TC2:high eff.
    fTrackCountingMap["Medium"] = 4.2; // use TC2:high eff.
    fTrackCountingMap["Tight"]  = 4.1;

    const int nptarray = 4;
    const int netaarray = 10;
    const int ncorrptarray = 5;
    const int ncorretaarray = 5;
    Double_t jetptbins[nptarray] = {30., 70., 120.,230.};
    Double_t jetetabins[netaarray] = {0.0,0.25,0.5,0.75,1.0,1.25,1.5,1.75,2.0,2.5};
    Double_t corrptbins[ncorrptarray] = {20.,40.,60.,80.,230.};
    Double_t corretabins[ncorrptarray] = {0.,0.5,1.,1.5,2.5};

    int nptbins = fJetPtAxis.GetNbins();;
    //const Double_t *jetptbins = (fJetPtAxis.GetXbins())->GetArray();
    int netabins = fJetEtaAxis.GetNbins();
    //const Double_t *jetetabins = (fJetEtaAxis.GetXbins())->GetArray();
    // int ncorrptbins = fCorrPtAxis.GetNbins();
    // const Double_t *corrptbins = (fCorrPtAxis.GetXbins())->GetArray();
    // int ncorretabins = fCorrEtaAxis.GetNbins();
    // const Double_t *corretabins = (fCorrEtaAxis.GetXbins())->GetArray();



    fstore->add( new TH2F("n_pT","MuTag pT vs pTrel",nptarray,jetptbins,50,0.,5.), "muon_in_jet" );
    fstore->add( new TH2F("p_pT","MuTag && CMBtag pT vs pTrel",nptarray,jetptbins,50,0.,5.), "muon_in_jet" );
    fstore->add( new TH2F("ntag_pT","opp tag: MuTag pT vs pTrel",nptarray,jetptbins,50,0.,5.), "muon_in_jet" );
    fstore->add( new TH2F("ptag_pT","opp tag MuTag && CMBtag pT vs pTrel",nptarray,jetptbins,50,0.,5.), "muon_in_jet" );
    fstore->add( new TH2F("nnoTag_pT","opp tag: MuTag pT vs pTrel",nptarray,jetptbins,50,0.,5.), "muon_in_jet" );
    fstore->add( new TH2F("pnoTag_pT","opp tag MuTag && CMBtag pT vs pTrel",nptarray,jetptbins,50,0.,5.), "muon_in_jet" );

    fstore->add( new TH2F("q_pT","other MuTag pT vs pTrel",nptarray,jetptbins,50,0.,5.), "muon_in_jet" );
    fstore->add( new TH2F("qtag_pT","other MuTag && Tagger pT vs pTrel",nptarray,jetptbins,50,0.,5.), "muon_in_jet" );


    fstore->add( new TH2F("n_eta","MuTag Eta vs pTrel",netaarray,jetetabins,50,0.,5.), "muon_in_jet" );
    fstore->add( new TH2F("p_eta","MuTag && CMBtag Eta vs pTrel",netaarray,jetetabins,50,0.,5.), "muon_in_jet" );
    fstore->add( new TH2F("ntag_eta","opp tag: MuTag Eta vs pTrel",netaarray,jetetabins,50,0.,5.), "muon_in_jet" );
    fstore->add( new TH2F("ptag_eta","opp tag MuTag && CMBtag Eta vs pTrel",netaarray,jetetabins,50,0.,5.), "muon_in_jet" );
    fstore->add( new TH2F("nnoTag_eta","opp tag: MuTag Eta vs pTrel",netaarray,jetetabins,50,0.,5.), "muon_in_jet" );
    fstore->add( new TH2F("pnoTag_eta","opp tag MuTag && CMBtag Eta vs pTrel",netaarray,jetetabins,50,0.,5.), "muon_in_jet" );

    fstore->add( new TH2F("q_eta","other MuTag pT vs pTrel",netaarray,jetetabins,50,0.,5.), "muon_in_jet" );
    fstore->add( new TH2F("qtag_eta","other MuTag && Tagger pT vs pTrel",netaarray,jetetabins,50,0.,5.), "muon_in_jet" );
//
    fstore->add( new TH2F("n_pT_b","b MuTag pT vs pTrel",nptarray,jetptbins,50,0.,5.), "MCTruth" );
    fstore->add( new TH2F("p_pT_b","b MuTag && CMBtag pT vs pTrel",nptarray,jetptbins,50,0.,5.), "MCTruth" );
    fstore->add( new TH2F("ntag_pT_b","b opp tag: MuTag pT vs pTrel",nptarray,jetptbins,50,0.,5.), "MCTruth" );
    fstore->add( new TH2F("ptag_pT_b","b opp tag MuTag && CMBtag pT vs pTrel",nptarray,jetptbins,50,0.,5.), "MCTruth" );
    fstore->add( new TH2F("nnoTag_pT_b","b opp tag: MuTag pT vs pTrel",nptarray,jetptbins,50,0.,5.), "MCTruth" );
    fstore->add( new TH2F("pnoTag_pT_b","b opp tag MuTag && CMBtag pT vs pTrel",nptarray,jetptbins,50,0.,5.), "MCTruth" );

    fstore->add( new TH2F("q_pT_b","other MuTag pT vs pTrel",nptarray,jetptbins,50,0.,5.), "MCTruth" );
    fstore->add( new TH2F("qtag_pT_b","other MuTag && Tagger pT vs pTrel",nptarray,jetptbins,50,0.,5.), "MCTruth" );

    fstore->add( new TH2F("n_eta_b","b MuTag Eta vs pTrel",netaarray,jetetabins,50,0.,5.), "MCTruth" );
    fstore->add( new TH2F("p_eta_b","b MuTag && CMBtag Eta vs pTrel",netaarray,jetetabins,50,0.,5.), "MCTruth" );
    fstore->add( new TH2F("ntag_eta_b","b opp tag: MuTag Eta vs pTrel",netaarray,jetetabins,50,0.,5.), "MCTruth" );
    fstore->add( new TH2F("ptag_eta_b","b opp tag MuTag && CMBtag Eta vs pTrel",netaarray,jetetabins,50,0.,5.), "MCTruth" );
    fstore->add( new TH2F("nnoTag_eta_b","b opp tag: MuTag Eta vs pTrel",netaarray,jetetabins,50,0.,5.), "MCTruth" );
    fstore->add( new TH2F("pnoTag_eta_b","b opp tag MuTag && CMBtag Eta vs pTrel",netaarray,jetetabins,50,0.,5.), "MCTruth" );

    fstore->add( new TH2F("q_eta_b","other MuTag pT vs pTrel",netaarray,jetetabins,50,0.,5.), "MCTruth" );
    fstore->add( new TH2F("qtag_eta_b","other MuTag && Tagger pT vs pTrel",netaarray,jetetabins,50,0.,5.), "MCTruth" );
//
    fstore->add( new TH2F("n_pT_cl","cl MuTag pT vs pTrel",nptarray,jetptbins,50,0.,5.), "MCTruth" );
    fstore->add( new TH2F("p_pT_cl","cl MuTag && CMBtag pT vs pTrel",nptarray,jetptbins,50,0.,5.), "MCTruth" );
    fstore->add( new TH2F("ntag_pT_cl","cl opp tag: MuTag pT vs pTrel",nptarray,jetptbins,50,0.,5.), "MCTruth" );
    fstore->add( new TH2F("ptag_pT_cl","cl opp tag MuTag && CMBtag pT vs pTrel",nptarray,jetptbins,50,0.,5.), "MCTruth" );
    fstore->add( new TH2F("nnoTag_pT_cl","cl opp tag: MuTag pT vs pTrel",nptarray,jetptbins,50,0.,5.), "MCTruth" );
    fstore->add( new TH2F("pnoTag_pT_cl","cl opp tag MuTag && CMBtag pT vs pTrel",nptarray,jetptbins,50,0.,5.), "MCTruth" );

    fstore->add( new TH2F("q_pT_cl","other MuTag pT vs pTrel",nptarray,jetptbins,50,0.,5.), "MCTruth" );
    fstore->add( new TH2F("qtag_pT_cl","other MuTag && Tagger pT vs pTrel",nptarray,jetptbins,50,0.,5.), "MCTruth" );

    fstore->add( new TH2F("n_eta_cl","cl MuTag Eta vs pTrel",netaarray,jetetabins,50,0.,5.), "MCTruth" );
    fstore->add( new TH2F("p_eta_cl","cl MuTag && CMBtag Eta vs pTrel",netaarray,jetetabins,50,0.,5.), "MCTruth" );
    fstore->add( new TH2F("ntag_eta_cl","cl opp tag: MuTag Eta vs pTrel",netaarray,jetetabins,50,0.,5.), "MCTruth" );
    fstore->add( new TH2F("ptag_eta_cl","cl opp tag MuTag && CMBtag Eta vs pTrel",netaarray,jetetabins,50,0.,5.), "MCTruth" );
    fstore->add( new TH2F("nnoTag_eta_cl","cl opp tag: MuTag Eta vs pTrel",netaarray,jetetabins,50,0.,5.), "MCTruth" );
    fstore->add( new TH2F("pnoTag_eta_cl","cl opp tag MuTag && CMBtag Eta vs pTrel",netaarray,jetetabins,50,0.,5.), "MCTruth" );

    fstore->add( new TH2F("q_eta_cl","other MuTag pT vs pTrel",netaarray,jetetabins,50,0.,5.), "MCTruth" );
    fstore->add( new TH2F("qtag_eta_cl","other MuTag && Tagger pT vs pTrel",netaarray,jetetabins,50,0.,5.), "MCTruth" );
//
    fstore->add( new TH2F("n_pT_c","c MuTag pT vs pTrel",nptarray,jetptbins,50,0.,5.), "MCTruth" );
    fstore->add( new TH2F("p_pT_c","c MuTag && CMBtag pT vs pTrel",nptarray,jetptbins,50,0.,5.), "MCTruth" );
    fstore->add( new TH2F("ntag_pT_c","c opp tag: MuTag pT vs pTrel",nptarray,jetptbins,50,0.,5.), "MCTruth" );
    fstore->add( new TH2F("ptag_pT_c","c opp tag MuTag && CMBtag pT vs pTrel",nptarray,jetptbins,50,0.,5.), "MCTruth" );
    fstore->add( new TH2F("nnoTag_pT_c","c opp tag: MuTag pT vs pTrel",nptarray,jetptbins,50,0.,5.), "MCTruth" );
    fstore->add( new TH2F("pnoTag_pT_c","c opp tag MuTag && CMBtag pT vs pTrel",nptarray,jetptbins,50,0.,5.), "MCTruth" );

    fstore->add( new TH2F("q_pT_c","other MuTag pT vs pTrel",nptarray,jetptbins,50,0.,5.), "MCTruth" );
    fstore->add( new TH2F("qtag_pT_c","other MuTag && Tagger pT vs pTrel",nptarray,jetptbins,50,0.,5.), "MCTruth" );

    fstore->add( new TH2F("n_eta_c","c MuTag Eta vs pTrel",netaarray,jetetabins,50,0.,5.), "MCTruth" );
    fstore->add( new TH2F("p_eta_c","c MuTag && CMBtag Eta vs pTrel",netaarray,jetetabins,50,0.,5.), "MCTruth" );
    fstore->add( new TH2F("ntag_eta_c","c opp tag: MuTag Eta vs pTrel",netaarray,jetetabins,50,0.,5.), "MCTruth" );
    fstore->add( new TH2F("ptag_eta_c","c opp tag MuTag && CMBtag Eta vs pTrel",netaarray,jetetabins,50,0.,5.), "MCTruth" );
    fstore->add( new TH2F("nnoTag_eta_c","c opp tag: MuTag Eta vs pTrel",netaarray,jetetabins,50,0.,5.), "MCTruth" );
    fstore->add( new TH2F("pnoTag_eta_c","c opp tag MuTag && CMBtag Eta vs pTrel",netaarray,jetetabins,50,0.,5.), "MCTruth" );

    fstore->add( new TH2F("q_eta_c","other MuTag pT vs pTrel",netaarray,jetetabins,50,0.,5.), "MCTruth" );
    fstore->add( new TH2F("qtag_eta_c","other MuTag && Tagger pT vs pTrel",netaarray,jetetabins,50,0.,5.), "MCTruth" );
//
    fstore->add( new TH2F("n_pT_l","l MuTag pT vs pTrel",nptarray,jetptbins,50,0.,5.), "MCTruth" );
    fstore->add( new TH2F("p_pT_l","l MuTag && CMBtag pT vs pTrel",nptarray,jetptbins,50,0.,5.), "MCTruth" );
    fstore->add( new TH2F("ntag_pT_l","l opp tag: MuTag pT vs pTrel",nptarray,jetptbins,50,0.,5.), "MCTruth" );
    fstore->add( new TH2F("ptag_pT_l","l opp tag MuTag && CMBtag pT vs pTrel",nptarray,jetptbins,50,0.,5.), "MCTruth" );
    fstore->add( new TH2F("nnoTag_pT_l","l opp tag: MuTag pT vs pTrel",nptarray,jetptbins,50,0.,5.), "MCTruth" );
    fstore->add( new TH2F("pnoTag_pT_l","l opp tag MuTag && CMBtag pT vs pTrel",nptarray,jetptbins,50,0.,5.), "MCTruth" );

    fstore->add( new TH2F("q_pT_l","other MuTag pT vs pTrel",nptarray,jetptbins,50,0.,5.), "MCTruth" );
    fstore->add( new TH2F("qtag_pT_l","other MuTag && Tagger pT vs pTrel",nptarray,jetptbins,50,0.,5.), "MCTruth" );

    fstore->add( new TH2F("n_eta_l","l MuTag Eta vs pTrel",netaarray,jetetabins,50,0.,5.), "MCTruth" );
    fstore->add( new TH2F("p_eta_l","l MuTag && CMBtag Eta vs pTrel",netaarray,jetetabins,50,0.,5.), "MCTruth" );
    fstore->add( new TH2F("ntag_eta_l","l opp tag: MuTag Eta vs pTrel",netaarray,jetetabins,50,0.,5.), "MCTruth" );
    fstore->add( new TH2F("ptag_eta_l","l opp tag MuTag && CMBtag Eta vs pTrel",netaarray,jetetabins,50,0.,5.), "MCTruth" );
    fstore->add( new TH2F("nnoTag_eta_l","l opp tag: MuTag Eta vs pTrel",netaarray,jetetabins,50,0.,5.), "MCTruth" );
    fstore->add( new TH2F("pnoTag_eta_l","l opp tag MuTag && CMBtag Eta vs pTrel",netaarray,jetetabins,50,0.,5.), "MCTruth" );

    fstore->add( new TH2F("q_eta_l","other MuTag pT vs pTrel",netaarray,jetetabins,50,0.,5.), "MCTruth" );
    fstore->add( new TH2F("qtag_eta_l","other MuTag && Tagger pT vs pTrel",netaarray,jetetabins,50,0.,5.), "MCTruth" );


    /*
    	fstore->add( new TH2F("ptVsEta","pt vs Eta",fJetPtAxis.GetNbins() , (fJetPtAxis.GetXbins())->GetArray(), fJetEtaAxis.GetNbins() , (fJetEtaAxis.GetXbins())->GetArray() ) );
    	fstore->add( new TH2F("ptVsEta_b","pt vs Eta",fJetPtAxis.GetNbins() , (fJetPtAxis.GetXbins())->GetArray(), fJetEtaAxis.GetNbins() , (fJetEtaAxis.GetXbins())->GetArray() ) );
    	fstore->add( new TH2F("taggedjet_ptVsEta","taggedjet pt vs Eta",fJetPtAxis.GetNbins() , (fJetPtAxis.GetXbins())->GetArray(), fJetEtaAxis.GetNbins() , (fJetEtaAxis.GetXbins())->GetArray() ) );
    	fstore->add( new TH2F("taggedjet_ptVsEta_b","taggedjet pt vs Eta",fJetPtAxis.GetNbins() , (fJetPtAxis.GetXbins())->GetArray(), fJetEtaAxis.GetNbins() , (fJetEtaAxis.GetXbins())->GetArray() ) );

    	fstore->add( new TH1D("alpha","alpha",nptarray,jetptbins), "MCTruth" );
    	fstore->add( new TH1D("beta","beta",nptarray,jetptbins), "MCTruth" );
    	fstore->add( new TH1D("kappa_cl","kappa_cl",nptarray,jetptbins), "MCTruth" );
    	fstore->add( new TH1D("kappa_b","kappa_b",nptarray,jetptbins), "MCTruth" );
    	fstore->add( new TH1D("delta","delta",nptarray,jetptbins), "MCTruth" );
    	fstore->add( new TH1D("gamma","gamma",nptarray,jetptbins), "MCTruth" );

    	fstore->add( new TH1D("alpha_eta","alpha_eta",netaarray,jetetabins), "MCTruth" );
    	fstore->add( new TH1D("beta_eta","beta_eta",netaarray,jetetabins), "MCTruth" );
    	fstore->add( new TH1D("kappa_eta_cl","kappa_eta_cl",netaarray,jetetabins), "MCTruth" );
    	fstore->add( new TH1D("kappa_eta_b","kappa_eta_b",netaarray,jetetabins), "MCTruth" );
    	fstore->add( new TH1D("delta_eta","delta_eta",netaarray,jetetabins), "MCTruth" );
    	fstore->add( new TH1D("gamma_eta","gamma_eta",netaarray,jetetabins), "MCTruth" );
    	*/
    fstore->add( new TH1D("jet_deltaR","#Delta R",60,0.,0.55) );
    fstore->add( new TH1D("jet_deltaR_b","#Delta R",60,0.,0.55) );
    fstore->add( new TH1D("jet_deltaR_c","#Delta R",60,0.,0.55) );
    fstore->add( new TH1D("jet_deltaR_udsg","#Delta R",60,0.,0.55) );
    // pending change colors
    fstore->add( new TH1D("jet_pTrel","p_{Trel} [GeV/c]" , 50, 0, 5 ) );
    fstore->add( new TH1D("jet_pTrel_b","p_{Trel} [GeV/c]" , 50, 0, 5 ) );
    fstore->add( new TH1D("jet_pTrel_c","p_{Trel} [GeV/c]" , 50, 0, 5 ) );
    fstore->add( new TH1D("jet_pTrel_udsg","p_{Trel} [GeV/c]" , 50, 0, 5 ) );

    fstore->add( new TH1F( "jet_pt", "jet pt", 30, 0, 150) );
    fstore->add( new TH1F( "muon_pt", "muon pt", 300, 0, 50) );
    fstore->add( new TH1F( "ptRel", "ptRel", 100, 0, 10) );


}

//______________________________________________________________________________________________________________________
//void PMHistograms::FillHistos(std::string type, TLorentzVector p4MuJet, double ptrel,
//			      int JetFlavor, std::map<std::string, bool> aMap)
void PMHistograms::FillHistos(std::string type, TLorentzVector p4MuJet, double ptrel,
                              int JetFlavor, bool tagged)
{
    //type == "n", "p"

    fstore->hist(type+"_pT")->Fill(p4MuJet.Pt(),ptrel);
    fstore->hist(type+"_eta")->Fill(TMath::Abs(p4MuJet.Eta()),ptrel);
    std::string flavor;
    if ( JetFlavor == 5 )
    {
        flavor = "b";
        fstore->hist(type+"_pT_"+flavor)->Fill(p4MuJet.Pt(),ptrel);
        fstore->hist(type+"_eta_"+flavor)->Fill(TMath::Abs(p4MuJet.Eta()),ptrel);
        if (tagged)
        {
            fstore->hist(type+"tag_pT_"+flavor)->Fill(p4MuJet.Pt(),ptrel);
            fstore->hist(type+"tag_eta_"+flavor)->Fill(TMath::Abs(p4MuJet.Eta()),ptrel);
        }
        else
        {

            fstore->hist(type+"noTag_pT_"+flavor)->Fill(p4MuJet.Pt(),ptrel);
            fstore->hist(type+"noTag_eta_"+flavor)->Fill(TMath::Abs(p4MuJet.Eta()),ptrel);
        }

    }
    if ( (JetFlavor>0 && JetFlavor<5) || JetFlavor == 21 )
    {
        flavor = "cl";
        fstore->hist(type+"_pT_"+flavor)->Fill(p4MuJet.Pt(),ptrel);
        fstore->hist(type+"_eta_"+flavor)->Fill(TMath::Abs(p4MuJet.Eta()),ptrel);
        if (tagged)
        {
            fstore->hist(type+"tag_pT_"+flavor)->Fill(p4MuJet.Pt(),ptrel);
            fstore->hist(type+"tag_eta_"+flavor)->Fill(TMath::Abs(p4MuJet.Eta()),ptrel);
        }
        else
        {

            fstore->hist(type+"noTag_pT_"+flavor)->Fill(p4MuJet.Pt(),ptrel);
            fstore->hist(type+"noTag_eta_"+flavor)->Fill(TMath::Abs(p4MuJet.Eta()),ptrel);
        }

    }
    if ( JetFlavor == 4 )
    {
        flavor = "c";
        fstore->hist(type+"_pT_"+flavor)->Fill(p4MuJet.Pt(),ptrel);
        fstore->hist(type+"_eta_"+flavor)->Fill(TMath::Abs(p4MuJet.Eta()),ptrel);
        if (tagged)
        {
            fstore->hist(type+"tag_pT_"+flavor)->Fill(p4MuJet.Pt(),ptrel);
            fstore->hist(type+"tag_eta_"+flavor)->Fill(TMath::Abs(p4MuJet.Eta()),ptrel);
        }
        else
        {

            fstore->hist(type+"noTag_pT_"+flavor)->Fill(p4MuJet.Pt(),ptrel);
            fstore->hist(type+"noTag_eta_"+flavor)->Fill(TMath::Abs(p4MuJet.Eta()),ptrel);
        }

    }
    if ( (JetFlavor>0 && JetFlavor<4) || JetFlavor == 21 )
    {
        flavor = "l";
        fstore->hist(type+"_pT_"+flavor)->Fill(p4MuJet.Pt(),ptrel);
        fstore->hist(type+"_eta_"+flavor)->Fill(TMath::Abs(p4MuJet.Eta()),ptrel);
        if (tagged)
        {
            fstore->hist(type+"tag_pT_"+flavor)->Fill(p4MuJet.Pt(),ptrel);
            fstore->hist(type+"tag_eta_"+flavor)->Fill(TMath::Abs(p4MuJet.Eta()),ptrel);
        }
        else
        {

            fstore->hist(type+"noTag_pT_"+flavor)->Fill(p4MuJet.Pt(),ptrel);
            fstore->hist(type+"noTag_eta_"+flavor)->Fill(TMath::Abs(p4MuJet.Eta()),ptrel);
        }

    }

    // avoid jets with flavor = 0

    if (tagged)
    {

        fstore->hist(type+"tag_pT")->Fill(p4MuJet.Pt(),ptrel);
        fstore->hist(type+"tag_eta")->Fill(TMath::Abs(p4MuJet.Eta()),ptrel);
    }
    else
    {

        fstore->hist(type+"noTag_pT")->Fill(p4MuJet.Pt(),ptrel);
        fstore->hist(type+"noTag_eta")->Fill(TMath::Abs(p4MuJet.Eta()),ptrel);
    }


}

/*

  for (std::map<std::string,bool>::const_iterator imap = aMap.begin(); imap != aMap.end(); ++imap )
    {

      if ( imap->second )
        {

	  if ( type == "n")
            {

	      TaggedMujetHistos->Fill2d(type+"tag_pT_"+(imap->first),p4MuJet.Pt(),ptrel,weight);
	      TaggedMujetHistos->Fill2d(type+"tag_eta_"+(imap->first),TMath::Abs(p4MuJet.Eta()),ptrel,weight);

	      if ( JetFlavor == 5 )
                {
		  std::string flavor = "b_";
		  TaggedMujetHistos_mc->Fill2d(type+"tag_pT_"+flavor+(imap->first),p4MuJet.Pt(),ptrel,weight);
		  TaggedMujetHistos_mc->Fill2d(type+"tag_eta_"+flavor+(imap->first),TMath::Abs(p4MuJet.Eta()),ptrel,weight);
                }
	      if ( (JetFlavor>0 && JetFlavor<5) || JetFlavor == 21 )
                {
		  std::string flavor = "cl_";
		  TaggedMujetHistos_mc->Fill2d(type+"tag_pT_"+flavor+(imap->first),p4MuJet.Pt(),ptrel,weight);
		  TaggedMujetHistos_mc->Fill2d(type+"tag_eta_"+flavor+(imap->first),TMath::Abs(p4MuJet.Eta()),ptrel,weight);
                }
	      if ( JetFlavor == 4 )
                {
		  std::string flavor = "c_";
		  TaggedMujetHistos_mc->Fill2d(type+"tag_pT_"+flavor+(imap->first),p4MuJet.Pt(),ptrel,weight);
		  TaggedMujetHistos_mc->Fill2d(type+"tag_eta_"+flavor+(imap->first),TMath::Abs(p4MuJet.Eta()),ptrel,weight);
                }
	      if ( (JetFlavor>0 && JetFlavor<4) || JetFlavor == 21 )
                {
		  std::string flavor = "l_";
		  TaggedMujetHistos_mc->Fill2d(type+"tag_pT_"+flavor+(imap->first),p4MuJet.Pt(),ptrel,weight);
		  TaggedMujetHistos_mc->Fill2d(type+"tag_eta_"+flavor+(imap->first),TMath::Abs(p4MuJet.Eta()),ptrel,weight);
                }
            }
	  else if ( type == "p")
            {


	      TaggedAwayjetHistos->Fill2d(type+"tag_pT_"+(imap->first),p4MuJet.Pt(),ptrel,weight);
	      TaggedAwayjetHistos->Fill2d(type+"tag_eta_"+(imap->first),TMath::Abs(p4MuJet.Eta()),ptrel,weight);

	      if ( JetFlavor == 5 )
                {
		  std::string flavor = "b_";
		  TaggedAwayjetHistos_mc->Fill2d(type+"tag_pT_"+flavor+(imap->first),p4MuJet.Pt(),ptrel,weight);
		  TaggedAwayjetHistos_mc->Fill2d(type+"tag_eta_"+flavor+(imap->first),TMath::Abs(p4MuJet.Eta()),ptrel,weight);
                }
	      if ( (JetFlavor>0 && JetFlavor<5) || JetFlavor == 21 )
                {
		  std::string flavor = "cl_";
		  TaggedAwayjetHistos_mc->Fill2d(type+"tag_pT_"+flavor+(imap->first),p4MuJet.Pt(),ptrel,weight);
		  TaggedAwayjetHistos_mc->Fill2d(type+"tag_eta_"+flavor+(imap->first),TMath::Abs(p4MuJet.Eta()),ptrel,weight);
                }
	      if ( JetFlavor == 4 )
                {
		  std::string flavor = "c_";
		  TaggedAwayjetHistos_mc->Fill2d(type+"tag_pT_"+flavor+(imap->first),p4MuJet.Pt(),ptrel,weight);
		  TaggedAwayjetHistos_mc->Fill2d(type+"tag_eta_"+flavor+(imap->first),TMath::Abs(p4MuJet.Eta()),ptrel,weight);
                }
	      if ( (JetFlavor>0 && JetFlavor<4) || JetFlavor == 21 )
                {
		  std::string flavor = "l_";
		  TaggedAwayjetHistos_mc->Fill2d(type+"tag_pT_"+flavor+(imap->first),p4MuJet.Pt(),ptrel,weight);
		  TaggedAwayjetHistos_mc->Fill2d(type+"tag_eta_"+flavor+(imap->first),TMath::Abs(p4MuJet.Eta()),ptrel,weight);
                }
            }
        }
    }
*/


