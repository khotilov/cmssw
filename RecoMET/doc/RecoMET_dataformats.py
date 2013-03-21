c = 7

title0="PFMET (in RECO and AOD)"

dict0 = {
'0':['pfMet','reco::PFMETCollection','MET in reconstructed particles in the particle flow algorithm']
}

title1="GenMET (in RECO and AOD)"

dict1 = {
'0':['genMetTrue','reco::GenMETCollection','MET in generated particles in simulation in their final states but excluding neutrinos, excited neutrinos, right-handed neutrinos, sneutrinos, neutralinos, gravitons, gravitinos'],
'1':['genMetCalo','reco::GenMETCollection','MET in generated particles in simulation in their final states but excluding neutrinos, excited neutrinos, right-handed neutrinos, sneutrinos, neutralinos, gravitons, gravitinos and also muons'],
'2':['genMetCaloAndNonPrompt','reco::GenMETCollection','MET in generated particles in simulation in their final states but excluding excited neutrinos, right-handed neutrinos, sneutrinos, neutralinos, gravitons, gravitinos and, additionally, excluding muons and neutrinos coming from the decay of gauge bosons and top quarks']
}

title2="CaloMET (in RECO and AOD)"

dict2 = {
'0':['met','reco::CaloMETCollection','MET in energy deposits in calorimeter towers in EB, EE, HB, HE, and HF'],
'1':['metNoHF','reco::CaloMETCollection','MET in energy deposits in calorimeter towers in EB, EE, HB, and HE'],
'2':['metHO','reco::CaloMETCollection','MET in energy deposits in calorimeter towers in EB, EE, HB, HE, HF, and HO'],
'3':['corMetGlobalMuons','reco::CaloMETCollection','MET in energy deposits in calorimeter towers in EB, EE, HB, HE, and HF with corrections for muons']
}

title3="TCMET (in RECO and AOD)"

dict3 = {
'0':['tcMet','reco::METCollection','Track Corrected MET using `met`, `muons`, `gsfElectrons` and `generalTracks`'],
'1':['htMetAK5','reco::METCollection','Raw Missing Transverse Energy calculated using anti-KT5 CaloJets'],
'2':['htMetAK7','reco::METCollection','Raw Missing Transverse Energy calculated using anti-KT7 CaloJets'],
'3':['htMetKT4','reco::METCollection','Raw Missing Transverse Energy calculated using FastKt4 CaloJets'],
'4':['htMetKT6','reco::METCollection','Raw Missing Transverse Energy calculated using FastKt6 CaloJets'],
'5':['htMetIC5','reco::METCollection','Raw Missing Transverse Energy calculated using IC5 CaloJets'],
'6':['tcMetWithPFclusters','reco::METCollection','Track Corrected MET using `particleFlowClusters`, `muons`, `gsfElectrons` and `generalTracks`']
}

title4="MET associations (in RECO and AOD)"

dict4 = {
'0':['muonMETValueMapProducer, muCorrData','reco::MuonMETCorrectionData','information on how muons were used to correct MET and what associated MIP deposits are used'],
'1':['muonTCMETValueMapProducer, muCorrData','reco::MuonMETCorrectionData','information on how muons were used to correct tcMET and what associated calo deposits are if the muon is treated as a pion']
}

title5="HCAL Noise (RECO and AOD)"

dict5 = {
'0':['hcalnoise','reco::HcalNoiseRBX',' '],
'1':['hcalnoise','HcalNoiseSummary',' '],
'2':['Beam Halo (RECO and AOD)',' ',' '],
'3':['BeamHaloSummary','reco::BeamHaloSummary',' ']
}

title6="CaloMET (in RECO only)"

dict6 = {
'0':['metNoHFHO','reco::CaloMETCollection','MET in energy deposits in calorimeter towers in EB, EE, HB, HE, and HO'],
'1':['metOpt','reco::CaloMETCollection','MET in energy deposits in calorimeter towers in EB, EE, HB, HE, and HF with optimized threshold parameters'],
'2':['metOptNoHF','reco::CaloMETCollection','MET in energy deposits in calorimeter towers in EB, EE, HB, and HE with optimized threshold parameters'],
'3':['metOptHO','reco::CaloMETCollection','MET in energy deposits in calorimeter towers in EB, EE, HB, HE, HF and HO with optimized threshold parameters'],
'4':['metOptNoHFHO','reco::CaloMETCollection','MET in energy deposits in calorimeter towers in EB, EE, HB, HE, and HO with optimized threshold parameters']
}



full_title = "RecoMET collections (in RECO and AOD)"

full = {
    '0':['met', 'recoCaloMETs', 'MET in energy deposits in calorimeter towers in EB, EE, HB, HE, and HF'] ,
    '1':['metNoHF', 'recoCaloMETs', 'MET in energy deposits in calorimeter towers in EB, EE, HB, and HE'] ,
    '2':['metHO', 'recoCaloMETs', 'MET in energy deposits in calorimeter towers in EB, EE, HB, HE, HF, and HO'] ,
    '3':['corMetGlobalMuons', 'recoCaloMETs', 'MET in energy deposits in calorimeter towers in EB, EE, HB, HE, and HF with corrections for muons'] ,
    '4':['tcMet', 'recoMETs', 'Track Corrected MET using `met`, `muons`, `gsfElectrons` and `generalTracks`'] ,
    '5':['tcMetWithPFclusters', 'recoMETs', 'Track Corrected MET using `particleFlowClusters`, `muons`, `gsfElectrons` and `generalTracks`'] ,
    '6':['pfMet', 'recoPFMETs', 'MET in reconstructed particles in the particle flow algorithm'] ,
    '7':['muonMETValueMapProducer', 'recoMuonMETCorrectionDataedmValueMap', 'Information on how muons were used to correct MET and what associated MIP deposits are used'] ,
    '8':['muonTCMETValueMapProducer', 'recoMuonMETCorrectionDataedmValueMap', 'Information on how muons were used to correct tcMET and what associated calo deposits are if the muon is treated as a pion'] ,
    '9':['hcalnoise', 'recoHcalNoiseRBXs', 'No documentation'] ,
    '10':['hcalnoise', 'HcalNoiseSummary', 'No documentation'] ,
    '11':['*', '*HaloData', 'No documentation'] ,
    '12':['BeamHaloSummary', '*BeamHaloSummary', 'No documentation'],
    '13':['genMetTrue', 'reco::GenMETCollection', 'MET in generated particles in simulation in their final states but excluding neutrinos, excited neutrinos, right-handed neutrinos, sneutrinos, neutralinos, gravitons, gravitinos'],
    '14':['genMetCalo', 'reco::GenMETCollection', 'MET in generated particles in simulation in their final states but excluding neutrinos, excited neutrinos, right-handed neutrinos, sneutrinos, neutralinos, gravitons, gravitinos and also muons'],
    '15':['genMetCaloAndNonPrompt', 'reco::GenMETCollection', 'MET in generated particles in simulation in their final states but excluding excited neutrinos, right-handed neutrinos, sneutrinos, neutralinos, gravitons, gravitinos and, additionally, excluding muons and neutrinos coming from the decay of gauge bosons and top quarks'],
    # Correction needed, because not matched with Event Content
    '16':['htMetAK5','reco::METCollection','Raw Missing Transverse Energy calculated using anti-KT5 CaloJets'],
    '17':['htMetAK7','reco::METCollection','Raw Missing Transverse Energy calculated using anti-KT7 CaloJets'],
    '18':['htMetKT4','reco::METCollection','Raw Missing Transverse Energy calculated using FastKt4 CaloJets'],
    '19':['htMetKT6','reco::METCollection','Raw Missing Transverse Energy calculated using FastKt6 CaloJets'],
    '20':['htMetIC5','reco::METCollection','Raw Missing Transverse Energy calculated using IC5 CaloJets'],
    '21':['metNoHFHO','reco::CaloMETCollection','MET in energy deposits in calorimeter towers in EB, EE, HB, HE, and HO'],
    '22':['metOpt','reco::CaloMETCollection','MET in energy deposits in calorimeter towers in EB, EE, HB, HE, and HF with optimized threshold parameters'],
    '23':['metOptNoHF','reco::CaloMETCollection','MET in energy deposits in calorimeter towers in EB, EE, HB, and HE with optimized threshold parameters'],
    '24':['metOptHO','reco::CaloMETCollection','MET in energy deposits in calorimeter towers in EB, EE, HB, HE, HF and HO with optimized threshold parameters'],
    '25':['metOptNoHFHO','reco::CaloMETCollection','MET in energy deposits in calorimeter towers in EB, EE, HB, HE, and HO with optimized threshold parameters'] 
}

reco_title = "RecoMET collections (in RECO only)"

reco = {
    '0':['met', 'recoCaloMETs', 'MET in energy deposits in calorimeter towers in EB, EE, HB, HE, and HF'] ,
    '1':['metNoHF', 'recoCaloMETs', 'MET in energy deposits in calorimeter towers in EB, EE, HB, and HE'] ,
    '2':['metHO', 'recoCaloMETs', 'MET in energy deposits in calorimeter towers in EB, EE, HB, HE, HF, and HO'] ,
    '3':['corMetGlobalMuons', 'recoCaloMETs', 'MET in energy deposits in calorimeter towers in EB, EE, HB, HE, and HF with corrections for muons'] ,
    '4':['tcMet', 'recoMETs', 'Track Corrected MET using `met`, `muons`, `gsfElectrons` and `generalTracks`'] ,
    '5':['tcMetWithPFclusters', 'recoMETs', 'Track Corrected MET using `particleFlowClusters`, `muons`, `gsfElectrons` and `generalTracks`'] ,
    '6':['pfMet', 'recoPFMETs', 'MET in reconstructed particles in the particle flow algorithm'] ,
    '7':['muonMETValueMapProducer', 'recoMuonMETCorrectionDataedmValueMap', 'Information on how muons were used to correct MET and what associated MIP deposits are used'] ,
    '8':['muonTCMETValueMapProducer', 'recoMuonMETCorrectionDataedmValueMap', 'Information on how muons were used to correct tcMET and what associated calo deposits are if the muon is treated as a pion'] ,
    '9':['hcalnoise', 'recoHcalNoiseRBXs', 'No documentation'] ,
    '10':['hcalnoise', 'HcalNoiseSummary', 'No documentation'] ,
    '11':['*', '*HaloData', 'No documentation'] ,
    '12':['BeamHaloSummary', '*BeamHaloSummary', 'No documentation'],
    '13':['genMetTrue', 'reco::GenMETCollection', 'MET in generated particles in simulation in their final states but excluding neutrinos, excited neutrinos, right-handed neutrinos, sneutrinos, neutralinos, gravitons, gravitinos'],
    '14':['genMetCalo', 'reco::GenMETCollection', 'MET in generated particles in simulation in their final states but excluding neutrinos, excited neutrinos, right-handed neutrinos, sneutrinos, neutralinos, gravitons, gravitinos and also muons'],
    '15':['genMetCaloAndNonPrompt', 'reco::GenMETCollection', 'MET in generated particles in simulation in their final states but excluding excited neutrinos, right-handed neutrinos, sneutrinos, neutralinos, gravitons, gravitinos and, additionally, excluding muons and neutrinos coming from the decay of gauge bosons and top quarks']
}

aod_title = "RecoMET collections (in AOD only)"

aod = {
    '0':['met', 'recoCaloMETs', 'MET in energy deposits in calorimeter towers in EB, EE, HB, HE, and HF'] ,
    '1':['metNoHF', 'recoCaloMETs', 'MET in energy deposits in calorimeter towers in EB, EE, HB, and HE'] ,
    '2':['metHO', 'recoCaloMETs', 'MET in energy deposits in calorimeter towers in EB, EE, HB, HE, HF, and HO'] ,
    '3':['corMetGlobalMuons', 'recoCaloMETs', 'MET in energy deposits in calorimeter towers in EB, EE, HB, HE, and HF with corrections for muons'] ,
    '4':['tcMet', 'recoMETs', 'Track Corrected MET using `met`, `muons`, `gsfElectrons` and `generalTracks`'] ,
    '5':['tcMetWithPFclusters', 'recoMETs', 'Track Corrected MET using `particleFlowClusters`, `muons`, `gsfElectrons` and `generalTracks`'] ,
    '6':['pfMet', 'recoPFMETs', 'MET in reconstructed particles in the particle flow algorithm'] ,
    '7':['muonMETValueMapProducer', 'recoMuonMETCorrectionDataedmValueMap', 'Information on how muons were used to correct MET and what associated MIP deposits are used'] ,
    '8':['muonTCMETValueMapProducer', 'recoMuonMETCorrectionDataedmValueMap', 'Information on how muons were used to correct tcMET and what associated calo deposits are if the muon is treated as a pion'] ,
    '9':['hcalnoise', 'HcalNoiseSummary', 'No documentation'] ,
    '10':['*', '*GlobalHaloData', 'No documentation'] ,
    '11':['BeamHaloSummary', '*BeamHaloSummary', 'No documentation'],
    '12':['genMetTrue', 'reco::GenMETCollection', 'MET in generated particles in simulation in their final states but excluding neutrinos, excited neutrinos, right-handed neutrinos, sneutrinos, neutralinos, gravitons, gravitinos'],
    '13':['genMetCalo', 'reco::GenMETCollection', 'MET in generated particles in simulation in their final states but excluding neutrinos, excited neutrinos, right-handed neutrinos, sneutrinos, neutralinos, gravitons, gravitinos and also muons'],
    '14':['genMetCaloAndNonPrompt', 'reco::GenMETCollection', 'MET in generated particles in simulation in their final states but excluding excited neutrinos, right-handed neutrinos, sneutrinos, neutralinos, gravitons, gravitinos and, additionally, excluding muons and neutrinos coming from the decay of gauge bosons and top quarks'] 
}

