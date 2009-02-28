import FWCore.ParameterSet.Config as cms

pfRecoTauProducer = cms.EDProducer("PFRecoTauProducer",
      PFCand_minPt         = cms.double(0.5),  #general pt cut on objects

      LeadPFCand_mintPt    = cms.double(5.0),  #cut on lead object (can be track, or gamma)

      #Pt cuts applied to PF constituents
      ChargedHadrCand_minPt = cms.double(1.0),
      GammaCand_minPt       = cms.double(1.5), 
      NeutrHadrCand_minPt   = cms.double(1.0),

      #SignalCone parameters
      TrackerSignalConeMetric = cms.string('DR'), ## * 
      TrackerSignalConeSizeFormula = cms.string('0.07'), ## **   
      TrackerSignalConeSize_min = cms.double(0.0),
      TrackerSignalConeSize_max = cms.double(0.6),

      ECALSignalConeMetric = cms.string('DR'), ## * 
      ECALSignalConeSizeFormula = cms.string('0.15'), ## **
      ECALSignalConeSize_min = cms.double(0.0),
      ECALSignalConeSize_max = cms.double(0.6),

      HCALSignalConeMetric = cms.string('DR'), ## *  
      HCALSignalConeSizeFormula = cms.string('0.10'), ## **
      HCALSignalConeSize_min = cms.double(0.0),
      HCALSignalConeSize_max = cms.double(0.6),

      #IsolationCone parameters
      TrackerIsolConeMetric = cms.string('DR'), ## * 
      TrackerIsolConeSizeFormula = cms.string('0.50'), ## ** 
      TrackerIsolConeSize_min = cms.double(0.0),
      TrackerIsolConeSize_max = cms.double(0.6),

      ECALIsolConeMetric = cms.string('DR'), ## *  
      ECALIsolConeSizeFormula = cms.string('0.50'), ## **         
      ECALIsolConeSize_min = cms.double(0.0),
      ECALIsolConeSize_max = cms.double(0.6),

      HCALIsolConeMetric = cms.string('DR'), ## * 
      HCALIsolConeSizeFormula = cms.string('0.50'), ## **
      HCALIsolConeSize_min = cms.double(0.0),
      HCALIsolConeSize_max = cms.double(0.6),

      # Cut on the number of tracker hits on isolation PF charged hadrons
      ChargedHadrCand_IsolAnnulus_minNhits = cms.uint32(8),

      #Electron rejection parameters
      ElectronPreIDProducer                = cms.InputTag("elecpreid"),
      EcalStripSumE_deltaPhiOverQ_minValue = cms.double(-0.1),
      EcalStripSumE_deltaPhiOverQ_maxValue = cms.double(0.5),
      EcalStripSumE_minClusEnergy          = cms.double(0.1),
      EcalStripSumE_deltaEta               = cms.double(0.03),
      ElecPreIDLeadTkMatch_maxDR           = cms.double(0.01),

      #Lead track matching cone parameters
      MatchingConeMetric      = cms.string('DR'),
      MatchingConeSizeFormula = cms.string('0.1'),
      MatchingConeSize_min    = cms.double(0.0), #min/max don't affect, as it is locked at 0.1
      MatchingConeSize_max    = cms.double(0.6), #just make sure they window 0.1!

      # PFTaus are seeded by TauTagInfos (basically a jet wrapper/quality filter)
      PFTauTagInfoProducer = cms.InputTag("pfRecoTauTagInfoProducer"),
      JetPtMin             = cms.double(0.0),
      #Filter lead charged hadron cand. by DZ to vertex?
      UseChargedHadrCandLeadChargedHadrCand_tksDZconstraint = cms.bool(True),
      ChargedHadrCandLeadChargedHadrCand_tksmaxDZ = cms.double(0.2),

      #Standard PV stuff
      PVProducer      = cms.InputTag('offlinePrimaryVertices'),
      smearedPVsigmaY = cms.double(0.0015),
      smearedPVsigmaX = cms.double(0.0015),
      smearedPVsigmaZ = cms.double(0.005),

      #FixedArea parameters
      AreaMetric_recoElements_maxabsEta = cms.double(2.5),

      DataType = cms.string("AOD"),            # Computring tier (modifies how some values are retrieved..)

      # The following parameters affect TRACKS ONLY.
      # i.e., they use tracks, from the jet tracks associator, not PFChargedHadronCandidates
      # Intended only as a cross check, don't use them unles syou need them!

      # Track - PV filtering
      LeadTrack_minPt               = cms.double(5.0),
      TrackLeadTrack_maxDZ          = cms.double(0.2),
      UseTrackLeadTrackDZconstraint = cms.bool(True),
      Track_minPt                   = cms.double(1.0),
      Track_IsolAnnulus_minNhits    = cms.uint32(3),

)
 # * possible metrics : "DR", "angle", "area";
    #   if the "area" metric is chosen, AreaMetric_recoElements_maxabsEta parameter is considered, the area of a cone is increased by increasing the angle of the cone;  
    #   functionnality to use a "DR" signal cone and an "area" isolation outer cone is not available;
    # ** may depend on E(energy) and/or PT(transverse momentum) of the initial PFJet, ex. : "3.0/E" or "3.0/ET" 
    #    if XXXConeSizeFormula>XXXConeSize_max then XXXConeSize_max is the considered cone size ; if XXXConeSizeFormula<XXXConeSize_min then XXXConeSize_min is the considered cone size;  
    # *** a PV is needed for computing a leading (charged hadron PFCand) rec. tk signed transverse impact parameter. 
    # For electron rejection variable

