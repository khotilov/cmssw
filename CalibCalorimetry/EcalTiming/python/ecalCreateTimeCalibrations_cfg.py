import FWCore.ParameterSet.Config as cms

process = cms.Process("EcalCreateTimeCalibrations")

# Global Tag -- for original timing calibrations
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.GlobalTag.globaltag = 'GR_R_43_V1::All'

process.GlobalTag.toGet = cms.VPSet(
  cms.PSet(record = cms.string("EcalTimeOffsetConstantRcd"),
           tag = cms.string("EcalTimeOffsetConstant_v01_offline"),
           connect = cms.untracked.string("frontier://FrontierProd/CMS_COND_43X_ECAL")
          )
)


process.MessageLogger = cms.Service("MessageLogger",
    cout = cms.untracked.PSet(threshold = cms.untracked.string('INFO')),
    categories = cms.untracked.vstring('*'),
    destinations = cms.untracked.vstring('cout')
)

process.TFileService = cms.Service("TFileService", 
    fileName = cms.string("ecalCreateTimeCalibs.root"),
    closeFileFast = cms.untracked.bool(True)
    )


process.createTimeCalibs = cms.EDAnalyzer("EcalCreateTimeCalibrations",
  OutputFileName = cms.string("file:converted1.root"),
  FileNameStart = cms.string("ecalCreateTimeCalibs"),
  SubtractDBcalibs = cms.bool(True),
  BxIncludeExclude = cms.string("-1"),
  OrbitIncludeExclude = cms.string("-1"),
  TriggerBitIncludeExclude = cms.string("-1"),
  TechTriggerBitIncludeExclude = cms.string("-1"),
  LumiIncludeExclude = cms.string("-1"),
  RunIncludeExclude = cms.string("-1"),
  AvgTimeMin = cms.double(-1),
  AvgTimeMax = cms.double(1),
  MinHitAmpEB = cms.double(15),
  MinHitAmpEE = cms.double(15),
  MaxSwissCross = cms.double(0.95),
  MaxHitTimeEB = cms.double(5),
  MinHitTimeEB = cms.double(-5),
  MaxHitTimeEE = cms.double(5),
  MinHitTimeEE = cms.double(-5),
  InputFileNames = cms.vstring(
                               #'file:/data/scooper/data/EcalTimeTree/2011A/EcalTimeTree_987654_100_1_9UZ.root',
                               #'file:/data/scooper/data/EcalTimeTree/2011A/EcalTimeTree_987654_100_2_f0s.root',
                               #'file:/data/scooper/data/EcalTimeTree/2011A/EcalTimeTree_987654_101_1_kyI.root',
                               #'file:/data/scooper/data/EcalTimeTree/2011A/EcalTimeTree_987654_101_2_HI2.root',
                               #'file:/data/scooper/data/EcalTimeTree/2011A/EcalTimeTree_987654_102_1_sOl.root',
                               #'file:/data/scooper/data/EcalTimeTree/2011A/EcalTimeTree_987654_102_2_dId.root',
                               #'file:/data/scooper/data/EcalTimeTree/2011A/EcalTimeTree_987654_103_2_0Mz.root',
                               #'file:/data/scooper/data/EcalTimeTree/2011A/EcalTimeTree_987654_104_1_yAK.root',
                               #'file:/data/scooper/data/EcalTimeTree/2011A/EcalTimeTree_987654_104_2_BTq.root',
                               #'file:/data/scooper/data/EcalTimeTree/2011A/EcalTimeTree_987654_105_2_oh1.root',
                               #'file:/data/scooper/data/EcalTimeTree/2011A/EcalTimeTree_987654_106_2_Ds5.root',
                               #'file:/data/scooper/data/EcalTimeTree/2011A/EcalTimeTree_987654_107_1_JfU.root',
                               #'file:/data/scooper/data/EcalTimeTree/2011A/EcalTimeTree_987654_107_2_c0Q.root',
                               #'file:/data/scooper/data/EcalTimeTree/2011A/EcalTimeTree_987654_109_1_74x.root',
                               #'file:/data/scooper/data/EcalTimeTree/2011A/EcalTimeTree_987654_109_2_Z8w.root',
                               #'file:/data/scooper/data/EcalTimeTree/2011A/EcalTimeTree_987654_110_2_JpT.root',
                               #'file:/data/scooper/data/EcalTimeTree/2011A/EcalTimeTree_987654_111_1_yeo.root',
                               #'file:/data/scooper/data/EcalTimeTree/2011A/EcalTimeTree_987654_11_1_mEY.root',
                               #'file:/data/scooper/data/EcalTimeTree/2011A/EcalTimeTree_987654_112_1_ocG.root',
                               #'file:/data/scooper/data/EcalTimeTree/2011A/EcalTimeTree_987654_11_2_ZNi.root',
                               #'file:/data/scooper/data/EcalTimeTree/2011A/EcalTimeTree_987654_113_2_U30.root',
                               #'file:/data/scooper/data/EcalTimeTree/2011A/EcalTimeTree_987654_114_1_gN0.root',
                               #'file:/data/scooper/data/EcalTimeTree/2011A/EcalTimeTree_987654_115_1_Hw2.root',
                               #'file:/data/scooper/data/EcalTimeTree/2011A/EcalTimeTree_987654_116_2_gf9.root',
                               #'file:/data/scooper/data/EcalTimeTree/2011A/EcalTimeTree_987654_117_2_34s.root',
                               #'file:/data/scooper/data/EcalTimeTree/2011A/EcalTimeTree_987654_118_2_QMw.root',
                               #'file:/data/scooper/data/EcalTimeTree/2011A/EcalTimeTree_987654_119_1_i7S.root',
                               #'file:/data/scooper/data/EcalTimeTree/2011A/EcalTimeTree_987654_119_2_SGD.root',
                               #'file:/data/scooper/data/EcalTimeTree/2011A/EcalTimeTree_987654_1_1_MfN.root',
                               #'file:/data/scooper/data/EcalTimeTree/2011A/EcalTimeTree_987654_120_1_lSb.root',
                               #'file:/data/scooper/data/EcalTimeTree/2011A/EcalTimeTree_987654_120_2_qqT.root',
                               #'file:/data/scooper/data/EcalTimeTree/2011A/EcalTimeTree_987654_121_2_LK0.root',
                               #'file:/data/scooper/data/EcalTimeTree/2011A/EcalTimeTree_987654_12_1_vTw.root',
                               #'file:/data/scooper/data/EcalTimeTree/2011A/EcalTimeTree_987654_122_2_2zJ.root',
                               #'file:/data/scooper/data/EcalTimeTree/2011A/EcalTimeTree_987654_123_1_1v9.root',
                               #'file:/data/scooper/data/EcalTimeTree/2011A/EcalTimeTree_987654_125_1_ayM.root',
                               #'file:/data/scooper/data/EcalTimeTree/2011A/EcalTimeTree_987654_126_1_L84.root',
                               #'file:/data/scooper/data/EcalTimeTree/2011A/EcalTimeTree_987654_126_2_Kc0.root',
                               #'file:/data/scooper/data/EcalTimeTree/2011A/EcalTimeTree_987654_127_1_ZCJ.root',
                               #'file:/data/scooper/data/EcalTimeTree/2011A/EcalTimeTree_987654_128_2_pay.root',
                               #'file:/data/scooper/data/EcalTimeTree/2011A/EcalTimeTree_987654_129_2_Wbi.root',
                               #'file:/data/scooper/data/EcalTimeTree/2011A/EcalTimeTree_987654_1_2_Hs0.root',
                               #'file:/data/scooper/data/EcalTimeTree/2011A/EcalTimeTree_987654_132_2_e9P.root',
                               #'file:/data/scooper/data/EcalTimeTree/2011A/EcalTimeTree_987654_134_2_GO2.root',
                               #'file:/data/scooper/data/EcalTimeTree/2011A/EcalTimeTree_987654_135_1_seo.root',
                               #'file:/data/scooper/data/EcalTimeTree/2011A/EcalTimeTree_987654_136_2_oSu.root',
                               #'file:/data/scooper/data/EcalTimeTree/2011A/EcalTimeTree_987654_137_2_7xj.root',
                               #'file:/data/scooper/data/EcalTimeTree/2011A/EcalTimeTree_987654_139_2_gj6.root',
                               #'file:/data/scooper/data/EcalTimeTree/2011A/EcalTimeTree_987654_141_1_29h.root',
                               #'file:/data/scooper/data/EcalTimeTree/2011A/EcalTimeTree_987654_142_2_dt2.root',
                               #'file:/data/scooper/data/EcalTimeTree/2011A/EcalTimeTree_987654_14_2_eXc.root',
                               #'file:/data/scooper/data/EcalTimeTree/2011A/EcalTimeTree_987654_143_2_Jd4.root',
                               #'file:/data/scooper/data/EcalTimeTree/2011A/EcalTimeTree_987654_145_1_wOf.root',
                               #'file:/data/scooper/data/EcalTimeTree/2011A/EcalTimeTree_987654_149_1_S5F.root',
                               #'file:/data/scooper/data/EcalTimeTree/2011A/EcalTimeTree_987654_150_1_TvM.root',
                               #'file:/data/scooper/data/EcalTimeTree/2011A/EcalTimeTree_987654_150_2_QQC.root',
                               #'file:/data/scooper/data/EcalTimeTree/2011A/EcalTimeTree_987654_151_2_y0n.root',
                               #'file:/data/scooper/data/EcalTimeTree/2011A/EcalTimeTree_987654_15_1_fMP.root',
                               #'file:/data/scooper/data/EcalTimeTree/2011A/EcalTimeTree_987654_152_1_cO9.root',
                               #'file:/data/scooper/data/EcalTimeTree/2011A/EcalTimeTree_987654_152_2_b25.root',
                               #'file:/data/scooper/data/EcalTimeTree/2011A/EcalTimeTree_987654_154_2_QMV.root',
                               #'file:/data/scooper/data/EcalTimeTree/2011A/EcalTimeTree_987654_155_1_gFA.root',
                               #'file:/data/scooper/data/EcalTimeTree/2011A/EcalTimeTree_987654_155_2_YWs.root',
                               #'file:/data/scooper/data/EcalTimeTree/2011A/EcalTimeTree_987654_156_1_CIN.root',
                               #'file:/data/scooper/data/EcalTimeTree/2011A/EcalTimeTree_987654_156_2_Bmq.root',
                               #'file:/data/scooper/data/EcalTimeTree/2011A/EcalTimeTree_987654_157_2_YeM.root',
                               #'file:/data/scooper/data/EcalTimeTree/2011A/EcalTimeTree_987654_16_1_IrK.root',
                               #'file:/data/scooper/data/EcalTimeTree/2011A/EcalTimeTree_987654_166_1_iIA.root',
                               #'file:/data/scooper/data/EcalTimeTree/2011A/EcalTimeTree_987654_17_1_13R.root',
                               #'file:/data/scooper/data/EcalTimeTree/2011A/EcalTimeTree_987654_18_2_zHv.root',
                               #'file:/data/scooper/data/EcalTimeTree/2011A/EcalTimeTree_987654_22_1_J9i.root',
                               #'file:/data/scooper/data/EcalTimeTree/2011A/EcalTimeTree_987654_2_2_2LB.root',
                               #'file:/data/scooper/data/EcalTimeTree/2011A/EcalTimeTree_987654_22_2_UvO.root',
                               #'file:/data/scooper/data/EcalTimeTree/2011A/EcalTimeTree_987654_23_2_0gu.root',
                               #'file:/data/scooper/data/EcalTimeTree/2011A/EcalTimeTree_987654_24_1_HfO.root',
                               #'file:/data/scooper/data/EcalTimeTree/2011A/EcalTimeTree_987654_24_2_pd2.root',
                               #'file:/data/scooper/data/EcalTimeTree/2011A/EcalTimeTree_987654_25_2_eZ8.root',
                               #'file:/data/scooper/data/EcalTimeTree/2011A/EcalTimeTree_987654_26_2_NAq.root',
                               #'file:/data/scooper/data/EcalTimeTree/2011A/EcalTimeTree_987654_28_2_ZxG.root',
                               #'file:/data/scooper/data/EcalTimeTree/2011A/EcalTimeTree_987654_29_2_JRG.root',
                               #'file:/data/scooper/data/EcalTimeTree/2011A/EcalTimeTree_987654_30_2_cpU.root',
                               #'file:/data/scooper/data/EcalTimeTree/2011A/EcalTimeTree_987654_31_1_Xu2.root',
                               #'file:/data/scooper/data/EcalTimeTree/2011A/EcalTimeTree_987654_32_2_9GI.root',
                               #'file:/data/scooper/data/EcalTimeTree/2011A/EcalTimeTree_987654_3_2_xTn.root',
                               #'file:/data/scooper/data/EcalTimeTree/2011A/EcalTimeTree_987654_33_2_Vfe.root',
                               #'file:/data/scooper/data/EcalTimeTree/2011A/EcalTimeTree_987654_35_1_TiC.root',
                               #'file:/data/scooper/data/EcalTimeTree/2011A/EcalTimeTree_987654_37_2_Ym2.root',
                               #'file:/data/scooper/data/EcalTimeTree/2011A/EcalTimeTree_987654_38_2_vxX.root',
                               #'file:/data/scooper/data/EcalTimeTree/2011A/EcalTimeTree_987654_39_1_aPZ.root',
                               #'file:/data/scooper/data/EcalTimeTree/2011A/EcalTimeTree_987654_39_2_nvy.root',
                               #'file:/data/scooper/data/EcalTimeTree/2011A/EcalTimeTree_987654_40_1_nlB.root',
                               #'file:/data/scooper/data/EcalTimeTree/2011A/EcalTimeTree_987654_41_1_AVv.root',
                               #'file:/data/scooper/data/EcalTimeTree/2011A/EcalTimeTree_987654_4_1_va8.root',
                               #'file:/data/scooper/data/EcalTimeTree/2011A/EcalTimeTree_987654_42_1_de4.root',
                               #'file:/data/scooper/data/EcalTimeTree/2011A/EcalTimeTree_987654_42_2_I24.root',
                               #'file:/data/scooper/data/EcalTimeTree/2011A/EcalTimeTree_987654_46_1_c5h.root',
                               #'file:/data/scooper/data/EcalTimeTree/2011A/EcalTimeTree_987654_46_2_0vv.root',
                               #'file:/data/scooper/data/EcalTimeTree/2011A/EcalTimeTree_987654_47_2_Kf0.root',
                               #'file:/data/scooper/data/EcalTimeTree/2011A/EcalTimeTree_987654_49_1_PuV.root',
                               #'file:/data/scooper/data/EcalTimeTree/2011A/EcalTimeTree_987654_49_2_imU.root',
                               #'file:/data/scooper/data/EcalTimeTree/2011A/EcalTimeTree_987654_54_1_JbF.root',
                               #'file:/data/scooper/data/EcalTimeTree/2011A/EcalTimeTree_987654_55_1_dvK.root',
                               #'file:/data/scooper/data/EcalTimeTree/2011A/EcalTimeTree_987654_56_1_UJK.root',
                               #'file:/data/scooper/data/EcalTimeTree/2011A/EcalTimeTree_987654_56_2_1y9.root',
                               #'file:/data/scooper/data/EcalTimeTree/2011A/EcalTimeTree_987654_57_1_q4n.root',
                               #'file:/data/scooper/data/EcalTimeTree/2011A/EcalTimeTree_987654_57_2_WNn.root',
                               #'file:/data/scooper/data/EcalTimeTree/2011A/EcalTimeTree_987654_58_2_HFK.root',
                               #'file:/data/scooper/data/EcalTimeTree/2011A/EcalTimeTree_987654_59_1_gqo.root',
                               #'file:/data/scooper/data/EcalTimeTree/2011A/EcalTimeTree_987654_60_1_zj0.root',
                               #'file:/data/scooper/data/EcalTimeTree/2011A/EcalTimeTree_987654_6_1_XKq.root',
                               #'file:/data/scooper/data/EcalTimeTree/2011A/EcalTimeTree_987654_63_1_BA1.root',
                               #'file:/data/scooper/data/EcalTimeTree/2011A/EcalTimeTree_987654_65_1_nbf.root',
                               #'file:/data/scooper/data/EcalTimeTree/2011A/EcalTimeTree_987654_65_2_Lp2.root',
                               #'file:/data/scooper/data/EcalTimeTree/2011A/EcalTimeTree_987654_70_1_XYV.root',
                               #'file:/data/scooper/data/EcalTimeTree/2011A/EcalTimeTree_987654_7_1_bEd.root',
                               #'file:/data/scooper/data/EcalTimeTree/2011A/EcalTimeTree_987654_7_2_zh5.root',
                               #'file:/data/scooper/data/EcalTimeTree/2011A/EcalTimeTree_987654_74_1_sjx.root',
                               #'file:/data/scooper/data/EcalTimeTree/2011A/EcalTimeTree_987654_74_2_xMC.root',
                               #'file:/data/scooper/data/EcalTimeTree/2011A/EcalTimeTree_987654_79_2_UlL.root',
                               'file:/data/scooper/data/EcalTimeTree/2011A/EcalTimeTree_987654_80_2_ijB.root',
                               'file:/data/scooper/data/EcalTimeTree/2011A/EcalTimeTree_987654_81_1_mnP.root',
                               'file:/data/scooper/data/EcalTimeTree/2011A/EcalTimeTree_987654_82_2_cAN.root',
                               'file:/data/scooper/data/EcalTimeTree/2011A/EcalTimeTree_987654_83_1_j3P.root',
                               'file:/data/scooper/data/EcalTimeTree/2011A/EcalTimeTree_987654_84_2_x10.root',
                               'file:/data/scooper/data/EcalTimeTree/2011A/EcalTimeTree_987654_86_1_zM0.root',
                               'file:/data/scooper/data/EcalTimeTree/2011A/EcalTimeTree_987654_86_2_xhe.root',
                               'file:/data/scooper/data/EcalTimeTree/2011A/EcalTimeTree_987654_87_1_3Ph.root',
                               'file:/data/scooper/data/EcalTimeTree/2011A/EcalTimeTree_987654_88_1_8A9.root',
                               'file:/data/scooper/data/EcalTimeTree/2011A/EcalTimeTree_987654_89_2_ny8.root',
                               'file:/data/scooper/data/EcalTimeTree/2011A/EcalTimeTree_987654_90_1_Wsc.root',
                               'file:/data/scooper/data/EcalTimeTree/2011A/EcalTimeTree_987654_91_2_yZA.root',
                               'file:/data/scooper/data/EcalTimeTree/2011A/EcalTimeTree_987654_93_2_k9W.root',
                               'file:/data/scooper/data/EcalTimeTree/2011A/EcalTimeTree_987654_94_1_Bsq.root',
                               'file:/data/scooper/data/EcalTimeTree/2011A/EcalTimeTree_987654_95_1_1VL.root',
                               'file:/data/scooper/data/EcalTimeTree/2011A/EcalTimeTree_987654_96_2_9Cv.root',
                               'file:/data/scooper/data/EcalTimeTree/2011A/EcalTimeTree_987654_97_2_MOR.root'

    )
)

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(1) )
process.source = cms.Source("EmptySource",
       numberEventsInRun = cms.untracked.uint32(1),
       firstRun = cms.untracked.uint32(888888)
)


process.p = cms.Path(process.createTimeCalibs)
