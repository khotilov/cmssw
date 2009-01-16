import FWCore.ParameterSet.Config as cms
import copy

from ElectroWeakAnalysis.EWKTau.tauForETauMod.tauForETauPFRecoTauProducer_cff import *
from ElectroWeakAnalysis.EWKTau.tauForETauMod.tauForETauPatConfig_cfi import *

PFTauForETauPat = cms.Sequence(
                                 pfRecoTauIsoDiscrForETau
                                *allLayer0PfTausForETau
                                *pfRecoTauLdgTrkFindForETau
                                *pfRecoTauLdgTrkPtCutForETau
                                *pfRecoTauTrkIsoDiscrForETau
                                *pfRecoTauEclIsoDiscrForETau
                                *pfRecoTauElecRejDiscrForETau
                                *pfRecoTauMuonRejDiscrForETau
                                *layer1PfTausForETau
                                *selectorPfTausForETauKin
                                *selectorPfTausForETauElecSep
                                *selectorPfTausForETauLdgTrk
                                *selectorPfTausForETauLdgTrkPt
                                *selectorPfTausForETauEclIsod
                                *selectorPfTausForETauTrkIsod
                                *selectorPfTausForETauElecRej
                                *selectorPfTausForETauProng
                                *selectorPfTausForOSSSChargeSep
                                *selectorPfTausForOSSSChargeSepAnti
                              )

PFTauForETauEffPat = cms.Sequence(
                                   pfRecoTauIsoDiscrForETauEff
                                  *allLayer0PfTausForETauEff
                                  *pfRecoTauLdgTrkFindForETauEff
                                  *pfRecoTauLdgTrkPtCutForETauEff
                                  *pfRecoTauTrkIsoDiscrForETauEff
                                  *pfRecoTauEclIsoDiscrForETauEff
                                  *pfRecoTauElecRejDiscrForETauEff
                                  *pfRecoTauMuonRejDiscrForETauEff
                                  *pfRecoTauIsoDiscrForETauEff
                                  *layer1PfTausForETauEff
                                  *selectorPfTausForETauKinEff
                                  *selectorPfTausForETauElecSepEff
                                  *selectorPfTausForETauLdgTrkEff
                                  *selectorPfTausForETauLdgTrkPtEff
                                  *selectorPfTausForETauEclIsodEff
                                  *selectorPfTausForETauTrkIsodEff
                                  *selectorPfTausForETauElecRejEff
                                  *selectorPfTausForETauProngEff
                                  *selectorPfTausForOSSSChargeSepEff
                                  *selectorPfTausForOSSSChargeSepAntiEff
                                 )






