import FWCore.ParameterSet.Config as cms

# Muon Digitization (CSC, DT, RPC electronics responce)
# CSC digitizer
#
from SimMuon.CSCDigitizer.muonCSCDigis_cfi import *
# DT digitizer
#
from SimMuon.DTDigitizer.muondtdigi_cfi import *
# RPC digitizer
# 
from SimMuon.RPCDigitizer.muonrpcdigi_cfi import *
muonDigi = cms.Sequence(muonCSCDigis+muonDTDigis+muonRPCDigis)

