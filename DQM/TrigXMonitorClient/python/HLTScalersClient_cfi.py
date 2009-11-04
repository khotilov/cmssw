import FWCore.ParameterSet.Config as cms

# HLT scalers client. wittich 8/08
# $Id: HLTScalersClient_cfi.py,v 1.1 2008/08/22 20:56:56 wittich Exp $
# $Log: HLTScalersClient_cfi.py,v $
# Revision 1.1  2008/08/22 20:56:56  wittich
# - add client for HLT Scalers
# - Move rate calculation to HLTScalersClient and slim down the
#   filter-farm part of HLTScalers
#
#
hltsClient = cms.EDFilter("HLTScalersClient",
  # no configuration yet
  dqmFolder = cms.untracked.string("HLT/HLTScalers_EvF")
)

