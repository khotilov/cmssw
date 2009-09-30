#! /usr/bin/env python

import sys
import subprocess
import FWCore.ParameterSet.Config as cms

config = sys.argv[1]

def extractBlock(config, blocks, target):
  #print 'configuration: %s' % config
  #print 'blocks:        %s' % ', '.join(blocks)
  #print 'target:        %s' % target
  #print
  commands = ','.join( block + '::outputCommands' for block in blocks )
  proc = subprocess.Popen(
    "edmConfigFromDB --configName %s --noedsources --nopaths --noes --nopsets --noservices --cff --blocks %s --format python | sed -e'/^streams/,/^)/d' -e'/^datasets/,/^)/d' > %s" % (config, commands, target),
    shell  = True,
    stdin  = None, 
    stdout = None,
    stderr = None,
  )
  proc.wait()

def extractBlocks(config):
  outputA    = ( 'hltOutputA', )
  outputALCA = ( 'hltOutputALCAPHISYM', 'hltOutputALCAPHISYMHCAL', 'hltOutputALCAP0', 'hltOutputRPCMON' )
  outputMON  = ( 'hltOutputDQM', 'hltOutputHLTDQM', 'hltOutputHLTMON', 'hltOutput8E29', 'hltOutput1E31', 'hltOutputHIon' )
  extractBlock(config, outputA,    'hltOutputA_cff.py')
  extractBlock(config, outputALCA, 'hltOutputALCA_cff.py')
  extractBlock(config, outputMON,  'hltOutputMON_cff.py')

def makePSet(statements):
  statements = list(statements)
  statements.sort()
  block = cms.PSet(
    outputCommands = cms.untracked.vstring( 'drop *_hlt*_*_*', )
  )
  block.outputCommands.extend( statements )
  return block


def buildPSet(blocks):
  statements = set()
  for block in blocks:
    statements.update( statement for statement in block if statement.find('drop') != 0 )
  return makePSet(statements)


def buildPSetWithoutRAWs(blocks):
  statements = set()
  for block in blocks:
    statements.update( statement for statement in block if statement.find('drop') != 0 and statement.find('keep FEDRawDataCollection') != 0)
  return makePSet(statements)


# extract the HLT layer event content
extractBlocks( config )
import hltOutputA_cff
import hltOutputALCA_cff
import hltOutputMON_cff

# hltDebugOutput
hltDebugOutputBlocks = (
  # the DQM, HLTDQM and HLTMON streams have the HLT debug outputs used online; hltOutput8E29 hltOutput1E31 hltOutputHIon have the optional extra offline event content
  hltOutputMON_cff.block_hltOutputDQM.outputCommands,
  hltOutputMON_cff.block_hltOutputHLTDQM.outputCommands,
  hltOutputMON_cff.block_hltOutputHLTMON.outputCommands,
  hltOutputMON_cff.block_hltOutput8E29.outputCommands,
  hltOutputMON_cff.block_hltOutput1E31.outputCommands,
  hltOutputMON_cff.block_hltOutputHIon.outputCommands,
)
hltDebugOutputContent = buildPSet(hltDebugOutputBlocks)


# hltDebugWithAlCaOutput
hltDebugWithAlCaOutputBlocks = (
  # the DQM, HLTDQM and HLTMON streams have the HLT debug outputs used online; hltOutput8E29 hltOutput1E31 hltOutputHIon have the optional extra offline event content
  hltOutputMON_cff.block_hltOutputDQM.outputCommands,
  hltOutputMON_cff.block_hltOutputHLTDQM.outputCommands,
  hltOutputMON_cff.block_hltOutputHLTMON.outputCommands,
  hltOutputMON_cff.block_hltOutput8E29.outputCommands,
  hltOutputMON_cff.block_hltOutput1E31.outputCommands,
  hltOutputMON_cff.block_hltOutputHIon.outputCommands,
  # the ALCA streams have the AlCa outputs
  hltOutputALCA_cff.block_hltOutputALCAPHISYM.outputCommands,
  hltOutputALCA_cff.block_hltOutputALCAP0.outputCommands,
  hltOutputALCA_cff.block_hltOutputRPCMON.outputCommands,
)
hltDebugWithAlCaOutputContent = buildPSet(hltDebugWithAlCaOutputBlocks)


# hltDefaultOutput
hltDefaultOutputBlocks = (
  # the A stream has the HLT default output, with FEDs - strip out the FEDRawDataCollection keep statements for hltDefaultOutput
  hltOutputA_cff.block_hltOutputA.outputCommands,
)
hltDefaultOutputContent         = buildPSetWithoutRAWs(hltDefaultOutputBlocks)
hltDefaultOutputWithFEDsContent = buildPSet(hltDefaultOutputBlocks)


# define the CMSSW default event content configurations

# RAW event content
HLTriggerRAW = cms.PSet(
    outputCommands = cms.vstring()
)
HLTriggerRAW.outputCommands.extend(hltDefaultOutputWithFEDsContent.outputCommands)

# RECO event content
HLTriggerRECO = cms.PSet(
    outputCommands = cms.vstring()
)
HLTriggerRECO.outputCommands.extend(hltDefaultOutputContent.outputCommands)

# AOD event content
HLTriggerAOD = cms.PSet(
    outputCommands = cms.vstring()
)
HLTriggerAOD.outputCommands.extend(hltDefaultOutputContent.outputCommands)

# HLTDEBUG RAW event content
HLTDebugRAW = cms.PSet(
    outputCommands = cms.vstring()
)
HLTDebugRAW.outputCommands.extend(hltDebugWithAlCaOutputContent.outputCommands)

# HLTDEBUG FEVT event content
HLTDebugFEVT = cms.PSet(
    outputCommands = cms.vstring()
)
HLTDebugFEVT.outputCommands.extend(hltDebugWithAlCaOutputContent.outputCommands)


# dump the expanded event content configurations to a python configuration fragment
dump = open('HLTrigger_EventContent_cff.py', 'w')
dump.write('''import FWCore.ParameterSet.Config as cms

# EventContent for HLT related products.

# This file exports the following five EventContent blocks:
#   HLTriggerRAW  HLTriggerRECO  HLTriggerAOD (without DEBUG products)
#   HLTDebugRAW   HLTDebugFEVT                (with    DEBUG products)
#
# as these are used in Configuration/EventContent
#
# All else is internal and should not be used directly by non-HLT users.
#
''')
dump.write('HLTriggerRAW = %s\n\n'  % HLTriggerRAW)
dump.write('HLTriggerRECO = %s\n\n' % HLTriggerRECO)
dump.write('HLTriggerAOD = %s\n\n'  % HLTriggerAOD)
dump.write('HLTDebugRAW = %s\n\n'   % HLTDebugRAW)
dump.write('HLTDebugFEVT = %s\n'    % HLTDebugFEVT)
dump.close()
