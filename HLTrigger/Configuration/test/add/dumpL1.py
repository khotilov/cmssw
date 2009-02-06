#! /usr/bin/env python

import sys

import FWCore.ParameterSet.Config as cms
from FWCore.ParameterSet.SequenceTypes import *
from table import process

for (name,path) in process.paths.iteritems():
  modules = []
  visitor = ModuleNodeVisitor(modules)
  path.visit(visitor)  
 
  seeds = [ module.L1SeedsLogicalExpression for module in modules if module.type_() == 'HLTLevel1GTSeed' ]
  if seeds:
    print '%-32s\t%s' % (name, seeds[0].value()) 
