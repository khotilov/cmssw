import FWCore.ParameterSet.Config as cms
import sys

from TauAnalysis.CandidateTools.tools.objSelConfigurator import *

#--------------------------------------------------------------------------------
# utility function for generation of a set of modules
# producing combinations of leptonic and hadronic decay products
# (pairs of pat::Electrons, pat::Muons and pat::Taus)
# for estimation of systematic uncertainties
#
# Author: Christian Veelken, UC Davis
#
#--------------------------------------------------------------------------------

class objProdConfigurator(cms._ParameterTypeBase):

    def __init__(self, objProd, systematics, pyModuleName = None):
        self.objProd = objProd
        self.systematics = systematics
        self.pyModuleName = pyModuleName,

    def _addModule(self, objProdItem, pyNameSpace, sysName, sysInputTags):        
        # create module
        moduleType = objProdItem.type_()
        module = cms.EDProducer(moduleType)

        # set module attributes
        # to default values
        for objProdAttrName in dir(objProdItem):
            objProdAttr = getattr(objProdItem, objProdAttrName)
            if isinstance(objProdAttr, cms._ParameterTypeBase) and not objProdAttrName in ["pluginName", "pluginType"]:
                setattr(module, objProdAttrName, objProdAttr)

        # set names of source collections
        # to objects shifted in energy/transverse momentum, theta, phi...
        for sysInputTagAttrName, sysInputTagValue in sysInputTags.items():
            setattr(module, sysInputTagAttrName, cms.InputTag(sysInputTagValue))
                
        moduleName = objSelConfigurator._composeModuleName(objSelConfigurator._getInstanceName(objProdItem, pyNameSpace), sysName)
        module.setLabel(moduleName)
               
        # register module in global python name-space
        pyModule = sys.modules[self.pyModuleName[0]]
        if pyModule is None:
            raise ValueError("'pyModuleName' Parameter invalid !!")
        setattr(pyModule, moduleName, module)

        # add module to sequence
        if self.sequence == None:
            self.sequence = module
        else:
            self.sequence *= module

    def configure(self, pyNameSpace = None):
        # configure set of modules
        # producing different combinations of leptonic and hadronic decay products
        # for estimation of systematic uncertainties

        # add original module (for production of central value) to sequence
        #self.sequence = None
        self.sequence = self.objProd

        for sysName, sysInputTags in self.systematics.items():
            self._addModule(self.objProd, pyNameSpace, sysName = sysName, sysInputTags = sysInputTags)

        return cms.Sequence(self.sequence)
