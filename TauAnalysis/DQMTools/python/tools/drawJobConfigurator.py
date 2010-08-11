import FWCore.ParameterSet.Config as cms
import copy

from TauAnalysis.DQMTools.tools.composeSubDirectoryName import composeSubDirectoryName

#--------------------------------------------------------------------------------
# utility function for generation of drawJob configurations
# for DQMHistPlotter tool
#
# Author: Christian Veelken, UC Davis
#
#--------------------------------------------------------------------------------

class drawJobConfigurator(cms._ParameterTypeBase):

    def __init__(self, template, dqmDirectory):

        self.setTemplate(template)
        self.setDQMdirectory(dqmDirectory)

        self.dqmSubDirectories = []
        self.plots = []
        
        self.drawJobs = cms.PSet()

    def setTemplate(self, template):
        self.template = template

    def setDQMdirectory(self, dqmDirectory):
        self.dqmDirectory = dqmDirectory

        if not self.dqmDirectory.endswith("/"):
            self.dqmDirectory += "/"

    def add(self, afterCut = None, beforeCut = None, plot = None, plots = None):
        # configure drawJob
        # and add to configuration object

        # check validity of parameters passed as function arguments
        if self.template is None:
            raise ValueError("Invalid 'template' Parameter !!")
        if plot is None and plots is None:
            raise ValueError("Invalid 'plot' and 'plots' Parameters !!")

        # check if need to call recursively
        # in case of multiple plots
        if plots is not None:
            for plot in plots:
                self.add(afterCut = afterCut, beforeCut = beforeCut, plot = plot)
            return
        
        dqmSubDirectory = composeSubDirectoryName(afterCut = afterCut, beforeCut = beforeCut)
        self.dqmSubDirectories.append(dqmSubDirectory)
        
        self.plots.append(plot)
        
    def configure(self):
        # return configuration object
        # for set of drawJobs

        for iPlot in range(len(self.plots)):

            drawJob = copy.deepcopy(self.template)

            dqmSubDirectory = self.dqmSubDirectories[iPlot]

            dqmDirectory = self.dqmDirectory
            if dqmSubDirectory != "":
                dqmDirectory += dqmSubDirectory
            if not dqmDirectory.endswith("/"):
                dqmDirectory += "/"

            plot = self.plots[iPlot]

            dqmMonitorElement = dqmDirectory + getattr(plot, "meName")

            setattr(drawJob.plots, "dqmMonitorElements", cms.vstring([ dqmMonitorElement, ]))
            if hasattr(plot, "PAR"):
                setattr(drawJob, "parameter", cms.vstring(getattr(plot, "PAR")))
            setattr(drawJob, "title", cms.string(getattr(plot, "title")))
            setattr(drawJob, "xAxis", cms.string(getattr(plot, "xAxis")))

            # add drawJob configuration to set of drawJobs
            setattr(self.drawJobs, getattr(plot, "name"), drawJob)

        return self.drawJobs

#--------------------------------------------------------------------------------
# auxiliary class stroring configuration parameters for one drawJob
#--------------------------------------------------------------------------------

class drawJobConfigEntry(cms._ParameterTypeBase):

    def __init__(self, meName, title, xAxis, name, PAR = None):
        self.meName = meName
        self.title = title
        self.xAxis = xAxis
        self.name = name
        if PAR is not None:
            self.PAR = PAR
