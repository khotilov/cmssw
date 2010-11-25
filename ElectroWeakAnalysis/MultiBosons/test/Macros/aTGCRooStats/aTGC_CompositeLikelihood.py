import sys
import os
from math import sqrt
from optparse import OptionParser

import ROOT

def main(options,args):

    par1Max = float(options.par1Max)
    par2Max = float(options.par2Max)

    output = ROOT.TFile.Open(options.workspaceName+".root","RECREATE")

    combinedWorkspace = ROOT.RooWorkspace(options.workspaceName)
    
    workspaces=[]
    NLLs = []

    par1 = ROOT.RooRealVar(options.couplingType+'_'+options.par1Name,options.par1Name,0,-par1Max,par1Max)
    par2 = ROOT.RooRealVar(options.couplingType+'_'+options.par2Name,options.par2Name,0,-par2Max,par2Max)
    
    for arg in args:
        thefile = ROOT.TFile.Open(arg)
        workspaces.append(thefile.Get(arg[:arg.find('.root')]))
        print 'Loaded workspace: ',arg[:arg.find('.root')]

    for ws in workspaces:
        #ws.var('err_x_gs').setConstant(False)
        #ws.var('err_x_gl').setConstant(False)
        #ws.var('err_x_gb').setConstant(False)

        #change workspace par1,par2 to point at the ones here
        getattr(combinedWorkspace,'import')(ws.function('nll_TopLevelPdf_aTGCDataUnitWeight_with_constr'),
                                            ROOT.RooFit.RenameAllNodes(ws.GetName()),
                                            ROOT.RooFit.RenameAllVariablesExcept(ws.GetName(),
                                                                                 par1.GetName()+','+par2.GetName())
                                            )
        
        NLLs.append('nll_TopLevelPdf_aTGCDataUnitWeight_with_constr_'+ws.GetName())
        
    print NLLs

    combinedNLL = 'sum::combinedNLL('

    for i in range(len(NLLs)):
        if i != len(NLLs)-1:
            combinedNLL += NLLs[i]+','
        else:
            combinedNLL += NLLs[i]
    
    combinedNLL+=')'

    print combinedNLL

    combinedWorkspace.factory(combinedNLL)

    combinedWorkspace.function('combinedNLL').Print()

    minuit = ROOT.RooMinuit(combinedWorkspace.function('combinedNLL'))

    minuit.setErrorLevel(.5) # force to NLL definition
    minuit.hesse()
    minuit.migrad()
    minuit.hesse()

    profileLL = combinedWorkspace.function('combinedNLL').createProfile(ROOT.RooArgSet(par1,par2))
    profileLL.getVal() # to cache the values of the constrained params
    
    level_68 = ROOT.TMath.ChisquareQuantile(.68,2)/2.0 # delta NLL for 68% confidence level for -log(LR)
    level_95 = ROOT.TMath.ChisquareQuantile(.95,2)/2.0 # delta NLL for 95% confidence level for -log(LR)

    print '68% CL Delta-NLL=',level_68
    print '95% CL Delta-NLL=',level_95

    fArgs = combinedWorkspace.function('combinedNLL').getParameters(ROOT.RooArgSet())
    
    profMinuit = profileLL.minuit()    
    profMinuit.setStrategy(2)
    profMinuit.setErrorLevel(.5) # force to likelihood error def.
    profMinuit.setPrintLevel(1)
                             
    profMinuit.migrad()
    profMinuit.hesse()
    #for ws in workspaces:        
    #    combinedWorkspace.var('err_x_gl_'+ws.GetName()).setConstant(True)
    #    combinedWorkspace.var('err_x_gs_'+ws.GetName()).setConstant(True)
    #    combinedWorkspace.var('err_x_gb_'+ws.GetName()).setConstant(True)    
    profMinuit.minos(ROOT.RooArgSet(fArgs.find(par1.GetName()),
                                    fArgs.find(par2.GetName())))    
    
    thePlot = profMinuit.contour(fArgs.find(par1.GetName()),
                                 fArgs.find(par2.GetName()),
                                 sqrt(2*level_68),sqrt(2*level_95))
    
    theCanvas = ROOT.TCanvas('contours','',500,500)
    
    thePlot.SetTitle('68% & 95% CL on the Best Fit Values of '+options.par1Name+' and '+options.par2Name)
    thePlot.Draw()
    
    theCanvas.Print(options.workspaceName+'_contour.root')

    output.cd()
    combinedWorkspace.Write()
    output.Close()
        
    #makePlots(theLHplot)    
    #really, that's all I had to do??
        
def makePlots(LLInterval,options):
    print "not done yet"
    theCanvas = ROOT.TCanvas("Likelihood Plot")

if __name__ == "__main__":
    parser = OptionParser(description="%prog : Calculates the composite likelihood of various aTGC measurements.",
                          usage="%prog <workspace1>.root <workspace2>.root ...  --couplingType=... --par1Max=# --par2Max=#")
    #    parser.add_option("--POI",dest="POIlist",help="The names of the parameters of interest")
    parser.add_option("--workspaceName",dest="workspaceName",help="name of the combined workspace")
    parser.add_option("--couplingType",dest="couplingType",help="ZZg or Zgg couplings?")
    parser.add_option("--par1Name",dest="par1Name",help="Name of parameter 1")
    parser.add_option("--par1Max",dest="par1Max",help="Bound on |par1|")
    parser.add_option("--par2Name",dest="par2Name",help="Name of parameter 2")
    parser.add_option("--par2Max",dest="par2Max",help="Bound on |par2|")
   
    (options,args) = parser.parse_args()

    miss_options = False

    #    if options.POIlist is None:
    #        print 'Need to specify --POI'
    #        miss_options=True
    if options.workspaceName is None:
        print "Need to specify --workspaceName"
        miss_options=True
    if options.couplingType is None:
        print 'Need to specify --couplingType (ZZg or Zgg)'
        miss_options=True
    if options.par1Max is None:
        print 'Need to specify --par1Max'
        miss_options=True
    if options.par2Max is None:
        print 'Need to specify --par2Max'
        miss_options=True
    if options.par1Name is None:
        print 'Need to specify --par1Name'
        miss_options=True
    if options.par2Name is None:
        print 'Need to specify --par2Name'
        miss_options=True

    if len(args) == 0:
        print 'You need to pass at least one root file!'
        miss_options=True
     
    if miss_options:
        exit(1)

    main(options,args)
