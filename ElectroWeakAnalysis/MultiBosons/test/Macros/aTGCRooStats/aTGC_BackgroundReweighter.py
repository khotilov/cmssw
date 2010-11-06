import sys
import os
from optparse import OptionParser

import ROOT

def main(options,args):

    out = ROOT.TFile.Open(options.output,'RECREATE')

    ROOT.gROOT.ProcessLine(
        "struct TreeContents {\
        Float_t photonEt;\
        Double_t weight;\
        }")

    treecontents = ROOT.TreeContents()

    inTreeContents = ROOT.TreeContents()
        
    outTree = ROOT.TTree(options.treeName,'The Background')
    outTree.Branch('photonEt',ROOT.AddressOf(treecontents,'photonEt'),'photonEt/F')
    outTree.Branch('weight',ROOT.AddressOf(treecontents,'weight'),'weight/D')

    for f in args:        
        print f
        currentFile = ROOT.TFile.Open(f)
        currentTree = currentFile.Get(options.treeName)
        currentTree.SetBranchAddress('photonEt',ROOT.AddressOf(inTreeContents,'photonEt'))        
        treecontents.weight = currentTree.GetWeight()*float(options.intLumi)/100.0

        for i in range(currentTree.GetEntries()):
            currentTree.GetEntry(i)
            treecontents.photonEt = inTreeContents.photonEt
            outTree.Fill()

        currentFile.Close()

    out.cd()
    outTree.Write()    
    out.Close()
    

if __name__ == "__main__":
    parser = OptionParser(description="%prog : A RooStats Implementation of Anomalous Triple Gauge Coupling Analysis.",
                          usage="%prog file1.root file2.root ... --output=<file> --treeName=<name>")
    parser.add_option("--output",dest="output",help="The name of your output file.")
    parser.add_option("--treeName",dest="treeName",help="The name of input TTrees.")
    parser.add_option("--intLumi",dest="intLumi",help="Integrated luminosity to scale to.")

    (options,args) = parser.parse_args()

    miss_options = False

    if options.output is None:
        print 'Need to specify --output'
        miss_options=True
    if options.treeName is None:
        print 'Need to specify --treeName'
        miss_options=True
    if options.intLumi is None:
        print 'Need to specify --intLumi'
        miss_options=True
    if len(args) == 0:
        print 'You need to pass at least one root file!'
        miss_options=True

    if miss_options:
        exit(1)

    main(options,args)
