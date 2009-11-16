from RecoTauTag.TauAnalysisTools.tools.ntauples import *
from array import array
from ROOT import *

cuts = {}

#Individual discriminators
cuts['byLeadTrackPt']               =  TauNtupleSelection("$ByLeadTrackPt")
cuts['byLeadPionPt']                =  TauNtupleSelection("$ByLeadPionPt")
cuts['byIsolation']                 =  TauNtupleSelection("$ByIsolation")
cuts['againstMuon']                 =  TauNtupleSelection("$AgainstMuon")
cuts['againstElectron']             =  TauNtupleSelection("$AgainstElectron")
cuts['charge']                      =  TauNtupleSelection("$charge == 1")
cuts['nTracks']                     =  TauNtupleSelection("$nTrks == 1 || $nTrks == 3")
cuts['byTaNCfrOne']                 =  TauNtupleSelection('$ByTaNCfrOne')
cuts['byTaNCfrHalf']                =  TauNtupleSelection('$ByTaNCfrHalf')
cuts['byTaNCfrQuarter']             =  TauNtupleSelection('$ByTaNCfrQuarter')
cuts['byTaNCfrTenth']               =  TauNtupleSelection('$ByTaNCfrTenth')

cuts['chargeAndTracks']             =  cuts['charge']*cuts['nTracks']

cuts['standardChain']               =  cuts['byLeadPionPt']*cuts['byIsolation']*cuts['chargeAndTracks']
cuts['standardChainNoMuon']         =  cuts['standardChain']*cuts['againstMuon']
cuts['standardChainNoElectron']     =  cuts['standardChain']*cuts['againstElectron']

cuts['tancChainOne']                =  cuts['chargeAndTracks']*cuts['byTaNCfrOne']
cuts['tancChainHalf']               =  cuts['chargeAndTracks']*cuts['byTaNCfrHalf']
cuts['tancChainQuarter']            =  cuts['chargeAndTracks']*cuts['byTaNCfrQuarter']
cuts['tancChainTenth']              =  cuts['chargeAndTracks']*cuts['byTaNCfrTenth']

cuts['tancChainOneNoMuon']          =  cuts['chargeAndTracks']*cuts['byTaNCfrOne']*cuts['againstMuon']
cuts['tancChainHalfNoMuon']         =  cuts['chargeAndTracks']*cuts['byTaNCfrHalf']*cuts['againstMuon']
cuts['tancChainQuarterNoMuon']      =  cuts['chargeAndTracks']*cuts['byTaNCfrQuarter']*cuts['againstMuon']
cuts['tancChainTenthNoMuon']        =  cuts['chargeAndTracks']*cuts['byTaNCfrTenth']*cuts['againstMuon']

cuts['tancChainOneNoElectron']      =  cuts['chargeAndTracks']*cuts['byTaNCfrOne']*cuts['againstElectron']
cuts['tancChainHalfNoElectron']     =  cuts['chargeAndTracks']*cuts['byTaNCfrHalf']*cuts['againstElectron']
cuts['tancChainQuarterNoElectron']  =  cuts['chargeAndTracks']*cuts['byTaNCfrQuarter']*cuts['againstElectron']
cuts['tancChainTenthNoElectron']    =  cuts['chargeAndTracks']*cuts['byTaNCfrTenth']*cuts['againstElectron']


def make_plots(files_and_weights, 
      output_file = "fakeWeights.root",
      x_expr="$jetPt", y_expr="$eta", z_expr="$jetWidth", 
      x_bins=[], y_bins=[], z_bins=[], selections=cuts):

   output = TFile(output_file, "RECREATE")

   x_bins = array('d', x_bins)
   y_bins = array('d', y_bins)
   z_bins = array('d', z_bins)

   def histo_maker(name, x_bins=x_bins, y_bins=y_bins, z_bins=z_bins):
      return TH3F(name, name, 
            len(x_bins)-1, x_bins,
            len(y_bins)-1, y_bins,
            len(z_bins)-1, z_bins)

   # Add the denominator
   selections["denominator"] = "1"

   # Build our histograms
   histograms = {}
   for cut_name, cut in selections.iteritems():
      histograms[cut_name] = {}
      histograms[cut_name]["histo"] = histo_maker(cut_name)
      histograms[cut_name]["selection"] = str(cut)
      # Make monitor histogram
      monitor_name = "%s_monitor" % cut_name
      histograms[cut_name]["monitor"] = TH1F(monitor_name, monitor_name, 100, -0.5, 99.5)

   for file, weight in files_and_weights:
      print "Processing file: ", file
      # Get events
      my_file = TFile(file, "READ")
      events = my_file.Get("Events")
      ntuples = TauNtupleManager(events)
      for name, histo_info in histograms.iteritems():
         print "Filling ", name
         output.cd()
         draw(events, ntuples.shrinkingConePFTau,
            # Yes, it actually does go z, y, x :/
            expr = '%s:%s:%s' % (z_expr, y_expr, x_expr),
            selection = histo_info["selection"],
            output_hist = "+%s" % name)
   
   # Fill monitor histograms
   for cut_name, cut in selections.iteritems():
      histo = histograms[cut_name]["histo"]
      monitor = histograms[cut_name]["monitor"]
      nBinsX = histo.GetNbinsX()+1
      nBinsY = histo.GetNbinsY()+1
      nBinsZ = histo.GetNbinsZ()+1
      for x in range(1, nBinsX):
         for y in range(1, nBinsY):
            for z in range(1, nBinsZ):
               monitor.Fill(histo.GetBinContent(x,y,z))

   output.Write()


if __name__ == "__main__":
   gROOT.SetBatch(True)
   make_plots( 
         [ ("taste_ntupled.root", 1.0) ],
         x_bins = [0, 10, 20, 50, 80, 120, 200],
         y_bins = [0, 0.5, 1.0, 1.5, 2.0, 2.5],
         z_bins = [0, 0.1, 0.2, 0.4, 0.6]
         )


