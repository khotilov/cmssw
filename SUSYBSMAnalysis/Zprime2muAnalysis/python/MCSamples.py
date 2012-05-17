#!/usr/bin/env python

import os
from SUSYBSMAnalysis.Zprime2muAnalysis.tools import big_warn, files_from_dbs
from SUSYBSMAnalysis.Zprime2muAnalysis.crabtools import dataset_from_publish_log

class sample(object):
    def __init__(self, name, nice_name, dataset, nevents, color, syst_frac, cross_section, k_factor=1, filenames=None, scheduler='condor', hlt_process_name='HLT', ana_dataset=None, is_zprime=False):
        self.name = name
        self.nice_name = nice_name
        self.dataset = dataset
        self.nevents = nevents
        self.color = color
        self.syst_frac = syst_frac
        self.cross_section = cross_section
        self.k_factor = k_factor
        self.filenames_ = filenames
        self.scheduler = scheduler
        self.hlt_process_name = hlt_process_name
        self.ana_dataset = ana_dataset
        self.is_zprime = is_zprime

    @property
    def partial_weight(self):
        return self.cross_section / float(self.nevents) * self.k_factor # the total weight is partial_weight * integrated_luminosity

    @property
    def filenames(self):
        # Return a list of filenames for running the histogrammer not
        # using crab.
        if self.filenames_ is not None:
            return self.filenames_
        return files_from_dbs(self.ana_dataset, ana02=True)

    def __getitem__(self, key):
        return getattr(self, key)

    def _dump(self, redump_existing=False):
        dst = os.path.join('/uscmst1b_scratch/lpc1/3DayLifetime/tucker', self.name)
        os.system('mkdir ' + dst)
        for fn in self.filenames:
            print fn
            if redump_existing or not os.path.isfile(os.path.join(dst, os.path.basename(fn))):
                os.system('dccp ~%s %s/' % (fn,dst))

class tupleonlysample(sample):
    def __init__(self, name, dataset, scheduler='condor', hlt_process_name='HLT'):
        super(tupleonlysample, self).__init__(name, 'dummy', dataset, 1, 1, 1, 1, scheduler=scheduler, hlt_process_name=hlt_process_name)

# https://twiki.cern.ch/twiki/bin/viewauth/CMS/StandardModelCrossSectionsat8TeV for xsecs (all below in pb)
# commented out samples are not available for Summer12/8 TeV
samples = [
    sample('zmumu',     '#gamma/Z #rightarrow #mu^{+}#mu^{-}',              '/DYToMuMu_M_20_TuneZ2star_8TeV_pythia6/Summer12-PU_S7_START50_V15-v1/AODSIM',             1963296, 432, 0.05, 1915.   - 9.414), # The subtraction here and below is because we drop the events with generated M in ranges in the other datasets, e.g. drop here M > 120.
    sample('dy120',     'DY120',                                            '/DYToMuMu_M_120_TuneZ2star_8TeV_pythia6/Summer12-PU_S7_START52_V9-v1/AODSIM',               49454, 433, 0.05, 9.414   - 1.176,    k_factor=1.3),
    sample('dy200',     'DY200',                                            '/DYToMuMu_M_200_TuneZ2star_8TeV_pythia6/Summer12-PU_S7_START52_V9-v1/AODSIM',               50400, 434, 0.05, 1.176   - 0.03557,  k_factor=1.3),
    sample('dy500',     'DY500',                                            '/DYToMuMu_M_500_TuneZ2star_8TeV_pythia6/Summer12-PU_S7_START52_V9-v1/AODSIM',               50560, 435, 0.05, 0.03557 - 0.00451,  k_factor=1.3),
    sample('dy800',     'DY800',                                            '/DYToMuMu_M_800_TuneZ2star_8TeV_pythia6/Summer12-PU_S7_START52_V9-v1/AODSIM',               50400,   3, 0.05, 0.00451 - 3.55E-4,  k_factor=1.3), # skip 1 TeV
##    sample('dy1000',    'DY1000',                                           '/DYToMuMu_M-1000_TuneZ2star_8TeV-pythia6/Summer12-PU_S7_START52_V9-v1/AODSIM',               XXXXX, 437, 0.05, 0.00149 - 3.55E-4,  k_factor=1.3),
    sample('dy1300',    'DY1300',                                           '/DYToMuMu_M-1300_TuneZ2star_8TeV-pythia6/Summer12-PU_S7_START52_V9-v1/AODSIM',              50560,   8, 0.05, 3.55E-4 - 9.28E-5,  k_factor=1.3),
    sample('dy1600',    'DY1600',                                           '/DYToMuMu_M-1600_TuneZ2star_8TeV-pythia6/Summer12-PU_S7_START52_V9-v1/AODSIM',              50112,  37, 0.05, 9.28E-5 - 0,        k_factor=1.3),
    sample('ttbar',     't#bar{t}',                                         '/TTJets_TuneZ2star_8TeV-madgraph-tauola/Summer12-PU_S7_START52_V9-v1/AODSIM',             6736135,   2, 0.15, 225.2),
    sample('tW',        'tW',                                               '/T_tW-channel-DR_TuneZ2star_8TeV-powheg-tauola/Summer12-PU_S7_START52_V9-v1/AODSIM',       497658,   1, 0.03, 11.2),
    sample('tbarW',     'tbarW',                                            '/Tbar_tW-channel-DR_TuneZ2star_8TeV-powheg-tauola/Summer12-PU_S7_START52_V9-v1/AODSIM',    493460,  12, 0.03, 11.2),
    sample('ww',        'WW',                                               '/WW_TuneZ2star_8TeV_pythia6_tauola/Summer12-PU_S7_START50_V15-v1/AODSIM',                10000431,   4, 0.035, 57),
#    sample('ww',        'WW',                                               '/WW_TuneZ2star_8TeV_pythia6_tauola/Summer12-PU_S7_START52_V9-v1/AODSIM',                 10000431,   4, 0.035, 57.1),
    sample('wz',        'WZ',                                               '/WZ_TuneZ2star_8TeV_pythia6_tauola/Summer12-PU_S7_START52_V9-v1/AODSIM',                  9996622,  30, 0.038, 32.3),
    sample('zz',        'ZZ',                                               '/ZZ_TuneZ2star_8TeV_pythia6_tauola/Summer12-PU_S7_START50_V15-v1/AODSIM',                 9799908,   6, 0.025,  8),
#    sample('zz',        'ZZ',                                               '/ZZ_TuneZ2star_8TeV_pythia6_tauola/Summer12-PU_S7_START52_V9-v1/AODSIM',                  9799908,   6, 0.025,  8.3),
    sample('ztautau',   'Z #rightarrow #tau^{+}#tau^{-}',                   '/DYToTauTau_M_20_TuneZ2star_8TeV_pythia6_tauola/Summer12-PU_S7_START52_V9-v1/AODSIM',     1987776,  46, 0.05, 1915.),
    sample('wjets',     'W+jets',                                           '/WJetsToLNu_TuneZ2Star_8TeV-madgraph-tarball/Summer12-PU_S7_START52_V9-v1/AODSIM',       18393090,   3, 0.05, 36257.),
    sample('inclmu15',  'QCD',                                              '/QCD_Pt_20_MuEnrichedPt_15_TuneZ2star_8TeV_pythia6/Summer12-PU_S7_START52_V9-v1/AODSIM',  7529312, 801, 0.1,  3.64E8 * 3.7E-4),
##    sample('zssm1000',  'Z\'_{SSM} (1 TeV) #rightarrow #mu^{+}#mu^{-}',     '/ZprimeSSMToMuMu_M-1000_TuneZ2_7TeV-pythia6/Summer11-PU_S4_START42_V11-v1/AODSIM',          20328,  38, 0.05,  0.089,   k_factor=1.3, is_zprime=True),
#    sample('zpsi750',   'Z\'_{#psi} (0.75 TeV) #rightarrow #mu^{+}#mu^{-}', '/ZprimePSIToMuMu_M-750_TuneZ2star_8TeV-pythia6/Summer12-PU_S7_START52_V9-v1/AODSIM',        25040,  48, 0.05,  0.14,    k_factor=1.3, is_zprime=True),
#    sample('zpsi1000',  'Z\'_{#psi} (1 TeV) #rightarrow #mu^{+}#mu^{-}',    '/ZprimePSIToMuMu_M-1000_TuneZ2star_8TeV-pythia6/Summer12-PU_S7_START52_V9-v1/AODSIM',       25040,  48, 0.05,  0.0369,  k_factor=1.3, is_zprime=True),
#    sample('zpsi1250',  'Z\'_{#psi} (1.25 TeV) #rightarrow #mu^{+}#mu^{-}', '/ZprimePSIToMuMu_M-1250_TuneZ2star_8TeV-pythia6/Summer12-PU_S7_START52_V9-v1/AODSIM',       25344,  48, 0.05,  0.0129,  k_factor=1.3, is_zprime=True),
#    sample('zpsi1500',  'Z\'_{#psi} (1.5 TeV) #rightarrow #mu^{+}#mu^{-}',  '/ZprimePSIToMuMu_M-1500_TuneZ2star_8TeV-pythia6/Summer12-PU_S7_START52_V9-v1/AODSIM',       25344,  48, 0.05,  0.00433, k_factor=1.3, is_zprime=True),
#    sample('zpsi1750',  'Z\'_{#psi} (1.75 TeV) #rightarrow #mu^{+}#mu^{-}', '/ZprimePSIToMuMu_M-1750_TuneZ2star_8TeV-pythia6/Summer12-PU_S7_START52_V9-v1/AODSIM',       25272,  48, 0.05,  0.00172, k_factor=1.3, is_zprime=True),
#    sample('zpsi2000',  'Z\'_{#psi} (2 TeV) #rightarrow #mu^{+}#mu^{-}',    '/ZprimePSIToMuMu_M-2000_TuneZ2star_8TeV-pythia6/Summer12-PU_S7_START52_V9-v1/AODSIM',       25092,  48, 0.05,  6.88E-4, k_factor=1.3, is_zprime=True),
#    sample('zpsi2250',  'Z\'_{#psi} (2.25 TeV) #rightarrow #mu^{+}#mu^{-}', '/ZprimePSIToMuMu_M-2250_TuneZ2star_8TeV-pythia6/Summer12-PU_S7_START52_V9-v1/AODSIM',       25104,  48, 0.05,  2.93E-4, k_factor=1.3, is_zprime=True),
#    sample('zpsi2500',  'Z\'_{#psi} (2.5 TeV) #rightarrow #mu^{+}#mu^{-}',  '/ZprimePSIToMuMu_M-2500_TuneZ2star_8TeV-pythia6/Summer12-PU_S7_START52_V9-v1/AODSIM',       25344,  48, 0.05,  1.27E-4, k_factor=1.3, is_zprime=True),
#    sample('zpsi2750',  'Z\'_{#psi} (2.75 TeV) #rightarrow #mu^{+}#mu^{-}', '/ZprimePSIToMuMu_M-2750_TuneZ2star_8TeV-pythia6/Summer12-PU_S7_START52_V9-v1/AODSIM',       25376,  48, 0.05,  5.55E-5, k_factor=1.3, is_zprime=True),
#    sample('zpsi3000',  'Z\'_{#psi} (3 TeV) #rightarrow #mu^{+}#mu^{-}',    '/ZprimePSIToMuMu_M-3000_TuneZ2star_8TeV-pythia6/Summer12-PU_S7_START52_V9-v1/AODSIM',       25040,  48, 0.05,  2.5E-5,  k_factor=1.3, is_zprime=True),
]

samples.reverse()

for sample in samples:
    exec '%s = sample' % sample.name
    sample.ana_dataset = '/%s/slava-datamc_%s-caca636ea661546409f4073c061b3e20/USER' % (sample.dataset.split('/')[1], sample.name)

ttbar.ana_dataset = '/TTJets_TuneZ2star_8TeV-madgraph-tauola/slava-datamc_ttbar-2e54d35c5572bb8d0ce3ffc532e82068/USER'
wz.ana_dataset     = '/WZ_TuneZ2star_8TeV_pythia6_tauola/slava-datamc_wz-2e54d35c5572bb8d0ce3ffc532e82068/USER'

#big_warn('nothing')

__all__ = ['samples'] + [s.name for s in samples]


if __name__ == '__main__':
    if False:
        from dbstools import dbsparents
        for s in samples:
            print s.dataset
            parents = dbsparents(s.dataset)
            for parent in parents:
                for line in os.popen('dbss rel %s' % parent):
                    if 'CMSSW' in line:
                        print parent, line,
            print

    if False:
        import os
        from dbstools import dbsparents
        for s in [ww,wz,zz]:
            print s.dataset
            parents = dbsparents(s.dataset)
            print parents
            os.system('dbsconfig %s > %s' % (parents[-1], s.name))

        os.system('dbss nevents %s' % x.replace('RECO','RAW'))
        os.system('dbss nevents %s' % x)

    if False:
        import os
        from dbstools import dbsparents
        for s in samples:
            print s.dataset
            def fuf(y):
                x = os.popen(y).read()
                for line in x.split('\n'):
                    try:
                        print int(line)
                    except ValueError:
                        pass
            fuf('dbss nevents %s' % s.dataset)
            fuf('dbss nevents %s' % s.dataset.replace('AODSIM','GEN-SIM-RECO'))

    if False:
        for s in samples:
            print s.name
            os.system('grep "total events" ~/nobackup/crab_dirs/384p3/publish_logs/publish.crab_datamc_%s' % s.name)
            os.system('grep "total events" ~/nobackup/crab_dirs/413p2/publish_logs/publish.crab_datamc_%s' % s.name)
            print

    if False:
        os.system('mkdir ~/scratch/wjets')
        for fn in wjets.filenames:
            assert fn.startswith('/store')
            fn = '/pnfs/cms/WAX/11' + fn
            cmd = 'dccp %s ~/scratch/wjets/' % fn
            print cmd
            os.system(cmd)

    if False:
        for s in samples:
            print s.name
            os.system('dbss site %s' % s.dataset)
            print

    if False:
        for s in samples:
            if s.ana_dataset is None:
                continue
            c = []
            for line in os.popen('dbss ana02 find file.numevents where dataset=%s' % s.ana_dataset):
                try:
                    n = int(line)
                except ValueError:
                    continue
                c.append(n)
            c.sort()
            print s.name, c
