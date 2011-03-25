#!/usr/bin/env python

import sys, os, datetime
from tuple_common import process, crab_cfg

process.source.fileNames = ['/store/data/Run2011A/SingleMu/AOD/PromptReco-v1/000/160/405/4642D954-D64F-E011-8280-003048F024DE.root']
process.source.fileNames = ['/store/data/Run2011A/SingleMu/AOD/PromptReco-v1/000/160/957/FCEB1657-5E55-E011-BCD7-000423D9890C.root']
process.maxEvents.input = 2000
process.GlobalTag.globaltag = 'GR_R_311_V2::All'

from SUSYBSMAnalysis.Zprime2muAnalysis.PATTools import removeMCUse
removeMCUse(process)

if __name__ == '__main__' and 'submit' in sys.argv:
    scheduler = 'condor'
    job_control_ex = '''
total_number_of_lumis = -1
lumis_per_job = %(lumis_per_job)s
%(lumi_mask)s
'''

    lumis_per_job = 200
    lumi_mask = ''

    if False:
        job_control_ex = '''
total_number_of_lumis = -1
number_of_jobs = 100
'''

    just_testing = 'testing' in sys.argv

    def submit(d):
        new_py = open('tuple_data.py').read()
        new_py += '\n\nprocess.GlobalTag.globaltag = "%(tag)s::All"\n' % d
        pset = 'psets/tuple_data_crab_%(name)s.py' % d
        open(pset, 'wt').write(new_py)

        job_control = job_control_ex % d
        for k,v in locals().iteritems():
            d[k] = v
        open('crab.cfg', 'wt').write(crab_cfg % d)
        if not just_testing:
            os.system('crab -create -submit all')
            os.system('rm -f crab.cfg tmp.json')

    run_limits = []
    for x in sys.argv:
        try:
            run_limits.append(int(x))
        except ValueError:
            pass

    if run_limits:
        if len(run_limits) != 2:
            raise RuntimeError('if any, must specify exactly two numeric arguments: min_run max_run')

        # Make up a fake lumi_mask that contains all lumis possible
        # for every run in the run range, since crab doesn't seem to
        # listen for a runselection parameter anymore.
        json = ['"%i": [[1,26296]]' % r for r in xrange(run_limits[0], run_limits[1] + 1)]
        open('tmp.json', 'wt').write('{' + ', '.join(json) + '}')
        lumi_mask = 'lumi_mask = tmp.json'

        name = 'promptB_' + datetime.datetime.today().strftime('%Y%m%d_%H%M%S')
        print name
        dataset = '/Mu/Run2010B-PromptReco-v2/RECO'
        tag = 'GR10_P_V10'
        submit(locals())
    else:
        x = [
            ('SingleMuRun2011A_testgrid1', '/SingleMu/Run2011A-PromptReco-v1/AOD', 'GR_R_311_V2'),
            ]
        for name, dataset, tag in x:
            scheduler = 'glite' if 'grid' in name else 'condor'
            submit(locals())

