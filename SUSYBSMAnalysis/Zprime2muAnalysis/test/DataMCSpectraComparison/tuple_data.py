#!/usr/bin/env python

import sys, os, datetime
from tuple_common import process, crab_cfg

just_testing = 'testing' in sys.argv

process.source.fileNames = ['/store/data/Run2010B/Mu/RECO/PromptReco-v2/000/147/043/B0286A14-17CE-DF11-A65D-001D09F292D1.root']
process.maxEvents.input = 100
tag = process.GlobalTag.globaltag = 'GR10_P_V10::All'

from SUSYBSMAnalysis.Zprime2muAnalysis.PATTools import removeMCUse
removeMCUse(process)

if __name__ == '__main__' and 'submit' in sys.argv:
    scheduler = 'condor'
    job_control_ex = '''
total_number_of_lumis = -1
lumis_per_job = 50
%(lumi_mask)s
'''

    run_limits = []
    for x in sys.argv:
        try:
            run_limits.append(int(x))
        except ValueError:
            pass
    if run_limits:
        if len(run_limits) != 2:
            raise RuntimeError('if any, must specify exactly two numeric arguments: min_run max_run')
        json = ['"%i": [[1,26296]]' % r for r in xrange(run_limits[0], run_limits[1] + 1)]
        open('tmp.json', 'wt').write('{' + ', '.join(json) + '}')
        lumi_mask = 'lumi_mask = tmp.json'
    else:
        lumi_mask = ''

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

    if run_limits:
        # Runs supplied in argv -- running on new data.
        name = 'promptB_' + datetime.datetime.today().strftime('%Y%m%d_%H%M%S')
        print name
        dataset = '/Mu/Run2010B-PromptReco-v2/RECO'
        submit(locals())
    else:
        x = [
            ('jul15',    '/Mu/Run2010A-Jul15thReReco-v1/RECO',   'GR_R_37X_V6D'),
            ('prompt',   '/Mu/Run2010A-PromptReco-v4/RECO',      'GR10_P_V7'),
            ]
        for name, dataset, tag in x:
            submit(locals())

        

