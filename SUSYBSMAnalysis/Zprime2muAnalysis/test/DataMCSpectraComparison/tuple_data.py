#!/usr/bin/env python

import sys, os, datetime
from tuple_common import process, crab_cfg

process.source.fileNames = ['/store/data/Run2011A/SingleMu/AOD/PromptReco-v5/000/172/163/A6E5AC06-18BC-E011-A9F4-003048D374F2.root']
process.source.lumisToProcess = cms.untracked.VLuminosityBlockRange('172163:115-172163:115')
process.GlobalTag.globaltag = 'GR_R_42_V13::All'

from SUSYBSMAnalysis.Zprime2muAnalysis.PATTools import removeMCUse
removeMCUse(process)

if __name__ == '__main__' and 'submit' in sys.argv:
    job_control_ex = '''
total_number_of_lumis = -1
lumis_per_job = %(lumis_per_job)s
%(lumi_mask)s
'''

    lumis_per_job = 200
    lumi_mask = ''

    create_only = 'create_only' in sys.argv
    just_testing = 'testing' in sys.argv
    scheduler = 'condor' if 'condor' in sys.argv else 'glite'

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
            if create_only:
                os.system('crab -create')
            else:
                os.system('crab -create -submit all')
            os.system('rm -f crab.cfg tmp.json')

    run_limits = []
    for x in sys.argv:
        try:
            run_limits.append(int(x))
        except ValueError:
            pass

    if run_limits:
        run1, run2 = run_limits
        if len(run_limits) != 2 or run1 > run2:
            raise RuntimeError('if any, must specify exactly two numeric arguments   min_run max_run  with max_run >= min_run')

        # Make up a fake lumi_mask that contains all lumis possible
        # for every run in the run range, since crab doesn't seem to
        # listen for a runselection parameter anymore.
        json = ['"%i": [[1,26296]]' % r for r in xrange(run_limits[0], run_limits[1] + 1)]
        open('tmp.json', 'wt').write('{' + ', '.join(json) + '}')
        lumi_mask = 'lumi_mask = tmp.json'

        name = 'SingleMu2011A_prompt_%i_%i_%s' % (run_limits[0], run_limits[1], datetime.datetime.today().strftime('%Y%m%d%H%M%S'))
        print name

        if run1 >= 170053 and run2 <= 172619:
            dataset = '/SingleMu/Run2011A-PromptReco-v5/AOD'
        elif run1 >= 172620:
            dataset = '/SingleMu/Run2011A-PromptReco-v6/AOD'
        else:
            raise ValueError("don't know how to do a run_limits production for run range [%i,%i]" % run_limits)
        
        tag = 'GR_R_42_V13'
        submit(locals())
    else:
        x = [
            ('SingleMuRun2011A_May10',   '/SingleMu/Run2011A-May10ReReco-v1/AOD', 'FT_R_42_V13A'),
            ('SingleMuRun2011A_Prompt4', '/SingleMu/Run2011A-PromptReco-v4/AOD',  'FT_R_42_V13A'),
            ]
        for name, dataset, tag in x:
            submit(locals())

