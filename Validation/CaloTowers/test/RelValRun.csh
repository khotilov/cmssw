#!/bin/tcsh

nohup cmsRun run_onRelVal_TTbar_MC_cfg.py      >& mcttb.txt &
nohup cmsRun run_onRelVal_TTbar_Startup_cfg.py >& stttb.txt &
nohup cmsRun run_onRelVal_QCD_MC_cfg.py        >& mcqcd.txt &
nohup cmsRun run_onRelVal_QCD_Startup_cfg.py   >& stqcd.txt &
