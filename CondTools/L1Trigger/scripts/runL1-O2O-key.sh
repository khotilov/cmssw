#!/bin/sh

source /nfshome0/cmssw2/scripts/setup.sh
eval `scramv1 run -sh`
export TNS_ADMIN=/nfshome0/popcondev/conddb

xflag=0
oflag=0
while getopts 'xoh' OPTION
  do
  case $OPTION in
      x) xflag=1
          ;;
      o) oflag=1
          ;;
      h) echo "Usage: [-x] [-o] tsckey tagbase"
	  echo "  -x: write to ORCON instead of sqlite file"
	  echo "  -o: overwrite keys"
	  exit
	  ;;
  esac
done
shift $(($OPTIND - 1))

tsckey=$1
tagbase=$2
shift 2

if [ ${oflag} -eq 1 ]
    then
    overwrite="overwriteKeys=1"
fi

if [ ${xflag} -eq 0 ]
    then
    echo "Writing to sqlite_file:l1config.db instead of ORCON."
    cmsRun $CMSSW_BASE/src/CondTools/L1Trigger/test/L1ConfigWritePayloadOnline_cfg.py tscKey=${tsckey} tagBase=${tagbase}_hlt outputDBConnect=sqlite_file:l1config.db ${overwrite} outputDBAuth=. print
    o2ocode=$?
    if [ ${o2ocode} -eq 0 ]
	then
	echo
	echo "`date` : checking O2O"
	if cmsRun $CMSSW_BASE/src/CondTools/L1Trigger/test/l1o2otestanalyzer_cfg.py tagBase=${tagbase}_hlt inputDBConnect=sqlite_file:l1config.db inputDBAuth=. printL1TriggerKeyList=1 | grep ${tsckey} ; then echo "L1TRIGGERKEY WRITTEN SUCCESSFULLY"
	else
	    echo "L1-O2O-ERROR: L1TRIGGERKEY WRITING FAILED" >&2
	    exit 199
	fi
    else
	if [ ${o2ocode} -eq 66 ]
	    then
	    echo "L1-O2O-ERROR: unable to connect to OMDS or ORCON.  Check that /nfshome0/centraltspro/secure/authentication.xml is up to date (OMDS)." >&2
        else
            if [ ${o2ocode} -eq 65 ]
                then
                echo "L1-O2O-ERROR: problem writing object to ORCON." >&2
            fi
        fi
	exit ${o2ocode}
    fi
else
    echo "Writing to cms_orcon_prod."
    cmsRun $CMSSW_BASE/src/CondTools/L1Trigger/test/L1ConfigWritePayloadOnline_cfg.py tscKey=${tsckey} tagBase=${tagbase}_hlt outputDBConnect=oracle://cms_orcon_prod/CMS_COND_31X_L1T outputDBAuth=/nfshome0/popcondev/conddb ${overwrite} print
    o2ocode=$?
    if [ ${o2ocode} -eq 0 ]
	then
	echo
	echo "`date` : checking O2O"
	if cmsRun $CMSSW_BASE/src/CondTools/L1Trigger/test/l1o2otestanalyzer_cfg.py tagBase=${tagbase}_hlt inputDBConnect=oracle://cms_orcon_prod/CMS_COND_31X_L1T inputDBAuth=/nfshome0/popcondev/conddb printL1TriggerKeyList=1 | grep ${tsckey} ; then echo "L1TRIGGERKEY WRITTEN SUCCESSFULLY"
	else
	    echo "L1-O2O-ERROR: L1TRIGGERKEY WRITING FAILED" >&2
	    exit 199
	fi
    else
	if [ ${o2ocode} -eq 66 ]
	    then
	    echo "L1-O2O-ERROR: unable to connect to OMDS or ORCON.  Check that /nfshome0/centraltspro/secure/authentication.xml is up to date (OMDS)." >&2
        else
            if [ ${o2ocode} -eq 65 ]
                then
                echo "L1-O2O-ERROR: problem writing object to ORCON." >&2
            fi
        fi
	exit ${o2ocode}
    fi
fi
exit
