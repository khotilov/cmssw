#!/bin/sh

ls /afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/data/OnlineData/original/00014xxxx/*/DQM_V00*_SiStrip_* > lsafs.dat
awk 'BEGIN { FS = "/" } ; { print $13 }' lsafs.dat > lsfiles.dat
rm lsafs.dat
sed 's/DQM_V0001_SiStrip_R000//' lsfiles.dat > lsruns
rm lsfiles.dat
echo "./extractErrorsvsTime /afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/data/OnlineData/original/00014xxxx/" `sed 's/\.root//' lsruns | wc -l` `sed 's/\.root//' lsruns | awk '{ORS=" "}{print $0}'`
# > lsruns.dat
rm lsruns

