#!/bin/bash

function die { echo Failure $1: status $2 ; exit $2 ; }

echo LOCAL_TMP_DIR = ${LOCAL_TMP_DIR}

pushd ${LOCAL_TMP_DIR}
  testConfig=create_run_only_file_cfg.py
  rm -f dqm_run_only.root
  echo ${testConfig} ------------------------------------------------------------
  cmsRun -p ${LOCAL_TEST_DIR}/${testConfig} || die "cmsRun ${testConfig}" $?

  testConfig=read_run_only_file_cfg.py
  echo ${testConfig} ------------------------------------------------------------
  cmsRun -p ${LOCAL_TEST_DIR}/${testConfig} || die "cmsRun ${testConfig}" $?

  checkFile=check_run_only_file.py
  echo ${checkFile} ------------------------------------------------------------
  python ${LOCAL_TEST_DIR}/${checkFile} || die "python ${checkFile}" $?

  testConfig=create_lumi_only_file_cfg.py
  rm -f dqm_lumi_only.root
  echo ${testConfig} ------------------------------------------------------------
  cmsRun -p ${LOCAL_TEST_DIR}/${testConfig} || die "cmsRun ${testConfig}" $?

  testConfig=read_lumi_only_file_cfg.py
  echo ${testConfig} ------------------------------------------------------------
  cmsRun -p ${LOCAL_TEST_DIR}/${testConfig} || die "cmsRun ${testConfig}" $?

  checkFile=check_lumi_only_file.py
  echo ${checkFile} ------------------------------------------------------------
  python ${LOCAL_TEST_DIR}/${checkFile} || die "python ${checkFile}" $?

  testConfig=create_run_lumi_file_cfg.py
  rm -f dqm_run_lumi.root
  echo ${testConfig} ------------------------------------------------------------
  cmsRun -p ${LOCAL_TEST_DIR}/${testConfig} || die "cmsRun ${testConfig}" $?

  testConfig=read_run_lumi_file_cfg.py
  echo ${testConfig} ------------------------------------------------------------
  cmsRun -p ${LOCAL_TEST_DIR}/${testConfig} || die "cmsRun ${testConfig}" $?

  checkFile=check_run_lumi_file.py
  echo ${checkFile} ------------------------------------------------------------
  python ${LOCAL_TEST_DIR}/${checkFile} dqm_run_lumi.root || die "python ${checkFile}" $?

  #read write
  testConfig=read_write_run_lumi_file_cfg.py
  rm -f dqm_run_lumi_copy.root
  echo ${testConfig} ------------------------------------------------------------
  cmsRun -p ${LOCAL_TEST_DIR}/${testConfig} || die "cmsRun ${testConfig}" $?

  checkFile=check_run_lumi_file.py
  echo ${checkFile} ------------------------------------------------------------
  python ${LOCAL_TEST_DIR}/${checkFile} dqm_run_lumi_copy.root || die "python ${checkFile}" $?

  #merging
  testConfig=create_file1_cfg.py
  rm -f dqm_file1.root
  echo ${testConfig} ------------------------------------------------------------
  cmsRun -p ${LOCAL_TEST_DIR}/${testConfig} || die "cmsRun ${testConfig}" $?

  testConfig=create_file2_cfg.py
  rm -f dqm_file2.root
  echo ${testConfig} ------------------------------------------------------------
  cmsRun -p ${LOCAL_TEST_DIR}/${testConfig} || die "cmsRun ${testConfig}" $?

  testConfig=read_file1_file2_cfg.py
  echo ${testConfig} ------------------------------------------------------------
  cmsRun -p ${LOCAL_TEST_DIR}/${testConfig} || die "cmsRun ${testConfig}" $?

  testConfig=create_file3_cfg.py
  rm -f dqm_file3.root
  echo ${testConfig} ------------------------------------------------------------
  cmsRun -p ${LOCAL_TEST_DIR}/${testConfig} || die "cmsRun ${testConfig}" $?

  testConfig=read_file1_file3_cfg.py
  echo ${testConfig} ------------------------------------------------------------
  cmsRun -p ${LOCAL_TEST_DIR}/${testConfig} || die "cmsRun ${testConfig}" $?

  testConfig=merge_file1_file2_cfg.py
  rm -f dqm_merged_file1_file2.root
  echo ${testConfig} ------------------------------------------------------------
  cmsRun -p ${LOCAL_TEST_DIR}/${testConfig} || die "cmsRun ${testConfig}" $?

  checkFile=check_merged_file1_file2.py
  echo ${checkFile} ------------------------------------------------------------
  python ${LOCAL_TEST_DIR}/${checkFile} || die "python ${checkFile}" $?

  testConfig=merge_file1_file3_file2_cfg.py
  rm -f dqm_merged_file1_file3_file2.root
  echo ${testConfig} ------------------------------------------------------------
  cmsRun -p ${LOCAL_TEST_DIR}/${testConfig} || die "cmsRun ${testConfig}" $?

  testConfig=read_merged_file1_file3_file2_cfg.py
  echo ${testConfig} ------------------------------------------------------------
  cmsRun -p ${LOCAL_TEST_DIR}/${testConfig} || die "cmsRun ${testConfig}" $?


popd

exit 0
