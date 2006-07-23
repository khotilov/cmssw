<?php
/* 
 * common.php
 * 
 * Common variables and functions for the application (basically a config
 * script)
 * $Id: common.php,v 1.3 2006/07/23 16:47:58 egeland Exp $
 */

function get_conn_params() {
  return array('user' => "cond01",
	       'pass' => "oracond01",
	       'sid'  => "ecalh4db");
}

function get_dqm_url($location, $run) {
  if ($location && $location == 'H4B') {
    $url = "http://pctorino1.cern.ch/html/";
  } else {
    $url = "http://lxcms201.cern.ch/html/";
  }
  return $url.str_pad($run, 9, '0', STR_PAD_LEFT);
}

function get_rootplot_path() {
  
}

function get_cache_dir() {

}

function get_datatype_array() {
  return array('Beam' => 'BEAM',
	       'Monitoring' => 'MON',
	       'DCU' => 'DCU',
	       'Laser' => 'LMF',
	       'DCS' => 'DCS');
}

function get_rootplot_handle($args) {
  putenv('ROOTSYS=/afs/cern.ch/cms/external/lcg/external/root/5.11.02/slc3_ia32_gcc323/root');
  putenv('LD_LIBRARY_PATH=/afs/cern.ch/cms/external/lcg/external/root/5.11.02/slc3_ia32_gcc323/root/lib:/afs/cern.ch/cms/external/lcg/external/Boost/1.33.1/slc3_ia32_gcc323/lib:$LD_LIBRARY_PATH');
  putenv('ROOTPLOT=CMSSW_0_8_0/bin/slc3_ia32_gcc323/cmsecal_rootplot');

  @system('rm rootplot_error.log');
  $handle = popen("\$ROOTPLOT $args > rootplot_error.log 2>&1", "w") or die('Failed to open rootplot program');

  if (! $handle ) {
    return 0;
  }

  flush();
  fflush($handle);
  if (get_rootplot_error()) {
    pclose($handle);
    return 0;
  }
  
  return $handle;
}

function get_rootplot_error() {
  $error_file = @fopen('rootplot_error.log', 'r');
  if (! $error_file) { 
    return 0;
  }

  $error_msg = "";
  while ($line = fgets($error_file)) {
    $error_msg .= $line;
  }
  fclose($error_file);
  return $error_msg;
}

function get_task_array() {
  return array('CI' => 'Channel Integrity Task',
	       'CS' => 'Cosmic Task',
	       'LS' => 'Laser Task',
	       'PD' => 'Pedestal Task',
	       'PO' => 'Pedestals Online Task',
	       'TP' => 'Test Pulse Task',
	       'BM' => 'Beam Task');
}

function get_task_outcome($list_bits, $outcome_bits) {
  if (!$list_bits && !$outcome_bits) { return false; }
  $tasks = get_task_array();
  
  $result = array();
  foreach(array_keys($tasks) as $i => $taskcode) {
    if ($list_bits & (1 << $i)) {
      $result[$taskcode] = $outcome_bits & (1 << $i);
    }
  }
  return $result;
}

function get_stylelinks() {
return "
<link rel='stylesheet' type='text/css' href='ecalconddb.css'/>
<!--[if IE]>
<link rel='stylesheet' type='text/css' href='fixie.css'/>
<![endif]-->
";
}

?>

