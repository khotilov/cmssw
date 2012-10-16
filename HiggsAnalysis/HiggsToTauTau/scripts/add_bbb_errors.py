#!/usr/bin/env python
'''

Add the bin-by-bin uncertainties to the H2Tau shape files.

For usage, see --help

Example:

    add_bbb_errors.py mt,et,em:7TeV,8TeV:01,03,05:QCD,W

will create a setup_bbb directory where the

    mt,et,em

channels, in 7,8 TeV, for the 1, 3, and 5 categories, have bbb systematics added
for the QCD and W systematics.

Author: Evan K. Friis, UW Madison

'''

from RecoLuminosity.LumiDB import argparse
import logging
import os
import shutil
import subprocess
import sys

def get_channel_name(finalstate, category):
    ''' Turn 'mt' + 00 -> muTau_0jet_low '''
    fs_map = {
        'mt' : 'muTau',
        'et' : 'eleTau',
        'em' : 'eMu',
        'tt' : 'tauTau',
    }
    cat_map = {
        '00' : '0jet_low',
        '01' : '0jet_high',
        '02' : 'boost_low',
        '03' : 'boost_high',
        '05' : 'vbf',
    }
    tt_cat_map = {
        '00' : 'boost',
        '01' : 'vbf',
    }
    if finalstate == 'tt':
        return '_'.join((fs_map[finalstate], tt_cat_map[category]))
    else:
        return '_'.join((fs_map[finalstate], cat_map[category]))


def get_shape_file(sourcedir, channel, period):
    # Ex: shape file for mt lives in
    #  setup/mt/htt_mt.inputs-sm-7TeV.root
    return os.path.join(sourcedir, channel,
                        'htt_' + channel + '.inputs-sm-' + period + '.root')

def get_card_config_files(sourcedir, channel, period, category):
    ''' Get the configuration files (cgs, unc.vals, etc)

    Returns a tuple of paths to the cgs, unc.conf, and unc.vals

    '''
    base_dir = os.path.join(sourcedir, channel)
    cgs = 'cgs-sm-%s-%s.conf' % (period, category)
    unc_c = 'unc-sm-%s-%s.conf' % (period, category)
    unc_v = 'unc-sm-%s-%s.vals' % (period, category)
    return (
        os.path.join(base_dir, cgs),
        os.path.join(base_dir, unc_c),
        os.path.join(base_dir, unc_v),
    )

def add_systematics(cat_name, process, systematics, unc_conf_file, unc_val_file):
    ''' Add the shape systematics in <systematics> to the unc. files

    <cat_name> specifies the category nice name - like muTau_0jet_high
    <unc_conf_file> and <unc_val_file> should be open file handles.

    '''
    # Write to the unc.conf file
    for systematic_name in systematics:
        unc_conf_file.write('%s shape\n' % systematic_name)
        unc_val_file.write(
            '%s %s %s 1.00\n' % (cat_name, process, systematic_name))

def create_systematics(channel, category, process, shape_file, threshold):
    ''' Create the bin-by-bin systematics in the shape file.

    Returns a tuple with (channel name, list of added systs)

    '''

    channel_name = get_channel_name(channel, category)

    root_path = os.path.join('/',channel_name, process)

    command = [
        'add_stat_shapes.py',
        shape_file, # input
        shape_file, # output (modded in place)
        '--filter',
        root_path, # the histogram we are bbb-ing
        '--prefix',
        # Make the prefix as short as possible, to avoid RooFit
        # string-to-long bug
        # This prefix is needed so the systematics names don't overlap
        # between ET/MT, boost/VBF, etc.
        'CMS_htt_%s%i' % (channel, int(category)),
        '--threshold',
        str(threshold),
    ]
    log.debug("Shape command:")
    log.debug(" ".join(command))
    # Run the command, get the list of new names (written to stdout)
    stdout = subprocess.Popen(
        command, stdout=subprocess.PIPE, stderr=sys.stderr).communicate()[0]
    added_systematics = []
    for line in stdout.split('\n'):
        if line and 'CMS_htt' in line:
            added_systematics.append(line.strip())
    return channel_name, added_systematics


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('commands', nargs='+',
                        help='Syntax: CHAN1,CHAN2:PERIDOD1,PER1:CAT0,CAT1:PROC1,PROC2')
    parser.add_argument('--threshold', type=float, default=0.05,
                        help='Minimum error for systematic creation,'
                        'default %(default)0.2f')

    parser.add_argument('--output-dir', dest='outputdir', default='setup_bbb',
                        help='Output directory for modified cards. Default: setup_bbb')

    parser.add_argument('--input-dir', dest='inputdir', default='setup',
                        help='Output directory for modified cars. Default: setup')

    parser.add_argument('-f', dest='force', action='store_true',
                        help='Force creation of new output dir')

    args = parser.parse_args()

    log = logging.getLogger('bin-by-bin')
    logging.basicConfig(stream=sys.stderr, level=logging.DEBUG)

    log.info("Parsing commands")

    # Flatten into a list of channels, periods, categories, and processes to do
    all_commands = []

    def bad_command(command):
        log.error("Could not parse command: %s", command)
        log.error("See syntax in --help")
        sys.exit(2)

    for command in args.commands:
        fields = command.split(':')
        if len(fields) != 4:
            bad_command(command)
        for channel in fields[0].split(','):
            for period in fields[1].split(','):
                for cat in fields[2].split(','):
                    for proc in fields[3].split(','):
                        all_commands.append((channel, period, cat, proc))

    log.info("Found %i shapes to mangle", len(all_commands))

    log.info("Copying initial cards from %s ==> %s", args.inputdir,
             args.outputdir)

    if os.path.exists(args.outputdir):
        if args.force:
            shutil.rmtree(args.outputdir)
        else:
            log.error("Output directory already exists. Use -f to override.")
            sys.exit(1)

    shutil.copytree(args.inputdir, args.outputdir)

    total_added_systematics = 0

    for command in all_commands:
        channel, period, cat, proc = command
        log.info("Mangling: %s", ' '.join(command))
        # Create the systematics
        shape_file = get_shape_file(args.outputdir, channel, period)
        nicename, systematics = create_systematics(
            channel, cat, proc, shape_file, args.threshold)
        log.info("Added systs for %i bins", len(systematics))
        total_added_systematics += len(systematics)
        cgs, unc_c, unc_v = get_card_config_files(
            args.outputdir, channel, period, cat)
        log.info("Adding systematics to files")
        with open(unc_c, 'a') as unc_c_file:
            with open(unc_v, 'a') as unc_v_file:
                add_systematics(nicename, proc, systematics,
                                unc_c_file, unc_v_file)

    log.info("Added %i new systematics!", total_added_systematics)
