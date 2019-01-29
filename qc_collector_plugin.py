#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
QC Collector plugin script.  Will pull a set of QC metrics from an Ion Torrent
run, as well as a set of associated VCF data, and generate a CSV file that can 
be imported into a QC tracking database.
"""
import sys
import os
import re
import json
import csv
import inspect
import argparse
import subprocess

from pprint import pprint as pp

# Global vars
loglevel = 'debug'
logfile = sys.stderr
plugin_params = {}
plugin_results = {}

def get_plugin_config():
    """
    Get config elems from startplugin.json and barcodes.json. Some of these data
    will be put into the results dict directly for later output.
    """
    global plugin_params
    global plugin_results

    parser = argparse.ArgumentParser(description = __doc__)
    parser.add_argument('start_plugin_file', metavar='startplugin.json')
    parser.add_argument('barcodes_file', metavar='barcodes.json')
    parser.add_argument('-V', '--version', dest='version')
    args = parser.parse_args()

    plugin_params['version'] = args.version

    start_plugin_json = read_json(args.start_plugin_file)

    plugin_params['plan_name'] = start_plugin_json['plan'].get('planName', None)
    instrument = start_plugin_json['expmeta'].get('instrument', None)

    plugin_params['instrument'] = plugin_results['instrument'] = instrument

    plugin_params['results_name'] = start_plugin_json['expmeta'].get(
        'results_name', None)
    plugin_params['analysis_dir'] = start_plugin_json['runinfo'].get(
        'analysis_dir', None)
    plugin_params['plugin_results'] = start_plugin_json['runinfo'].get(
        'results_dir', None)
    plugin_params['pluginconfig']['ir_ip'] = start_plugin_json['runinfo']['plugin']['pluginconfig'].get('ip_address', None)
    plugin_params['pluginconfig']['ir_token'] = start_plugin_json['runinfo']['plugin']['pluginconfig'].get('api_token', None)
    plugin_params['plugin_root'] = start_plugin_json['runinfo']['plugin'].get(
        'plugin_dir', None)

    # Get some results from some of the data collected so far.
    regex = re.compile('Auto.*?%s-([0-9]{3})-%s_[0-9]+$' %
        (plugin_params['instrument'], plugin_params['plan_name']))
    plugin_results['chip_data']['run_num'] = regex.search(plugin_params['results_name']).group(1)

    plugin_results['chip_data']['analysis_date'] = start_plugin_json['expmeta'].get(
        'analysis_date', None)

    plugin_results['chip_data']['pluginconfig'] = start_plugin_json['pluginconfig']


    pp(plugin_params)
    print()
    pp(plugin_results)
    sys.exit()
    '''
    TODO: Here would be a good place to proc the barcodes file to start filling
    out those data.
    '''

def get_chip_metrics(basecaller_file, serialized_json):
    """
    Read the BaseCaller.json and serialized_<expt>.json files, and extract
    the pertinent info.
    """

    # Read the important data from serialized JSON first. Will need later.
    ser_data = read_json(serialized_json)
    for block in ser_data:
        if block['model'] == 'rundb.analysismetrics':
            plugin_results['unfiltered_lib_reads'] = block['fields']['lib']
            plugin_results['filtered_lib_reads']   = block['fields']['libFinal']
            plugin_results['bead_load_pct']        = block['fields']['loading']
        elif block['model'] == 'rundb.libmetrics':
            plugin_results['key_signal']    = block['fields']['aveKeyCounts']
            plugin_results['mean_raw_acc']  = block['fields']['raw_accuracy']

    # Now get the metrics from the BC data file.
    bc_data = read_json(basecaller_file)['Filtering']['LibraryReport']
    for elem in bc_data:
        plugin_results[elem] = '%0.1f' % float(
            int(bc_data[elem]) / int(plugin_results['unfiltered_lib_reads']) * 100
        )

    # Rename the "final_library_reads" to usable reads to match the final
    # output, and to better represent what we have here.  We have the 
    # actual final library reads number elsewhere.
    plugin_results['usable_lib_pct'] = plugin_results.pop('final_library_reads')
    
def get_sample_stats(sample_list):
    """
    Get the basic sample level stats.  Will add more data from the VCFs later.

    TODO: Need to get paths to the <sample>_rawlib.ionstats_alignment.json files.
    """

    results = {}

    # Proces all samples in the list. Store DNA and RNA separately.
    for sample in sample_list:
        ionstats_file = file_dir + '/%s_rawlib.ionstats_alignment.json' % sample
        if os.path.exists(ionstats_file):
            results[sample_list[sample]] = read_ionstats_file(ionstats_file, 'sample')
        else:
            writelog('error', 'No ionstats_alignment file for %s!\n' % sample_list[sample])
            sys.exit(1)

    # Add in the test fragment data
    results['tf_1'] = read_ionstats_file('%s/TFStats.json' % file_dir, 'tf')

    # TODO: Do we want to update the main dict?
    return results

def read_ionstats_file(ifile, sample_type):
    jdata = read_json(ifile)
    if sample_type == 'sample':
        return {'mean_read_length' : jdata['full']['mean_read_length']}
    elif sample_type == 'tf':
        return {'TF_1_50AQ17' : jdata['TF_1']['Percent 50Q17']}

def get_vcf_metrics(ip, token):
    """
    Get the corresponding VCF files for each sample, and from them, extract
    the QC metric data we need.  Will require `ir_api_retrieve.py`, 
    `extract_ir_data.sh`, and `get_metrics_from_vcf.py` helper scripts.
    """
    vcfs = get_vcfs(ip, token)
    results = get_metrics_from_vcf(vcfs, dna_only=False)
    return results

def get_vcfs(ip, token, analysis_date):

    """
    Use ir_api_retrieve to get a set of VCF files that can be used to generate
    QC metrics.

    TODO:
        Get analysis date
        Get paths of helper scripts.
    """
    cmd = ['ir_api_retrieve.py', '--ip', ip, '--token', token, '--date_range',
            analysis_date]
    writelog('info', "Retrieving VCF files for date %s...\n" % analysis_date)
    p = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    stdout, stderr = p.communicate()

    if p.returncode != 0:
        writelog('error', "There was a problem running `ir_api_retrieve.py`. "
            "Can not continue!\n")
        writelog(stderr)
        sys.exit(p.returncode)
    else:
        sys.stderr.write('Retrieved VCFs successfully. Getting RNA VCFs...\n')
        # Use the ir_api_retrieve companion script `extract_ir_data.sh` to unzip
        # and collect vcfs.  
        p = subprocess.run('extract_ir_data.sh', stderr=subprocess.DEVNULL,
            stdout=subprocess.DEVNULL)

        # TODO: Fix the path here.
        vcfdir = os.path.join(os.getcwd(), 'vcfs')
        return [os.path.join(vcfdir, x) for x in os.listdir(vcfdir)]

def get_metrics_from_vcf(vcfs, dna_only):

    """
    Use get_metrics_from_vcf.py to obtain the QC metrics from a VCF. Return a 
    list (or tuple, or maybe even dict if it makes sense) of metrics.

    # TODO:
        get paths of helper scripts.
    """
    cmd = ['get_metrics_from_vcf.py'] 
    if dna_only:
        cmd.append('--dna_only')
    cmd.extend(vcfs)

    p = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    stdout, stderr = p.communicate()

    if p.returncode != 0:
        writelog('error', "There was a problem running `get_metrics_from_vcf.py`"
            ". Can not continue!\n")
        writelog(stderr.decode('ascii'))
        sys.exit(p.returncode)
    
    outdata = stdout.decode('ascii').split('\n')
    header = outdata[0].split()

    results = defaultdict(dict)
    for line in outdata[1:]:
        fields = line.split()
        # We can get a blank line at the end; skip it!
        if (len(fields) < 1):
            continue
        for i,e in enumerate(header):
            results[fields[0]].update({e : fields[i]})
    return results

def read_json(jfile):
    with open(jfile) as fh:
        return json.load(fh)

def writelog(flag, msg):
    """
    Simple level based logging function. 
    """
    now = datetime.datetime.now().strftime('%c')
    log_levels = {
        'info'  : (0, "INFO:"),
        'warn'  : (1, 'WARN:'),
        'error' : (2, 'ERROR:'),
        'debug' : (3, "DEBUG:"),
    }
    log_threshold = log_levels[loglevel][0] 

    if flag is not None:
        flag = flag.lower() # Make sure it's lowercase to avoid key error
        tier, annot = log_levels.get(level, None)
        annot = 'something wrong' if annot == None else annot
        logstr = '{:26s} {:6s} {}\n'.format(now, annot, msg)

        if tier <= log_threshold:
            logfile.write(logstr)
            logfile.flush()

def __exit__(msg=None):
    sys.stderr.write('\n\033[38;5;196mExited at line: {}, with message: '
        '{}\033[00m\n'.format(inspect.stack()[1][2], msg))
    sys.exit()

def plugin_main():
    get_plugin_config()
    # pp(plugin_params)

    '''
    TODO: Fix this for plugin
    basecaller_json = os.path.join(plugin_params['analysis_dir'], 
        'basecaller_results', 'BaseCaller.json')
    serialized_json = os.path.join(plugin_params['analysis_dir'], 
        'serialized_%s.json' % plugin_params['results_name'])
    '''
    basecaller_json = os.path.join(os.path.dirname(__file__), 'resource_files',
        'BaseCaller.json')
    serialized_json = os.path.join(os.path.dirname(__file__), 'resource_files',
        'serialized_%s.json' % plugin_params['results_name'])

    get_chip_metrics(basecaller_json, serialized_json)
    pp(plugin_results)
    __exit__()

    '''
    get_sample_stats()
    get_vcf_metrics()
    '''

if __name__ == "__main__":
    exit(plugin_main())
