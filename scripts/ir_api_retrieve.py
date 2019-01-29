#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
# 4/1/2015 - D Sims
################################################################################
"""
Starting with a run name from IR or a batch list of run names in a file, grab the
filtered and unfiltered variant call results from the IR API. In addition, we can
now download the RNA BAM file for a run, based on the analysis ID. This is helpful
because the mapped RNA BAM is only found on the IR server, and we often need the
mapped RNA BAM to visualize and verify fusion results.

This script requires an external config file with the API Token in order to access 
the server, which can be generated using the associated config_gen.py script.  
"""
import sys
import os
import io
import re
import argparse
import json
import requests
import zipfile
import datetime
import progressbar

from termcolor import colored,cprint
from pprint import pprint as pp

version = '6.1.012619' 
config_file = os.path.dirname(
        os.path.realpath(__file__)) + '/config/ir_api_retrieve_config.json'


class Config(object):
    def __init__(self,config_file):
        self.config_file = config_file
        self.config_data = Config.read_config(self.config_file)

    def __repr__(self):
        return '%s:%s' % (self.__class__,self.__dict__)

    def __getitem__(self,key):
        return self.config_data[key]

    def __iter__(self):
        return self.config_data.itervalues()

    @classmethod
    def read_config(cls,config_file):
        '''Read in a config file of params to use in this program'''
        try:
            with open(config_file) as fh:
                data = json.load(fh)
        except IOError:
            sys.stderr.write("ERROR: No configuration file found. Do you need "
                "to run the config_gen.py script first?\n")
            sys.exit(1)
        return data


def get_args():
    parser = argparse.ArgumentParser(description = __doc__)
    parser.add_argument(
        'analysis_id', 
        nargs='?', 
        help='Analysis ID to retrieve if not using a batchfile'
    )
    parser.add_argument(
        '-H', '--Host', 
        metavar="<server_hostname>",
        help="Hostname of server from which to gather data. Use '?' to print "
            "out all valid hosts as defined in your config file."
    )
    parser.add_argument(
        '-b','--batch', 
        metavar='<batch_file>',
        help='Batch file of experiment names to retrieve. Batchfile must '
            'contain an analysis ID, one per line.'
    )
    parser.add_argument(
        '-i', '--ip', 
        metavar='<ip_address>',
        help='IP address if not entered into the config file yet.'
    )
    parser.add_argument(
        '-t','--token', 
        metavar='<ir_token>',
        help='API token if not entered into the config file yet.'
    )
    parser.add_argument(
        '-m','--method', 
        metavar='<api_method_call>', 
        default='getvcf', 
        help='Method call / entry point to API.'
    )
    parser.add_argument(
        '-r', '--rna', 
        action='store_true', 
        help='Download the RNA BAM file instead of the VCF data.'
    )
    parser.add_argument(
        '-d', '--dna',
        action = 'store_true',
        help='Download the DNA BAM file instead of the VCF data. This can take a '
            'very long time, and is only recommended if you have no other way to '
            'retrieve a DNA BAM file.  Getting directly from TS is preferred!'
    )
    parser.add_argument(
        '--date_range', 
        metavar='<YYYY-MM-dd,YYYY-MM-dd>', 
        help='Range of dates in the format of "start,end" where each date is in '
        'the format YYYY-MM-dd. This will be the range which will be used to '
        'pull out results. One can input only 1 date if the range is only going '
        'to be one day. Note that this method is a bit slow as every run seems '
        'to be checked for inclusion, as well, as the fact that anything within '
        'the range will be downloaded. So, choose this method only if it is '
        'faster than just simply copy / pasting a discrete list of IDs you want!'
    )
    parser.add_argument(
        '-v', '--version', 
        action='version', 
        version='%(prog)s - v' + version
    )
    cli_args = parser.parse_args()

    if not cli_args.Host:
        if cli_args.ip and not cli_args.token:
            sys.stderr.write("ERROR: You must enter a custom token with the '-t'"
                " option if you are using a custom IP.\n")
            sys.exit(1)
        elif cli_args.token and not cli_args.ip:
            sys.stderr.write("ERROR: You must enter an IP with the '-i' option if "
                "you are using a custom token.\n")
            sys.exit(1)
        elif not cli_args.ip and not cli_args.token:
            sys.stderr.write("ERROR: You must either enter a host name or a custom "
                "IP and token!\n")
            sys.exit(1)
    return cli_args

def __validate_date(date):
    try:
        datetime.datetime.strptime(date, '%Y-%M-%d')
    except ValueError:
        sys.stderr.write("ERROR: the date '%s' is not in a valid format. "
            "You must use YYYY-MM-dd.\n")
        sys.exit(1)

def get_host(hostname, hostdata=None):
    """
    Return the IP address and API token of the server from which we wish to 
    retrieve data
    """
    if hostname == '?':
        sys.stderr.write("Current list of valid list of hosts are: \n")
        for host in hostdata:
            sys.stderr.write("\t{}\n".format(host))
        sys.exit()
    else:
        try:
            ip = hostdata[hostname]['ip']
            token = hostdata[hostname]['token']
        except KeyError:
            sys.stderr.write("ERROR: '{}' is not a valid IR Server "
                "name.\n".format(hostname))
            get_host('?',hostdata)
            sys.exit(1)
    return ip, token

def format_url(ip):
    pieces = ip.lstrip('https://').split('.')
    if len(pieces) != 4: 
        sys.stderr.write("ERROR: the IP address you entered, '{}', does not "
            "appear to be valid!\n".format(ip))
        sys.exit(1)
    if all(0<=int(p)<256 for p in pieces):
        return 'https://{}'.format('.'.join(pieces))
    else:
        sys.stderr.write("ERROR: the IP address you entered, '{}', does not "
            "appear to be valid!\n".format(ip))
        sys.exit(1)

def jdump(json_data):
    print(json.dumps(json_data, indent=4, sort_keys=True))

def proc_batchfile(batchfile):
    with open(batchfile) as fh:
        return [line.rstrip() for line in fh if line != '\n']

def make_bam_datalink(na_type, run_summary, session, header):
    """
    Have to get the DNA or RNA BAM file name, which is going to be stored in an 
    RRS file that contains the sample name. 

    The RRS file will contain the information we need to find the BAM file. 
    In the case of the DNA BAM, we can just follow the path outlined.  In the
    case of the RNA BAM, we need the one from the /outputs/RNACountsActor-00
    dir, but we don't know the barcode without reading the RRS file, and we 
    need to add the '_merged' string to the name.
    """
    data_dir = os.path.dirname(run_summary['data_links']['unfiltered_variants'])
    try:
        ir_sample_name = run_summary['samples'][na_type]
    except KeyError:
        sys.stderr.write("ERROR: No {} sample information for for run: {}! "
            "Was there an RNA run with this set?\n".format(
                na_type, run_summary['name'])
        )
        return None

    rrs_file = ir_sample_name + '.rrs'

    response = session.get(data_dir + '/' + rrs_file, headers=header,
        verify=False)
    z = zipfile.ZipFile(io.BytesIO(response.content))
    data = z.read(rrs_file).decode('ascii')
    elems = data.split()
    
    if na_type == 'RNA':
        rna_bam = os.path.basename(elems[-1]).replace('.bam', '', 1) + '_merged.bam'
        api_path = data_dir + '/outputs/RNACountsActor-00/' + rna_bam 
        return api_path
    elif na_type == 'DNA':
        dna_bam = elems.pop().rstrip('\n')
        api_path = '{}={}'.format(data_dir.split('=')[0], dna_bam)
        return api_path

def api_call(url, query, header, batch_type, get_rna, get_dna, name=None):
    requests.packages.urllib3.disable_warnings()
    s = requests.Session()

    # If we want to get RNA files, we need to get the officially entered RNA 
    # name, which means we need to get a summary call.
    if get_rna or get_dna:
        url = url.replace('getvcf', 'analysis')

    request = s.get(url, headers=header, params=query, verify=False)
    try:
        request.raise_for_status()
    except requests.exceptions.HTTPError as error:
        cprint('\n\n\t{}'.format(error), 'red', attrs=['bold'], file=sys.stderr)
        if batch_type == 'range':
            cprint('\tThere may be no data available for the range input. Check '
                'the date range and try again.\n','red', 
                attrs=['bold'], file=sys.stderr)
        else:
            cprint('\tSkipping analysis id: %s. Check ID for this run and try '
                'again.\n' % query['name'], 'red', attrs=['bold'], file=sys.stderr)
        return None

    json_data = request.json()
    datatype = 'VCF data'

    if batch_type == 'range':
        sys.stdout.write('Done!\n')
        sys.stdout.write('Total number to retrieve: %s.\n' % len(json_data))
        sys.stdout.flush()
        return [x['name'] for x in json_data]

    for analysis_set in json_data:
        if get_rna:
            data_link = make_bam_datalink('RNA', analysis_set, s, header)
            datatype = 'RNA BAM file'
            if not data_link:
                return None
        elif get_dna:
            data_link = make_bam_datalink('DNA', analysis_set, s, header)
            datatype = 'DNA BAM file'
            if not data_link:
                return None
        else:
            data_link = analysis_set['data_links']

        zip_name = name + '_download.zip'
        with open(zip_name, 'wb') as zip_fh:
            response = s.get(data_link, headers=header, verify=False, 
                stream=True)
            
            total_size = response.headers.get('content-length', None)
            prog_bar2(response, total_size, zip_fh)
    sys.stderr.write('Done!\n\n')

def prog_bar2(response, size, fh):
    """
    Using the ProgressBar2 library, generate a progress bar to help determine
    the download speed, time, etc.  Not very useful for VCF files as they come
    down almost instantly.  But, for BAM files, helpful to get an idea of how 
    long it'll take to get here.  For DNA BAM files, we don't know the size, so
    just output the amount downloaded, speed, etc.  For RNA BAM, can get a whole
    set of data.
    """
    wrote = 0

    if size is None:
        # Have a DNA BAM and don't know the actual size.
        widgets = [
            'Downloaded: ',
            progressbar.DataSize(),
            '; (',
            progressbar.FileTransferSpeed(), ' | ', progressbar.Timer(), ')',
        ]
        pbar = progressbar.ProgressBar(
            widgets=widgets, 
            maxval=progressbar.UnknownLength
        ).start()
    else:
        widgets = [
            progressbar.Percentage(),
            progressbar.Bar(marker="=", left=" [", right="] "),
            progressbar.DataSize(), 
            ' of ', 
            progressbar.DataSize(variable='max_value'), 
            " (", progressbar.FileTransferSpeed(), ", ", progressbar.ETA(), " )",
        ]
        size = int(size)
        pbar = progressbar.ProgressBar(widgets=widgets, maxval=size, term_width=100).start()

    for buf in response.iter_content(1024):
        if buf:
            fh.write(buf)
            wrote += len(buf)
            pbar.update(wrote)
    pbar.finish()

def main():
    cli_args = get_args()
    if cli_args.rna:
        datatype = 'RNA BAM file'
    elif cli_args.dna:
        datatype = 'DNA BAM file'
    else:
        datatype = 'VCF data'
    
    server = cli_args.Host if cli_args.Host else cli_args.ip

    if cli_args.ip and cli_args.token:
        server_url = format_url(cli_args.ip) 
        api_token = cli_args.token
    else:
        program_config = Config.read_config(config_file)
        server_url, api_token = get_host(cli_args.Host, program_config['hosts'])

    server_url += '/api/v1/'

    analysis_ids=[]
    if cli_args.batch:
        analysis_ids = proc_batchfile(cli_args.batch)
    elif cli_args.analysis_id:
        analysis_ids.append(cli_args.analysis_id)
    elif not cli_args.date_range:
        sys.stderr.write("ERROR: No analysis ID or batch file loaded!\n")
        sys.exit(1)

    header = {
        'Authorization':api_token,
        'content-type':'application/x-www-form-urlencoded'
    }
    method=cli_args.method
    url = server_url + method
    #print('::DEBUG:: formated base url: {}'.format(url))

    if cli_args.date_range:
        # Allow for one to just put one date to look for data on that date alone
        start, end = (cli_args.date_range.split(',') + [None]*2)[:2]
        if end is None:
            end = start
        __validate_date(start)
        __validate_date(end)

        sys.stdout.write('Getting list of results from IR {} for dates from {} '
            'to {}...'.format(server, start, end))
        sys.stdout.flush()
        query = {
            'format' : 'json', 
            'start_date' : start, 
            'end_date' : end, 
            'exclude' : 'filteredvariants' 
        }

        analysis_ids = api_call(url, query, header, 'range', cli_args.rna, 
            cli_args.dna)
    
    sys.stdout.write('Getting data from IR {} (total runs: {}).\n\n'.format(
        server, len(analysis_ids)))
    sys.stdout.flush()
    count = 0
    for expt in analysis_ids:
        count += 1
        sys.stdout.write('[{}/{}]  Retrieving {} for analysis ID: {}...\n'.format(
            count, len(analysis_ids), datatype, expt)
        )
        sys.stdout.flush()
        query = {
            'format'  : 'json', 
            'name'    : expt, 
            'exclude' : 'filteredvariants'
        }
        api_call(url, query, header, 'single', cli_args.rna, cli_args.dna, expt)

if __name__ == '__main__':
    try:
        main()
    except KeyboardInterrupt:
        sys.exit(1)
