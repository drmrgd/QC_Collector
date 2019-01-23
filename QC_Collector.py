#!/usr/bin/python3
# -*- coding: utf-8 -*-
import os
import sys
import subprocess

from ion.plugin import *

class QC_Collector(IonPlugin):
    """
    Plugin to collect QC metrics from an Ion Torrent run, and generate a CSV
    file that can be loading into tracking databases downstream.
    """
    version = '0.2.012219'
    major_block = False
    runtypes = [RunType.FULLCHIP, RunType.THUMB, RunType.COMPOSITE]
    runlevels = [RunLevel.DEFAULT]

    def launch(self, data=None):
        sys.stderr.write('Writing output to "data.csv"\n')

        with open('data.csv', 'w') as outfh:
            for var in os.environ:
                outfh.write('{},{}\n'.format(var, os.environ[var]))

        rundate = os.environ.get('TSP_RUN_DATE', None)
        runid = os.environ.get('TSP_RUNID', None)
        sys.stderr.write('run date: %s\n' % rundate)
        sys.stderr.write("run id: %s\n" % runid)

        sys.exit()
    #def launch(self, data=None):
        #cmd = ['qc_collector_plugin.py', '-v', self.version, '-d', 
                #'startplugin.json', 'barcodes.json']
        #plugin = subprocess.Popen(cmd, stdout=subprocess.PIPE, 
                #stderr=subprocess.PIPE)
        #plugin.communicate()
        #sys.exit(plugin.poll())

if __name__ == '__main__':
    PluginCLI()
