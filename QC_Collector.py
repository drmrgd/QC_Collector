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
    version = '0.6.013119'
    major_block = False
    runtypes = [RunType.FULLCHIP, RunType.THUMB, RunType.COMPOSITE]
    runlevels = [RunLevel.DEFAULT]
    depends = ['coverageAnalysis']

    def launch(self, data=None):
        cmd = ['qc_collector_plugin.py', '-V', self.version, 'startplugin.json', 
            'barcodes.json']
        plugin = subprocess.Popen(cmd, stdout=subprocess.PIPE, 
                stderr=subprocess.PIPE)
        plugin.communicate()
        sys.exit(plugin.poll())

if __name__ == '__main__':
    PluginCLI()
