#-*-Python-*-
# Created by likai at 29 Jun 2020  14:48

"""
This script <fill in purpose>

defaultVars parameters
----------------------
:param kw1: kw1 can be passed to this script as <path to script>.run(kw1='hello')
"""

import subprocess
import os,sys
import namelist
import os.path
from glob import glob

height=root['INPUTS']['reped_input']['height']

heightt = int(height*500) + 5
for ii in range(heightt):
    varpeddir = str(ii).zfill(6)

    root['OUTPUTS']['TEQ'][varpeddir]['dskgato'] = OMFITascii(root['OUTPUTS']['TEQ'][varpeddir]['tfile'].filename)
    root['OUTPUTS']['TEQ'][varpeddir]['dskgato'].saveas(os.path.dirname(root['OUTPUTS']['TEQ'][varpeddir]['dskgato'].filename)+os.sep+'dskgato')
    inputs = [(root['OUTPUTS']['TEQ'][varpeddir]['dskgato'],'dskgato'),(root['INPUTS']['elite.in'],'elite.in')]

    outputs = ['elite.log']

    executable = "/project/CFETR/bin/pbsMonitor -jn 1 -cn 1 -jq batch -wt 06:00:00 -exe eliteng.x elite 5 30 5"

    OMFITx.executable(root,inputs,outputs,executable=executable)
    root['OUTPUTS']['TEQ'][varpeddir]['elite.log']=OMFITascii('elite.log')
