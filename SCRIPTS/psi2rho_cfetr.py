#! /bin/env python
import sys, os
import numpy as np
import tokdat
tokdat.tok = 'CFETR'

if len(sys.argv) !=5:
	print('Usage: psi2rho.py shot time efit_dir datafile.\n \
      efit_dir can be efit_cfetr or the directory of gfile.\n \
      datafile has one column of psi value.')
	exit()

MDS_SERVER = '202.127.204.12'
shot = int(sys.argv[1])
time = float(sys.argv[2])
efit_tree = sys.argv[3]
File = sys.argv[4]

from pmds import mdsconnect
from efit_eqdsk import get_gdat, read_gfiles, eqdsk_to_1t
from profiles_mapping import get_mapping

if efit_tree == 'efitrt_east' or efit_tree =='efit_east':
    time = time/1000.
    mdsconnect(MDS_SERVER)
    # Only the keys for profile mapping
    gdat_keys =  ['gtime','rmaxis','zmaxis','cpasma','ssibry', 'ssimag','psirz',\
          'fpol','rhovn','qpsi','nbdry','bdry','lim']
    gdat = get_gdat(shot, tmin=time, tmax=time, efit=efit_tree, gnames=gdat_keys, open_tree=True)
else:
    gdat = read_gfiles(shot,tmin=time,tmax=time, efdir=efit_tree)
gdat = eqdsk_to_1t(gdat, time)

psi = np.loadtxt(File)
mapping = get_mapping(1,1,efittree='efit_east', efdat=gdat)
rhob = mapping['rhob']
rho = rhob(psi)
#for ii in rho: print(ii)

np.savetxt('rhonorm.txt', rho, fmt="%10.6f")
print('\nNormal termination')
