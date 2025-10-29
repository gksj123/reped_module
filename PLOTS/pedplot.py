#-*-Python-*-
# Created by likai at 29 Jun 2020  14:49

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
import numpy as np
from pylab import *
import matplotlib.pyplot as plt
from matplotlib.pyplot import plot,savefig
from scipy import optimize
from sympy import *
import time

for jj in range(100000,100001,1):
    varpeddir = str(jj).zfill(6)
    os.chdir(varpeddir)
    psibar1 = []
    ne1 = []
    te1 = []
    pressure1 = []
    jpar1 = []
    with open('ntpj.txt','r') as ntpjfile:
        for i in ntpjfile.readlines()[0:]:
            (tmp) = i.split()
            psibar1.append(tmp[0])
            ne1.append(tmp[1])
            te1.append(tmp[2])
            pressure1.append(tmp[3])
            jpar1.append(tmp[4])
    ntpjfile.close()
    psibar = map(eval, psibar1)
    ne = map(eval, ne1)
    te = map(eval, te1)
    pressure = map(eval, pressure1)
    jpar = map(eval, jpar1)  # MA/m^2
    plt.rcParams['xtick.direction'] = 'in'
    plt.rcParams['ytick.direction'] = 'in'
    plt.ion()
    plt.figure(1)
    ax = plt.subplot(1,1,1)
    ax.spines['bottom'].set_linewidth(1.5)
    ax.spines['left'].set_linewidth(1.5)
    ax.spines['top'].set_linewidth(1.5)
    ax.spines['right'].set_linewidth(1.5)
    plt.plot(psibar, ne, linestyle='-', color='red', linewidth=3, label='ne')
    xlim(0., 1)
    ylim(0.)
    xlabel('$\psi$',fontsize=16)
    ylabel('n$_{e}$ (10$^{19}$/m$^{3}$)',fontsize=16)
    plt.xticks(fontsize = 15)
    plt.yticks(fontsize = 15)
    plt.legend(loc='best',fontsize=20)
    plt.savefig('./ne.png')
#   plt.savefig('./ne.eps')
    time.sleep(5)
    plt.close(1)
    plt.figure(2)
    ax = plt.subplot(1,1,1)
    ax.spines['bottom'].set_linewidth(1.5)
    ax.spines['left'].set_linewidth(1.5)
    ax.spines['top'].set_linewidth(1.5)
    ax.spines['right'].set_linewidth(1.5)
    plt.plot(psibar, te, linestyle='-', color='blue', linewidth=3, label='Te')
    xlim(0., 1)
    ylim(0.)
    xlabel('$\psi$',fontsize=16)
    ylabel('T$_{e}$ (KeV)',fontsize=16)
    plt.xticks(fontsize = 15)
    plt.yticks(fontsize = 15)
    plt.legend(loc='best',fontsize=20)
    plt.savefig('./te.png')
#   plt.savefig('./te.eps')
    time.sleep(5)
    plt.close(2)
    plt.figure(3)
    ax = plt.subplot(1,1,1)
    ax.spines['bottom'].set_linewidth(1.5)
    ax.spines['left'].set_linewidth(1.5)
    ax.spines['top'].set_linewidth(1.5)
    ax.spines['right'].set_linewidth(1.5)
    plt.plot(psibar, pressure, linestyle='-', color='red', linewidth=3, label='Pressure')
    xlim(0., 1)
    ylim(0.)
    xlabel('$\psi$',fontsize=16)
    ylabel('Pressure (kPa)',fontsize=16)
    plt.xticks(fontsize = 15)
    plt.yticks(fontsize = 15)
    plt.legend(loc='best',fontsize=20)
    plt.savefig('./pressure.png')
#   plt.savefig('./pressure.eps')
    time.sleep(5)
    plt.close(3)
    plt.figure(4)
    ax = plt.subplot(1,1,1)
    ax.spines['bottom'].set_linewidth(1.5)
    ax.spines['left'].set_linewidth(1.5)
    ax.spines['top'].set_linewidth(1.5)
    ax.spines['right'].set_linewidth(1.5)
    plt.plot(psibar, jpar, linestyle='-', color='blue', linewidth=3, label='Current density')
    xlim(0.26, 1)
    ylim(0.)
    xlabel('$\psi$',fontsize=16)
    ylabel('Current density (MA/m$^{2}$)',fontsize=16)
    plt.xticks(fontsize = 15)
    plt.yticks(fontsize = 15)
    plt.legend(loc='best',fontsize=20)
    plt.savefig('./current.png')
#   plt.savefig('./current.eps')
    time.sleep(5)
    plt.close(4)
    os.chdir('../')
