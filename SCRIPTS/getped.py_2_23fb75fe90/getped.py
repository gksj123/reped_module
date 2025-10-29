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

import sympy as sp
#from sympy.abc import x,y,a,b,c,d
import time

height=root['INPUTS']['reped_input']['height']
nesepp=root['INPUTS']['reped_input']['nesep']
tesep=root['INPUTS']['reped_input']['tesep']
zefff=root['INPUTS']['reped_input']['zeffped']
betann=root['OUTPUTS']['reped_output']['betan']
nebarr=root['OUTPUTS']['reped_output']['neavg']
nepeak=root['OUTPUTS']['reped_output']['nepeak']
Ip=root['OUTPUTS']['reped_output']['Ip']
polline=root['OUTPUTS']['reped_output']['polline']
zionn=1.0
zimpp=6.0
an00=1.0
alphan0=0.9
alphan1=1.8
at11=2.7
alphat0=1.2
alphat1=1.4
c11=0.7
c22=0.9

nnn = root['SETTINGS']['PHYSICS']['width']  #width
nnnn = root['SETTINGS']['PHYSICS']['height'] #height
pppp = 0
height0 = int(height*500) + 5



press1 = []
width = []
Height = []
for ii in range(height0):
    varpeddir = str(ii).zfill(6)



    if "pressure.txt" in root['OUTPUTS']['TEQ'][varpeddir]:    #root['OUTPUTS']['TEQ'][varpeddir].has_key('pressure.txt'):
          width0 = root['OUTPUTS']['TEQ'][varpeddir]['width']
          height0 = root['OUTPUTS']['TEQ'][varpeddir]['height']
          press0 = root['OUTPUTS']['TEQ'][varpeddir]['pressure.txt']['data']['col2'][(200-int(width0/0.005))]

          press1.append(press0)
          width.append(width0)
          Height.append(height0)

print (press1,width,Height)



def f_1(x, A, B):
    return A*x + B

plt.rcParams['xtick.direction'] = 'in'
plt.rcParams['ytick.direction'] = 'in'
wid = width
at = Height
pres = press1
x0 = [n for n in wid if n>0]
y0 = [n for n in pres if n>0]
x000 = [n for n in wid if n>0]
y000 = [n for n in pres if n>0]
plt.figure()
ax = plt.subplot(1,1,1)
ax.spines['bottom'].set_linewidth(1.5)
ax.spines['left'].set_linewidth(1.5)
ax.spines['top'].set_linewidth(1.5)
ax.spines['right'].set_linewidth(1.5)
plt.plot(x0, y0, 'sb', markersize=8)
A1, B1 = optimize.curve_fit(f_1, x000, y000)[0]
print (A1, B1)
x1 = np.arange(0.02, 0.11, 0.01)
y1 = A1*x1 + B1
plt.plot(x1, y1, 'b-', linewidth=3, label='PBM stability boundary')

x00 = [n for n in at if n>0]
y00 = [n for n in pres if n>0]
#plt.plot(x0, y0, 'b-o', linewidth=3)
A2, B2 = optimize.curve_fit(f_1, x00, y00)[0]
print (A2, B2)
x11 = np.arange(0.02, 0.11, 0.01)
y11 = A2*x11 + B2

cc = 0.1 #0.07-0.1
Ipp = Ip*1.e6 # plasma current, in Ampere
ll = polline # length of the poloidal perimeter of the plasma, in m
mu0= 4.* 3.14159 * 1.e-7
x2 = np.linspace(0., 0.08, 51)
betap_ped = (x2 / cc) **2
y2 = 0.5 * 1.e-3 * betap_ped * mu0 * Ipp**2 / ll**2
plt.plot(x2, y2, linestyle='--', color='green', linewidth=3.5, label='KBM stability boundary')

#--------------------------------before is OK---------------------------------------------


x=sp.Symbol('x')
y=sp.Symbol('y')





kkk=[-A1*x-B1+y, - 0.5 * 1.e-3 * ((x / cc) **2) * mu0 * Ipp**2 / ll**2+y]
tt = sp.solve(kkk,[x,y])
width = tt[1][0]
hight = tt[1][1]
pedestal = [width,hight,0,0,0]
print ("width=",width,"hight=",hight)
plt.plot(width,hight,'ro',markersize=9)

xlim(0.02, 0.08)
ylim(0, 120)
xlabel('Pedestal width ($\psi$)',fontsize=16)
ylabel('Pedestal height (p$_{ped}$, kPa)',fontsize=16)
plt.text(0.03,5.1,'[width=%f,height=%f]'%(width,hight),fontsize=16)
plt.xticks(fontsize = 12)
plt.yticks(fontsize = 12)
plt.legend(loc='best',fontsize=15)
#plt.savefig('./eped.png')
#plt.savefig('./eped.eps')
plt.show()

at000 = 0
at000 = sp.solve(A2*x+B2-hight,x)
at0000 = at000[0]
print (at0000)
