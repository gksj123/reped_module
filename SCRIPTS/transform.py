#-*-Python-*-
# Created by wmq at 26 Mar 2021  23:03

"""
This script <fill in purpose>

defaultVars parameters
----------------------
:param kw1: kw1 can be passed to this script as <path to script>.run(kw1='hello')
"""

import numpy as np


data=np.array(root['OUTPUTS']['TEQ']['tot_result']['prfiles_ntpj']['data']['col1'])


data.transpose()
data=np.transpose(data)
numpy.savetxt("result.txt", data)
root['OUTPUTS']['psi']=OMFITascii("result.txt")

shot = root['INPUTS']['gEQDSK']['CASE'][3][-6:]

time = root['INPUTS']['gEQDSK']['CASE'][4][:-2]
print('python psi2rho_cfetr.py '+shot+' '+time+' ./ result.txt')
run='python psi2rho_cfetr.py '+shot+' '+time+' ./ result.txt'

inputs = [root['INPUTS']['gEQDSK'],root['SCRIPTS']['psirho'],root['OUTPUTS']['psi']]

outputs = ['rhonorm.txt']
executable = 'source /scratch/wmq/myPython/env_set'
executable +='\n'+ run
OMFITx.executable(root,inputs=inputs,outputs=outputs,executable=executable)

root['OUTPUTS']['rho1']=OMFITasciitable('rhonorm.txt')


bb = root['OUTPUTS']['rho1']['data']['0.000000']
aa = [0]


for i in range(200):

     aa.append(bb[i])



root['OUTPUTS']['profiles']=OMFITasciitable(root['OUTPUTS']['TEQ']['tot_result']['prfiles_ntpj'].filename)
root['OUTPUTS']['profiles']['data']['col1']=aa

root['OUTPUTS']['profiles']['header']=str(['rho','ne(*1e13cm^-3)','Te(keV)','Pressure(kPa)','currrent density(MA/m^2)'])

rya=root['OUTPUTS']['profiles']['data']['col1']
rho=root['OUTPUTS']['TEQ']['tot_result']['prfiles_ntpj']['data']['col1']
ne=root['OUTPUTS']['profiles']['data']['col2']
te=root['OUTPUTS']['profiles']['data']['col3']
press=root['OUTPUTS']['profiles']['data']['col4']
cur=root['OUTPUTS']['profiles']['data']['col5']

mask=[i for i in range(201)]

nee=rho*0
tee=rho*0
presss=rho*0
curr=rho*0

nee[mask]=interpolate.interp1d(rya,ne,kind='slinear')(rho[mask])

tee[mask]=interpolate.interp1d(rya,te,kind='slinear')(rho[mask])

presss[mask]=interpolate.interp1d(rya,press,kind='slinear')(rho[mask])

curr[mask]=interpolate.interp1d(rya,cur,kind='slinear')(rho[mask])

root['OUTPUTS']['profiles']['data']['col1']=rho
root['OUTPUTS']['profiles']['data']['col2']=nee
root['OUTPUTS']['profiles']['data']['col3']=tee
root['OUTPUTS']['profiles']['data']['col4']=presss
root['OUTPUTS']['profiles']['data']['col5']=curr
