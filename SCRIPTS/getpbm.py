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
nesepp=root['INPUTS']['reped_input']['nesep']
tesep=root['INPUTS']['reped_input']['tesep']
zefff=root['INPUTS']['reped_input']['zeffped']
betann=root['OUTPUTS']['reped_output']['betan']
nebarr=root['OUTPUTS']['reped_output']['neavg']
nepeak=root['OUTPUTS']['reped_output']['nepeak']
zionn=1.0
zimpp=10.0
an00=1.0
alphan0=0.9
alphan1=1.8
at11=20
alphat0=1.2
alphat1=1.4
c11=0.7
c22=0.9


nnn = root['SETTINGS']['PHYSICS']['width']  #width
nnnn = root['SETTINGS']['PHYSICS']['height'] #height
pppp = 0
height0 = int(height*100) + 1
height1 = int(height*200) + 2
height2 = int(height*300) + 3
height3 = int(height*400) + 4
height4 = int(height*500) + 5
heightnn = [pppp,height0,height1,height2,height3,height4]
print(heightnn)

for kk in range(5):
    for ii in range(heightnn[kk],heightnn[kk+1],1):
           varpeddir = str(ii).zfill(6)

           at00 = 0   #height
           elitein = root['OUTPUTS']['TEQ'][varpeddir]['elite.log']

           elite = open(str(elitein))
           try:
               elfile = elite.read()
               if 'Stability check finds mode is **UNSTABLE**' in elfile:
                        print('unstable')
                        print(varpeddir)
                        tt = ii-1
                        if tt < heightnn[kk]:
                               break   #break smallest for cycle
                        else:
                               varpeddir = str(tt).zfill(6)
                               teqfile = os.path.basename(str(root['OUTPUTS']['TEQ'][varpeddir]['tfile']))
                               A = []
                               A = teqfile[9:13]
                               at00 = float(A)/1000
                               print(at00)
                               print(tt)
                               widd = nnn[kk]

                               varpedbas1 = OMFITascii('varpedbas1')
                               open(varpedbas1.filename,'w').write(
"""
teq_inv(0,0)
lsrf = 0; msrf=201; map=201; nht=60
teq_inv(0,0)
read jbootstrap.bas
integer gradb = 1
real betan = %f
real nebar = %f
real zion = %f
real zimp = %f
integer ii
real width = %f
real psi_mid = 1.0 - width/2.
real psi_ped = psi_mid - width/2.
real ne_peak = %f
real c_nesep = %f
real an0 = %f
real alpha_n0 = %f
real alpha_n1 = %f
real ne_ped =  an0*( tanh(2.* (1.-psi_mid)/width) - tanh(2.* (psibar-psi_mid)/width))
real ne_sep = c_nesep /(1.-c_nesep) * ne_ped(1)
ne_ped = ne_ped + ne_sep
real an1 = (ne_peak-1.) * ne_ped(1)
real ne = ne_ped
integer peak_ok = 0
while (peak_ok == 0)
    do ii=1,msrf
        if ( psibar(ii) < psi_ped ) then
            ne(ii) = ne_ped(ii) + an1 * (1. - (psibar(ii)/psi_ped)**alpha_n0 )**alpha_n1
        endif
    enddo
    real ne_ave = sum(ne * vpsrf * abs(dpsi00) /msrf) / volume
    if ( abs( ne( int( msrf*0.2*0.2 ) ) /ne_ave - ne_peak ) < 0.01 ) then
        peak_ok = 1
    else
        an1 = an1 + (ne_peak - ne( int( msrf*0.2*0.2 ) )/ne_ave)
    endif
        an1
endwhile
real raxis_o = squeeze( rlsd(1,:,1) )
real raxis_i = squeeze( rlsd(map,:,1) )
real dr_i(msrf-1), dr_o(msrf-1)
do ii=1, msrf-1
    dr_i(ii) = raxis_i(ii+1) - raxis_i(ii)
    dr_o(ii) = raxis_o(ii+1) - raxis_o(ii)
enddo
real ne_bar = sum( ne(1:msrf-1) * (abs(dr_i) + abs(dr_o)) ) / rbore / 2.
ne = ne * (nebar/ne_bar)
ne_ave = ne_ave * (nebar/ne_bar)
ne_ped = ne_ped * (nebar/ne_bar)
real te_sep = %f
real at0 = %f
real te_ped =  te_sep + at0*( tanh(2.* (1.-psi_mid)/width) - tanh(2.* (psibar-psi_mid)/width))
real pres = ne * te_ped * 1602 * 2.
psave = pres * 10
teq_inv(0,0)
real betan_ped = ctroy
betan_ped
real at1 = %f
real alpha_t0 = %f
real alpha_t1 = %f
real te = te_ped
integer betan_ok = 0
while ( betan_ok == 0)
    do ii=1,msrf
        if ( psibar(ii) < psi_ped ) then
            te(ii) = te_ped(ii) + at1 * (1. - (psibar(ii)/psi_ped)**alpha_t0 )**alpha_t1
        endif
    enddo
    pres = ne * te * 1602 * 2.
    psave = pres * 10
    teq_inv(0,0)
    if ( abs(betan-ctroy) < 0.01) then
        betan_ok = 1
    else
        at1 =  (betan-betan_ped)/(ctroy-betan_ped) * at1
    endif
    at1; ctroy
endwhile
real zeff(msrf)
zeff = %f
real te_ave = sum(te * vpsrf * abs(dpsi00) /msrf) / volume
real nu_eff = 0.1 * zeff(1) * ne_ave * rcntr/100. / (te_ave*te_ave)
real pbeam(msrf)
pbeam  = 0.
real tee = te * 1000.
real tii = tee
real xxx = jbootstrap(pres, ne, tee, tii, zeff, pbeam, zion, zimp, gradb)
real jbs = xxx(,1)
real jbs_teqpar = jbs * 1.e-4 * frsrf(1)/rcntr /bsqrf
real jpar0 = jparsave
real placur0 = placur
real placur1 = 0.
real peak0 = max(jbs_teqpar(msrf*0.9 : msrf))
real peak1 = 0.
real c1=%f, c2=%f
while ( abs((placur0-placur1)/placur0) > 0.001 | abs(peak0-peak1)/peak0 > 0.001)
    jparsave = c1 * jpar0 + c2 *jbs_teqpar
    teq_inv(3,0)
    placur
    placur1 = placur
    peak1 = max(jparsave(msrf*0.9 : msrf))
    c1 = c1 - (placur1-placur0)/placur0
    c2 = c2 - (peak1 - peak0)/peak0
endwhile
real press = pres/1000
integer pressure = basopen("pressure.txt","w")
integer l
do l = 1, msrf
    pressure << format(psibar(l),12,5,1) << format(press(l),12,5,1)
enddo
basclose(pressure)
quit
"""%(betann, nebarr, zionn, zimpp, widd, nepeak, nesepp, an00, alphan0, alphan1, tesep, at00,
     at11, alphat0, alphat1, zefff, c11, c22))

                               root['OUTPUTS']['varpedbas1'] = varpedbas1
                               inputs = [(root['OUTPUTS']['SAV']['cfetr_inv.sav'],'cfetr_inv.sav'),(root['INPUTS']['jbootstrap'],'jbootstrap.bas'),(root['OUTPUTS']['varpedbas1'],'varpedbas1')]

                               at00 = round(at00,2)
                               pppp='t'+varpeddir+'.'+str(int(at00*1000)).zfill(5)+'_inv_teq'

                               print(pppp)
                               outputs = [pppp,'pressure.txt']

                               executable = "/project/CFETR/bin/pbsMonitor -jn 1 -cn 1 -jq batch -wt 06:00:00 -exe caltrans cfetr_inv.sav  varpedbas1"
        #executable=str(root['SETTINGS']['SETUP']['executable'])+' '+'east_inv.sav '+os.path.split(varpedbas.filename)[1]+' | tee CORSICA_log.txt'
                               OMFITx.executable(root,inputs,outputs,executable=executable)




                               root['OUTPUTS']['TEQ'][varpeddir]['newtfile']=OMFITascii(pppp)

                               root['OUTPUTS']['TEQ'][varpeddir]['pressure.txt']=OMFITasciitable('pressure.txt')
                               root['OUTPUTS']['TEQ'][varpeddir]['width']=widd
                               root['OUTPUTS']['TEQ'][varpeddir]['height']=at00

           finally:
               elite.close()
           if at00 > 0:

               break
           else:

               continue
