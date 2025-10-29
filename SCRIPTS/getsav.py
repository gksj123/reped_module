#-*-Python-*-
# Created by likai at 29 Jun 2020  10:35

"""
This script <fill in purpose>

defaultVars parameters
----------------------
:param kw1: kw1 can be passed to this script as <path to script>.run(kw1='hello')
"""

gEQDSK=os.path.basename(str(root['INPUTS']['gEQDSK']))

basename=gEQDSK[1:]
    ## 'g'+'*.*'

gfile = re.sub('\.0*','_',basename)+'_inv.sav'



bas_script=OMFITascii('bas_script')
open(bas_script.filename,'w').write(
"""


lsrf = 0; msrf=201; map=201; nht=60
run
saveq("cfetr_inv.sav")

integer polline = basopen("polline.txt","w")
integer l
do l = 1, msrf
    polline << format(alsrf(l),12,5,1)
enddo
basclose(polline)
quit
""")







inputs = [root['OUTPUTS']['SAV'][gfile],
          bas_script,
	  root['INPUTS']['reped_input']]

outputs = ['polline.txt',
           'CORSICA_log.txt',
           re.sub('\.0*','_',basename)+'.sav',re.sub('\.0*','_',basename)+'_inv.sav',
           'cfetr_inv.sav']

OMFITx.executable(root,inputs=inputs,outputs=outputs,
                  executable=str(root['SETTINGS']['SETUP']['executable'])+' '+'*_inv.sav '+os.path.split(bas_script.filename)[1]+' | tee CORSICA_log.txt')

root.addBranchPath("['OUTPUTS']['SAV']")
root['OUTPUTS']['SAV']['CORSICA_log']=OMFITascii('CORSICA_log.txt')

for file in set(outputs).intersection(set(glob.glob('*'))):
    if os.path.splitext(file)[1]=='.sav':
        root['OUTPUTS']['SAV'][file]=OMFITpdb(file)
    else:
        root['OUTPUTS']['SAV'][file]=OMFITasciitable(file)

pol = []
with open('polline.txt','r') as polfile:
    for i in polfile.readlines()[0:]:
        (tmpp) = i.split()
        pol.append(tmpp[0])
polfile.close()
polline = float(pol[200])/100
root['OUTPUTS']['reped_output']['polline']=polline
root['OUTPUTS']['reped_output']['betan']=root['INPUTS']['aEQDSK']['betan']
root['OUTPUTS']['reped_output']['Ip']=root['INPUTS']['gEQDSK']['CURRENT']/1.e6
root['OUTPUTS']['reped_output']['a']=root['INPUTS']['aEQDSK']['aout']/100
root['OUTPUTS']['reped_output']['neavg']=root['INPUTS']['reped_input']['ne_avg']/1.e19
root['OUTPUTS']['reped_output']['nepeak']=root['INPUTS']['reped_input']['ne0']/root['INPUTS']['reped_input']['ne_avg']
