#-*-Python-*-
# Created by likai at 29 Jun 2020  10:35

"""
This script <fill in purpose>

defaultVars parameters
----------------------
:param kw1: kw1 can be passed to this script as <path to script>.run(kw1='hello')
"""

gEQDSK=os.path.basename(str(root['INPUTS']['gEQDSK']))
print(gEQDSK)
basename=gEQDSK[1:]
print(basename)    ## 'g'+'*.*'

bas_script=OMFITascii('bas_script')
open(bas_script.filename,'w').write(
"""
read cfetr.bas
gfiles
cfetr(gfilelist(1),0)

quit
""")







inputs = [root['INPUTS']['gEQDSK'],
          root['INPUTS']['aEQDSK'],
          bas_script]

outputs = ['CORSICA_log.txt',
           re.sub('\.0*','_',basename)+'.sav',re.sub('\.0*','_',basename)+'_inv.sav']

OMFITx.executable(root,inputs=inputs,outputs=outputs,
                  executable=str(root['SETTINGS']['SETUP']['executable'])+' '+os.path.split(bas_script.filename)[1]+' | tee CORSICA_log.txt')

root.addBranchPath("['OUTPUTS']['SAV']")
root['OUTPUTS']['SAV']['CORSICA_log']=OMFITascii('CORSICA_log.txt')

for file in set(outputs).intersection(set(glob.glob('*'))):
    if os.path.splitext(file)[1]=='.sav':
        root['OUTPUTS']['SAV'][file]=OMFITpdb(file)
    else:
        root['OUTPUTS']['SAV'][file]=OMFITasciitable(file)
