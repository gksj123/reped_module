#-*-Python-*-
# Created by likai at 29 Jun 2020  10:37

"""
This script <for the REPED model>

defaultVars parameters
----------------------
"""

OMFITx.TitleGUI('Run REPED GUI')

OMFITx.FilePicker("scratch['geqdsk']",'gFile',default='')
OMFITx.FilePicker("scratch['aeqdsk']",'aFile',default='')

OMFITx.Entry("root['INPUTS']['reped_input']['ne0']",'ne0 (m^-3)',help='The electron density at rho=0')
OMFITx.Entry("root['INPUTS']['reped_input']['ne_avg']",'ne_avg (m^-3)',help='The line-averaged electron density')
OMFITx.Entry("root['INPUTS']['reped_input']['height']",'Height',help='The parameter to constran the pedestal height,0.05-0.1',default=0.05)
OMFITx.Entry("root['INPUTS']['reped_input']['nesep']",'nesep',help='The ratio of electron density at separatrix to electron density at pedestal top,0.1-0.4',default=0.25)
OMFITx.Entry("root['INPUTS']['reped_input']['tesep']",'Tesep (keV)',help='The electron temperature at the separatrix,0.05-0.1keV',default=0.075)
OMFITx.Entry("root['INPUTS']['reped_input']['zeffped']",'Zeffped',help='The effective charge number at the pedestal top',default=2.0)

OMFITx.Entry("root['SETTINGS']['PHYSICS']['width']",'width set',help='set width')
OMFITx.Entry("root['SETTINGS']['PHYSICS']['height']",'height set',help='set height')


OMFITx.Separator()
OMFITx.Button('Run REPED',"root['SCRIPTS']['run_REPED'].run")
