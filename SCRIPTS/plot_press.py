press = OMFIT['new_kineticEFIT']['REPED']['INPUTS']['gEQDSK1']['PRES']/1000
rho = OMFIT['new_kineticEFIT']['REPED']['INPUTS']['gEQDSK1']['RHOVN']

press_REPED = OMFIT['new_kineticEFIT']['REPED']['OUTPUTS']['profiles']['data']['col4']
rho_REPED = OMFIT['new_kineticEFIT']['REPED']['OUTPUTS']['profiles']['data']['col1']

xlim(0.8,1.005)
ylim(0,120)


plt.plot(rho,press,'--',label='pressure_exp',linewidth=1.0,c='black')
plt.plot(rho_REPED,press_REPED,'-',label='mtanh',linewidth=1.0,c='red')
xlabel(r'$\rho$',fontsize=16)
ylabel('pressure (kPa)',fontsize=16)
plt.xticks(fontsize=13)
plt.yticks(fontsize=13)
plt.legend(loc='lower left',fontsize=16)
