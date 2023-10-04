import numpy as np
import matplotlib.pyplot as plt

# Rod
Nrod_list = np.array([10, 100, 400, 1000, 2304, 5184, 9216, 16384])
rod_fcm_time_list = np.array([7.526577e-03, 1.233129e-02, 2.188203e-02,\
                           4.718907e-02, 1.061803e-01, 2.155749e-01,\
                           3.752160e-01, 6.713320e-01])
rod_total_time_list = np.array([7.963457e-03, 1.540168e-02, 3.568006e-02,\
                             8.603893e-02, 2.140005e-01, 4.673536e-01,\
                             8.318663e-01, 1.590558e+00 ])
rod_Niter_list = np.array([6.6, 14.5, 24, 27, 31.5, 34, 39, 49])

fig1 = plt.figure(figsize=(4.8, 3.2))
ax11 = fig1.add_subplot(1,1,1)
ax12 = ax11.twinx()
ax11.plot(Nrod_list, rod_fcm_time_list, marker='+', linestyle='dotted', c='black', label='Fast FCM')
ax11.plot(Nrod_list, rod_total_time_list, marker='+', linestyle='-', c='black', label='Total')
ax12.plot(Nrod_list, rod_Niter_list , marker='+', c='blue', linestyle='dashed', label='No. of iterations')
ax12.tick_params(axis='y', colors='blue') 
ax11.set_xlabel(r'$N_{rod}$')
ax11.set_ylabel(r'Wall time per iteration/s')
ax11.set_ylim(0)
ax11.set_xlim(0)
ax12.set_xlim(0)
ax12.set_ylim(0, 70)
ax11.legend()
ax12.legend()



# Fil
Nfil_list = np.array([1, 16, 64, 256, 1024, 2704, 4096 ])
fil_fcm_time_list = np.array([2.384288e-02, 2.455505e-02, 1.493097e-02, 1.763091e-02, 
                              4.634576e-02, 1.231594e-01, 1.932845e-01])
fil_total_time_list = np.array([2.389592e-02, 2.495763e-02, 1.707681e-02, 2.719358e-02, 
                                9.131879e-02, 2.813344e-01, 4.684991e-01 ])
fil_Niter_list = np.array([4, 5, 7, 8, 
                           11, 17, 25])
fig2 = plt.figure(figsize=(4.8, 3.2))
ax21 = fig2.add_subplot(1,1,1)
ax22 = ax21.twinx()
ax21.plot(Nfil_list, fil_fcm_time_list, marker='+', linestyle='dotted', c='black', label='Fast FCM')
ax21.plot(Nfil_list, fil_total_time_list, marker='+', linestyle='-', c='black', label='Total')
ax22.plot(Nfil_list, fil_Niter_list , marker='+', c='blue', linestyle='dashed', label='No. of iterations')
ax22.tick_params(axis='y', colors='blue') 
ax21.set_xlabel(r'$N_{fil}$')
ax21.set_ylabel(r'Wall time per iteration/s')
ax21.set_ylim(0)
ax21.set_xlim(0)
ax22.set_xlim(0)
ax22.set_ylim(0, 30)
ax21.legend()
ax22.legend()

fig1.savefig('fig/fastfcm_scaling_rod.pdf', format='pdf', bbox_inches='tight')
fig2.savefig('fig/fastfcm_scaling_fil.pdf', format='pdf', bbox_inches='tight')

plt.show()
