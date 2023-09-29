import numpy as np
import matplotlib.pyplot as plt

Nrod_list = np.array([10, 100, 400, 1000, 2304, 5184, 9216, 16384])
fcm_time_list = np.array([7.526577e-03, 1.233129e-02, 2.188203e-02,\
                           4.718907e-02, 1.061803e-01, 2.155749e-01,\
                           3.752160e-01, 6.713320e-01])
total_time_list = np.array([7.963457e-03, 1.540168e-02, 3.568006e-02,\
                             8.603893e-02, 2.140005e-01, 4.673536e-01,\
                             8.318663e-01, 1.590558e+00 ])
Niter_list = np.array([6.6, 14.5, 24, 27, 31.5, 34, 39, 49])

ax = plt.figure().add_subplot(1,1,1)
ax2 = ax.twinx()
ax.plot(Nrod_list, fcm_time_list, marker='+', linestyle='dotted', c='black', label='Fast FCM')
ax.plot(Nrod_list, total_time_list, marker='+', linestyle='-', c='black', label='Total')
ax2.plot(Nrod_list, Niter_list, marker='+', c='blue', label='No. of iterations')
ax2.tick_params(axis='y', colors='blue') 
ax.set_xlabel(r'$N_{rod}$')
ax.set_ylabel(r'Wall time per iteration/s')
ax.set_ylim(0)
ax.set_xlim(0)
ax2.set_xlim(0)
ax2.set_ylim(0, 70)
ax.legend()
ax2.legend()

plt.savefig('fig/fastfcm_scaling_rod.pdf', format='pdf')
plt.show()
