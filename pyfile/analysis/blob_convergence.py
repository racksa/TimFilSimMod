import numpy as np
import matplotlib.pyplot as plt

AR=12.65
L = 19*2.6
R = 0.5*AR*L

dissipation_list = np.array([6153.945039746058, 6155.715800350974, 6154.048970565967, 
                             6166.055095573409, 6156.84580080739, 6209.91026945243, 
                             6256.858922291375, 6343.197801544943, 6464.95522387998,
                             6615.488752320283, 6760.217836996505, 6867.205154016618, 
                             6933.9848578267065, 6966.285556960928,])

v_list = np.array([0.9527672397629782, 0.9521426595065936, 0.9533044139949042, 
                   0.9484213083968568, 0.9514570905604488, 0.9355456568687464, 
                   0.9195028224280672, 0.8972491963934501, 0.870002751022489, 
                   0.8459229101541402, 0.8317494836413607, 0.8256686488844301, 
                   0.8235959141938096, 0.8234286658255042,])

blob_list = np.array([int(10*4**(0.5*i)+2) for i in range(len(v_list))])
density_list = blob_list*L**2/R**2

v_accurate_value = v_list[-1]
v_error_list = np.abs(v_list - v_accurate_value)/v_accurate_value

dis_accurate_value = dissipation_list[-1]
dis_error_list = np.abs(dissipation_list - dis_accurate_value)/dis_accurate_value


print(blob_list)
print(v_list)
print(density_list)

fig = plt.figure()
ax = fig.add_subplot(1,1,1)
fig2 = plt.figure()
ax2 = fig2.add_subplot(1,1,1)
# fig3 = plt.figure()
# ax3 = fig3.add_subplot(1,1,1)

ax.plot(blob_list[:-1], v_error_list[:-1], c='r', marker='+')
ax.set_xlabel(r'$N_b$')
ax.set_ylabel(r'% V error')
ax.set_xlim(0)
ax.set_yscale('log')

ax2.plot(blob_list[:-1], dis_error_list[:-1], c='r', marker='+')
ax2.set_xlabel(r'$N_b$')
ax2.set_ylabel(r'% dis. error')
ax2.set_xlim(0)
ax2.set_yscale('log')


fig.savefig('fig/v_convergence.pdf', format='pdf')
fig2.savefig('fig/dis_convergence.pdf', format='pdf')
plt.show()