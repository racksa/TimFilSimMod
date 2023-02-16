import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.lines as mlines
import matplotlib.colors as mcolors
import random
from batchelor1982 import *

# Fourier modes = 0.5 0
a = 20.533623615045787
c = 4.*np.pi*a**2/1400
num_blob_list = np.array([12, 24, 48, 96, 
                          192, 384, 768, 1536])
predicted_a_list = np.sqrt(num_blob_list*c/4/np.pi)
fourier_mode_list = predicted_a_list*0.5/20.52
print('fourier modes', fourier_mode_list)

a_list = np.array([1.9022991125314643, 2.690261716530194, 3.804604749247682, 5.380542640881668,
                   7.609213419601921, 10.761086958962553, 15.218399297665027, 21.522126532652486]) 
ra_list = np.array([2.01, 2.02, 2.03,
                    2.05, 2.08, 2.12, 2.15, 
                    2.2, 2.3, 2.4, 2.5,
                    2.6, 2.8, 3.0, 3.5,
                    4.0, 4.5
                    ])
a_list = np.reshape(a_list, (-1, 1))
r_list = (a_list * ra_list)
print('r_list', r_list)
# 0.04632174 0.06550883 0.09264348 0.13101766 
# 0.18528695 0.26203532 0.37057391 0.52407065

v_prl_list = np.array([[4.2791227227e-01, 4.2753576811e-01, 4.2716469508e-01,
                        4.2641803912e-01, 4.2530818224e-01, 4.2384664454e-01, 4.2276431223e-01,
                        4.2098650780e-01, 4.1752108143e-01, 4.1417896636e-01, 4.1094725796e-01,
                        4.0782327055e-01, 4.0187379304e-01, 3.9627969538e-01, 3.8365634359e-01,
                        3.7276392384e-01, 3.6339710679e-01
                        ],
                       [3.1108870930e-01, 3.1082014800e-01, 3.1055246589e-01,
                        3.1001973602e-01, 3.0922725709e-01, 3.0818290575e-01, 3.0740878287e-01,
                        3.0614963674e-01, 3.0366797947e-01, 3.0126586581e-01, 2.9894176614e-01,
                        2.9669182371e-01, 2.9239895971e-01, 2.8835980440e-01, 2.7924896736e-01,
                        2.7140282849e-01, 2.6467350257e-01
                        ],
                       [2.2514252397e-01, 2.2494676386e-01, 2.2475171391e-01,
                        2.2436371135e-01, 2.2378681559e-01, 2.2302657842e-01, 2.2246299659e-01,
                        2.2153615483e-01, 2.1972897915e-01, 2.1798122669e-01, 2.1629014897e-01,
                        2.1465312700e-01, 2.1153113037e-01, 2.0859679935e-01, 2.0199040354e-01,
                        1.9633087864e-01, 1.9149688356e-01
                        ],
                       [1.6202180842e-01, 1.6187935075e-01, 1.6173744613e-01,
                        1.6145520355e-01, 1.6103562765e-01, 1.6048310283e-01, 1.6007388177e-01,
                        1.5940139142e-01, 1.5808978014e-01, 1.5682219113e-01, 1.5559646579e-01,
                        1.5441020861e-01, 1.5214918865e-01, 1.5002615964e-01, 1.4526031310e-01,
                        1.4118899241e-01, 1.3772372319e-01
                        ],
                       [1.1593401025e-01, 1.1583164020e-01, 1.1572962500e-01,
                        1.1552664998e-01, 1.1522480144e-01, 1.1482716854e-01, 1.1453259611e-01,
                        1.1404843763e-01, 1.1310476921e-01, 1.1219268736e-01, 1.1131072103e-01,
                        1.1045746995e-01, 1.0883187909e-01, 1.0730648423e-01, 1.0388779078e-01,
                        1.0097539057e-01, 9.8502700676e-02],
                       [8.2608756791e-02, 8.2535571876e-02, 8.2462630653e-02,
                        8.2317471667e-02, 8.2101514467e-02, 8.1816874080e-02, 8.1605865073e-02,
                        8.1266387675e-02, 8.0591654378e-02, 7.9939486495e-02, 7.9308858720e-02,
                        7.8698791716e-02, 7.7536675531e-02, 7.6446576592e-02, 7.3999492628e-02,
                        7.1926335830e-02, 7.0166912424e-02],
                       [5.8724269863e-02, 5.8671999876e-02, 5.8619916599e-02,
                        5.8516220702e-02, 5.8362814035e-02, 5.8160283747e-02, 5.8010247470e-02,
                        5.7763605216e-02, 5.7281493106e-02, 5.6816482276e-02, 5.6366887943e-02,
                        5.5932038448e-02, 5.5104259542e-02, 5.4327462524e-02, 5.2589675796e-02,
                        5.1109380476e-02, 4.9865525017e-02],
                       [4.1664469450e-02, 4.1627519891e-02, 4.1590705215e-02,
                        4.1517477854e-02, 4.1408623013e-02, 4.1265279337e-02, 4.1159074641e-02,
                        4.0984449957e-02, 4.0643719401e-02, 4.0308450957e-02, 3.9990271960e-02,
                        3.9682178830e-02, 3.9094863963e-02, 3.8544115999e-02, 3.7310976094e-02,
                        3.6252309492e-02, 3.5368497666e-02]])
v_isolate_list = np.array([2.6872499589e-01, 1.9746582184e-01, 1.4382615990e-01, 1.0384768970e-01,
                           7.4492869800e-02, 5.3170590322e-02, 3.7835511441e-02, 2.6866000007e-02])
                          
########################################################################
# Plot Error
########################################################################
def random_color(seed):
    random.seed(3*seed)
    cr = random.random()
    random.seed(3*seed+1)
    cg = random.random()
    random.seed(3*seed+2)
    cb = random.random()
    return (cr, cg, cb)
    
vw_prl_list = v_prl_list/np.reshape(v_isolate_list, (-1,1))
vw_prl_batchelor_list = A11(ra_list, 1) + A12(ra_list,1)
vw_prl_error_list = abs(vw_prl_list - vw_prl_batchelor_list) / vw_prl_batchelor_list

ax2 = plt.figure().add_subplot(1,1,1)
for i in range(len(num_blob_list)):
    ax2.plot(ra_list, vw_prl_error_list[i], marker='+', color=random_color(i), label=r'$N_{blob}=$' + str(num_blob_list[i]))
ax2.set_yscale('log')
ax2.set_xlabel(r'$r/a$')
ax2.set_ylabel(r'$\frac{|V/W - V/W_{exact}|}{V/W_{exact}}$')
ax2.set_title('Fixed number density on surface')
ax2.set_xlim(2, 4.5)
ax2.set_ylim(1e-4, 1e-1)
plt.savefig('resolution_test_fixed_number_density.eps', format='eps')
ax2.legend(ncol=2)
plt.show()