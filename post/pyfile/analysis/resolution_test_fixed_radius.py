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

ra_list = np.array([2.01, 2.02, 2.03,
                    2.05, 2.08, 2.12, 2.15, 
                    2.2, 2.3, 2.4, 2.5,
                    2.6, 2.8, 3.0, 3.5,
                    4.0, 4.5
                    ])
a_list = np.array([a for i in range(len(num_blob_list))])
a_list = np.reshape(a_list, (-1, 1))
r_list = (a_list * ra_list)
print('r_list', r_list)

v_prl_list = np.array([[
                        ],
                       [
                        ],
                       [
                        ],
                       [
                        ],
                       [],
                       [],
                       [],
                       []])

v_isolate = 2.8145485528e-02
                          
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
    
vw_prl_list = v_prl_list/v_isolate
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