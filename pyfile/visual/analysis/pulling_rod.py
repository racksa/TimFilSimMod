import numpy as np
import matplotlib.pyplot as plt


blob_sym_list = np.array([22, 30, 42, \
    50, 62, 70, 82])
v_hor_sym_list = np.array([1.9598670072e-01, 1.6083591412e-01, 1.2843079173e-01, \
    1.1385833459e-01, 9.7811970955e-02, 8.9642491850e-02, 7.9883812079e-02])
v_ver_sym_list = np.array([2.5969245307e-01, 2.2026104495e-01, 1.8175742518e-01, \
    1.6368815209e-01, 1.4319528085e-01, 1.3251698722e-01, 1.1952499763e-01])

blob_asym_list = np.array([20, 32, 40, \
    52, 60, 72, 80])
v_hor_asym_list = np.array([2.1070878503e-01, 1.5553951082e-01, 1.3374338650e-01, \
    1.1140790113e-01, 1.0063186726e-01, 8.8207575298e-02, 8.1657447005e-02])
v_ver_asym_list = np.array([2.7248169070e-01, 2.1251516653e-01, 1.8704118787e-01, \
    1.5980288883e-01, 1.4618844119e-01, 1.3012790300e-01, 1.2148755043e-01 ])

vh_ratio_sym_list = v_hor_sym_list/v_ver_sym_list
vh_ratio_asym_list = v_hor_asym_list/v_ver_asym_list


########################################################################
# Plot
########################################################################
ax = plt.figure().add_subplot(1,1,1)
ax.scatter(blob_sym_list, vh_ratio_sym_list, facecolors='none', edgecolors='r', label='Symmetric')
ax.scatter(blob_asym_list, vh_ratio_asym_list, facecolors='none', edgecolors='b', label='Asymmetric')
ax.set_xlabel(r'$Blob\ number$')
ax.set_ylabel(r'$V_{vertical} / V_{horizontal}$')
# ax.set_title('Rod polarisation speed graph')
ax.legend()

plt.savefig('fig/rod_polarisation.eps', format='eps')
plt.show()