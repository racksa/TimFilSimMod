import numpy as np
import matplotlib.pyplot as plt

# Nblob = 1400
a = 20.52
v_isolate = 2.8145485528e-02
r_list = np.array([43, 47, 50, 55,\
    60, 65, 70, 75, 80,\
    85, 90, 95, 100, 105,\
    120, 135, 150, 165, 180,\
    200, 220, 240, 260, 280
    ])
v_ver_list = np.array([4.3330605504e-02, 4.2620063376e-02, 4.2116847789e-02, 4.1330910102e-02,\
    4.0604915553e-02, 3.9933188305e-02, 3.9311993202e-02, 3.8738259007e-02, 3.8209018617e-02,\
    3.7706797263e-02, 3.7258035973e-02, 3.6844123261e-02, 3.6462003039e-02, 3.6108788790e-02,\
    3.5196939310e-02, 3.4462238311e-02, 3.3860622583e-02, 3.3360424574e-02, 3.2938786281e-02,\
    3.2470567208e-02, 3.2084481573e-02, 3.1760965101e-02, 3.1486128594e-02, 3.1249868122e-02])

v_hor_list = np.array([3.9304687476e-02, 3.8359559074e-02, 3.7700578582e-02, 3.6755490273e-02,\
    3.5958657052e-02, 3.5296898002e-02, 3.4737697762e-02, 3.4259881445e-02, 3.3847293250e-02,\
    3.3487614337e-02, 3.3171372551e-02, 3.2891196575e-02, 3.2641270680e-02, 3.2416940368e-02,\
    3.1863043569e-02, 3.1437938325e-02, 3.1100687974e-02, 3.0826456227e-02, 3.0599053635e-02,\
    3.0350038381e-02, 3.0147127153e-02, 2.9978565764e-02, 2.9836288538e-02, 2.9714576130e-02])

########################################################################
# Cooley & Oneil
########################################################################
ra_cooleyOneil_ver_list = np.array([2, 2.01, 2.05, 2.1, 2.25,\
    2.5, 2.75, 3, 3.5, 4,\
    5, 7])
vw_cooleyOneil_ver_list = np.array([1.5500, 1.5487, 1.5431, 1.5363, 1.5167,\
    1.4861, 1.4580, 1.4320, 1.3861, 1.3472,\
    1.2866, 1.2100])
ra_cooleyOneil_hor_list = np.array([2, 2.001, 2.005, 2.01, 2.1, 3.0, 7.0])
vw_cooleyOneil_hor_list = np.array([1.381, 1.4004, 1.4027, 1.4032, 1.3918, 1.2668, 1.1086])

########################################################################
# Batchelor (1982)
########################################################################
def A11(p, lam):
    return 1 - 60*lam**2/(1+lam)**4/p**4 - 192*lam**3*(5-22*lam**2+3*lam**4)/(1+lam)**8/p**8

def A12(p, lam):
    return 1.5/p - 2*(1+lam**2)/((1+lam)**2*p**3) + 1200*lam**3/(1+lam)**6/p**7

def B11(p, lam):
    return 1 - 68*lam**5/(1+lam)**6/p**6 - 32*lam**3*(10-9*lam**2+9*lam**4)/(1+lam)**8/p**8

def B12(p, lam):
    return 0.75/p + (1+lam**2)/(1+lam)**2/p**3

########################################################################
# Plot
########################################################################
vverw_list = v_ver_list/v_isolate
vhorw_list = v_hor_list/v_isolate

ra_batchelor_ver_list = np.linspace(2, 14, 100)
vw_batchelor_ver_list = A11(ra_batchelor_ver_list, 1) + A12(ra_batchelor_ver_list, 1)
ra_batchelor_hor_list = np.linspace(2, 14, 100)
vw_batchelor_hor_list = B11(ra_batchelor_hor_list, 1) + B12(ra_batchelor_hor_list, 1)

ax = plt.figure().add_subplot(1,1,1)
ax.scatter(r_list/a, vverw_list, facecolors='none', edgecolors='r', label='Vertical alignment')
ax.scatter(r_list/a, vhorw_list, facecolors='none', edgecolors='b', label='Horizontal alignment')
# ax.scatter(ra_cooleyOneil_ver_list, vw_cooleyOneil_ver_list, c='r', label='Vertical alignment cooleyOneil')
# ax.scatter(ra_cooleyOneil_hor_list, vw_cooleyOneil_hor_list, c='b', label='Horizontal alignment cooleyOneil')
ax.plot(ra_batchelor_ver_list, vw_batchelor_ver_list, c='r', label='Vertical alignment Batchelor(1982)')
ax.plot(ra_batchelor_hor_list, vw_batchelor_hor_list, c='b', label='Horizontal alignment Batchelor(1982)')
ax.set_xlabel(r'$r/a$')
ax.set_ylabel(r'V/W')
ax.set_title('Two sphere settling speed')
ax.set_xlim(2, 14)
ax.set_ylim(1, 1.7)
ax.legend()

plt.savefig('two_sphere.eps', format='eps')
plt.show()







