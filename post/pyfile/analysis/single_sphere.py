import numpy as np
import matplotlib.pyplot as plt

# Nblob = 1400
a = 20.52
v_rpy = 2.8145485528e-02
L_list = np.array([43, 46, 52, 58, 66,\
    76, 84, 92, 100, 110,\
    128, 156, 180, 256, 400,\
    512, 700, 1024, 2048, 3072])
v_list = np.array([1.0353384131e-03, 1.5042830716e-03, 2.8763976390e-03, 4.3915972638e-03, 6.3767955099e-03,\
    8.5738858409e-03, 1.0099684246e-02, 1.1431189235e-02, 1.2596349965e-02, 1.3863575235e-02,\
    1.5706014730e-02, 1.7810764363e-02, 1.9131256985e-02, 2.1746806895e-02, 2.4029385356e-02,\
    2.4925992848e-02, 2.5786169908e-02, 2.6536708964e-02, 2.7343904923e-02, 2.7612200528e-02])
phi_list = 4./3.*np.pi*a**3/L_list**3


def hasimoto(c_13):
    return 1 - 1.7601*c_13 + c_13**3 - 1.5593*c_13**6

def sangani(c_13):
    return 1 - 1.7601*c_13 + c_13**3 - 1.5593*c_13**6 + 3.9799*c_13**8 - 3.0734*c_13**10

c_13_list = phi_list**(1./3)
vw_list = v_list/v_rpy
hasimoto_x_list = np.linspace(0, 0.8, 100)
hasimoto_y_list = hasimoto(hasimoto_x_list)
sangani_x_list = np.linspace(0, 0.8, 100)
sangani_y_list = sangani(sangani_x_list)

ax = plt.figure().add_subplot(1,1,1)
ax.plot(sangani_x_list, sangani_y_list, c='r', label='Sangani & Acrivos(1982)')
ax.scatter(c_13_list, vw_list, facecolors='none', edgecolors='r', label='Rigid sphere')
ax.set_xlabel(r'$c^{1/3}$')
ax.set_ylabel(r'V/W')
# ax.set_title('Single sphere settling speed')
ax.set_xlim(0, 0.8)
ax.legend()

plt.savefig('fig/single_sphere.eps', format='eps')
plt.show()




















#