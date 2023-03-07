import numpy as np
import matplotlib.pyplot as plt
import scipy.integrate as spi

L = np.pi
aspect_ratio = 5
num_modes = 43
P = 2*L

a = L
ar = aspect_ratio

def spheroid(x):
    return np.sqrt(ar**2 * (a**2 - x**2))

def spheroid_r(theta):
    return np.sqrt(ar**2/(np.cos(theta)**2 + (ar*np.sin(theta))**2))

def DFT(x):
    # N = len(x)
    # n = np.arange(N)
    # k = n.reshape((N, 1))
    # e = np.exp(-2j * np.pi * k * n / N)/N

    N = len(x)
    n = np.arange(N)
    k = np.linspace(0, 2*np.pi-(2*np.pi)/N, N).reshape((N, 1))
    e = np.exp(-1j*k*n) / N

    X = np.dot(e, x)
    
    return X

def IDFT(k):
    # N = len(k)
    # n = np.arange(N)
    # x = n.reshape((N, 1))
    # e = np.exp(2j * np.pi * x * n / N)

    N = len(k)
    n = np.arange(N)
    x = np.linspace(0, 2*np.pi-(2*np.pi)/N, N).reshape((N, 1))
    e = np.exp(1j*x*n)

    # cosine, sine = e.real, e.imag
    # an, bn = k.real, -k.imag

    X = np.dot(e, k)

    return X

# spheroid_x_array = np.linspace(-L, L, num_modes)
# spheroid_y_array = spheroid(spheroid_x_array)
# K = DFT(spheroid_y_array)
# spheroid_y_array_IDFT = IDFT(K)

spheroid_theta_array = np.linspace(0, 2*np.pi, num_modes)
spheroid_r_array = spheroid_r(spheroid_theta_array)
Kr = DFT(spheroid_r_array)
spheroid_r_array_IDFT = IDFT(Kr)

file1 = open('spheroid.fourier_modes', 'w')
file1.close()
file1 = open('spheroid.fourier_modes', 'a')
file1.write(str(num_modes) + ' ')
for i in range(num_modes):
    file1.write(str(Kr[i].real) + ' ' + str(-Kr[i].imag) + ' ')
file1.close()

########################################################################
# Plot
########################################################################
# frameon=False
ax = plt.figure().add_subplot(1,1,1)
# ax.plot(spheroid_x_array, spheroid_y_array, c='black', label='Prolate func.')
# ax.scatter(spheroid_x_array, spheroid_y_array_IDFT, facecolors='none', edgecolors='b', label='DFT')

# ax.plot(spheroid_theta_array, spheroid_r_array, c='black', label=r'Prolate $r(\theta)$.')
# ax.scatter(spheroid_theta_array, spheroid_r_array_IDFT, facecolors='none', edgecolors='b', label='DFT')

ax.plot(np.cos(spheroid_theta_array)*spheroid_r_array, np.sin(spheroid_theta_array)*spheroid_r_array, c='black', label=r'Prolate $r(\theta)$.')
# ax.scatter(np.cos(spheroid_theta_array)*spheroid_r_array_IDFT, np.sin(spheroid_theta_array)*spheroid_r_array_IDFT, facecolors='none', edgecolors='b', label='DFT')


# ax.set_xlabel(r'$Blob\ number$')
# ax.set_ylabel(r'$V_{vertical} / V_{horizontal}$')
# ax.set_title('Rod polarisation speed graph')
# ax.legend()
ax.set_aspect('equal')
ax.axis('off')

plt.savefig('fig/spheroid_shape.tif', format='tif')
plt.show()