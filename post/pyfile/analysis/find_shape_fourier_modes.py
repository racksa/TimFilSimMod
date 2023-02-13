import numpy as np
import matplotlib.pyplot as plt
import scipy.integrate as spi

L = np.pi
aspect_ratio = 0.1
num_modes = 201
P = 2*L

a = L
ar = aspect_ratio

def spheroid(x):
    return np.sqrt(ar**2 * (a**2 - x**2))


spheroid_x_array = np.linspace(-L, L, num_modes)
spheroid_y_array = spheroid(spheroid_x_array)

def compute_real_fourier_coeffs(func, N):
    result = []
    for n in range(N+1):
        an = (2./P) * spi.quad(lambda t: func(t) * np.cos(2 * np.pi * n * t / P), 0, P)[0]
        bn = (2./P) * spi.quad(lambda t: func(t) * np.sin(2 * np.pi * n * t / P), 0, P)[0]
        result.append((an, bn))
    return np.array(result)

#function that computes the real form Fourier series using an and bn coefficients
def fit_func_by_fourier_series_with_real_coeffs(t, AB):
    result = 0.
    A = AB[:,0]
    B = AB[:,1]
    for n in range(0, len(AB)):
        if n > 0:
            result +=  A[n] * np.cos(2. * np.pi * n * t / P) + B[n] * np.sin(2. * np.pi * n * t / P)
        else:
            result +=  A[0]/2.
    return result

def DFT(x):
    N = len(x)
    n = np.arange(N)
    k = n.reshape((N, 1))

    cosine = np.cos(2*np.pi*k*n/P)
    sine = np.sin(2*np.pi*k*n/P)

    a_coeff = np.dot(cosine, x)/ num_modes
    b_coeff = np.dot(sine, x)  / num_modes

    return a_coeff, b_coeff

def IDFT(an, bn):
    N = len(an)
    n = np.arange(N)
    x = n.reshape((N, 1))

    cosine = np.cos(2*np.pi*x*n/P)
    sine = np.sin(2*np.pi*x*n/P)

    X = np.dot(an[1:], cosine[1:]) + np.dot(bn[1:], sine[1:])
    return X

# def DFT(x):
#     N = len(x)
#     n = np.arange(N)
#     k = n.reshape((N, 1))
#     e = np.exp(-2j * np.pi * k * n / N)
#     X = np.dot(e, x)
    
#     return X

# def IDFT(k):
#     N = len(k)
#     n = np.arange(N)
#     x = n.reshape((N, 1))
#     e = np.exp(2j * np.pi * x * n / N)/N
#     X = np.dot(e, k)
    
#     return X

an, bn = DFT(spheroid_y_array)
# print(an, bn)
spheroid_y_array_IDFT = IDFT(an, bn)

# coeffs = compute_real_fourier_coeffs(spheroid, num_modes)
# spheroid_y_array_IDFT = fit_func_by_fourier_series_with_real_coeffs(spheroid_x_array, coeffs)

########################################################################
# Plot
########################################################################
ax = plt.figure().add_subplot(1,1,1)
ax.plot(spheroid_x_array, spheroid_y_array, c='r', label='Prolate func.')
ax.scatter(spheroid_x_array, spheroid_y_array_IDFT, facecolors='none', edgecolors='b', label='DFT')
ax.set_xlabel(r'$Blob\ number$')
ax.set_ylabel(r'$V_{vertical} / V_{horizontal}$')
ax.set_title('Rod polarisation speed graph')
ax.legend()
ax.set_aspect('equal')

plt.savefig('rod_polarisation.eps', format='eps')
plt.show()