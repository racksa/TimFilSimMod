import numpy as np
import matplotlib.pyplot as plt
import svd_test_data
import random

m = 400
n = 200
x_array = np.linspace(-10, 10, m)
t_array = np.linspace(0, 4*np.pi, n)
X1 = np.zeros((m, n))
X2 = np.zeros((m, n))

def func1(x, t):
    return np.sin(x-0.5*t)
    return 1./np.cosh(x-t) * (np.exp(1j*2.3*t))

def func2(x, t):
    # return 2./np.cosh(x)*np.tanh(x)*np.exp(1j*2.8*t)
    return np.cos(3*t)

for i in range(m):
    for j in range(n):
        X1[i][j] = func1(x_array[i], t_array[j]) + 4*random.random()
        X2[i][j] = func2(x_array[i], t_array[j])

X = X1 #+ X2

U, sigma, V = np.linalg.svd(X, full_matrices=False)
Sigma = np.diag(sigma)

r = 15
V = V.conj().T
U = U[:, :r]
Sigma = Sigma[:r, :r]
V = V[:, :r]

fig = plt.figure()
ax = fig.add_subplot(1,1,1)
fig2 = plt.figure()
ax2 = fig2.add_subplot(1,1,1)
fig3 = plt.figure()
ax3 = fig3.add_subplot(1,1,1)
fig4 = plt.figure()
ax4 = fig4.add_subplot(1,1,1)
# fig5 = plt.figure()
# ax5 = fig5.add_subplot(1,1,1)

ax.imshow(X)
ax.set_xlabel('t')
ax.set_ylabel('x')

ax2.plot(X[:,0])
ax2.set_xlabel('x')
ax2.set_ylabel('f(x, t=0)')


ax3.scatter(np.arange(0, len(sigma), 1), sigma, marker='*')
ax3.set_xlabel('Mode')
ax3.set_ylabel('Eigenvalue')

for i in range(3):
    ax4.plot(U[:,i], label=f'Mode {i}')
ax4.set_xlabel('x')
ax4.set_ylabel('f(x, t=0)')
ax4.legend()

# ax5.plot(U[:,1], c='red')

fig.savefig(f'fig/pod_signal.pdf', bbox_inches = 'tight', format='pdf')
fig2.savefig(f'fig/pod_initial_snapshot.pdf', bbox_inches = 'tight', format='pdf')
fig3.savefig(f'fig/pod_eigenvalues.pdf', bbox_inches = 'tight', format='pdf')
fig4.savefig(f'fig/pod_modes.pdf', bbox_inches = 'tight', format='pdf')
plt.show()

plt.show()






















#