import numpy as np
import matplotlib.pyplot as plt

def box(x, box_size):
    return x - np.floor(x/box_size)*box_size

mrow = 1000000000
data = np.loadtxt("data/JFNK_sims/20240214_test_solution_d_fluctuation/ciliate_159fil_9000blob_8.00R_0.0600torsion_true_states.dat", max_rows=mrow)
data = data[:, 2:]

a = data[1:301]
b = data[301:]
nfil = int(a.shape[1]/2)
length = a.shape[0]
aux = np.ones(nfil)*np.pi

print(a.shape, b.shape)
error_norm_array = np.zeros(a.shape[0])

for i in range(length):
    diff = a[i] - b[i]

    a[i][:nfil] = box(a[i][:nfil], 2*np.pi)

    error_norm = np.linalg.norm(diff[:nfil])
    error_norm_array[i] = error_norm
    error = np.linalg.norm(diff[:nfil]) / np.linalg.norm(a[i][:nfil])
    print(i, error_norm, error, np.linalg.norm(a[i][:nfil]))


fig = plt.figure()
ax = fig.add_subplot(1,1,1)

ax.plot(np.arange(length), error_norm_array)

plt.show()