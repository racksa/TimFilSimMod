import numpy as np
import matplotlib.pyplot as plt
import util

prop_cycle = plt.rcParams['axes.prop_cycle']
colors = prop_cycle.by_key()['color']

fig = plt.figure()
ax = fig.add_subplot()
fig2 = plt.figure()
ax2 = fig2.add_subplot(projection='3d')


R = 1
theta0 = 0.02*np.pi
N0 = 16
a = 2
dtheta = 0.02*np.pi
dx = 2*R*np.pi*np.sin(theta0)/N0
dy = R*np.sin(dtheta)

AR_min = 2*R/dy
N0_max = 11*(2*np.pi*np.sin(theta0)/R*AR_min)
print(f'AR > {AR_min} to minimise collision')
print(f'N0 < {N0_max} to minimise collision')

theta_list = list()
N_list = list()


while (theta0<0.5*np.pi):
    theta_list.append(theta0)
    N_list.append(N0)

    # sin_next_theta = np.sin(theta0)*(1 + a/N0)
    # if(sin_next_theta <= 1):
    #     theta0 = np.arcsin(sin_next_theta)
    #     N0 += 2
    # else:
    #     break

    theta0 += dtheta
    l0 = 2*np.pi*R*np.sin(theta0)
    N0 = int(2*np.pi*np.sin(theta0)/dx)

num_azimuth = len(theta_list)
num_fil = np.sum(N_list)

print(f'num fil = {num_fil}')

# C-style
index = 0
phi0 = 0
dphi = 0
fil_pos = np.zeros((num_fil, 3))
for i in range(num_azimuth):
    phi0 += 0.5*dphi
    dphi = 2*np.pi/N_list[i]
    start_index = index
    for j in range(N_list[i]):
        fil_pos[index] = util.spherical_to_cartesian(R, theta_list[i], phi0+dphi*j)
        index += 1
    ax2.plot(fil_pos[start_index:index,0], fil_pos[start_index:index, 1], fil_pos[start_index:index, 2], c=colors[i%len(colors)])


# # Python style
# N_cumsum = np.zeros(len(N_list)+1, dtype=int)
# N_cumsum[1:] = np.cumsum(N_list)
# dphi_array = 2*np.pi/np.array(N_list)
# phi0_array = np.cumsum(0.5*dphi_array)
# fil_pos = np.zeros((num_fil, 3))
# for i in range(num_azimuth):
#     for j in range(N_list[i]):
#         index = N_cumsum[i] + j
#         fil_pos[index] = util.spherical_to_cartesian(R, theta_list[i], phi0_array[i]+dphi_array[i]*j)
#     ax2.plot(fil_pos[N_cumsum[i]:N_cumsum[i+1],0],\
#                 fil_pos[N_cumsum[i]:N_cumsum[i+1], 1],\
#                     fil_pos[N_cumsum[i]:N_cumsum[i+1], 2],\
#                         c=colors[i%len(colors)])

ax.plot(theta_list)
ax2.set_box_aspect((np.ptp(fil_pos[:,0]), np.ptp(fil_pos[:,1]), np.ptp(fil_pos[:,2]))) 
ax2.scatter(fil_pos[:,0], fil_pos[:,1], fil_pos[:,2])



plt.show()























#