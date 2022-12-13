import numpy as np
import math

def rot_mat(quaternion):
    ret = np.zeros((3,3))
    ret[0, 0] = 1.0
    ret[1, 1] = 1.0
    ret[2, 2] = 1.0

    temp = 2.0*quaternion[1]*quaternion[1]
    ret[1, 1] -= temp;
    ret[2, 2] -= temp

    temp = 2.0*quaternion[2]*quaternion[2]
    ret[0, 0] -= temp
    ret[2, 2] -= temp

    temp = 2.0*quaternion[3]*quaternion[3]
    ret[0, 0] -= temp
    ret[1, 1] -= temp

    temp = 2.0*quaternion[1]*quaternion[2]
    ret[1, 0] = temp
    ret[0, 1] = temp

    temp = 2.0*quaternion[1]*quaternion[3]
    ret[2, 0] = temp
    ret[0, 2] = temp

    temp = 2.0*quaternion[2]*quaternion[3]
    ret[1, 2] = temp
    ret[2, 1] = temp

    temp = 2.0*quaternion[0]*quaternion[3]
    ret[1, 0] += temp;
    ret[0, 1] -= temp;

    temp = 2.0*quaternion[0]*quaternion[2]
    ret[2, 0] -= temp;
    ret[0, 2] += temp;

    temp = 2.0*quaternion[0]*quaternion[1];
    ret[2, 1] += temp;
    ret[1, 2] -= temp;

    return ret


def blob_point_from_data(body_states, blob_references):
    blob_pos = np.matmul(rot_mat(body_states[3:7]), blob_references)
    # print('blobref, blobpos', blob_references, blob_pos)

    x = body_states[0] + blob_pos[0]
    y = body_states[1] + blob_pos[1]
    z = body_states[2] + blob_pos[2]

    return x, y, z


def set_axes_equal(ax):
    '''Make axes of 3D plot have equal scale so that spheres appear as spheres,
    cubes as cubes, etc..  This is one possible solution to Matplotlib's
    ax.set_aspect('equal') and ax.axis('equal') not working for 3D.

    Input
      ax: a matplotlib axis, e.g., as output from plt.gca().
    '''

    x_limits = ax.get_xlim3d()
    y_limits = ax.get_ylim3d()
    z_limits = ax.get_zlim3d()

    x_range = abs(x_limits[1] - x_limits[0])
    x_middle = np.mean(x_limits)
    y_range = abs(y_limits[1] - y_limits[0])
    y_middle = np.mean(y_limits)
    z_range = abs(z_limits[1] - z_limits[0])
    z_middle = np.mean(z_limits)

    # The plot bounding box is a sphere in the sense of the infinity
    # norm, hence I call half the max range the plot radius.
    plot_radius = 0.5*max([x_range, y_range, z_range])

    ax.set_xlim3d([x_middle - plot_radius, x_middle + plot_radius])
    ax.set_ylim3d([y_middle - plot_radius, y_middle + plot_radius])
    ax.set_zlim3d([z_middle - plot_radius, z_middle + plot_radius])


def box(x, box_size):
    return x - math.floor(x/box_size)*box_size;

def two_points_at_boundary(two_points_x, two_points_y, two_points_z, rod_length):
    points = np.array((two_points_x, two_points_y, two_points_z))
    vector = points[:,1] - points[:,0]
    if np.linalg.norm(vector) > rod_length:
        return [0,0], [0,0], [0,0]
        # two_points_x, two_points_y, two_points_z = [], [], []
    return two_points_x, two_points_y, two_points_z




#