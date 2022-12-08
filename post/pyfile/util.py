import numpy as np

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

    x = body_states[0] + blob_pos[0]
    y = body_states[1] + blob_pos[1]
    z = body_states[2] + blob_pos[2]

    return x, y, z










#