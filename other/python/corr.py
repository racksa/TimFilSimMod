import numpy as np

deg = 1
original = np.array([1, 3, 6, 10, 15, 21, 16, 12, 8, 4])


# diff2 = np.zeros(np.shape(original))

diff = np.diff(original, prepend=original[-1:])

corr = np.abs(diff[:-1]) + np.abs(diff[1:]) # |a_i-a_{i-1}| + |a_{i+1}-a_i|



print('original', original)
print('diff', diff)
print('diff2', corr)