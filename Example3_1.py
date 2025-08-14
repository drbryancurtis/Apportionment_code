import numpy as np

#
# Verification of the calculations in Example 3.1
#

# Define the matrix A
A = np.array([
    [0, 1, 0, 0, 0],
    [0, 0, 1, 0, 0],
    [0, 0, 0, 0, 0],
    [0, 0, 0, 0, 1],
    [0, 0, 0, 0, 0]
], dtype=complex)

# Define the values for d_j
d1 = np.exp(0 * np.pi/3 * 1j)
d2 = np.exp(1 * np.pi/3 * 1j)
d3 = np.exp(2 * np.pi/3 * 1j)
d4 = np.exp(3 * np.pi/3 * 1j)
d5 = -np.exp(2 * np.pi/3 * 1j)

# Define the matrix M
M = np.array([
    [1, 0, 0, 0, d5],
    [d1, 1, 0, 0, 0],
    [0, d2, 1, 0, 0],
    [0, 0, d3, 1, 0],
    [0, 0, 0, d4, 1]
], dtype=complex)

# Compute the inverse of M
M_inv = np.linalg.inv(M)

# Compute MAM_inv
B= M @ A @ M_inv

# Making the displayed matrices more readable
np.set_printoptions(precision=3, suppress=True)

print('Matrix A = \n' + str(A))
print('Matrix M = \n' + str(M))
print('Inverse of M is \n' + str(M_inv))
print('MAM^{-1} = \n' + str(B))
print('Modulus of MAM^{-1} is \n' + str(np.abs(B)))