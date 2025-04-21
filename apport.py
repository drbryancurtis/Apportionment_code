import numpy as np
from sympy import Matrix
from cmath import sqrt

# Input: 2x2 ndarray and optional apportionment constant
# Output: apportioned matrix, apportioning matrix
def apport2(mat: np.ndarray, k = 1):
    # Check inputs are valid
    if not isinstance(mat, np.ndarray):
        raise TypeError("Input must be a numpy array")
    if mat.shape != (2, 2):
        raise ValueError("Input array must have shape (2,2)")
    if k <= 0:
        raise ValueError("k must be positive")
    
    evals, evecs = np.linalg.eig(mat)

    # check if all zeros matrix
    if not np.any(mat):
        matM = np.array([[1,0],[0,1]])
        matMinv = np.array([[1,0],[0,1]])

    # check if nilpotent
    elif evals[0] == 0 and evals[1] == 0: 
        M = Matrix(mat) # convert to matrix - jordan_form() comes from sympy 
        P = M.jordan_form()[0]

        matP = np.array(P).astype(np.float64)
        matPinv = np.linalg.inv(matP)

        matM = np.matmul(np.array([[np.sqrt(k), 1],[np.sqrt(k), 1/np.sqrt(k) + 1]]), matP)
        matMinv = np.matmul(matPinv, np.array([[1/np.sqrt(k) + 1, -1],[-np.sqrt(k), np.sqrt(k)]]))

    # check if rank 1 (not nilpotent)
    elif evals[0] == 0 or evals[1] == 0:
        rho = max(abs(evals))
        r = k/rho

        if k < rho/2:
            print('cannot apportion for k =', k)
            print('instead using k =', rho/2)
            k = rho/2

        matP = evecs
        matPinv = np.linalg.inv(evecs)

        matM = np.matmul(np.array([[1, 1/2 - np.sqrt(r**2 - 0.25)*1j], [1, -1/2 + np.sqrt(r**2 - 0.25)*1j]]), matPinv)
        matMinv =  np.matmul(matP, np.array([[1/2 + np.sqrt(r**2 - 0.25)*1j, 1/2 - np.sqrt(r**2 - 0.25)*1j], [1, -1]]))

    # mat is rank 2
    else:
        if evals[0] == evals[1]:
            raise ValueError('matrix is not apportionable')

        g = (evals[0] + evals[1])/(evals[0] - evals[1])
        rho = max(abs(evals))

        if g == 0:
            if k < 1/np.sqrt(2):
                print('cannot apportion for k =', k)
                print('instead using k =', 1/np.sqrt(2))
                k = 1/np.sqrt(2)

            w = np.sqrt(k/(2*rho) + 0.25) + np.sqrt(k/(2*rho) - 0.25)*1j

        elif (g**2).real >= abs(g)**4 or abs(g) > 1:
            raise ValueError('matrix is not apportionable')

        elif abs(g) == 1:
            w = 0

        else:
            w = g * np.sqrt((1 - abs(g)**4)/(2*(abs(g)**4 - (g**2).real))) * 1j

        matP = evecs
        matPinv = np.linalg.inv(evecs)

        b = sqrt(w**2 - 1)/2
        matM = np.matmul(np.array([[1, b],[(w - 1)/(2 * b), (w + 1)/2]]),matPinv)
        matMinv= np.matmul(matP,np.array([[(w + 1)/2, -b],[(1 - w)/(2 * b), 1]]))

    matB = np.matmul(np.matmul(matM, mat), matMinv)

    return matB, matM

#
# Testing
#

mat = np.array([[1,1],[0,-1]])
apport, matM = apport2(mat, 0.2)
print(apport)
print(abs(apport))
print(matM)
