import numpy as np
from sympy import Matrix
from cmath import sqrt

# Input: 2x2 ndarray and optional apportionment constant
# Output: apportioned matrix, apportioning matrix
def apport2x2(mat: np.ndarray, k = 1):
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
        print('k = 0 is the only apportionment constant')
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
            k = rho/2
            print('k must be in the interval ['+ k +', inf)')
            print('instead using k =', k)

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
                k = rho/np.sqrt(2)
                print('k must be in the interval ['+ str(k) +', inf)')
                print('instead using k =', k)

            w = sqrt(k**2/(2*rho**2) + 0.25) + sqrt(k**2/(2*rho**2) - 0.25)*1j

        elif (g**2).real >= abs(g)**4 or abs(g) > 1:
            raise ValueError('matrix is not apportionable')

        elif abs(g) == 1:
            k = abs(evals[0]+evals[1])/2
            print('k =', k, 'is the only apportionment constant')
            w = 0

        else:
            k = abs(evals[0]+evals[1])/2 * np.sqrt(1 + (1 + abs(g)**4)/(2*(abs(g)**4 - (g**2).real)))
            print('k =', k, 'is the only apportionment constant')
            w = g * np.sqrt((1 - abs(g)**4)/(2*(abs(g)**4 - (g**2).real))) * 1j

        matP = evecs
        matPinv = np.linalg.inv(evecs)

        b = sqrt(w**2 - 1)/2
        matM = np.matmul(np.array([[1, b],[(w - 1)/(2 * b), (w + 1)/2]]),matPinv)
        matMinv= np.matmul(matP,np.array([[(w + 1)/2, -b],[(1 - w)/(2 * b), 1]]))

    matB = np.matmul(np.matmul(matM, mat), matMinv)

    return matB, matM

#
# Example
#
if __name__ == '__main__':

    # Making the displayed matrices more readable
    np.set_printoptions(precision=2, suppress=True)

    # Starting matrix and desired apportionment constant
    mat = np.array([[1,2],[0,-1]])
    k = 0.001

    print('Starting matrix: \n' + str(mat))
    print('Attempting to apportion matrix with apportionment constant k = ' + str(k))

    # Using apportionment function
    apport, matM = apport2x2(mat, k)

    print('Apportioned matrix: \n' + str(apport))
    print('Modulus of apportioned matrix: \n' + str(abs(apport)))
    print('Apportioning matrix: \n' + str(matM))
