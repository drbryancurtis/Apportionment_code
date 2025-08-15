import numpy as np
from sympy import Matrix
from cmath import sqrt

#
# Function for apportioning 2x2 matrices
#

def apport2x2(mat: np.ndarray, k: float = 1):
    """
    Args:
        mat (np.ndarray): The 2x2 complex matrix to be apportioned.
        k (float): The desired positive apportionment constant. Defaults to 1.

    Returns:
        tuple: (apportioned_matrix (np.ndarray), apportioning_matrix (np.ndarray))
               Returns (None, None) if the matrix is not apportionable.
               If mat is not apportionable for input k but is apportionable, then k will be changed to an apportionment constant

    Raises:
        TypeError: If mat is not a numpy array.
        ValueError: If mat is not 2x2, k is negative
    """
    # Check inputs are valid
    if not isinstance(mat, np.ndarray):
        raise TypeError("Input must be a numpy array.")
    if mat.shape != (2, 2):
        raise ValueError("Input array must have shape (2,2).")
    if k < 0:
        raise ValueError("k must be nonnegative.")

    # Convert matrix to complex dtype
    mat = mat.astype(complex)

    # Case 1: mat is all zeros matrix
    if not np.any(mat):
        print('Input is the zero matrix; k = 0 is the only apportionment constant.')
        return np.zeros((2,2)), np.eye(2)
  
    evals, evecs = np.linalg.eig(mat)

    # Case 2: mat is nilpotent
    if evals[0] == 0 and evals[1] == 0: 
        M = Matrix(mat) # convert to matrix - jordan_form() comes from sympy 
        P = M.jordan_form()[0]

        matP = np.array(P).astype(complex)
        matPinv = np.linalg.inv(matP)

        matM = np.matmul(np.array([[np.sqrt(k), 1],[np.sqrt(k), 1/np.sqrt(k) + 1]]), matP)
        matMinv = np.matmul(matPinv, np.array([[1/np.sqrt(k) + 1, -1],[-np.sqrt(k), np.sqrt(k)]]))

    # Case 3: mat is rank 1 and not nilpotent
    elif evals[0] == 0 or evals[1] == 0:
        rho = max(abs(evals))
        r = k/rho

        if k < rho/2:
            print('Cannot apportion for k =', k)
            k = rho/2
            print('k must be in the interval ['+ k +', inf)')
            print('instead using k =', k)

        matP = evecs
        matPinv = np.linalg.inv(evecs)

        matM = np.matmul(np.array([[1, 1/2 - np.sqrt(r**2 - 0.25)*1j], [1, -1/2 + np.sqrt(r**2 - 0.25)*1j]]), matPinv)
        matMinv = np.matmul(matP, np.array([[1/2 + np.sqrt(r**2 - 0.25)*1j, 1/2 - np.sqrt(r**2 - 0.25)*1j], [1, -1]]))

    # Case 4: mat is nonsingular
    # Follows Theorem 5.7
    else:
        if evals[0] == evals[1]:
            print('Matrix is not apportionable.')
            return None, None

        g = (evals[0] + evals[1])/(evals[0] - evals[1])
        rho = max(abs(evals))

        if g == 0:
            if k < 1/np.sqrt(2):
                print('Cannot apportion for k =', k)
                k = rho/np.sqrt(2)
                print('k must be in the interval ['+ str(k) +', inf)')
                print('instead using k =', k)

            w = sqrt(k**2/(2*rho**2) + 0.25) + sqrt(k**2/(2*rho**2) - 0.25)*1j

        elif (g**2).real >= abs(g)**4 or abs(g) > 1:
            print('Matrix is not apportionable.')
            return None, None

        elif abs(g) == 1:
            k = abs(evals[0]+evals[1])/2
            print('Using k =', k, 'as this is the only apportionment constant.')
            w = 0

        else:
            k = abs(evals[0]+evals[1])/2 * np.sqrt(1 + (1 + abs(g)**4)/(2*(abs(g)**4 - (g**2).real)))
            print('Using k =', k, 'as this is the only apportionment constant.')
            w = g * np.sqrt((1 - abs(g)**4)/(2*(abs(g)**4 - (g**2).real))) * 1j

        matP = evecs
        matPinv = np.linalg.inv(evecs)

        b = sqrt(w**2 - 1)/2
        matM = np.matmul(np.array([[1, b],[(w - 1)/(2 * b), (w + 1)/2]]),matPinv)
        matMinv= np.matmul(matP,np.array([[(w + 1)/2, -b],[(1 - w)/(2 * b), 1]]))

    matB = np.matmul(np.matmul(matM, mat), matMinv)

    return matB, matM


# Main execution block for examples
if __name__ == '__main__':
    # Making the displayed matrices more readable
    np.set_printoptions(precision=3, suppress=True)

    # Example 1:
    print('--- Exampe 1: Poor choice of k')
    mat = np.array([[1,2],[0,-1]])
    k = 0.001

    print('Starting matrix: \n' + str(mat))
    print('Attempting to apportion matrix with apportionment constant k = ' + str(k))

    apport, matM = apport2x2(mat, k)

    if not apport is None:
        print('Apportioned matrix: \n' + str(apport))
        print('Modulus of apportioned matrix: \n' + str(abs(apport)))
        print('Apportioning matrix: \n' + str(matM))

    # Example 2:
    print('--- Example 2: Not apportionable')
    mat = np.array([[1,0],[0,1]])
    k = 0.001

    print('Starting matrix: \n' + str(mat))
    print('Attempting to apportion matrix with apportionment constant k = ' + str(k))

    apport, matM = apport2x2(mat, k)

    if not apport is None:
        print('Apportioned matrix: \n' + str(apport))
        print('Modulus of apportioned matrix: \n' + str(abs(apport)))
        print('Apportioning matrix: \n' + str(matM))