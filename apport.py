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
    if np.allclose(mat, np.zeros((2, 2)), atol=1e-10):
        print('Input is the zero matrix; k = 0 is the only apportionment constant.')
        return np.zeros((2,2)), np.eye(2)
  
    evals, evecs = np.linalg.eig(mat)

    # Case 2: mat is nilpotent
    # See Theorem 3.2
    if np.isclose(evals[0], 0, atol=1e-10) and np.isclose(evals[1], 0, atol=1e-10):
        M = Matrix(mat)
        P = M.jordan_form()[0]

        matP = np.array(P).astype(complex)

        matM = np.array([[1, 1],[np.exp(2 * np.pi/3 * 1j), 1]]) @ np.diag([1,k * np.sqrt(3)]) @ matP
        matMinv = np.linalg.inv(matM)

    # Case 3: mat is rank 1 and not nilpotent
    elif np.isclose(evals[0], 0, atol=1e-10) or np.isclose(evals[1], 0, atol=1e-10):
        rho = max(np.abs(evals))
        r = k/rho

        if k < rho/2:
            print('Cannot apportion for k =', str(k))
            k = rho/2
            print('k must be in the interval ['+ str(k) +', inf)')
            print('instead using k =', str(k))

        matP = evecs
        matPinv = np.linalg.inv(evecs)

        matM = np.matmul(np.array([[1, 1/2 - np.sqrt(r**2 - 0.25)*1j], [1, -1/2 + np.sqrt(r**2 - 0.25)*1j]]), matPinv)
        matMinv = np.matmul(matP, np.array([[1/2 + np.sqrt(r**2 - 0.25)*1j, 1/2 - np.sqrt(r**2 - 0.25)*1j], [1, -1]]))

    # Case 4: mat is nonsingular
    # See Theorem 5.7
    else:
        # Ensure eigenvalues are distinct
        if np.isclose(evals[0], evals[1], atol=1e-10):
            print('Matrix is not apportionable.')
            return None, None

        g = (evals[0] + evals[1])/(evals[0] - evals[1])
        rho = max(np.abs(evals))

        # Case g = 0
        if np.isclose(g, 0, atol=1e-10):
            if k < 1/np.sqrt(2):
                print('Cannot apportion for k =', str(k))
                k = rho/np.sqrt(2)
                print('k must be in the interval ['+ str(k) +', inf)')
                print('instead using k =', str(k))

            w = sqrt(k**2/(2*rho**2) + 0.25) + sqrt(k**2/(2*rho**2) - 0.25)*1j

        # Case |g| > 1 or Re(g^2) >= |g|^4
        elif (g**2).real >= np.abs(g)**4 or np.abs(g) > 1:
            print('Matrix is not apportionable.')
            return None, None

        # Case |g| = 1
        elif np.isclose(np.abs(g), 1, atol=1e-10):
            k = np.abs(evals[0]+evals[1])/2
            print('Using k =', str(k), 'as this is the only apportionment constant.')
            w = 0

        else:
            k = np.abs(evals[0]+evals[1])/2 * np.sqrt(1 + (1 + np.abs(g)**4)/(2*(np.abs(g)**4 - (g**2).real)))
            print('Using k =', k, 'as this is the only apportionment constant.')
            w = g * np.sqrt((1 - np.abs(g)**4)/(2*(np.abs(g)**4 - (g**2).real))) * 1j

        matP = evecs
        matPinv = np.linalg.inv(evecs)

        b = sqrt(w**2 - 1)/2
        matM = np.matmul(np.array([[1, b],[(w - 1)/(2 * b), (w + 1)/2]]),matPinv)
        matMinv= np.matmul(matP,np.array([[(w + 1)/2, -b],[(1 - w)/(2 * b), 1]]))

    matB = matM @ mat @ matMinv

    return matB, matM


# Main execution block for examples
if __name__ == '__main__':
    # Making the displayed matrices more readable
    np.set_printoptions(precision=3, suppress=True)

    # Example 1:
    print('--- Exampe 1: Nonsingular and poor choice of k')
    mat = np.array([[1,2],[0,-1]])
    k = 0.001

    print('Starting matrix: \n' + str(mat))
    print('Attempting to apportion matrix with apportionment constant k = ' + str(k))

    apport, matM = apport2x2(mat, k)

    if not apport is None:
        print('Apportioned matrix: \n' + str(apport))
        print('Modulus of apportioned matrix: \n' + str(np.abs(apport)))
        print('Apportioning matrix: \n' + str(matM))

    # Example 2:
    print('--- Example 2: Not apportionable')
    mat = np.array([[1,0],[0,1]])
    k = 1

    print('Starting matrix: \n' + str(mat))
    print('Attempting to apportion matrix with apportionment constant k = ' + str(k))

    apport, matM = apport2x2(mat, k)

    if not apport is None:
        print('Apportioned matrix: \n' + str(apport))
        print('Modulus of apportioned matrix: \n' + str(np.abs(apport)))
        print('Apportioning matrix: \n' + str(matM))

    # Example 3:
    print('--- Example 3: Nilpotent matrix')
    mat = np.array([[1j,1],[1,-1j]])
    k = 5

    print('Starting matrix: \n' + str(mat))
    print('Attempting to apportion matrix with apportionment constant k = ' + str(k))

    apport, matM = apport2x2(mat, k)

    if not apport is None:
        print('Apportioned matrix: \n' + str(apport))
        print('Modulus of apportioned matrix: \n' + str(np.abs(apport)))
        print('Apportioning matrix: \n' + str(matM))
        
    # Example 4:
    print('--- Example 4: Rank 1 and not nilpotent')
    mat = np.array([[2j,0],[2,0]])
    k = 2

    print('Starting matrix: \n' + str(mat))
    print('Attempting to apportion matrix with apportionment constant k = ' + str(k))

    apport, matM = apport2x2(mat, k)

    if not apport is None:
        print('Apportioned matrix: \n' + str(apport))
        print('Modulus of apportioned matrix: \n' + str(np.abs(apport)))
        print('Apportioning matrix: \n' + str(matM))