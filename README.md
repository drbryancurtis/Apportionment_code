# Matrix Apportionment

This repository contains Python scripts for exploring **matrix apportionment**, a mathematical process that transforms a matrix into a new matrix with a specific property: the sum of the magnitudes of each row and column is a constant value.

## What is Matrix Apportionment?

Matrix apportionment transforms a given matrix **A** into a new matrix **B** with the property that the sum of the magnitudes of each row and column is equal to a constant, **k**. This transformation is achieved by finding a non-singular matrix **M** such that:

$B = MAM^{-1}$

The primary objective of the code is to find a suitable matrix **M** that satisfies this apportionment property for a chosen value of **k**. The scripts in this repository demonstrate this for various types of matrices.

---

## Repository Contents

* `apport.py`: This is the core script for the repository. The function `apport2x2(mat, k)` takes a 2x2 NumPy array and an optional constant **k** as input. It then returns the apportioned matrix **B** and the apportioning matrix **M**. The script can handle various matrix types, including those that are rank 0, nilpotent, rank 1 (non-nilpotent), or rank 2. An example of how to use the function is included at the end of the script.

* `Example3_1.py`: This script verifies the calculations from a specific example in a related research paper. It demonstrates the process of manually calculating and verifying the apportionment for a larger 5x5 matrix. The code defines the matrices **A** and **M**, then computes the apportioned matrix **B** using the formula $B = MAM^{-1}$.

* `figure.py`: This script generates a static plot to visualize the admissible eigenvalues for a 2x2 matrix. The plot shows the region of complex numbers that can be a second eigenvalue given a specific first eigenvalue. This visualization helps in understanding the conditions under which a matrix can be apportioned.

* `interactive_evals.py`: This script provides an interactive version of the `figure.py` plot. It uses sliders to allow the user to change the values of a known eigenvalue (real and imaginary parts). As the user adjusts the sliders, the plot dynamically updates to show the corresponding admissible region for the second eigenvalue. This allows for a more intuitive exploration of the mathematical constraints.

---

## Requirements

To run the scripts in this repository, you'll need **Python 3.3.1** and the following libraries. You can install them using pip:

```bash
pip install numpy sympy matplotlib