# Matrix Apportionment

This repository contains Python scripts for exploring **matrix apportionment**. Matrix apportionment transforms a given matrix $A$ into a new matrix $B$ with the property that the modulus of every entry is equal to a nonnegative real number. This transformation is achieved by finding a non-singular matrix $M$ such that:

$B = MAM^{-1}$

---

## Repository Contents

* `apport.py`: This script contains the function `apport2x2(mat, k)` which takes a $2\times 2$ NumPy array and an optional constant $k$ as input. It returns the apportioned matrix $B$ and the apportioning matrix $M$. An example of how to use the function is included at the end of the script.

* `Example3_1.py`: This script verifies the calculations from a specific example in a related research paper.

* `figure.py`: This script generates a static plot to visualize the admissible eigenvalues for a $2\times 2$ matrix. The plot shows the region of complex numbers that can be an eigenvalue given the other eigenvalue is 1.

* `interactive_evals.py`: This script provides an interactive version of the `figure.py` plot. It uses sliders to allow the user to change the values of a known eigenvalue (real and imaginary parts). As the user adjusts the sliders, the plot dynamically updates to show the corresponding admissible region for the second eigenvalue.

---

## Requirements

These scripts were written using the following:
* Python 3.13.3
* numpy 2.2.4
* sympy 1.13.3
* matplotlib 3.10.1