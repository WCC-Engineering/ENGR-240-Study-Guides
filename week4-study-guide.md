# Week 4 Study Guide: Linear Systems of Equations

## Overview

This study guide covers numerical methods for solving systems of linear equations, a fundamental problem in engineering computation. The activities, worksheets, and programming assignments for Week 4 introduce matrix notation, Gauss elimination, LU factorization, built-in NumPy/SciPy functions for solving linear systems, and the concept of condition number.

## 1. Matrix Notation for Linear Systems

### Key Concepts

- **Systems of Linear Equations**: A set of equations where each variable appears at most to the first power
- **Matrix Notation**: A compact way to represent systems using matrices and vectors
- **Standard Form**: Ax = b, where A is the coefficient matrix, x is the vector of unknowns, and b is the right-hand side vector
- **Solution Types**:
  - Unique solution: Exactly one value for each unknown
  - Infinite solutions: Underdetermined system
  - No solution: Overdetermined or inconsistent system

### Matrix Representation of Linear Systems

A system of linear equations:

$$
\begin{align}
a_{11}x_1 + a_{12}x_2 + \ldots + a_{1n}x_n &= b_1\\
a_{21}x_1 + a_{22}x_2 + \ldots + a_{2n}x_n &= b_2\\
\vdots\\
a_{m1}x_1 + a_{m2}x_2 + \ldots + a_{mn}x_n &= b_m
\end{align}
$$

Can be represented in matrix form as:

$$
\begin{bmatrix}
a_{11} & a_{12} & \cdots & a_{1n} \\
a_{21} & a_{22} & \cdots & a_{2n} \\
\vdots & \vdots & \ddots & \vdots \\
a_{m1} & a_{m2} & \cdots & a_{mn}
\end{bmatrix}
\begin{bmatrix}
x_1 \\
x_2 \\
\vdots \\
x_n
\end{bmatrix}
=
\begin{bmatrix}
b_1 \\
b_2 \\
\vdots \\
b_m
\end{bmatrix}
$$

Or simply:

$$\textbf{A}\textbf{x} = \textbf{b}$$

### Properties of Linear Systems

- **Square System**: Number of equations equals number of unknowns (m = n)
- **Overdetermined System**: More equations than unknowns (m > n)
- **Underdetermined System**: Fewer equations than unknowns (m < n)
- **Homogeneous System**: Right-hand side is zero (b = 0)
- **Non-homogeneous System**: Right-hand side is non-zero (b ≠ 0)

### Existence and Uniqueness of Solutions

- A square system (n equations, n unknowns) has a unique solution if and only if the determinant of the coefficient matrix is non-zero (det(A) ≠ 0)
- An equivalent condition is that the matrix A must be invertible (non-singular)
- When det(A) = 0, the matrix is singular and either:
  - The system has no solution (inconsistent)
  - The system has infinitely many solutions (dependent equations)

### Engineering Applications of Linear Systems

- Structural analysis (force distributions in trusses)
- Circuit analysis (Kirchhoff's laws)
- Heat transfer calculations
- Material balance equations in chemical engineering
- Traffic flow modeling
- Mechanical equilibrium problems

### Resources
- **Davishahl Numerical Methods Videos**
  - [Intro to Linear Systems](https://youtu.be/yj0diLlHwgs?si=oYAkEv14GQwlWqEJ)
  - [Intro to Matrix Algebra](https://youtu.be/rFNg9MQPau8?si=7__kTM7Tdnapmi9f)

- **Berkeley Numerical Methods**
   - [Basics of Linear Algebra](https://pythonnumericalmethods.studentorg.berkeley.edu/notebooks/chapter14.01-Basics-of-Linear-Algebra.html)
   - [Systems of Linear Equations](https://pythonnumericalmethods.studentorg.berkeley.edu/notebooks/chapter14.03-Systems-of-Linear-Equations.html)
     
- **Sullivan Numerical Methods Textbook**:
  - [Linear Algebra](https://numericalmethodssullivan.github.io/ch-linearalgebra.html#intro-to-numerical-linear-algebra)

## 2. Gauss Elimination

### Key Concepts

- **Principle**: Systematic elimination of variables to convert the system to an upper triangular form
- **Forward Elimination**: The process of obtaining an upper triangular matrix
- **Back Substitution**: The process of finding solutions starting from the last equation
- **Pivoting**: Techniques to improve numerical stability
  - Partial pivoting: Row exchanges
  - Complete pivoting: Row and column exchanges
- **Elementary Row Operations**:
  - Multiply a row by a non-zero constant
  - Add a multiple of one row to another
  - Interchange two rows

### Gauss Elimination Algorithm

1. **Forward Elimination Phase**:
   - For each pivot position (diagonal element):
     - Select the pivot (optionally implement pivoting strategies)
     - Eliminate all elements below the pivot by subtracting multiples of the pivot row

2. **Back Substitution Phase**:
   - Start with the last equation (now with only one unknown)
   - Solve for that unknown
   - Substitute the value into the previous equation and solve
   - Continue until all unknowns are found

### Mathematical Representation

For a 3×3 system:

$$
\begin{bmatrix}
a_{11} & a_{12} & a_{13} \\
a_{21} & a_{22} & a_{23} \\
a_{31} & a_{32} & a_{33}
\end{bmatrix}
\begin{bmatrix}
x_1 \\
x_2 \\
x_3
\end{bmatrix}
=
\begin{bmatrix}
b_1 \\
b_2 \\
b_3
\end{bmatrix}
$$

**Forward Elimination**:

First, eliminate $x_1$ from the second and third equations:
- Compute multipliers: $m_{21} = a_{21}/a_{11}$ and $m_{31} = a_{31}/a_{11}$
- Update second row: $[a_{21}~a_{22}~a_{23}~b_2] = [a_{21}~a_{22}~a_{23}~b_2] - m_{21}[a_{11}~a_{12}~a_{13}~b_1]$
- Update third row: $[a_{31}~a_{32}~a_{33}~b_3] = [a_{31}~a_{32}~a_{33}~b_3] - m_{31}[a_{11}~a_{12}~a_{13}~b_1]$

Then, eliminate $x_2$ from the third equation:
- Compute multiplier: $m_{32} = a_{32}/a_{22}$
- Update third row: $[a_{32}~a_{33}~b_3] = [a_{32}~a_{33}~b_3] - m_{32}[a_{22}~a_{23}~b_2]$

Resulting in upper triangular form:

$$
\begin{bmatrix}
a_{11} & a_{12} & a_{13} \\
0 & a_{22}' & a_{23}' \\
0 & 0 & a_{33}'
\end{bmatrix}
\begin{bmatrix}
x_1 \\
x_2 \\
x_3
\end{bmatrix}
=
\begin{bmatrix}
b_1 \\
b_2' \\
b_3'
\end{bmatrix}
$$

**Back Substitution**:

Solve equations in reverse order:
- $x_3 = b_3'/a_{33}'$
- $x_2 = (b_2' - a_{23}'x_3)/a_{22}'$
- $x_1 = (b_1 - a_{12}x_2 - a_{13}x_3)/a_{11}$

### Operation Counting

For an n×n system, Gauss elimination requires:
- Forward elimination: approximately $\frac{n^3}{3} + \frac{n^2}{2} - \frac{5n}{6}$ operations
- Back substitution: approximately $n^2$ operations
- Total operations: $O(n^3)$ (cubic complexity)

This operation count demonstrates why solving large linear systems is computationally intensive.

### Pivoting for Numerical Stability

When the pivot element is zero or very small, numerical errors can be magnified. Pivoting strategies include:

**Partial Pivoting**:
- Search for the largest absolute value in the current column below the diagonal
- Swap rows to place this element in the pivot position
- Reduces numerical errors and prevents division by zero

**Complete Pivoting**:
- Search for the largest absolute value in the remaining submatrix
- Swap both rows and columns to place this element in the pivot position
- Maximum numerical stability but more computationally expensive

### Resources

- **Davishahl Numerical Methods Videos**:
  - [Gauss Elimination](https://youtu.be/n7XYMsypvA0?si=50-qAFvBM8v9MZ5Z)
  - [Gauss Elimination Efficiency](https://youtu.be/5MgBn5WyaUU?si=M8JP7CC3JlxBYrLt)
 
- **Berkeley Numerical Methods**
  - [Solutions to Systems of Linear Equations] (https://pythonnumericalmethods.studentorg.berkeley.edu/notebooks/chapter14.04-Solutions-to-Systems-of-Linear-Equations.html)

## 3. LU Factorization

### Key Concepts

- **Principle**: Decompose the coefficient matrix A into a product of lower (L) and upper (U) triangular matrices
- **Form**: A = LU
- **Advantages**:
  - Efficient for solving multiple systems with the same coefficient matrix but different right-hand sides
  - Simpler operations for computing determinants and inverses
  - Reduced computational cost for repeated solutions
- **Variants**:
  - LU decomposition without pivoting
  - LU decomposition with partial pivoting (LUP factorization)
  - Cholesky decomposition for symmetric positive-definite matrices

### LU Factorization Process

The LU factorization decomposes a matrix A into:

$$A = LU$$

Where:
- L is a lower triangular matrix with ones on the diagonal
- U is an upper triangular matrix

For a 3×3 system:

$$
\begin{bmatrix}
a_{11} & a_{12} & a_{13} \\
a_{21} & a_{22} & a_{23} \\
a_{31} & a_{32} & a_{33}
\end{bmatrix}
=
\begin{bmatrix}
1 & 0 & 0 \\
l_{21} & 1 & 0 \\
l_{31} & l_{32} & 1
\end{bmatrix}
\begin{bmatrix}
u_{11} & u_{12} & u_{13} \\
0 & u_{22} & u_{23} \\
0 & 0 & u_{33}
\end{bmatrix}
$$

### Computing LU Factorization

The LU factorization can be computed using several methods:

1. **Doolittle Algorithm**:
   - U is the same as the upper triangular matrix obtained from Gauss elimination
   - L contains the multipliers used during elimination

2. **Crout Algorithm**:
   - Different algorithm leading to the same result
   - More numerically stable in some cases

3. **LU with Partial Pivoting**:
   - Adds row interchanges with a permutation matrix P
   - Results in PA = LU

### Solving Systems with LU Factorization

Once A = LU factorization is available, solving Ax = b becomes:
1. LUx = b
2. Let Ux = y, then Ly = b

The solution process:
1. Forward substitution: Solve Ly = b for y
2. Back substitution: Solve Ux = y for x

This is particularly efficient when solving multiple systems with the same A but different b:
- Factorize A = LU once: $O(n^3)$ operations
- Solve Ly = b and Ux = y for each new b: $O(n^2)$ operations per system

### Advantages of LU Factorization

1. **Computational Efficiency**: For multiple right-hand sides (b vectors)
   - Single factorization: $O(n^3)$ operations
   - Each additional solution: $O(n^2)$ operations

2. **Matrix Properties**:
   - Determinant calculation: det(A) = det(L) × det(U) = product of diagonal elements of U
   - Inverse calculation: More efficient through LU factorization

3. **Stability Analysis**:
   - Error analysis is more straightforward
   - Pivoting strategies can be incorporated

### Resources

- **Davishahl Numerical Methods Videos**:
  - [LU Factorization](https://youtu.be/8_0Xo9pwY14?si=vFa5c8KVIQCBKUOf) 

## 4. NumPy/SciPy Functions for Linear Systems

### Key Concepts

- **NumPy Functions**: np.linalg module provides core linear algebra operations
- **SciPy Functions**: scipy.linalg module offers additional specialized functions
- **Key Operations**:
  - Direct solving of linear systems
  - Computing matrix inverse
  - LU decomposition
  - Analysis of matrix properties (determinant, rank, etc.)
- **Advantages**: Optimized implementations, vectorized operations, numerical stability

### NumPy Functions for Linear Systems

NumPy's `numpy.linalg` module provides essential functions for linear algebra operations:

```python
import numpy as np

# Define coefficient matrix and right-hand side
A = np.array([[4, 2, 3], [3, 1, -2], [2, -3, 1]])
b = np.array([7, -1, 2])

# Solve the system using numpy.linalg.solve
x = np.linalg.solve(A, b)
print("Solution:", x)

# Check the solution
print("Verification:", np.allclose(np.dot(A, x), b))

# Compute the inverse of A
A_inv = np.linalg.inv(A)
print("Matrix inverse:\n", A_inv)

# Compute the determinant
det_A = np.linalg.det(A)
print("Determinant:", det_A)

# Compute the rank
rank_A = np.linalg.matrix_rank(A)
print("Rank:", rank_A)
```

### SciPy Functions for Linear Systems

SciPy's `scipy.linalg` module extends NumPy's capabilities with more specialized functions:

```python
import numpy as np
import scipy.linalg as spla

# Define coefficient matrix and right-hand side
A = np.array([[4, 2, 3], [3, 1, -2], [2, -3, 1]])
b = np.array([7, -1, 2])

# LU decomposition with SciPy
P, L, U = spla.lu(A)
print("L matrix:\n", L)
print("U matrix:\n", U)
print("P matrix:\n", P)

# Verify A = P.T @ L @ U
print("Verification:", np.allclose(A, P.T @ L @ U))

# Solve using LU factorization
# Forward substitution
y = spla.solve_triangular(L, P @ b, lower=True)
# Back substitution
x = spla.solve_triangular(U, y, lower=False)
print("Solution via LU:", x)

# More efficient diagonal solving for triangular systems
y = spla.solve_triangular(L, P @ b, lower=True, unit_diagonal=True)
```

### SciPy's Specialized Solvers

SciPy provides specialized solvers for different matrix structures:

{% raw %}
```python
import numpy as np
import scipy.linalg as spla

# For symmetric positive definite matrices
A_spd = np.array([[4, 1, 1], [1, 3, 1], [1, 1, 3]])
b = np.array([6, 5, 5])

# Cholesky decomposition
L = spla.cholesky(A_spd, lower=True)
print("Cholesky factor:\n", L)

# Solve using Cholesky decomposition
x = spla.cho_solve((L, True), b)
print("Solution via Cholesky:", x)

# For banded matrices
# Tridiagonal system with diagonals (-1, 0, 1)
A_banded = np.array([[2, -1, 0, 0], 
                     [-1, 2, -1, 0],
                     [0, -1, 2, -1],
                     [0, 0, -1, 2]])
b = np.array([1, 0, 0, 1])

# Solve banded system
ab = np.array([[-1, -1, -1, 0], [2, 2, 2, 2], [-1, -1, -1, 0]])  # Banded storage
x = spla.solve_banded((1, 1), ab, b)
print("Solution for banded system:", x)
```
{% endraw %}

### Efficiency and Best Practices

1. **Choose the Right Function**:
   - Use `np.linalg.solve` instead of computing `A_inv @ b` for better numerical stability
   - For specialized matrix structures, use the corresponding SciPy functions
   - When solving with multiple right-hand sides, factorize once and reuse

2. **Sparse Matrices**:
   - For large, sparse systems, use `scipy.sparse` and its solvers
   - Tremendous savings in memory and computation time

3. **Error Checking**:
   - Verify solutions with `np.allclose(np.dot(A, x), b)`
   - Check matrix properties before applying specialized algorithms

### Resources

- **Berkeley Numerical Methods**
  - [Solve Systems of Linear Equations in Python](https://pythonnumericalmethods.studentorg.berkeley.edu/notebooks/chapter14.05-Solve-Systems-of-Linear-Equations-in-Python.html)
    
- **NumPy and SciPy Documentation**:
  - [NumPy Linear Algebra](https://numpy.org/doc/stable/reference/routines.linalg.html)
  - [SciPy Linear Algebra](https://docs.scipy.org/doc/scipy/reference/linalg.html)

## 5. Matrix Condition Number

### Key Concepts

- **Definition**: A measure of how sensitive the solution of a linear system is to small changes in the input data
- **Formula**: The condition number of a matrix A is κ(A) = ||A|| · ||A⁻¹||
- **Interpretation**:
  - κ(A) ≈ 1: Well-conditioned matrix (stable solutions)
  - κ(A) >> 1: Ill-conditioned matrix (small input changes can cause large output changes)
  - κ(A) = ∞: Singular matrix (no unique solution)
- **Relation to numerical error**: Error in solution ≤ κ(A) × (error in input data)

### Understanding Condition Number

The condition number quantifies how errors in the input data affect the solution of a linear system. A high condition number indicates that small changes in the coefficient matrix A or the right-hand side b can lead to large changes in the solution x.

For a system Ax = b:

- If relative error in b is εᵦ, then relative error in x can be up to κ(A) × εᵦ
- If relative error in A is εₐ, then relative error in x can be up to κ(A) × εₐ

### Computing the Condition Number

In practice, the condition number is computed using the ratio of the largest to the smallest singular value:

```python
import numpy as np

# Define a matrix
A = np.array([[4, 2, 3], [3, 1, -2], [2, -3, 1]])

# Compute condition number using 2-norm
cond_num = np.linalg.cond(A)
print("Condition number:", cond_num)

# Compute singular values
singular_values = np.linalg.svd(A, compute_uv=False)
print("Singular values:", singular_values)
print("Ratio of largest to smallest singular value:", singular_values[0]/singular_values[-1])
```

### Effects of Ill-Conditioning

Examples of ill-conditioned matrices:

1. **Nearly Dependent Rows/Columns**:
   ```python
   # Nearly linearly dependent rows
   A_ill = np.array([[1, 2], [0.999, 1.999]])
   print("Condition number:", np.linalg.cond(A_ill))
   ```

2. **Hilbert Matrix** (notoriously ill-conditioned):
   ```python
   # Create a 5x5 Hilbert matrix
   n = 5
   H = np.zeros((n, n))
   for i in range(n):
       for j in range(n):
           H[i, j] = 1 / (i + j + 1)
   print("Hilbert matrix condition number:", np.linalg.cond(H))
   ```

### Improving Conditioning

Techniques to handle ill-conditioned systems:

1. **Scaling**: Normalize rows or columns to similar magnitudes
   ```python
   # Row scaling example
   D = np.diag(1 / np.sqrt(np.sum(A**2, axis=1)))
   A_scaled = D @ A
   print("Original condition number:", np.linalg.cond(A))
   print("Scaled condition number:", np.linalg.cond(A_scaled))
   ```

2. **Regularization**: Add a small perturbation to make the matrix better conditioned
   ```python
   # Tikhonov regularization
   lambda_param = 0.01
   A_reg = A.T @ A + lambda_param * np.eye(A.shape[1])
   b_reg = A.T @ b
   x_reg = np.linalg.solve(A_reg, b_reg)
   ```

3. **Singular Value Decomposition (SVD)**: Solve using SVD and truncate small singular values
   ```python
   U, s, Vh = np.linalg.svd(A, full_matrices=False)
   # Truncate small singular values
   threshold = 1e-10 * s[0]
   s_inv = np.array([1/s_val if s_val > threshold else 0 for s_val in s])
   x_svd = Vh.T @ (s_inv * (U.T @ b))
   ```

### Resources

- **Davishahl Numerical Methods Videos**:
  - [Floating Point Precision](https://youtu.be/iDRkUGTHNXw?si=ewDWuqiNGS1wyjg9)
  - [Gauss Elimination Accuracy](https://youtu.be/qOlr__85-o0?si=grkefh43G-fNgdV4)

## 6. Best Practices for Linear System Solving

### Key Considerations

- **Algorithm Selection**:
  - For standard systems: NumPy's `linalg.solve`
  - For multiple right-hand sides: LU factorization
  - For symmetric positive-definite matrices: Cholesky decomposition
  - For large sparse systems: Iterative methods from `scipy.sparse.linalg`

- **Numerical Stability**:
  - Always check the condition number before solving
  - Use pivoting strategies for direct methods
  - Consider scaling for poorly conditioned systems
  - Validate solutions by substituting back

- **Computational Efficiency**:
  - Exploit matrix structure when possible
  - Avoid explicit inverse calculation
  - Reuse factorizations for multiple right-hand sides
  - Consider sparse storage for large systems

### Common Pitfalls and How to Avoid Them

1. **Using Inverse for Solving**:
   - Avoid `np.linalg.inv(A) @ b`
   - Instead use `np.linalg.solve(A, b)` for better numerical stability and efficiency

2. **Ignoring Ill-Conditioning**:
   - Always check condition number: `np.linalg.cond(A)`
   - Be cautious when condition number > 10⁶
   - Use regularization techniques for ill-conditioned systems

3. **Overlooking Matrix Structure**:
   - Take advantage of specialized solvers for banded, triangular, or symmetric matrices
   - Use sparse matrices when appropriate

4. **Inefficient Implementation**:
   - Avoid Python loops for matrix operations
   - Use vectorized operations
   - For very large systems, consider compiled solutions (Fortran/C++ via SciPy)

## Tips for Week 4 Assignments

1. **First, check the condition number** of any system you're working with to understand potential numerical issues.

2. **Implement at least one solution algorithm from scratch** (Gauss elimination or LU factorization) to understand the underlying process.

3. **Compare your implementations with NumPy/SciPy functions** to validate correctness and benchmark performance.

4. **For multiple right-hand sides**, implement LU factorization to demonstrate the efficiency advantage.

5. **Experiment with ill-conditioned matrices** to observe how perturbations in input affect the solution.

6. **Use visualization** to help understand the systems you're solving, especially for engineering applications.

7. **Document your code thoroughly** with comments explaining the mathematical steps.

8. **Implement error checking** to validate solutions by substituting back into the original equations.
