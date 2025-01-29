# FEA py
Simple Finite Element Analysis (FEA) code written in Python,
it is for 2D thin plate of arbitrary shape.
Subplot left to right is FEA with linear and quadratic triangular elements respectively.
This figure is created with ![test.py](test.py).
Currently only support these 2 elements.
However, additional elements can be easily added by specifying their element stiffness matrix (K) and shape function derivative matrix (B) implementations.
Gaussian quadrature is used for the intergral approximation in K for quadratic triangular element,
values for the sampling points and weights are taken from my lecture notes.

When implementing K,
it is important to keep parameters accepted by K the same across different element types,
as it is called the same way within the evaluate function.

# Dependencies
- NumPy
- Matplotlib

# References
- Textbook of Finite Element Analysis, by P. Seshu, ch 5.4

![Figure_1.png](Figure_1.png)
