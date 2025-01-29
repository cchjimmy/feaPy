# FEA py
Simple Finite Element Analysis (FEA) code written in Python,
it is for 2D thin plate of arbitrary shape.
Subplot left to right is FEA with linear and quadratic triangular elements respectively.
This figure is created with ![test.py](test.py).
Currently only support these 2 elements.
However, additional elements can be easily added by specifying their element stiffness matrix (K) and shape function derivative matrix (B) implmentations.
![Gaussian quadrature](https://en.wikipedia.org/wiki/Gaussian_quadrature) is used for the intergral approximation in K for quadratic triangular element,
values for the sampling points and weights are taken from my lecture notes.

# Dependencies
- ![NumPy](https://numpy.org/)
- ![Matplotlib](https://matplotlib.org/)

# References
- ![Textbook of Finite Element Analysis, by P. Seshu](https://soaneemrana.com/onewebmedia/TEXT%20BOOKOF%20FINITE%20ELEMENT%20ANALYSIS%20BY%20P.%20SESHU%20(1).pdf), ch 5.4

![Figure_1.png](Figure_1.png)
