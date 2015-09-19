# Example 3 - Solving the Laplace equation

This examples shows how to solve the Laplace equation in `deal.II`. Recall the Laplace equation:
```
-div(grad(u(x))) = f(x)
```
where x lives in some domain. We choose `f(x)=1` and complete the equation with zero Dirichlet conditions on the boundary.

This example uses the `GridFactory` class from the example 2 to create `Grid`s. There are a number of abstract classes defined to help solve general problems:
	
1. `FESpace`: this class holds objects representing a finite element space. It holds a `Grid` object to represent the domain and a `FiniteElement` object to represent the function space. It handles the degrees of freedom using its `DoFHandler` object.

2. `FEProblem`: this class holds objects that are needed to actually define and solve the problem in the application. It holds objects correponding to the LHS matrix and the RHS and solutions vectors and also a reference to an `FESpace` object.

To solve the Laplace equation in this example, we use continuous, piecewise polynomials for the finite element space. This is enabled in the code by the `FE_Q_Space` class which is derived from `FESpace`. The assembly of the problem is defined in the `LaplaceProblem` class.

The application is defined in the `Application` class. It has a `run()` method which scripts the whle process of solving the Laplace equation.

