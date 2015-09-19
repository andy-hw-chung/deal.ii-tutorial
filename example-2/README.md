# Example 2 - Creating grids

This examples shows how to create triangulations in `deal.II`. In this application, we have wrapped the `Triangulation` class into a user-defined `Grid` class. The `Grid` class comes equipped with functions to perform helpful operations on the underlying `Triangulation` object it contains. If the `Triangulation` is set to have a corresponding `Manifold` object then the `Manifold` object needs to live as long as the `Triangulation` is still using it. Thus, also included in the `Grid` class is a smart pointer to a `Manifold` object to store the pointer provided in the `Grid::setManifold` method.

To create the `Grids`, there is the `GridFactory` class. This has methods to create a cube/square, sphere/circle or ring shaped domain. When this application runs, it creates one of each grid and outputs them in `vtu` format.
