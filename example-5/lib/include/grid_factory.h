//=================================
// include guard
#ifndef GRID_FACTORY_H
#define GRID_FACTORY_H
//=================================

#include <deal.II/base/point.h>
#include <grid.h>

/**************
*
*  Grid Factory
*
**************/

/**
 * A class that generates Grid objects on request. The space dimension of the
 * grids
 * is controlled by the template parameter.
 */
template <int dim>
class GridFactory {
public:
    Grid<dim> cubeGrid(const double& left, const double& right);

    /** Return a spherical shell grid */
    Grid<dim> sphericalShellGrid(const Point<dim>& centre, const double& innerRadius, const double& outerRadius);

    /** Return a spherical grid */
    Grid<dim> sphereGrid(const Point<dim>& centre, const double& radius);
};


#endif

