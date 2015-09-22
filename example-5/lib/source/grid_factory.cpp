#include <deal.II/base/point.h>
#include <deal.II/grid/grid_generator.h>
#include <grid.h>
#include <grid_factory.h>


using namespace dealii;

template <int dim>
Grid<dim> GridFactory<dim>::cubeGrid(const double& left, const double& right)
{
    Grid<dim> grid;
    // use the deal.ii GridGenerator class to fill the triangulation with a hyper
    // cube mesh
    GridGenerator::hyper_cube(grid.triangulation, left, right);
    return grid;
}


template <int dim>
Grid<dim> GridFactory<dim>::sphereGrid(const Point<dim>& centre,
    const double& radius)
{
    Grid<dim> grid;

    // use the deal.ii GridGenerator class to fill the triangulation with a hyper
    // shell mesh
    GridGenerator::hyper_ball(grid.triangulation, centre, radius);

    grid.setManifold(new SphericalManifold<dim>(centre), true);

    return grid;
}

template <int dim>
Grid<dim> GridFactory<dim>::sphericalShellGrid(const Point<dim>& centre,
    const double& innerRadius,
    const double& outerRadius)
{
    Grid<dim> grid;

    // use the deal.ii GridGenerator class to fill the triangulation with a hyper
    // shell mesh
    GridGenerator::hyper_shell(grid.triangulation, centre, innerRadius,
        outerRadius, 6);

    grid.setManifold(new SphericalManifold<dim>(centre));

    return grid;
}


template class GridFactory<2>;
template class GridFactory<3>;

