/**
 * Include all the necessary deal.II headers
 */
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/manifold_lib.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/grid_tools.h>
#include <iostream>
#include <fstream>
#include <limits>
#include <memory>
#include <string>

/** Use the deal.ii namespace for convenience */
using namespace dealii;

/** Print the sentence "`name` is `t`" */
template <typename T>
void print_is(const std::string name, const T t) {
  std::cout << name << " is " << t << std::endl;
}


/**************
*
*  Grid
*
**************/

/** A class containing a deal.ii triangulation and some asscociated utility
 * functions */
template <int dim>
class Grid {
public:
    /** The deal.ii traingulation object */
    Triangulation<dim> triangulation;

    /** Set the pointer to the manifold to the given one */
    void setManifold(Manifold<dim>* m, bool onlyOnBoundary = false);

    /** Output the given grid  */
    void output(const std::string& fileName);

    /** refine the mesh. Each cell is halved in size in each refinement iteration
   */
    void refineGlobal(const int& nRefinements);

    /** the minumum dimeter of all cells in the mesh */
    double h();

private:
    /** A pointer to a manifold object used by the triangulation */
    std::unique_ptr<Manifold<dim> > manifold;
};

template <int dim>
void Grid<dim>::setManifold(Manifold<dim>* m, bool onlyOnBoundary)
{
    // set the manifold id of each cell and face of the triangulation to the same
    // number to indicate
    // each cell and face share the same manifold description
    if (onlyOnBoundary)
        triangulation.set_all_manifold_ids_on_boundary(0);
    else
        triangulation.set_all_manifold_ids(0);

    // set each cell and face of the triangulation with manifold if 0 to have the
    // given manifold description
    triangulation.set_manifold(0, *m);
    // save the pointer to the manifold
    manifold = std::unique_ptr<Manifold<dim> > (m);
}

template <int dim>
void Grid<dim>::output(const std::string& fileName)
{
    std::ofstream out(fileName);
    GridOut gridOut;
    gridOut.write_vtu(triangulation, out);
    std::cout << "Grid written to " << fileName << std::endl;
    out.close();
}

template <int dim>
void Grid<dim>::refineGlobal(const int& nRefinements)
{
    triangulation.refine_global(nRefinements);
}

template <int dim>
double Grid<dim>::h()
{
    double min_diameter = std::numeric_limits<double>::max();
    // get a triangulation cell iterator
    typedef typename Triangulation<dim>::active_cell_iterator cellIterator;
    // loop through each cell of the triangulation
    for (cellIterator cell = triangulation.begin_active();
         cell != triangulation.end(); ++cell) {
        const double diameter = cell->diameter();
        if (diameter < min_diameter)
            min_diameter = diameter;
    }
    return min_diameter;
}

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
    Grid<dim> sphericalShellGrid(const Point<dim>& centre,
        const double& innerRadius,
        const double& outerRadius);

    /** Return a spherical grid */
    Grid<dim> sphereGrid(const Point<dim>& centre, const double& radius);
};

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




int main ()
{
  // Set the dimension of the grids here
  const int dim = 2;

  // Set the number of refinements to make
  const int nRefinements = 3;

  // Get a grid factory
  GridFactory<dim> gridFactory;

  // get a square grid
  const double left = 0.;
  const double right = 1.; 
  Grid<dim> cubeGrid = gridFactory.cubeGrid(left, right);
  cubeGrid.refineGlobal(nRefinements);
  cubeGrid.output("square-grid.vtu");

     print_is("Cube grid: h", cubeGrid.h());


  // create a cube triangulation object
  const Point<dim> centre;
  const double innerRadius = .5,
               outerRadius = 1.;


  Grid<dim> sphereGrid = gridFactory.sphereGrid(centre, outerRadius);
  sphereGrid.refineGlobal(nRefinements);
  sphereGrid.output("sphere-grid.vtu");
  
     print_is("Sphere grid: h", sphereGrid.h());


  Grid<dim> shellGrid = gridFactory.sphericalShellGrid(centre, innerRadius, outerRadius);
  shellGrid.refineGlobal(nRefinements);
  shellGrid.output("shell-grid.vtu");

     print_is("Shell grid: h", shellGrid.h());




  return 0;
}
