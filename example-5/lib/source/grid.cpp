#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/grid/manifold_lib.h>
#include <fstream>
#include <limits>
#include <memory>
#include <string>
#include <grid.h>

using namespace dealii;

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

template class Grid<2>;
template class Grid<3>;
