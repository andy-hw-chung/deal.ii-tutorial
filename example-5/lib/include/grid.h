//=================================
// include guard
#ifndef GRID_H
#define GRID_H
//=================================

#include <deal.II/grid/tria.h>
#include <deal.II/grid/manifold_lib.h>
#include <memory>
#include <string>

using namespace dealii;

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


#endif

