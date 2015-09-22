//=================================
// include guard
#ifndef FE_SPACE_H
#define FE_SPACE_H
//=================================

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/fe/fe.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/lac/sparsity_pattern.h>
#include <memory>
#include <grid.h>

using namespace dealii;

/**************
*
*  FESpace
*
**************/

/**
  * A class which encapsulates a finite element space wih the given grid.
 * Handles the distibution of
  * the degrees of freedom and the corresponding sparsity pattern. Call init()
 * before using this class.
 */
template <int dim>
class FESpace {
public:
    FESpace(Grid<dim>* grid, FiniteElement<dim>* fe);

    ~FESpace();

    /** Initialise the dof handler and generte the sparsity pattern */
    void init();

    DoFHandler<dim> dof_handler;
    std::unique_ptr<FiniteElement<dim> > fe;
    SparsityPattern sparsity_pattern;

private:
    Grid<dim>* grid;
};


/**************
*
*  FE_Q_Space
*
**************/

/**
 *  A class deriving from FESpace representing the finite element space of
 * continuous, piecewise polynomials of degree p in each coordinate direction, where p is
 * given asscociated
 *  the second parameter in the constructor
 */
template <int dim>
class FE_Q_Space : public FESpace<dim> {
public:
    FE_Q_Space(Grid<dim>* grid, const int& polynomial_order);
};



#endif

