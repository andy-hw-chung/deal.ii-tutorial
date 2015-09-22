#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_renumbering.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/fe/fe.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/lac/sparsity_pattern.h>
#include <memory>
#include <fe_space.h>
#include <grid.h>


using namespace dealii;

/**************
*
*  FESpace
*
**************/

template <int dim>
FESpace<dim>::FESpace(Grid<dim>* grid, FiniteElement<dim>* fe):
	dof_handler(grid->triangulation), // associate the dof handler with the triangulation in the grid object
        fe(fe),
        grid(grid)
{
}

/** Make sure the dof handler does not still hold a pointer to the fe object */
template <int dim>
FESpace<dim>::~FESpace()
{
    this->dof_handler.clear();
}

template <int dim>
void FESpace<dim>::init()
{
    // distribute the degrees of freedoms according to the finite element
    // description supplied
    dof_handler.distribute_dofs(*fe);

    DoFRenumbering::Cuthill_McKee (dof_handler);

    // create the sparsity pattern
    DynamicSparsityPattern dsp(dof_handler.n_dofs());
    DoFTools::make_sparsity_pattern(dof_handler, dsp);

    sparsity_pattern.copy_from(dsp);
}


/**************
*
*  FE_Q_Space
*
**************/

template <int dim>
FE_Q_Space<dim>::FE_Q_Space(Grid<dim>* grid, const int& polynomial_order)
    : FESpace<dim>(grid, new FE_Q<dim>(polynomial_order))
{
}



template class FESpace<2>;
template class FESpace<3>;
template class FE_Q_Space<2>;
template class FE_Q_Space<3>;
