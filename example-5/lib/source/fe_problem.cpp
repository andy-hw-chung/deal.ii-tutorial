#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/vector.h>
#include <string>
#include <fe_problem.h>
#include <fe_space.h>
#include <grid.h>

using namespace dealii;

/**************
*
*  FEProblem
*
**************/

template <int dim>
FEProblem<dim>::FEProblem(FESpace<dim>* fe_space):
fe_space(fe_space)
{}


template <int dim>
void FEProblem<dim>::init()
{
    fe_space->init();
    system_matrix.reinit(fe_space->sparsity_pattern);
    solution.reinit(fe_space->dof_handler.n_dofs());
    system_rhs.reinit(fe_space->dof_handler.n_dofs());
}

template class FEProblem<2>;
template class FEProblem<3>;
