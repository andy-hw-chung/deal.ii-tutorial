//=================================
// include guard
#ifndef LAPLACE_PROBLEM_H
#define LAPLACE_PROBLEM_H
//=================================

#include <fe_problem.h>
#include <fe_space.h>
#include <parameters.h>

using namespace dealii;

/**************
*
*  LaplaceProblem
*
**************/

template <int dim>
class LaplaceProblem : public FEProblem<dim> {
public:
    LaplaceProblem(FESpace<dim>* fe_space, Parameters *parameters);
    using FEProblem<dim>::FEProblem;

    /** Assemble the system corresponding to the laplace equation
     *
     *  -laplace(u(x)) = f(x)
     * 
     * in the domain with f(x) = 1 and u(x) = 0 on the boundary
     */
    void assemble_system();

    void solve();

    void output_results(const std::string &filename) const;

    private:
    Parameters *parameters;
};



#endif

