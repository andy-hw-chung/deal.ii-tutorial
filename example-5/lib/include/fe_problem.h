//=================================
// include guard
#ifndef FE_PROBLEM_H
#define FE_PROBLEM_H
//=================================

#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/vector.h>
#include <string>
#include <fe_space.h>
#include <grid.h>

using namespace dealii;

/**************
*
*  FEProblem
*
**************/

/**
 * A class with the tools needed to solve an FE problem.
  */
template <int dim>
class FEProblem {
public:
    explicit FEProblem(FESpace<dim>* fe_space);

    virtual ~FEProblem() {};

    void init();

protected:
    // let these functions be pverridden by child classes

    /** Assemble the system matrix, solution and system rhs */
    virtual void assemble_system() = 0;

    /** Solve the matrix system */
    virtual void solve() = 0;

    /** Output the solution in vtu form with the given filename*/
    virtual void output_results(const std::string &filename) const = 0;

    // deal.ii objects needed to hold the solution and rhs
    FESpace<dim>* fe_space;
    SparseMatrix<double> system_matrix;
    Vector<double> solution;
    Vector<double> system_rhs;
};


#endif

