#include <deal.II/base/function.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/precondition.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/vector.h>
#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/matrix_tools.h>
#include <deal.II/numerics/vector_tools.h>

#include <fstream>

#include <fe_problem.h>
#include <fe_space.h>
#include <laplace_problem.h>
#include <parameters.h>

using namespace dealii;

template <int dim>
LaplaceProblem<dim>::LaplaceProblem(FESpace<dim>* fe_space, Parameters *parameters):
	FEProblem<dim>::FEProblem(fe_space),
	parameters(parameters)
{}


template <int dim>
void LaplaceProblem<dim>::assemble_system()
{
    // get a quadrature formula object to perform quadrature with	 
    QGauss<dim> quadrature_formula(parameters -> n_gauss_quad);

    // the following object holds the values and gradients of the shape functions in each cell
    FEValues<dim> fe_values(*(this->fe_space->fe), quadrature_formula,
        update_values | update_gradients | update_JxW_values);

    // get the number of degrees of freedom per cell
    const unsigned int dofs_per_cell = this->fe_space->fe->dofs_per_cell;
    // gett the number of quadrature points for the quadrature formula
    const unsigned int n_q_points = quadrature_formula.size();

    // gt objects to hold the local contribution of each cell to the global system
    // matrix and the global rhs vector - override each of them at each iteration
    FullMatrix<double> cell_matrix(dofs_per_cell, dofs_per_cell);
    Vector<double> cell_rhs(dofs_per_cell);

    // an object representing a map from the local dof indices to the global dof indices
    std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);

    // get iterators to iterate over cells
    typename DoFHandler<dim>::active_cell_iterator
        cell = this->fe_space->dof_handler.begin_active(),
        endc = this->fe_space->dof_handler.end();

    // loop through each cell
    for (; cell != endc; ++cell) {
    	 // reinitialise the values of the shape functions for this cell
        fe_values.reinit(cell);
        // reinitialise the local contibutions
        cell_matrix = 0;
        cell_rhs = 0;
   
        /**
          * remember: we are building the matrix system corresponding to:
          * 
          *   \int_{Omega} \phi_i\phi_j = \int_{Omega} f(x)\phi_i    
          * 
          * with f(x) = 1.
          */    	

        // loop through each quadrature point
        for (unsigned int q_index = 0; q_index < n_q_points; ++q_index) {
            // loop through the i-th dof
            for (unsigned int i = 0; i < dofs_per_cell; ++i) {
            	
              // add the local contribution the rhs
              cell_rhs(i) += (fe_values.shape_value(i, q_index) * 1. * fe_values.JxW(q_index));
                // loop through the j-th dof

                for (unsigned int j = 0; j < dofs_per_cell; ++j) {
                    // add the local contrbution to the system mtrix
                    cell_matrix(i, j) += (fe_values.shape_grad(i, q_index) * fe_values.shape_grad(j, q_index) * fe_values.JxW(q_index));
                } // end for dof j   
        } // end for dof i
      } // end for quad point

        // get the map from local to global indices
        cell->get_dof_indices(local_dof_indices);
        // add the local contributions to the global objects
        for (unsigned int i = 0; i < dofs_per_cell; ++i) {
          // add the local rhs to the global rhs
          this->system_rhs(local_dof_indices[i]) += cell_rhs(i);

            for (unsigned int j = 0; j < dofs_per_cell; ++j) {
            	 // add the local matrix to the global matrix
                this->system_matrix.add(local_dof_indices[i], local_dof_indices[j],
                    cell_matrix(i, j));
             } // end for dof j
        } // end for dof i
    } // end for cells

    // apply zero boundary conditions
    std::map<types::global_dof_index, double> boundary_values;

    VectorTools::interpolate_boundary_values(
        this->fe_space->dof_handler, 0, ZeroFunction<dim>(), boundary_values);

    MatrixTools::apply_boundary_values(boundary_values, this->system_matrix,
        this->solution, this->system_rhs);
}

template <int dim>
void LaplaceProblem<dim>::solve()
{
    // set the iterative solver with the maximum number of iterations and the solver tolerance
    SolverControl solver_control(parameters -> matrix_solver_max_iterations, parameters -> matrix_solver_tolerance);
   // use a CG solver for a symmetrical matrix
    SolverCG<> solver(solver_control);
   // solve the system and put the solution into the solution vector. We do no preconditioning in this example -
   // using `PreconditionIdentity()` uses the identity operation 
    solver.solve(this->system_matrix, this->solution, this->system_rhs,
        PreconditionIdentity());
}

template <int dim>
void LaplaceProblem<dim>::output_results(const std::string &filename) const
{
    // get a dataout object which works with deal.ii objects to output solutions
    DataOut<dim> data_out;
    // attach the dof handler to the data out object
    data_out.attach_dof_handler(this -> fe_space -> dof_handler);
    // attach the solution vector
    data_out.add_data_vector(this -> solution, "solution");
    // process the objects we just gave it
    data_out.build_patches();
    // get an output file stream 
    std::ofstream output(parameters -> output_directory + "/" + filename);
    // write vtu using the file stream
    data_out.write_vtu(output);
    // close the filestream
    output.close();
}

template class LaplaceProblem<2>;
template class LaplaceProblem<3>;
