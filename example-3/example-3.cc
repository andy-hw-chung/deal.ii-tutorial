/**
 * Include all the necessary deal.II headers
 */
#include <deal.II/base/function.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/smartpointer.h>
#include <deal.II/base/timer.h>
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_renumbering.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/manifold_lib.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/precondition.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/vector.h>
#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/matrix_tools.h>
#include <deal.II/numerics/vector_tools.h>
#include <iostream>
#include <fstream>
#include <limits>
#include <memory>
#include <string>

/** Use the deal.ii namespace for convenience */
using namespace dealii;

/** Print the sentence "`name` is `t`" */
template <typename T>
void print_is(const std::string name, const T t)
{
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

template <int dim>
FESpace<dim>::FESpace(Grid<dim>* grid, FiniteElement<dim>* fe)
    : dof_handler(grid->triangulation)
    , // associate the dof handler with the
    // triangulation in the grid object
    fe(fe)
    , grid(grid)
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

template <int dim>
FE_Q_Space<dim>::FE_Q_Space(Grid<dim>* grid, const int& polynomial_order)
    : FESpace<dim>(grid, new FE_Q<dim>(polynomial_order))
{
}

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

/**************
*
*  LaplaceProblem
*
**************/

template <int dim>
class LaplaceProblem : public FEProblem<dim> {
public:
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
};

template <int dim>
void LaplaceProblem<dim>::assemble_system()
{
    // get a quadrature formula object to perform quadrature with	 
    QGauss<dim> quadrature_formula(2);

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
    SolverControl solver_control(1000, 1e-12);
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
    std::ofstream output(filename);
    // write vtu using the file stream
    data_out.write_vtu(output);
    // close the filestream
    output.close();
}

/**************
*
*  Application
*
**************/

template <int dim>
class Application {
public:
    void run();
};

template <int dim>
void Application<dim>::run()
{
    // Set the number of refinements to make
    const int nRefinements = 3;

    // order of the fe_q
    const int polynomial_order = 1;

    // Get a grid factory
    GridFactory<dim> gridFactory;

    // get a grid
    // Grid<dim> grid = gridFactory.cubeGrid(-1., 1.);
    Grid<dim> grid = gridFactory.sphereGrid(Point<dim>(), 1. );
    grid.refineGlobal(nRefinements);
    grid.output("grid.vtu");

    print_is("Grid: h", grid.h());

    // create the fe space
    FE_Q_Space<dim> fe_q_space(&grid, polynomial_order);

    // create a laplace problem
    LaplaceProblem<dim> problem(&fe_q_space);

    problem.init();

    // assemble the system
    problem.assemble_system();
    // solve it
    problem.solve();
    // output the solution
    problem.output_results("solution.vtu");
}

/**************
*
*  main
*
**************/

int main()
{
    deallog.depth_console(0);

    // get a timer to time the application
    Timer  timer;

    // Set the dimension of the grids here
    const int dim = 3;

    Application<dim> application;

    // start the timer before running the application
    timer.start();
    application.run();
    // stop the timer once the application has finished
    timer.stop();

    std::cout << "Application complete: CPU time " << timer() << "s, Wall time " << timer.wall_time() << "s" << std::endl;

    return 0;
}
