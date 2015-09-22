/**
 * Include all the necessary deal.II headers
 */
#include <deal.II/base/function.h>
#include <deal.II/base/parameter_handler.h>
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
#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>

/** Use the deal.ii namespace for convenience */
using namespace dealii;

/** Print the sentence "`name` is `t`" */
template <typename T>
void print_is(const std::string name, const T t)
{
    std::cout << name << " is " << t << std::endl;
}

enum class Shape {
cube,
sphere,
spherical_shell
};


void make_output_directory(const std::string output_directory)
{
 // check if output directory already exists - check with user before deleting
    if ( access( output_directory.c_str(), 0 ) == 0 )
   {
      struct stat status;
      stat( output_directory.c_str(), &status );

      if ( status.st_mode & S_IFDIR )
      {
        std::cout << "The directory \"" << output_directory
                   << "\" already exists. OK to delete directory and its contents? [y/n]" << std::endl;
		std::string go_ahead;		   
		 while (go_ahead != "y")
		 {
            std::cin >> go_ahead;
			if (go_ahead == "n")
			{
            std::printf("Goodbye!\n");	
            abort();			
			}
	    }
		
		std::string tmp= "rm -rf "+ output_directory;
		system(tmp.c_str());
		mkdir(output_directory.c_str(),0777);
      }
   }
   else
   {
      std::cout << "Creating output directory \"" << output_directory << "\"" << std::endl;
	  mkdir(output_directory.c_str(),0777);
   } 

}



/**************
*
*  Parameters
*
**************/


class Parameters {
	
public:
   Parameters(const std::string filename);
   void init();

   std::string output_directory;
   int         grid_dimension,
               grid_n_refinements;
   Shape       grid_shape;
   double      grid_cube_left,
               grid_cube_right,
               grid_sphere_radius,
               grid_shell_inner_radius,
               grid_shell_outer_radius;
   int         fe_order;
   int         n_gauss_quad,
               matrix_solver_max_iterations;
   double      matrix_solver_tolerance;

private:
  std::string filename;
  ParameterHandler  parameter_handler;
  void declare_parameters();
  void get_parameters();
  Shape stringtoShape(const std::string shape_name);
	
};


Parameters::Parameters(const std::string filename):
  filename(filename)
{}

void Parameters::init() {
  declare_parameters();
  parameter_handler.read_input(filename);
  get_parameters();	

}

void Parameters::declare_parameters() {
	
   parameter_handler.declare_entry ("output directory", "test",  Patterns::Anything(), "The name of the directry to output to");

   /*
    * declare grid parameters
    */
    parameter_handler.enter_subsection ("grid");
    {
        parameter_handler.declare_entry ("dimension", "2", Patterns::Integer(), "The dimension of the grid");
        parameter_handler.declare_entry ("shape", "cube", Patterns::Selection("cube|sphere|spherical shell"), "The dimension of the grid");
        parameter_handler.declare_entry ("number of refinements", "0", Patterns::Integer(), "The number of refinement levels to perform on the grid");

    parameter_handler.enter_subsection ("cube");
    {
    	   parameter_handler.declare_entry ("left point",  "0.", Patterns::Double(), "left point of cube interval");
    	   parameter_handler.declare_entry ("right point", "1.", Patterns::Double(), "right point of cube interval");
    } // end section cube
    parameter_handler.leave_subsection();
    

    parameter_handler.enter_subsection ("sphere");
    {
    	   parameter_handler.declare_entry ("radius",  "1.", Patterns::Double(), "radius of sphere");
    } // end section sphere
    parameter_handler.leave_subsection();

    parameter_handler.enter_subsection ("spherical shell");
    {
    	   parameter_handler.declare_entry ("inner radius",  ".5", Patterns::Double(), "inner radius of spherical shell");
    	   parameter_handler.declare_entry ("outer radius",  "1.", Patterns::Double(), "outer radius of spherical shell");
    } // end section spherical shell
    parameter_handler.leave_subsection();

  } // end section grid
    parameter_handler.leave_subsection();


   /**
    *  declare grid parameters
    */
    parameter_handler.enter_subsection ("finite elements");
  {
      parameter_handler.declare_entry ("polynomial order", "1", Patterns::Integer(), "The polynomial of the piecewise polynomial shape functions");
  	
  } // end section finite elements
  parameter_handler.leave_subsection();
  
   /**
    *  declare numerics parameters
    */
    parameter_handler.enter_subsection ("numerics");
  {
      parameter_handler.declare_entry ("number of gaussian quadrate points", "2", Patterns::Integer(),
                       "The number of quadrature points for gaussian quadrature");

    parameter_handler.enter_subsection ("matrix solver");
    {
      parameter_handler.declare_entry ("max number of iterations", "1000", Patterns::Integer(), "The maximum number of iterations for the matrix solver");
      parameter_handler.declare_entry ("tolerance", "1.e-10", Patterns::Double(), "The tolerance for the matrix solver");
    } // end section matrix solver
    parameter_handler.leave_subsection();
 	
  } // end section numerics
  parameter_handler.leave_subsection();

    	
}


Shape Parameters::stringtoShape(const std::string shape_name) {
	
  if(shape_name == "cube")
    return Shape::cube;
  else if(shape_name == "sphere")
    return Shape::cube;
  else if(shape_name == "spherical shell")
    return Shape::spherical_shell;
  else
   std::cout << "Cannot determine the grid shape from name: " << shape_name << std::endl;
   abort();
}

void Parameters::get_parameters() {
	
   output_directory = parameter_handler.get("output directory") ;
   
   parameter_handler.enter_subsection ("grid"); {
     grid_dimension     =  parameter_handler.get_integer("dimension") ;
     grid_shape         =  stringtoShape(parameter_handler.get("shape")) ;
     grid_n_refinements = parameter_handler.get_integer("number of refinements") ;

   parameter_handler.enter_subsection ("cube"); {
     grid_cube_left  = parameter_handler.get_double("left point") ;
     grid_cube_right = parameter_handler.get_double("right point") ;
   } // end section cube
   parameter_handler.leave_subsection();

   parameter_handler.enter_subsection ("sphere"); {
     grid_sphere_radius  = parameter_handler.get_double("radius") ;
   } // end section sphere
   parameter_handler.leave_subsection();

   parameter_handler.enter_subsection ("spherical shell"); {
     grid_shell_inner_radius  = parameter_handler.get_double("inner radius") ;
     grid_shell_outer_radius  = parameter_handler.get_double("outer radius") ;
   } // end section spherical shell
   parameter_handler.leave_subsection();

  } // end section grid
  parameter_handler.leave_subsection();

  parameter_handler.enter_subsection ("finite elements");
  {
      fe_order = parameter_handler.get_integer("polynomial order");
  	
  } // end section finite elements
  parameter_handler.leave_subsection();
  
   /**
    *  declare numerics parameters
    */
    parameter_handler.enter_subsection ("numerics");
  {
      n_gauss_quad = parameter_handler.get_integer("number of gaussian quadrate points");

    parameter_handler.enter_subsection ("matrix solver");
    {
      matrix_solver_max_iterations = parameter_handler.get_integer("max number of iterations");
      matrix_solver_tolerance = parameter_handler.get_double("tolerance");
    } // end section matrix solver
    parameter_handler.leave_subsection();
 	
  } // end section numerics
  parameter_handler.leave_subsection();
	
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

/**************
*
*  Application
*
**************/

template <int dim>
class Application {
public:
    Application(Parameters *parameters);
    void run();

private:
	Parameters *parameters;
	Grid<dim> get_grid();
};


template <int dim>
Application<dim>::Application(Parameters *parameters):
	parameters(parameters)
{}

template <int dim>
Grid<dim> Application<dim>::get_grid() {
   
  // Get a grid factory
  GridFactory<dim> gridFactory;	
  Shape shape = parameters -> grid_shape;

  if(shape == Shape::cube)
    return gridFactory.cubeGrid(parameters -> grid_cube_left, parameters -> grid_cube_right);

  else if(shape == Shape::sphere)
    return gridFactory.sphereGrid(Point<dim>(), parameters -> grid_sphere_radius);

  else if(shape == Shape::spherical_shell)
    return gridFactory.sphericalShellGrid(Point<dim>(), parameters -> grid_shell_inner_radius, parameters -> grid_shell_outer_radius);

  else
    abort();
	
}
	
template <int dim>
void Application<dim>::run()
{
    // get a grid
    // Grid<dim> grid = gridFactory.cubeGrid(-1., 1.);
    Grid<dim> grid = get_grid();
 
    grid.refineGlobal(parameters -> grid_n_refinements);
    grid.output(parameters -> output_directory + "/" + "grid.vtu");

    print_is("Grid: h", grid.h());

    // create the fe space
    FE_Q_Space<dim> fe_q_space(&grid, parameters -> fe_order);

    // create a laplace problem
    LaplaceProblem<dim> problem(&fe_q_space, parameters);

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
template <int dim>
void runApplication(Application<dim> &application) {
   // get a timer to time the application
   Timer timer;
   timer.start();
   application.run();
    // stop the timer once the application has finished
   timer.stop();

    std::cout << "Application complete: CPU time " << timer() << "s, Wall time " << timer.wall_time() << "s" << std::endl;

}


int main(int argc, char* argv[])
{
    deallog.depth_console(0);

    // get the parameter file from the command line arguments
    std::string parameters_file;
    std::string default_parameters_file = "parameters.param";
    
    if(argc < 1) {
    	std::cout  << "No parameter file passed through the command line. Using default parameter file: " << default_parameters_file << std::endl;
    	parameters_file = default_parameters_file;
    }
    else if(argc > 1) {
    	std::cout << "Too many parameters passed: received " << argc -1 << " instead of 1" << std::endl;
        abort();
    }
    else {
    	std::string file = argv[1];
    	std::cout << "Using the parameter file: " << file << std::endl;
        parameters_file = file;
    }

    Parameters parameters(parameters_file);
    parameters.init();

    // make the output directory
    make_output_directory(parameters.output_directory);    

    // Set the dimension of the grids here
    const int dim = parameters.grid_dimension;

    if(dim == 2) {
     Application<2> application(&parameters);
     runApplication(application);
    }
    else if(dim == 3) {
     Application<3> application(&parameters);
     runApplication(application);
    }




    return 0;
}
