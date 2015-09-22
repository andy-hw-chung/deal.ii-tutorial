#include <application.h>
#include <fe_space.h>
#include <grid.h>
#include <grid_factory.h>
#include <parameters.h>
#include <laplace_problem.h>
#include <shape.h>
#include <utility.h>


/**************
*
*  Application
*
**************/

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

template class Application<2>;
template class Application<3>;
