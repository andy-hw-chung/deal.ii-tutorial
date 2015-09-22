#include <deal.II/base/parameter_handler.h>
#include <string>
#include <parameters.h>
#include <shape.h>

/**************
*
*  Parameters
*
**************/


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

