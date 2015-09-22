//=================================
// include guard
#ifndef PARAMETERS_H
#define PARAMETERS_H
//=================================

#include <deal.II/base/parameter_handler.h>
#include <string>
#include <shape.h>
using namespace dealii;

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



#endif
