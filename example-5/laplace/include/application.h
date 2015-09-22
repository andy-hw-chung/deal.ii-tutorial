//=================================
// include guard
#ifndef APPLICATION_H
#define APPLICATION_H
//=================================

#include <grid.h>
#include <parameters.h>

using namespace dealii;

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


#endif
