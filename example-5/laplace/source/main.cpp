#include <deal.II/base/logstream.h>
#include <deal.II/base/timer.h>
#include <application.h>
#include <parameters.h>
#include <utility.h>

/** Use the deal.ii namespace for convenience */
using namespace dealii;

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

template class Application<2>;
template class Application<3>;
