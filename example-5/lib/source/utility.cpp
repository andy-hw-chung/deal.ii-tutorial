#include <sys/stat.h>
#include <sys/types.h>
#include <iostream>
#include <string>
#include <unistd.h>
#include <utility.h>

/**
 * Create the output directory with the given name. If it already exists, ask the user whether it
 * is ok to overwrite the contents of the directoy first - note that if this happens, this will
 * block the application until an answer is given.  
 * 
 */
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
