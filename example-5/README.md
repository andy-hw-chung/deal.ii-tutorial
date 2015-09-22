# Example 5 - Separate compilation

This example shows how to split code up into more manageable parts. It uses the same code as in `example 4` and splits it up into separate units. The generic parts which can be re-used in other applications are placed in the `lib` directory and application specific code is placed in the `laplace` directory. In each of these directories, the code is further split into header files placed in the `include` directories and the source code in the `source` directory. 

The project structure is defined in the `CMakeLists.txt` files. The top level CMakeList defines the subdirectories
```
ADD_SUBDIRECTORY(lib)
ADD_SUBDIRECTORY(laplace)
```

The code in the `lib` directory is compiled and built into a library object. The `CMakeList` in the `lib` directory reads as follows:
```
INCLUDE_DIRECTORIES(include)

ADD_LIBRARY(libandy
  source/fe_problem.cpp
  source/fe_space.cpp
  source/grid.cpp
  source/grid_factory.cpp
  source/utility.cpp
  )


DEAL_II_SETUP_TARGET(libandy)
```
This points to the header files in the `include` directory for the includes in the source code. The source code making up the library is defined and then a target is configured with the name `libandy`. The resulting library will be `libandy.a`.

The `CMakeList` in the `laplace` directory reads:
```
INCLUDE_DIRECTORIES(
  include
  ${CMAKE_SOURCE_DIR}/lib/include
  )


ADD_EXECUTABLE(laplace source/main.cpp source/parameters.cpp source/application.cpp source/laplace_problem.cpp )
DEAL_II_SETUP_TARGET(laplace)
TARGET_LINK_LIBRARIES(laplace libandy)
```
This includes all the `include` directories for the headers. The source code making up the target application is then defined. Finally we link our user-defined library `libandy` to the application. 



