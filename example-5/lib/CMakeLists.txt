# lib/CMakeLists.txt
INCLUDE_DIRECTORIES(include)

ADD_LIBRARY(libandy
  source/fe_problem.cpp
  source/fe_space.cpp
  source/grid.cpp
  source/grid_factory.cpp
  source/utility.cpp
  )


DEAL_II_SETUP_TARGET(libandy)

