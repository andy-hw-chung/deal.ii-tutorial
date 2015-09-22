# laplace/CMakeLists.txt
# cleanup some files

INCLUDE_DIRECTORIES(
  include
  ${CMAKE_SOURCE_DIR}/lib/include
  )


ADD_EXECUTABLE(laplace source/main.cpp source/parameters.cpp source/application.cpp source/laplace_problem.cpp )
DEAL_II_SETUP_TARGET(laplace)
TARGET_LINK_LIBRARIES(laplace libandy)



