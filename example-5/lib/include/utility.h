//=================================
// include guard
#ifndef UTILITY_H
#define UTILITY_H
//=================================


#include <iostream>

/**
 * Print the sentence "`name` is `t`" 
 */
template <typename T>
inline void print_is(const std::string name, const T t)
{
    std::cout << name << " is " << t << std::endl;
}


void make_output_directory(const std::string output_directory);

#endif
