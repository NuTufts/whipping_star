#include <iostream>

#include "GridDefinition.h"

int main( int nargs, char** argv )
{

  GridDefinition grid_def;
  grid_def.define_proper_bounds_course();

  std::cout << "dm2 bins: ";
  for ( auto& dm2 : grid_def.dm2_bin_v )
    std::cout << dm2 << " ";
  std::cout << std::endl;

  std::cout << "Ue4 bins: ";
  for ( auto& ue4 : grid_def.Ue4_bin_v )
    std::cout << ue4 << " ";
  std::cout << std::endl;

  std::cout << "Um4 bins: ";
  for ( auto& um4 : grid_def.Umu4_bin_v )
    std::cout << um4 << " ";
  std::cout << std::endl;
  
  return 0;
  
}
