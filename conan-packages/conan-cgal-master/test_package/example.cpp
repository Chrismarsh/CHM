#include "exact.h"
#include "array_convex_hull_2.h"
#include "vector_convex_hull_2.h"
#include <iostream>


int main()
{
  exact();
  array_convex_hull_2();
  vector_convex_hull_2();
  std::cout << "Pass" << std::endl;
  return 0;
}