#ifndef MHTESTFUNCTIONS_H
#define MHTESTFUNCTIONS_H
#include <array>
namespace MHTestFunctions {
  double Rosenbrock (double *x);
  double Himmelblau (double *x);
  std::array<double, 2> RosenbrockDerivative (double *x);
} // namespace MHTestFunctions
#endif // MHTESTFUNCTIONS_H
