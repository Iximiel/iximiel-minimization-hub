#include "mhtestfunctions.hpp"

namespace MHTestFunctions {

double Rosenbrock(double *x) {
  return 1.0 * (x[0] - 1) * (x[0] - 1) +
         100.0 * (x[1] - x[0] * x[0]) * (x[1] - x[0] * x[0]);
}

double Himmelblau(double *x) {
  double a = (x[0] * x[0] + x[1] - 11);
  double b = (x[0] + x[1] * x[1] - 7);
  return a * a + b * b;
}

} // namespace MHTestFunctions
