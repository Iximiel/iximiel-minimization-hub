#ifndef MHTESTFUNCTIONS_H
#define MHTESTFUNCTIONS_H
#include <array>
namespace MHTestFunctions {
  double Rosenbrock (double *x);
  std::array<double, 2> RosenbrockDerivative (double *x);

  double Himmelblau (double *x);
  std::array<double, 2> HimmelblauDerivative (double *x);

  double Parabola2D (double *x);
  std::array<double, 2> Parabola2DDerivative (double *x);
} // namespace MHTestFunctions
#endif // MHTESTFUNCTIONS_H
