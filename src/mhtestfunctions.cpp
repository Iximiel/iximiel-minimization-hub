#include "mhtestfunctions.hpp"

namespace MHTestFunctions {

  double Rosenbrock (double *x) {
    // minimum @ (1,1)
    return 1.0 * (x[0] - 1) * (x[0] - 1) +
           100.0 * (x[1] - x[0] * x[0]) * (x[1] - x[0] * x[0]);
  }

  std::array<double, 2> RosenbrockDerivative (double *x) {
    return {2.0 * (200.0 * x[0] * x[0] * x[0] - 200.0 * x[0] * x[1] + x[0] - 1),
            200.0 * (x[1] - x[0] * x[0])};
  }

  double Himmelblau (double *x) {
    // minimima
    // @ (3.0, 2.0)
    // @ (-2.805118, 3.131312)
    // @ (-3.779310, -3.283186)
    // @ (3.584428, -1.848126)
    double a = (x[0] * x[0] + x[1] - 11);
    double b = (x[0] + x[1] * x[1] - 7);
    return a * a + b * b;
  }

  std::array<double, 2> HimmelblauDerivative (double *x) {
    double a = (x[0] * x[0] + x[1] - 11);
    double b = (x[0] + x[1] * x[1] - 7);
    return {2.0 * (2.0 * x[0] * a + b), 2.0 * (a + 2 * x[1] * b)};
  }

  double Parabola2D (double *x) {
    // minimum @ (0,0)
    return x[0] * x[0] + x[1] * x[1];
  }
  std::array<double, 2> Parabola2DDerivative (double *x) {
    return {2.0 * x[0], 2.0 * x[1]};
  }

} // namespace MHTestFunctions
