#ifndef MHTESTFUNCTIONS_H
#define MHTESTFUNCTIONS_H
#include <array>
#include <numeric>
namespace MHTestFunctions {
  double Rosenbrock (double *x);
  std::array<double, 2> RosenbrockDerivative (double *x);

  double Himmelblau (double *x);
  std::array<double, 2> HimmelblauDerivative (double *x);

  double Parabola2D (double *x);
  std::array<double, 2> Parabola2DDerivative (double *x);

  template <unsigned dim> double Sphere (double *x) {
    return std::inner_product (x, x + dim, x, 0.0);
  }

  template <unsigned dim> std::array<double, dim> SphereDerivative (double *x) {
    std::array<double, dim> ret;
    for (size_t i = 0; i < dim; ++i) {
      ret[i] = 2.0 * x[i];
    }
    return ret;
  }
} // namespace MHTestFunctions
#endif // MHTESTFUNCTIONS_H
