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

  template <unsigned Nat> double LennardJones12_6 (double *x) {
    constexpr unsigned dim = 3 * Nat;
    constexpr double sigma = 1.0;
    constexpr double epsilon = 1.0;
    constexpr double A = 4.0 * epsilon * sigma * sigma * sigma * sigma * sigma *
                         sigma * sigma * sigma * sigma * sigma * sigma * sigma;
    constexpr double B =
        4.0 * epsilon * sigma * sigma * sigma * sigma * sigma * sigma;
    double dx, dy, dz, d;
    double E = 0;
    for (unsigned i = 0; i < dim; i += 3) {
      for (unsigned j = i + 3; j < dim; j += 3) {
        dx = x[i] - x[j];
        dy = x[i + 1] - x[j + 1];
        dz = x[i + 2] - x[j + 2];
        d = dx * dx + dy * dy + dz * dz;
        d = d * d * d; // r^6
        E += A / (d * d) - B / d;
      }
    }
    // 2.0*E because we are using a triangular matrix
    return 2.0 * E;
  }

  template <unsigned Nat>
  std::array<double, 3 * Nat> LennardJones12_6Derivative (double *x) {
    constexpr unsigned dim = 3 * Nat;
    constexpr double sigma = 1.0;
    constexpr double epsilon = 1.0;
    constexpr double A = 4.0 * epsilon * sigma * sigma * sigma * sigma * sigma *
                         sigma * sigma * sigma * sigma * sigma * sigma * sigma;
    constexpr double B =
        4.0 * epsilon * sigma * sigma * sigma * sigma * sigma * sigma;
    std::array<double, dim> derivative;
    std::fill (derivative.begin (), derivative.end (), 0.0);
    double dx, dy, dz, d2, d6, E;

    for (unsigned i = 0; i < dim; i += 3) {
      for (unsigned j = i + 3; j < dim; j += 3) {
        dx = x[i] - x[j];
        dy = x[i + 1] - x[j + 1];
        dz = x[i + 2] - x[j + 2];
        d2 = dx * dx + dy * dy + dz * dz;
        d6 = d2 * d2 * d2; // r^6

        E = (12.0 * A / (d6 * d6) - 6.0 * B / d6) / d2;
        // not extracting the square root of 2:
        // force is (dx,dy,dz)/d*F,
        // where F is (12.0*A / (d6 * d6) - 6.0*B / d6)/d
        //=> (dx,dy,dz)/d*F = (dx,dy,dz)*(12.0*A / (d6 * d6) - 6.0*B / d6)/d2

        derivative[i] += E * dx;
        derivative[i + 1] += E * dy;
        derivative[i + 2] += E * dz;

        derivative[j] -= E * dx;
        derivative[j + 1] -= E * dy;
        derivative[j + 2] -= E * dz;
      }
    }

    return derivative;
  }
} // namespace MHTestFunctions
#endif // MHTESTFUNCTIONS_H
