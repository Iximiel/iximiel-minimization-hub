#ifndef MHLINESEARCH_H
#define MHLINESEARCH_H
#include "../mhutilities.hpp"
#include <algorithm>
#include <array>
#include <cmath>
#include <functional>
#include <iostream>
#include <limits>
#include <numeric>
using std::cout;
using std::endl;
namespace MHMethods {
  // Quasi-newtown linesearchBacktracking, following Numerical recipes
  template <typename T, unsigned dim>
  std::array<T, dim>
  linesearchBacktracking (std::array<T, dim> initialPosition, const T fStart,
                          const std::array<T, dim> &derivative,
                          const T maxStepLenght,
                          std::function<T (T *)> function) {
    using myArray = std::array<T, dim>;
    constexpr T Alpha = 1e-4;
    constexpr T xTol = std::numeric_limits<T>::epsilon ();
    myArray direction = [] (const myArray &x, const myArray &g) {
      myArray t;
      for (unsigned i = 0; i < dim; ++i) {
        t[i] = x[i] - g[i];
      }
      return t;
    }(initialPosition, derivative);

    /// TODO: the accumulation of the sum can be done in the very same loop of
    /// the slope calculation
    {
      T sum = sqrt (std::inner_product (direction.begin (), direction.end (),
                                        direction.begin (), T (0.0)));
      cout << sum << endl;
      if (sum > maxStepLenght) {
        sum = maxStepLenght / sum;
        std::transform (direction.begin (), direction.end (),
                        direction.begin (), [=] (T i) { return i * sum; });
      }
    }
    /// TODO: the accumulation of the sum can be done in the very same loop of
    /// the slope calculation
    T slope = std::inner_product (direction.begin (), direction.end (),
                                  derivative.begin (), T (0.0));
    cout << slope << endl;
    if (slope >= 0.0) {
      throw "Roundoff problem in linesearchBacktracking: slope is not negative";
    }

    T lambdaMIN =
        xTol /
        std::inner_product (
            direction.begin (), direction.end (), initialPosition.begin (),
            T (0.0),
            [] (T init, T result) { return init > result ? init : result; },
            [] (T dir, T xold) {
              return std::abs (dir / (std::abs (xold) > 1.0 ? xold : 1.0));
            });
    T lambda = 1.0;
    T tmpLambda, rhs1, rhs2, lambdaPrec, fPrec, a, b, disc;
    myArray xnew;
    for (;;) {
      for (unsigned i = 0; i < dim; ++i) {
        xnew[i] = initialPosition[i] + lambda * direction[i];
      }
      T f = function (xnew.data ());
      if (lambda < lambdaMIN) {
        // this should throw or something
        swap (initialPosition, xnew);
        return xnew;
      } else if (f <= (fStart + Alpha * lambda * slope)) {
        return xnew;
      } else {
        if (lambda == 1.0) {
          tmpLambda = -slope / (2.0 * (f - fStart - slope));
        } else {
          rhs1 = f - fStart - lambda * slope;
          rhs2 = fPrec - fStart - lambdaPrec * slope;
          a = (rhs1 / (lambda * lambda) - rhs2 / (lambdaPrec * lambdaPrec)) /
              (lambda - lambdaPrec);
          b = (lambda * rhs2 / (lambdaPrec * lambdaPrec) -
               lambdaPrec * rhs1 / (lambda * lambda)) /
              (lambda - lambdaPrec);
          if (a == 0.0) {
            tmpLambda = -slope / (2.0 * b);
          } else {
            disc = b * b - 3.0 * a * slope;
            if (disc < 0.0) {
              tmpLambda = 0.5 * lambda;
            } else if (b <= 0.0) {
              tmpLambda = (sqrt (disc) - b) / (3.0 * a);
            } else
              tmpLambda = -slope / (b + sqrt (disc));
          }
          if (tmpLambda > 0.5 * lambda) {
            tmpLambda = 0.5 * lambda;
          }
        } // not first step
      }
      lambdaPrec = lambda;
      fPrec = f;
      lambda = std::max (tmpLambda, 0.1 * lambda);
    }
  }
} // namespace MHMethods
#endif // MHLINESEARCH_H1