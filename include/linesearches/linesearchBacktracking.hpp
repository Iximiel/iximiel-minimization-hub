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

  template <typename Iterator0, typename Iterator1, typename Iterator2,
            typename T>
  inline void multipyByFactorAndSum (Iterator0 x0, Iterator0 x0End,
                                     Iterator1 direction, Iterator2 result,
                                     const T lambda) {

    for (; x0 != x0End; ++x0, (void)++direction, (void)++result)
      *result = *x0 + lambda * (*direction);
  }

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
      // cout << sum << endl;
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
    // cout << slope << endl;
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
    T tmpLambda, rhs, rhsPrec, lambdaPrec, fPrec, a, b, disc, f;
    myArray xnew;
#define USENEW
#ifdef USENEW
    { // first step
      multipyByFactorAndSum (initialPosition.begin (), initialPosition.end (),
                             direction.begin (), xnew.begin (), lambda);
      /*for (unsigned i = 0; i < dim; ++i) {
        xnew[i] = initialPosition[i] + lambda * direction[i];
      }*/
      f = function (xnew.data ());
      if (lambda < lambdaMIN) {
        // this should throw or something
        swap (initialPosition, xnew);
        return xnew;
      } else if (f <= (fStart + Alpha * lambda * slope)) {
        return xnew;
      } else {
        tmpLambda = -slope / (2.0 * (f - fStart - slope));
      }
      rhs = f - fStart - lambda * slope;
    }
    for (;;) {
      rhsPrec = rhs;
      lambdaPrec = lambda;
      fPrec = f;
      lambda = std::max (tmpLambda, 0.1 * lambda);

      multipyByFactorAndSum (initialPosition.begin (), initialPosition.end (),
                             direction.begin (), xnew.begin (), lambda);
      /*for (unsigned i = 0; i < dim; ++i) {
        xnew[i] = initialPosition[i] + lambda * direction[i];
      }*/
      f = function (xnew.data ());
      if (lambda < lambdaMIN) {
        // this should throw or something
        swap (initialPosition, xnew);
        return xnew;
      } else if (f <= (fStart + Alpha * lambda * slope)) {
        return xnew;
      } else {
        rhs = f - fStart - lambda * slope;
        // rhsPrec = fPrec - fStart - lambdaPrec * slope;
        a = (rhs / (lambda * lambda) - rhsPrec / (lambdaPrec * lambdaPrec)) /
            (lambda - lambdaPrec);
        b = (lambda * rhsPrec / (lambdaPrec * lambdaPrec) -
             lambdaPrec * rhs / (lambda * lambda)) /
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
      }
    }
#else
    for (;;) {
      multipyByFactorAndSum (initialPosition.begin (), initialPosition.end (),
                             direction.begin (), xnew.begin (), lambda);
      // for (unsigned i = 0; i < dim; ++i) {xnew[i] = initialPosition[i] +
      lambda *direction[i];
    }
    f = function (xnew.data ());
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
        rhs = f - fStart - lambda * slope;
        rhsPrec = fPrec - fStart - lambdaPrec * slope;
        a = (rhs / (lambda * lambda) - rhsPrec / (lambdaPrec * lambdaPrec)) /
            (lambda - lambdaPrec);
        b = (lambda * rhsPrec / (lambdaPrec * lambdaPrec) -
             lambdaPrec * rhs / (lambda * lambda)) /
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
    // rhsPrec=rhs;
    lambdaPrec = lambda;
    fPrec = f;
    lambda = std::max (tmpLambda, 0.1 * lambda);
  }
#endif
  } // namespace MHMethods
} // namespace MHMethods
#endif // MHLINESEARCH_H1