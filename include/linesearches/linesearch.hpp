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
  // Quasi-newtown linesearch, following Numerical recipes
  template <typename T, unsigned dim>
  void linesearch (
      std::array<T, dim> initialPosition, std::array<T, dim> direction,
      const T maxStepLenght,
      std::function<MHUtilities::derivativeAndValue<T, dim> (T *)> function) {
    using myArray = std::array<T, dim>;
    MHUtilities::derivativeAndValue<T, dim> fval =
        function (initialPosition.data ());
    constexpr T Alpha = 1e-4;
    constexpr T xTol = std::numeric_limits<T>::epsilon ();
    /// TODO: the accumulation of the sum can be done in the very same loop of
    /// the slope calculation
    T sum = sqrt (std::inner_product (direction.begin (), direction.end (),
                                      direction.begin (), T (0.0)));
    cout << sum << endl;
    if (sum > maxStepLenght) {
      sum = maxStepLenght / sum;
      std::transform (direction.begin (), direction.end (), direction.begin (),
                      [=] (T i) { return i * sum; });
    }
    /// TODO: the accumulation of the sum can be done in the very same loop of
    /// the slope calculation
    T slope = std::inner_product (direction.begin (), direction.end (),
                                  fval.derivative.begin (), T (0.0));
    cout << slope << endl;
    if (slope >= 0.0) {
      throw "Roundoff problem in linesearch: slope is not negative";
    }
    /*{
      myArray temp;
      T slope = std::inner_product (direction.begin (), direction.end (),
                                    initialPosition.begin (), T (0.0), ());
      template <class InputIterator1, class InputIterator2, class T>
      T inner_product (InputIterator1 first1, InputIterator1 last1,
                       InputIterator2 first2, T init) {
        while (first1 != last1) {
          init = init + (*first1) * (*first2);
          // or: init = binary_op1 (init, binary_op2(*first1,*first2));
          ++first1;
          ++first2;
        }
        return init;
      }
    }*/
  }
} // namespace MHMethods
#endif // MHLINESEARCH_H