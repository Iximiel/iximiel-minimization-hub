#ifndef MHUTILITIES_H
#define MHUTILITIES_H
#include <array>
namespace MHUtilities {
  template <typename T, unsigned dim> struct derivativeAndValue {
    T value;
    std::array<T, dim> derivative;
  };
} // namespace MHUtilities
#endif // MHUTILITIES_H