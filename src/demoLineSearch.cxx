#include "linesearches/linesearchBacktracking.hpp"
#include "mhtestfunctions.hpp"
#include "mhutilities.hpp"
#include <iostream>
using std::cout;
using std::endl;
constexpr unsigned dim = 2;
MHUtilities::derivativeAndValue<double, dim> func (double *c) {
  return {MHTestFunctions::Rosenbrock (c),
          MHTestFunctions::RosenbrockDerivative (c)};
}

template <typename T, unsigned dim>
void test (std::function<T (T *x)> func,
           std::function<std::array<T, dim> (T *x)> dev,
           std::array<double, dim> x = {1, 1}) {

  std::array<double, dim> derivative = dev (x.data ());
  cout << "Init:" << func (x.data ()) << endl;
  auto t = MHMethods::linesearchBacktracking<double, dim> (
      x, func (x.data ()), derivative, 1.0, func);
  cout << "res:" << func (t.data ()) << endl;
  for (size_t i = 0; i < dim; i++) {
    cout << x[i] << " -> " << t[i] << endl;
  }
}

int main (int, char **) {
  std::cout << "Rosenbrock{5,5}\n:";
  test<double, 2> (MHTestFunctions::Rosenbrock,
                   MHTestFunctions::RosenbrockDerivative, {5, 5});
  std::cout << "Parabola{0.1,0.1}\n:";
  test<double, 2> (MHTestFunctions::Parabola2D,
                   MHTestFunctions::Parabola2DDerivative, {1, 1});
  test<double, 2> (MHTestFunctions::Parabola2D,
                   MHTestFunctions::Parabola2DDerivative, {0.1, 0.1});
  return 0;
}
