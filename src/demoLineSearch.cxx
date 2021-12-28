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

int main (int, char **) {
  std::array<double, dim> x = {5, 5};
  std::array<double, dim> derivative =
      MHTestFunctions::RosenbrockDerivative (x.data ());
  cout << "Init:" << MHTestFunctions::Rosenbrock (x.data ()) << endl;
  auto t = MHMethods::linesearchBacktracking<double, dim> (
      x, MHTestFunctions::Rosenbrock (x.data ()), derivative, 1.0,
      MHTestFunctions::Rosenbrock);
  cout << "res:" << MHTestFunctions::Rosenbrock (t.data ()) << endl;
  for (size_t i = 0; i < dim; i++) {
    cout << x[i] << " -> " << t[i] << endl;
  }
  return 0;
}