#include "linesearches/linesearch.hpp"
#include "mhtestfunctions.hpp"
#include "mhutilities.hpp"

MHUtilities::derivativeAndValue<double, 2> func (double *c) {
  return {MHTestFunctions::Rosenbrock (c),
          MHTestFunctions::RosenbrockDerivative (c)};
}

int main (int, char **) {
  std::array<double, 2> x = {5, 5};
  std::array<double, 2> dir = func (x.data ()).derivative;
  dir[0] = x[0] - dir[0];
  dir[1] = x[1] - dir[1];
  MHMethods::linesearch<double, 2> (x, dir, 1.0, func);
  return 0;
}