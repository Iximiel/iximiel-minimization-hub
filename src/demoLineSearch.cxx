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

template <typename T> class counter {
public:
  counter (std::function<T (T *x)> func) : f_ (func) {}
  T operator() (T *x) {
    ++count;
    return f_ (x);
  }
  static unsigned count;

private:
  std::function<T (T *x)> f_;
};
template <typename T> unsigned counter<T>::count = 0;

template <typename T, unsigned dim>
std::array<T, dim> test (std::string name, std::function<T (T *x)> func,
                         std::function<std::array<T, dim> (T *x)> dev,
                         std::array<double, dim> x = {1, 1}) {
  std::cout << name << " {";
  for (unsigned i = 0; i < dim; ++i) {
    std::cout << ((i == 0) ? "" : ", ") << x[i];
  }
  std::cout << "}. ";
  counter<T>::count = 0;
  counter<T> funcCounter{func};
  std::array<T, dim> derivative = dev (x.data ());
  T init = funcCounter (x.data ());
  cout << "Initial value:" << init << endl;
  auto t = MHMethods::linesearchBacktracking<double, dim> (x, init, derivative,
                                                           1.0, funcCounter);
  std::cout << "->{";
  for (unsigned i = 0; i < dim; ++i) {
    std::cout << ((i == 0) ? "" : ", ") << t[i];
  }
  std::cout << "}. Res : " << func (t.data ()) << " in "
            << counter<double>::count << " steps." << endl
            << endl;
  return t;
}

int main (int, char **) {

  test<double, 2> ("Rosenbrock", MHTestFunctions::Rosenbrock,
                   MHTestFunctions::RosenbrockDerivative, {5.0, 5.0});
  {
    std::array<double, 2> t{3.5, -2.0};
    for (int i = 0; i < 4; ++i) {
      t = test<double, 2> ("Himmelblau", MHTestFunctions::Himmelblau,
                           MHTestFunctions::HimmelblauDerivative, t);
    }
  }
  test<double, 2> ("Parabola", MHTestFunctions::Parabola2D,
                   MHTestFunctions::Parabola2DDerivative, {1.0, 1.0});

  test<double, 2> ("Parabola", MHTestFunctions::Parabola2D,
                   MHTestFunctions::Parabola2DDerivative, {2.0, 0.0});

  test<double, 2> ("Parabola", MHTestFunctions::Parabola2D,
                   MHTestFunctions::Parabola2DDerivative, {0.1, 0.1});

  return 0;
}
