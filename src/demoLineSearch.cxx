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
void test (std::string name, std::function<T (T *x)> func,
           std::function<std::array<T, dim> (T *x)> dev,
           std::array<double, dim> x = {1, 1}) {
  std::cout << name << " {";
  for (unsigned i = 0; i < dim; ++i) {
    std::cout << ((i == 0) ? "" : ", ") << x[i];
  }
  std::cout << "}:\n";
  counter<T>::count = 0;
  counter<T> funcCounter{func};
  std::array<T, dim> derivative = dev (x.data ());
  T init = funcCounter (x.data ());
  cout << "Initial value:" << init << endl;
  auto t = MHMethods::linesearchBacktracking<double, dim> (x, init, derivative,
                                                           1.0, funcCounter);
  std::cout << "count:" << counter<double>::count << "\n";
  cout << "res:" << func (t.data ()) << endl;
  for (size_t i = 0; i < dim; i++) {
    cout << x[i] << " -> " << t[i] << endl;
  }
}

int main (int, char **) {

  test<double, 2> ("Rosenbrock", MHTestFunctions::Rosenbrock,
                   MHTestFunctions::RosenbrockDerivative, {5, 5});

  test<double, 2> ("Parabola", MHTestFunctions::Parabola2D,
                   MHTestFunctions::Parabola2DDerivative, {1, 1});

  test<double, 2> ("Parabola", MHTestFunctions::Parabola2D,
                   MHTestFunctions::Parabola2DDerivative, {0.1, 0.1});

  return 0;
}
