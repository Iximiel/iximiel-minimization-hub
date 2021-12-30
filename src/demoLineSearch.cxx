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

template <typename T, int Nat> void exampleSolverForLJ (std::string title) {
  constexpr unsigned N = 3 * Nat;
  std::array<T, N> initialPosition;
  for (int i = 0; i < Nat; ++i) {
    initialPosition[i * 3] = i;
    initialPosition[i * 3 + 1] = i;
    initialPosition[i * 3 + 2] = i;
  }
  counter<T>::count = 0;
  counter<T> funcCounter{MHTestFunctions::LennardJones12_6<Nat>};
  T init = funcCounter (initialPosition.data ());
  auto derivative = MHTestFunctions::LennardJones12_6Derivative<Nat> (
      initialPosition.data ());
  auto t = MHMethods::linesearchBacktracking<double, N> (
      initialPosition, init, derivative, 1.0, funcCounter);
  std::cout << "Res : " << MHTestFunctions::LennardJones12_6<Nat> (t.data ())
            << " in " << counter<double>::count << " steps." << endl
            << endl;

  for (int i = 0; i < N; i += 3) {
    std::cout << i / 3 << ": " << t[i] << " " << t[i + 1] << " " << t[i + 2]
              << '\n';
  }
  std::cout << std::flush;
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

  constexpr unsigned spheredim = 8;
  test<double, spheredim> ("Sphere8D", MHTestFunctions::Sphere<spheredim>,
                           MHTestFunctions::SphereDerivative<spheredim>,
                           {1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0});

  exampleSolverForLJ<double, 3> ("LJ 3 atoms");
  exampleSolverForLJ<double, 6> ("LJ 6 atoms");
  exampleSolverForLJ<double, 9> ("LJ 9 atoms");
  return 0;
}
