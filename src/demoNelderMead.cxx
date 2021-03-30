#include "simplexNelderMeadMinimizator.h"
#include <iostream>

template <typename T, int N>
void exampleSolverFromStarting (std::string title, std::function<T(T *)> function) {
    typename simplexNelderMeadMethodMinimization::simplex<T, N>::vertex
        startingvertex;
    for (int i = 0; i < N; ++i) {
        startingvertex[i] = 1.1;
    }
    auto t = simplexNelderMeadMethodMinimization::
        minimizerNelderMeadFromStartingVertex<double, N>(2.5, startingvertex,
                                                         function);
    std::cout << title << " minimum found: "
    << t.getValue() << std::endl
              << " at\n";

    for (int i = 0; i < N; ++i) {
        std::cout << i << ": " << t[i] << std::endl;
    }
}

int main(int, char **) {
  auto Rosenbrock = [](double *x) -> double {
    return 1.0 * (x[0] - 1) * (x[0] - 1) +
           100.0 * (x[1] - x[0] * x[0]) * (x[1] - x[0] * x[0]);
  };
  auto Himmelblau = [](double *x) -> double {
    double a = (x[0] * x[0] + x[1] - 11);
    double b = (x[0] + x[1] * x[1] - 7);
    return a * a + b * b;
  };

  exampleSolverFromStarting<double,2> ("Rosenbrock",Rosenbrock);
  exampleSolverFromStarting<double,2> ("Himmelblau",Himmelblau);
  return 0;
}
