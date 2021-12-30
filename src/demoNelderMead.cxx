#include "mhtestfunctions.hpp"
#include "simplexNelderMeadMinimizator.hpp"

#include <iostream>

template <typename T, int N>
void exampleSolverFromStarting (std::string title,
                                std::function<T (T *)> function) {
  typename MHMethods::simplex<T, N>::vertex startingvertex;
  for (int i = 0; i < N; ++i) {
    startingvertex[i] = -1.1;
  }
  auto t = MHMethods::minimizerNelderMeadFromStartingVertex<double, N> (
      2.5, startingvertex, 50000, function);
  std::cout << title << " minimum found: " << t.bestVertex.getValue () << " in "
            << t.functionEvaluations << " evaluations" << std::endl
            << " at\n";

  for (int i = 0; i < N; ++i) {
    std::cout << i << ": " << t.bestVertex[i] << std::endl;
  }
}

int main (int, char **) {
  exampleSolverFromStarting<double, 2> ("Rosenbrock",
                                        MHTestFunctions::Rosenbrock);
  exampleSolverFromStarting<double, 2> ("Himmelblau",
                                        MHTestFunctions::Himmelblau);
  constexpr unsigned spheredim = 10;
  exampleSolverFromStarting<double, spheredim> (
      "Sphere8D", MHTestFunctions::Sphere<spheredim>);
  return 0;
}
