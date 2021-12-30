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

template <typename T, int Nat> void exampleSolverForLJ (std::string title) {
  constexpr unsigned N = 3 * Nat;
  typename MHMethods::simplex<T, N>::vertex startingvertex;
  for (int i = 0; i < Nat; ++i) {
    startingvertex[i * 3] = i;
    startingvertex[i * 3 + 1] = i;
    startingvertex[i * 3 + 2] = i;
  }
  auto t = MHMethods::minimizerNelderMeadFromStartingVertex<double, N> (
      2.5, startingvertex, 50000, MHTestFunctions::LennardJones12_6<Nat>);
  std::cout << title << " minimum found: " << t.bestVertex.getValue () << " in "
            << t.functionEvaluations << " evaluations" << std::endl
            << " at\n";

  for (int i = 0; i < N; i += 3) {
    std::cout << i / 3 << ": " << t.bestVertex[i] << " " << t.bestVertex[i + 1]
              << " " << t.bestVertex[i + 2] << '\n';
  }
  std::cout << std::flush;
}

int main (int, char **) {
  exampleSolverFromStarting<double, 2> ("Rosenbrock",
                                        MHTestFunctions::Rosenbrock);
  exampleSolverFromStarting<double, 2> ("Himmelblau",
                                        MHTestFunctions::Himmelblau);
  exampleSolverForLJ<double, 3> ("LJ 3 atoms");
  exampleSolverForLJ<double, 6> ("LJ 6 atoms");
  exampleSolverForLJ<double, 9> ("LJ 9 atoms");
  constexpr unsigned spheredim = 10;
  exampleSolverFromStarting<double, spheredim> (
      "Sphere8D", MHTestFunctions::Sphere<spheredim>);
  return 0;
}
