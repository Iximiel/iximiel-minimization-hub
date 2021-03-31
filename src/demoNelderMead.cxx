#include "simplexNelderMeadMinimizator.hpp"
#include "mhtestfunctions.hpp"

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
    exampleSolverFromStarting<double,2> ("Rosenbrock", MHTestFunctions::Rosenbrock);
    exampleSolverFromStarting<double,2> ("Himmelblau", MHTestFunctions::Himmelblau);
  return 0;
}
