#include <iostream>
#include "simplexMinimizator.h"

int main(int, char **) {
    constexpr int N=2;
    simplexMethodMinimization::simplex<double,N>::vertex startingvertex;
    for(int i=0;i<N;++i) {
        startingvertex[i]=1.0;
    }
    auto t = simplexMethodMinimization::simplexMinimizerFromStartingVertex<double,N>(2.5,startingvertex,[](double* x)->double{
        //return (x[0]-1)*(x[0]-1)+x[1]*x[1];
        //return (x[0]-1)*(x[0]-1) + 100.0*(x[1]-x[0]*x[0])*(x[1]-x[0]*x[0]);
        double a=(x[0]*x[0]+x[1]-11);
        double b =(x[1]*x[1]+x[0]-7);
        return a*a+b*b;
    });
    std::cout << t.getValue() << std::endl;
    for(int i=0;i<N;++i) {
        std::cout <<i << ": "<< t[i] << std::endl;
    }
    return 0;
}
