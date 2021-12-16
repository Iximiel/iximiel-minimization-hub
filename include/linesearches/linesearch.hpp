#ifndef MHLINESEARCH_H
#define MHLINESEARCH_H
#include <array>
#include <functional>
namespace MHMethods {
//Quasi-newtown linesearch, following Numerical recipes
template <typename T, unsigned dim>
void 
linesearch(std::array<T, dim> initial, 
                    std::function<T(T *)> function);
}
#endif // MHLINESEARCH_H