#ifndef SIMPLEXMINIMIZATOR_H
#define SIMPLEXMINIMIZATOR_H
#include <algorithm>
#include <array>
#include <iostream>
#include <functional>

namespace simplexMethodMinimization {
constexpr unsigned maxFunctionEvaluations = 5000;

template <typename T, unsigned N> struct simplex {
  struct vertex {
    vertex() = default;
    vertex(const vertex &other) : value_(other.value_) {
      for (unsigned i = 0; i < N; ++i) {
        coordinates_[i] = other.coordinates_[i];
      }
    }
    vertex(vertex &&other)
        : value_(std::move(other.value_)), coordinates_(other.coordinates_) {}
    void evaluation(std::function<T (T *)> function) {
      value_ = function(coordinates_.data());
      // std::cout << value <<std::endl;
    }
    vertex &operator=(const vertex &other) {
      if (this != &other) {
        for (unsigned i = 0; i < N; ++i) {
          coordinates_[i] = other.coordinates_[i];
        }
        value_ = other.value_;
      }
      return *this;
    }
    vertex &operator=(vertex &&other) {
      if (this != &other) {
        std::swap(coordinates_, other.coordinates_);
        value_ = std::move(other.value_);
      }
      return *this;
    }
    // sorting
    bool operator>(const vertex &other) { return value_ > other.value_; }
    bool operator<(const vertex &other) { return value_ < other.value_; }
    T operator[](const unsigned &i) const { return coordinates_[i]; }
    T &operator[](const unsigned &i) { return coordinates_[i]; }
    T getValue() const { return value_; }

  protected:
    T value_{0.0};
    std::array<T, N> coordinates_;
  };
  simplex() = default;
  simplex(const simplex& s):vertices_(s.vertices_){}
  simplex(const simplex&& s):vertices_(s.vertices_){}
  simplex& operator=(const simplex& s){
      if(this!=&s) {
          vertices_=s.vertices_;
      }
      return *this;
  }
  simplex& operator=(const simplex&& s){
      if(this!=&s) {
          vertices_.swap(s.vertices_);
      }
      return *this;
  }
  simplex(T intitialDisplacement, const vertex &startingvertex) {
    vertices_[0] = startingvertex;
    for (unsigned i = 1; i < N + 1; ++i) {
      vertices_[i] = startingvertex;
      vertices_[i][i - 1] += intitialDisplacement;
    }
  }
  inline void centroidCalculation() {
    for (unsigned i = 0; i < N; ++i) {
      centroid_[i] = vertices_[0][i];
      // Numerical Recipes uses also the last point in the centroid, wikipedia
      // not
      for (unsigned j = 1; j < N /*+1*/; ++j) {
        centroid_[i] += vertices_[j][i];
      }
    }
  }
  vertex extrapolateOnTheWorst(const T &factor, std::function<T (T *)> function) const {
    vertex newTry;
    T f1 = (1 - factor) / N;
    T f2 = f1 - factor;
    for (unsigned i = 0; i < N + 1; ++i) {
      newTry[i] = centroid_[i] * f1 - vertices_[N][i] * f2;
    }
    newTry.evaluation(function);
    return newTry;
  }
  std::array<vertex, N + 1> vertices_;
  std::array<T, N> centroid_;
};

template <typename T, unsigned dim>
typename simplex<T, dim>::vertex
simplexMinimizer(simplex<T, dim> s,
                 std::function<T (T *)> function) {
  constexpr T notDen0 = 1e-10;

  using Vertex = typename simplex<T, dim>::vertex;

  T rtol;
  T ftol = 1e-5;
  unsigned functionEvaluations = 0;
  while (functionEvaluations < maxFunctionEvaluations) {
    // order
    std::sort(s.vertices_.begin(), s.vertices_.end());
    s.centroidCalculation();
    rtol = 2.0 *
           std::abs(s.vertices_[dim].getValue() - s.vertices_[0].getValue()) /
           (std::abs(s.vertices_[dim].getValue()) +
            std::abs(s.vertices_[0].getValue()) + notDen0);
    if (rtol < ftol) {
      break;
    }
    // reflection
    Vertex reflection = s.extrapolateOnTheWorst(-1.0, function);
    ++functionEvaluations;
    if (reflection < s.vertices_[0]) { // expand if is a good move
      s.vertices_[dim] = reflection;
      Vertex expansion = s.extrapolateOnTheWorst(2.0, function);
      ++functionEvaluations;
      if (expansion < reflection) {
        s.vertices_[dim] = expansion;
      }
    } else if (reflection < s.vertices_[dim - 1]) { // reflection
      s.vertices_[dim] = reflection;
    } else { // contraction or if contraction fails shrink
      Vertex contraction = s.extrapolateOnTheWorst(0.5, function);
      ++functionEvaluations;
      if (contraction < s.vertices_[dim]) { // contraction
        s.vertices_[dim] = contraction;
      } else { // shrink
        for (unsigned j = 1; j < dim + 1; ++j) {
          for (unsigned i = 0; i < dim; ++i) {
            s.vertices_[j][i] = 0.5 * (s.vertices_[j][i] - s.vertices_[0][i]);
            // x0 + alpha *(xi-x0)
          }
          s.vertices_[j].evaluation(function);
          ++functionEvaluations;
        }
      }
    }
  }
  std::cout << functionEvaluations << " / " << maxFunctionEvaluations
            << std::endl;
  return s.vertices_[0];
}

template <typename T, unsigned dim>
typename simplex<T, dim>::vertex
simplexMinimizerFromStartingVertex(T intitialDisplacement,
                 typename simplex<T, dim>::vertex startingvertex,
                 std::function<T (T *)> function) {
    using Simplex = simplex<T, dim>;

    Simplex s(intitialDisplacement, startingvertex);
    for (unsigned i = 0; i < dim + 1; ++i) {
        s.vertices_[i].evaluation(function);
    }
    return simplexMinimizer(s,function);
}
} // namespace simplexMethodMinimization
#endif // SIMPLEXMINIMIZATOR_H
