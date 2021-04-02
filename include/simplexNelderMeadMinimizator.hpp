#ifndef SIMPLEXMINIMIZATOR_H
#define SIMPLEXMINIMIZATOR_H
#include <algorithm>
#include <array>
#include <functional>
#include <iostream>

namespace MHMethods {

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

    void evaluation(std::function<T(T *)> function) {
      value_ = function(coordinates_.data());
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

    // accessing
    T operator[](const unsigned &i) const { return coordinates_[i]; }
    T &operator[](const unsigned &i) { return coordinates_[i]; }

    T getValue() const { return value_; }

  protected:
    T value_{0.0};
    std::array<T, N> coordinates_;
  };

  simplex() = default;

  simplex(const simplex &s) : vertices_(s.vertices_) {}

  simplex(const simplex &&s) : vertices_(s.vertices_) {}

  simplex &operator=(const simplex &s) {
    if (this != &s) {
      vertices_ = s.vertices_;
    }
    return *this;
  }

  simplex &operator=(const simplex &&s) {
    if (this != &s) {
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

  inline void updateCentroid() {
    for (unsigned i = 0; i < N; ++i) {
      centroid_[i] = vertices_[0][i];
      // Numerical Recipes uses also the last point in the centroid
      for (unsigned j = 1; j < N /*+1*/; ++j) {
        centroid_[i] += vertices_[j][i];
      }
      centroid_[i] /= N;
    }
  }

  void orderVertices() { std::sort(vertices_.begin(), vertices_.end()); }

  vertex newVertexFromWorstAndCentroid(const T &factor,
                                       std::function<T(T *)> function) const {
    vertex newTry;
    for (unsigned i = 0; i < N; ++i) {
      newTry[i] = (centroid_[i] - vertices_[N][i]) * factor + centroid_[i];
    }
    newTry.evaluation(function);
    return newTry;
  }

  // accessing
  vertex operator[](const unsigned &i) const { return vertices_[i]; }
  vertex &operator[](const unsigned &i) { return vertices_[i]; }

protected:
  std::array<vertex, N + 1> vertices_;
  std::array<T, N> centroid_;
};

template <typename T, unsigned dim> struct NelderMeadReturnType {
  typename simplex<T, dim>::vertex bestVertex;
  unsigned functionEvaluations;
};

template <typename T, unsigned dim>
NelderMeadReturnType<T, dim>
minimizerNelderMead(simplex<T, dim> s, const unsigned maxFunctionEvaluations,
                    std::function<T(T *)> function) {
  constexpr T notDen0 = 1e-10;

  using Vertex = typename simplex<T, dim>::vertex;

  T rtol;
  T ftol = 1e-5;
  unsigned functionEvaluations = 0;
  while (functionEvaluations < maxFunctionEvaluations) {
    // order
    s.orderVertices();
    s.updateCentroid();
    rtol = 2.0 * std::abs(s[dim].getValue() - s[0].getValue()) /
           (std::abs(s[dim].getValue()) + std::abs(s[0].getValue()) + notDen0);
    if (rtol < ftol) {
      break;
    }
    // reflection
    Vertex reflection = s.newVertexFromWorstAndCentroid(1.0, function);
    ++functionEvaluations;
    if (reflection < s[0]) { // expand if is a good move
      s[dim] = reflection;
      Vertex expansion = s.newVertexFromWorstAndCentroid(-2.0, function);
      ++functionEvaluations;
      if (expansion < reflection) {
        s[dim] = expansion;
      }
    } else if (reflection < s[dim - 1]) { // reflection
      s[dim] = reflection;
    } else { // contraction or if contraction fails shrink
      Vertex contraction = s.newVertexFromWorstAndCentroid(-0.5, function);
      ++functionEvaluations;
      if (contraction < s[dim]) { // contraction
        s[dim] = contraction;
      } else { // shrink
        for (unsigned j = 1; j < dim + 1; ++j) {
          for (unsigned i = 0; i < dim; ++i) {
            s[j][i] = 0.5 * (s[j][i] - s[0][i]);
            // x0 + alpha *(xi-x0)
          }
          s[j].evaluation(function);
          ++functionEvaluations;
        }
      }
    }
  }
  // std::cout << functionEvaluations << " / " << maxFunctionEvaluations <<
  // std::endl;
  return {s[0], functionEvaluations};
}

template <typename T, unsigned dim>
NelderMeadReturnType<T, dim> minimizerNelderMeadFromStartingVertex(
    T intitialDisplacement, typename simplex<T, dim>::vertex startingvertex,
    const unsigned maxFunctionEvaluations, std::function<T(T *)> function) {
  using Simplex = simplex<T, dim>;

  Simplex s(intitialDisplacement, startingvertex);
  for (unsigned i = 0; i < dim + 1; ++i) {
    s[i].evaluation(function);
  }
  return minimizerNelderMead(s, maxFunctionEvaluations, function);
}
} // namespace MHMethods
#endif // SIMPLEXMINIMIZATOR_H
