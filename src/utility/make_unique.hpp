#ifndef MAKE_UNIQUE_HPP
#define MAKE_UNIQUE_HPP
/*
 * Bare bones implementation of make_unique for use with c++11
 * - from https://stackoverflow.com/questions/10149840/c-arrays-and-make-unique
 * - works with arrays
 * - should be able to remove directly in favor of the std::make_unique in C++14 and beyond
 */
#include <type_traits>
#include <utility>
#include <memory>
namespace std {

  template <class T, class ...Args>
  typename std::enable_if
  <
    !std::is_array<T>::value,
    std::unique_ptr<T>
    >::type
  make_unique(Args&& ...args)
  {
    return std::unique_ptr<T>(new T(std::forward<Args>(args)...));
  }

  template <class T>
  typename std::enable_if
  <
    std::is_array<T>::value,
    std::unique_ptr<T>
    >::type
  make_unique(std::size_t n)
  {
    typedef typename std::remove_extent<T>::type RT;
    return std::unique_ptr<T>(new RT[n]);
  }
}
#endif // MAKE_UNIQUE_HPP
