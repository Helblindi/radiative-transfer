/*
 * Copyright 2020-2021 The Texas A&M University System, an agency of the State of
 * Texas, and Comet contributors
 * 
 * Redistribution and use in source and binary forms, with or without modification,
 * are permitted provided that the following conditions are met:
 * 
 * 1. Redistributions of source code must retain the above copyright notice, this
 *    list of conditions and the following disclaimer.
 * 
 * 2. Redistributions in binary form must reproduce the above copyright notice,
 *    this list of conditions and the following disclaimer in the documentation
 *    and/or other materials provided with the distribution.
 * 
 * 3. Neither the name of the copyright holder nor the names of its contributors
 *    may be used to endorse or promote products derived from this software
 *    without specific prior written permission.
 * 
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
 * ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
 * WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 * DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
 * FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
 * DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
 * SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
 * CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
 * OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 * OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/

#pragma once

#include <cassert>
#include <stdexcept>
#include <atomic>

/*! \file NDArray.h
 *  \brief A shared multi-dimensional array of elements.
 */

namespace ndarray
{

/*  ndarray::array is a shared, contiguous container that provides support for fixed
 *  size multi-dimensional arrays. Several arrays may own the same data, and memory
 *  is deallocated when the last remaining array owning the data is destroyed. All
 *  array copies are shallow. std::copy or similar should be used to obtain value
 *  copies of the data.
 * 
 *  All elements are stored contiguously in memory. N-dimensional address tuples
 *  are linearized to provide access to elements. Storage size is fixed between 
 *  calls to resize().
 * 
 *  Template parameters:
 *  T - Type of elements
 *  N - Rank
 * 
 *  Member functions:
 *  (constructor) - constructs array
 *  (destructor)  - destructs array
 *  operator=     - assigns the array
 *  swap          - swaps the contents with another array
 *  size          - returns number of array elements
 *  begin         - returns an iterator to the beginning
 *  cbegin        - returns a const iterator to the beginning
 *  end           - returns an iterator to the end
 *  cend          - returns a const iterator to the end
 *  data          - direct access to the array data
 *  empty         - checks if array is empty
 *  rank          - returns array rank
 *  resize        - resizes the array
 *  operator()    - access the specifed element
 * 
 *  Example:
 *  #include "NDArray.h"
 * 
 *  #include <algorithm>
 *  #include <iostream>
 * 
 *  int main()
 *  {
 *    // Create an array of rank three with size 2x3x5 that will contain doubles
 *    ndarray::array<double, 3> ar(2, 3, 5);
 * 
 *    // Fill ar with a constant value
 *    std::fill(ar.begin(), ar.end(), 5.0);
 * 
 *    // Print out ar
 *    for(double e : ar)
 *      std::cout << e << " ";
 *    std::cout << "\n";
 * 
 *    // Resize ar 
 *    ar.resize(1, 3, 3);
 *   
 *    // Fill ar and add 1 to the last value
 *    std::fill(ar.begin(), ar.end(), 1.0);
 *    ar(0, 2, 2) += 1.0;
 * 
 *    // Print out ar
 *    for(double e: ar)
 *      std::cout << e << " ";
 *    std::cout << "\n";
 *  }
 * 
 *  Output:
 *  5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5
 *  1 1 1 1 1 1 1 1 2
 */

template<typename T, size_t N>
class array
{
private:
  T *m_base;
  size_t *m_mpcs;
  size_t m_size;
  std::atomic<size_t> *m_counter;

  template<bool...>
  struct bool_pack{};

  template<class... U>
  using conjunction = std::is_same<bool_pack<true,U::value...>, bool_pack<U::value..., true>>;

  template<typename... U>
  using AllIntegral = typename conjunction<std::is_integral<U>...>::type;

public:
  /*! \brief Creates an array with the specified number of elements in each dimentsion.
   *  \param args Comma separated list of the number of elememnts in each dimension.
   *  \throw std::bad_alloc if memory allocation fails.
   *
   *  This constructor creates an array with the specified size and initializes the refernce count to one.
   */
  template<typename ... Args>
  explicit array(Args const&...args)
    : m_base(nullptr), m_mpcs(nullptr), m_size(0), m_counter(nullptr)
  {
    static_assert(N > 0, "ndarray::array: Size must be >= 1");
    static_assert(AllIntegral<Args...>::value, "ndarray::array: All parameters must be of integral type");
    static_assert(sizeof...(args) == N, "ndarray::array: Extents must match array size");
    size_t dims[] { static_cast<size_t>(args)... };

    m_mpcs = new size_t[N];
    m_size = 1;
    for(size_t i = 0;i < N;++i)
    {
      assert(dims[i] > 0);
      m_size *= dims[i];
      m_mpcs[i] = 1;
      for(size_t j = i+1;j < N;++j)
        m_mpcs[i] *= dims[j];
    }

    m_base = new T[m_size];
    m_counter = new std::atomic<size_t>();
    m_counter->store(1);
  }

  /*! \brief Creates an empty array.
   *  \throw std::bad_alloc if memory allocation fails.
   *
   *  This constructor creates an empty array and initializes the refernce count to one.
   */
  array()
    : m_base(nullptr), m_mpcs(nullptr), m_size(0), m_counter(nullptr)
  {
    m_mpcs = new size_t[N];
    m_counter = new std::atomic<size_t>();
    m_counter->store(1);
  }

  /*! \brief Deletes the array.
   *
   *  The destructor substracts one from the reference count and deletes 
   *  the underlying array data if the reference count is zero.
   */
  ~array()
  {
    --(*m_counter);
    if(*m_counter == 0)
    {
      delete [] m_base;
      m_base = nullptr;
      delete [] m_mpcs;
      m_mpcs = nullptr;
      delete m_counter;
      m_counter = nullptr;
    }  
  }

  /*! Swap the contents of this array with another array.
   *  \param other The array to swap with.
   */
  void swap(array<T, N> &other)
  {
    std::swap(m_base, other.m_base);
    std::swap(m_mpcs, other.m_mpcs);
    std::swap(m_size, other.m_size);
    std::swap(m_counter, other.m_counter);
  }

  /*! Copy construct from another array.
   *  \param other The array to copy.
   */
  explicit array(array<T, N> const &other)
    : m_base(other.m_base), m_mpcs(other.m_mpcs), m_size(other.m_size), m_counter(other.m_counter)
  {
    ++(*m_counter);
  }

  /*! Assign from another array.
   *  \param other The array to copy.
   */
  array<T, N> & operator=(array<T, N> const &other)
  {
    array<T, N>(other).swap(*this);
    return *this;
  }

  array(array<T, N> &&) = delete;
  array<T, N> & operator=(array<T, N> &&) = delete;

  /*! Returns the number of elements in the array.
   *  \return The number of elements in the array.
   */
  size_t size() const noexcept { return m_size; }

  /*! Returns true if the array has no elements.
   *  \return True if the array contains no elements; false, otherwise.
   */
  bool empty() const noexcept { return m_size == 0; }

  /*! Returns an iterator pointing to the beginning of the array.
   *  \return An iterator pointing the the beginning of the array.
   */
  T * begin() const noexcept { return m_base; }

  /*! Returns a constant iterator pointing to the beginning of the array.
   *  \return A constant iterator pointing to the beginning of the array.
   */
  const T * cbegin() const noexcept { return m_base; }

  /*! Returns an iterator pointing to one element past the end of the array.
   *  \return An iterator pointing to one element past the end of the array.
   */
  T * end() const noexcept { return m_base + m_size; }

  /*! Returns a constant iterator pointing to one element past the end of the array. 
   *  \return A constant iterator pointing to one element past the end of the array.
   */
  const T * cend() const noexcept { return m_base + m_size; }

  /*! Returns a pointer to the underlying array data.
   *  \return A pointer to the underlying array data.
   */
  T * data() const noexcept { return m_base; }

  /*! Returns the rank of the array.
   *  \return The rank of the array.
   */
  size_t rank() const noexcept { return N; }

  /*! \brief Resizes the array.
   *  \param args Comma separated list of the number of elememnts in each dimension.
   *  \throw std::bad_alloc if memory allocation fails.
   *
   *  This method resizes the array to the specified number of elements. If the current
   *  size is equal to the new size, no memory allocation occurs.
   */
  template<typename ... Args>
  void resize(Args const &...args)
  {
    static_assert(N > 0, "ndarray::array: Size must be >= 1");
    static_assert(AllIntegral<Args...>::value, "ndarray::array: All parameters must be of integral type");
    static_assert(sizeof...(args) == N, "ndarray::array: Extents must match array size");
    size_t dims[] { static_cast<size_t>(args)... };

    size_t curr_size = m_size;
    m_size = 1;
    for(size_t i = 0;i < N;++i)
    {
      m_size *= dims[i];
      m_mpcs[i] = 1;
      for(size_t j = i+1;j < N;++j)
        m_mpcs[i] *= dims[j];
    }
    if(m_size != curr_size)
    {
      if(m_base)
      {
        delete [] m_base;
        m_base = nullptr;
      }
      m_base = new T[m_size];
    }
  }

  /*! Accesses the specified element.
   *  \param args The indices of the desired element.
   *  \return Read/write reference to the element.
   */
  template<typename... Args>
  T & operator()(Args... args) noexcept
  {
    static_assert(AllIntegral<Args...>::value, "ndarray::array: All parameters must be of integral type");
    static_assert(sizeof...(args) == N, "ndarray::array: Extents must match array size");
    size_t indices[] { static_cast<size_t>(args)... };
    T *address = m_base + indices[N-1];
    for(size_t i = 0;i < N-1;++i)
      address += m_mpcs[i]*indices[i];
    return *(address);
  }

  /*! Accesses the specified element.
   *  \param args The indices of the desired element.
   *  \return Read reference to the element.
   */
  template<typename... Args>
  T const & operator()(Args... args) const noexcept
  {
    static_assert(AllIntegral<Args...>::value, "ndarray::array: All parameters must be of integral type");
    static_assert(sizeof...(args) == N, "ndarray::array: Extents must match array size");
    size_t indices[] { static_cast<size_t>(args)... };
    T *address = m_base + indices[N-1];
    for(size_t i = 0;i < N-1;++i)
      address += m_mpcs[i]*indices[i];
    return *(address);
  }
};

}

