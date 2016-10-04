#ifndef ARRAY_HPP
#define ARRAY_HPP

#include <stdexcept>
#include <iostream>

namespace mathx {

/**
* @brief This class is a variable size, random-access data structure for type T
*/
template<class T>
class array{
private:
  /**
  * Base dynamic primative array
  */
  T* container;

  /**
  * Current size of the array
  */
  int mysize;

  /**
  * Current capacity of the array
  */
  int mycapacity;

  /**
  * Function to increase capacity
  */
  void grow();

  /**
  * Function to decrease capacity
  * @param the new capacity
  */
  void shrink(int);
public:
  /**
  * Default constructor initializing everything to 0
  */
  array<T>() : container(new T[0]), mysize(0), mycapacity(0){};

  /**
  * Constructor initializing container to capacity c
  * @param c - an int value to set the initial capacity
  */
  array<T>(int c) : container(new T[c]), mysize(0), mycapacity(c){};

  /**
  * Constructor accepting an initializer list
  * @details Allows assignment of array like array<int> arr = {1,2,3,4};
  * @param c - an initializer list (i.e {0,1,2,3})
  */
  array<T>(std::initializer_list<T> c){
    container = new T[c.size()];
    mycapacity = mysize = c.size();
    std::copy(c.begin(), c.end(), container);
  }

  /**
  * Method to add element to end of array
  */
  void push(T el);

  /**
  * Gets an element from the array
  * @param i - index of element. i < mysize
  */
  T get(int i){ if(i < mysize) return container[i]; else throw std::runtime_error("index out of bounds"); };

  /**
  * Method to pop element from end of array
  */
  T pop();

  /**
  * Method to destroy all elements in array
  */
  void clear();

  /**
  * Method to get the size of the array
  */
  int size() { return mysize; };

  /**
  * Method to get the capacity of the array
  */
  int capacity() { return mycapacity; };

  /**
  * Overload of array index operators. Same as get(i)
  * @param i - index of element. i < mysize
  */
  T& operator[](std::size_t i){ return container[i]; };

  /**
  * Overload of array index operators. Same as get(i)
  * @param i - index of element. i < mysize
  */
  T& operator[](std::size_t i) const { return container[i]; };

  /**
  * @brief An iterator class for array
  * @detail This class enables use of array in
  * functions expecting iterators such as for each loops
  */
  class iterator {
  public:
    /**
    * Constructor that takes a pointer to the iteratee
    */
    iterator(T* ptr) : ptr(ptr){};

    /**
    * Overload of ++ operator to move iterator through array
    */
    iterator operator++() { iterator i(ptr); ++ptr; return i; }

    /**
    * Overload of != operator for comparisons of iterators
    */
    bool operator!=(const iterator& other) { return ptr != other.ptr; }

    /**
    * Overload of * operator for access to data pointed to by ptr
    */
    const T& operator*() const { return *ptr; }
  private:
    /**
    * Stores the reference to the iteraree
    */
    T* ptr;
  };

  /**
  * Points to the first element of the array
  */
  iterator begin() { return iterator(container); };

  /**
  * Points to the "past-the-end element" of the array
  */
  iterator end() { return iterator(container + mysize); };
};

// PRIVATE METHODS
template<typename T>
/**
* Implementation of the private method array::grow()
*/
void array<T>::grow(){
  // Set the new capacity to double the current. 0 -> 2, 2 -> 4, ...
  // Doubling the array every grow leads to amortized O(n) grow operations
  array<T>::mycapacity = array<T>::mycapacity == 0 ? 2 : 2 * array<T>::mycapacity;

  // Initialize a temporary primative
  // array to shuffle old array for resize
  T* tmp = new T[array<T>::mycapacity];

  // Copy the old array to the temporary array
  std::copy(array<T>::container, array<T>::container + array<T>::mysize, tmp);

  // Delete old array
  delete[] array<T>::container;

  // Assign container to be the new array with
  // larger capacity
  array<T>::container = tmp;
};

/**
* Implementation of the private method array::shrink()
* @param capacity the new capacity to shrink to
*/
template<typename T>
void array<T>::shrink(int capacity){
  // Set capacity to new capacity
  array<T>::mycapacity = capacity;

  // Initialize temporary array
  T* tmp = new T[array<T>::mycapacity];

  // Copy the old array to smaller temp array
  std::copy(array<T>::container, array<T>::container + array<T>::mysize, tmp);

  // Destroy old array
  delete[] array<T>::container;

  // Assign container to be the new array with
  // smaller capacity
  array<T>::container = tmp;
};

// PUBLIC METHODS
/**
* Implementation of public method push
* @param el - element to append to end of array
*/
template<typename T>
void array<T>::push(T el){
  // If there is not room
  // in the array grow it
  if(array<T>::mysize >= array<T>::mycapacity){
    array<T>::grow();
  }

  // Add element to end of array
  array<T>::container[array<T>::mysize] = el;

  // Increment size
  array<T>::mysize++;
};

/**
* Implementation of public method pop
*/
template<typename T>
T array<T>::pop(){
  // Shrink array
  array<T>::shrink(array<T>::mysize - 1 == 0 ? 0 : array<T>::mysize);

  // Return the last element decrementing size
  return array<T>::container[--mysize];
};

/**
* Implementation of public method clear
*/
template<typename T>
void array<T>::clear(){
  // Set size to 0
  array<T>::mysize = 0;

  // Set capacity to 0
  array<T>::mycapacity = 0;

  // Destroy elements
  delete[] array<T>::container;
};

}

#endif
