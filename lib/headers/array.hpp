#ifndef ARRAY_HPP
#define ARRAY_HPP

#include <sstream>
#include <stdexcept>
#include <iostream>

namespace mathx {

/**
* @brief This class is a variable size, random-access data structure for type T
* @details As working with raw pointers can be dangerous and hard to debug, this
* class serves as a wrapper of a pointer container of type T. Included in this
* class are helper methods to accomplish tasks such as adding to vectors, doting
* two vectors, and dynamically extending the container via a push method.
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
  int my_size;

  /**
  * Current capacity of the array
  */
  int my_capacity;

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
  array<T>() : container(new T[0]), my_size(0), my_capacity(0){};

  /**
  * Constructor initializing container to capacity c
  * @param c - an int value to set the initial capacity
  */
  array<T>(int c) : container(new T[c]), my_size(0), my_capacity(c){};

  /**
  * Constructor initializing container to capacity c with value v
  * @param c - an int value to set the initial capacity
  * @param v - a value to set all elements to
  */
  array<T>(int c, T v): container(new T[c]), my_size(c), my_capacity(c){
    for(int i = 0; i < c; i++){
      container[i] = v;
    }
  };

  /**
  * Constructor accepting an initializer list
  * @details Allows assignment of array like array<int> arr = {1,2,3,4};
  * @param c - an initializer list (i.e {0,1,2,3})
  */
  array<T>(std::initializer_list<T> c){
    container = new T[c.size()];
    my_capacity = my_size = c.size();
    std::copy(c.begin(), c.end(), container);
  }

  /**
  * Copy constructor
  */
  array<T>(const array<T> &a):container(new T[a.my_capacity]), my_size(a.my_size), my_capacity(a.my_capacity){
    for(int i = 0; i < my_size; i++)
      container[i] = a[i];
  };

  /**
  * assignment operator
  */
  array<T>& operator=(const array<T>& rhs){
    if(this == &rhs)
      return *this;
    my_size = rhs.my_size;
    my_capacity = rhs.my_capacity;
    container = new T[my_capacity];
    for(int i = 0; i < my_size; i++)
      container[i] = rhs[i];
    return *this;
  }

  /**
  * Destructor
  */
  ~array<T>(){
    delete[] container;
  }

  /**
  * Method to add element to end of array
  */
  void push(T el);

  /**
  * Gets an element from the array
  * @param i - index of element. i < mysize
  */
  T get(int i){ if(i < my_size) return container[i]; else throw std::runtime_error("index out of bounds"); };

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
  int size() { return my_size; };

  /**
  * Method to get the capacity of the array
  */
  int capacity() { return my_capacity; };

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
  * Overload of plus operator to add two vectors
  * @param rhs - another array to add to this one
  */
  array<T> operator+(const array<T> rhs){
    array<T> n(my_size);
    for(int i = 0; i < my_size; i++)
      n.push(container[i] + rhs.container[i]);

    return n;
  }

  /**
  * Overload of minus operator to subtract two vectors
  * @param rhs - another array to subtract from this one
  */
  array<T> operator-(const array<T> rhs){
    array<T> n(my_size);
    for(int i = 0; i < my_size; i++)
      n.push(container[i] - rhs.container[i]);

    return n;
  }

  /**
  * Overload of mult operator for dot product
  * @param rhs - another array to dot
  */
  T operator*(const array<T>& rhs){
    T product = 0;
    for (int i = 0; i < my_size; i++) {
      product += container[i] * rhs[i];
    }

    return product;
  }

  /**
  * Overload of mult operator for scalar
  * @param rhs - value to mult this array by
  */
  array<T> operator*(const T& rhs){
    array<T> n = *this;
    for(int i = 0; i < mysize; i++){
      n[i] *= rhs;
    }

    return n;
  }

  /**
  * Prints a string representation of array
  */
  std::string to_string(){
    std::stringstream ss;
    ss << "[ ";
    for(int i = 0; i < mysize; i++){
      ss << container[i] << " ";
    }
    ss << " ]^T";

    std::string output = ss.str();

    return output;
  }

  /**
  * @brief An iterator class for array
  * @details This class enables use of array in
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
/**
* Implementation of the private method array::grow()
*/
template<typename T>
void array<T>::grow(){
  // Set the new capacity to double the current. 0 -> 2, 2 -> 4, ...
  // Doubling the array every grow leads to amortized O(n) grow operations
  array<T>::my_capacity = array<T>::my_capacity == 0 ? 2 : 2 * array<T>::my_capacity;

  // Initialize a temporary primative
  // array to shuffle old array for resize
  T* tmp = new T[array<T>::my_capacity];

  // Copy the old array to the temporary array
  std::copy(array<T>::container, array<T>::container + array<T>::my_size, tmp);

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
  array<T>::my_capacity = capacity;

  // Initialize temporary array
  T* tmp = new T[array<T>::my_capacity];

  // Copy the old array to smaller temp array
  std::copy(array<T>::container, array<T>::container + array<T>::my_size, tmp);

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
  if(array<T>::my_size >= array<T>::my_capacity){
    array<T>::grow();
  }

  // Add element to end of array
  array<T>::container[array<T>::my_size] = el;

  // Increment size
  array<T>::my_size++;
};

/**
* Implementation of public method pop
*/
template<typename T>
T array<T>::pop(){
  // Shrink array
  array<T>::shrink(array<T>::my_size - 1 == 0 ? 0 : array<T>::my_size);

  // Return the last element decrementing size
  return array<T>::container[--my_size];
};

/**
* Implementation of public method clear
*/
template<typename T>
void array<T>::clear(){
  // Set size to 0
  array<T>::my_size = 0;

  // Set capacity to 0
  array<T>::my_capacity = 0;

  // Destroy elements
  delete[] array<T>::container;
};

}

#endif
