#ifndef ARRAY_HPP
#define ARRAY_HPP

#include <stdexcept>
#include <iostream>

namespace mathx {

template<class T>
class array{
private:
  T* container;
  int mysize;
  int mycapacity;
  void grow();
  void shrink();
public:
  array<T>() : container(new T[0]), mysize(0), mycapacity(0){};
  array<T>(int c) : container(new T[c]), mysize(0), mycapacity(c){};
  array<T>(std::initializer_list<T> c){
    container = new T[c.size()];
    mycapacity = mysize = c.size();
    std::copy(c.begin(), c.end(), container);
  }
  void push(T el);
  T get(int i){ if(i < mysize) return container[i]; else throw std::runtime_error("You suck!"); };
  T pop();
  void clear();
  int size() { return mysize; };
  int capacity() { return mycapacity; };
  T& operator[](std::size_t i){ return container[i]; };
  T& operator[](std::size_t i) const { return container[i]; };
  class iterator {
  public:
    iterator(T* ptr) : ptr(ptr){};
    iterator operator++() { iterator i(ptr); ++ptr; return i; }
    bool operator!=(const iterator& other) { return ptr != other.ptr; }
    const T& operator*() const { return *ptr; }
  private:
    T* ptr;
  };
  iterator begin() { return iterator(container); };
  iterator end() { return iterator(container + mysize); };
};

// PRIVATE METHODS
template<typename T>
void array<T>::grow(){
  array<T>::mycapacity = array<T>::mycapacity == 0 ? 2 : 2 * array<T>::mycapacity;
  T* tmp = new T[array<T>::mycapacity];
  std::copy(array<T>::container, array<T>::container + array<T>::mysize, tmp);
  delete[] array<T>::container;
  array<T>::container = tmp;
};

template<typename T>
void array<T>::shrink(){
  T* tmp = new T[array<T>::mycapacity];
  std::copy(array<T>::container, array<T>::container + array<T>::mysize, tmp);
  delete[] array<T>::container;
  array<T>::container = tmp;
};

// PUBLIC METHODS
template<typename T>
void array<T>::push(T el){
  if(array<T>::mysize >= array<T>::mycapacity){
    array<T>::grow();
  }
  array<T>::container[array<T>::mysize] = el;
  array<T>::mysize++;
};

template<typename T>
T array<T>::pop(){
  array<T>::mycapacity = array<T>::mysize - 1 == 0 ? 0 : array<T>::mysize;
  array<T>::shrink();
  return array<T>::container[--mysize];
};

template<typename T>
void array<T>::clear(){
  array<T>::mysize = 0;
  array<T>::mycapacity = 0;
  delete[] array<T>::container;
};

}

#endif
