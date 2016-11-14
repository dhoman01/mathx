#ifndef MATRIX_HPP
#define MATRIX_HPP

#include "goodrand.hpp"
#include <cfloat>
#include <iostream>
#include <iomanip>
#include <sstream>

namespace mathx {

/**
* @brief This class is a variable size, random-access data structure for matrices
*/
template<class T>
class matrix {
private:
  /**
  * Base dynamic 2D primative array
  */
  T** container;

  /**
  * Number of rows in the matrix
  */
  int col;

  /**
  * Number of columns in the matrix
  */
  int row;
public:
  /**
  * Default constructor
  */
  matrix<T>(): container(new T*[0]), col(0), row(0){};

  /**
  * Constructor initializing container to r rows and c columns. If rand is true the matrix will be initialized with random numbers
  * @param r - int value to set the number of rows
  * @param c - int value to set the number of columns
  * @param rand = true - bool to initialize with random numbers
  */
  matrix<T>(int r, int c, bool rand = false): container(new T*[r]), col(c), row(r){
    for(int i = 0; i < r; i++){
      container[i] = new T[c];
    }
    if(rand){
      for(int i = 0; i < row; i++)
        for(int j = 0; j < col; j++)
          container[i][j] = (i == j ? 10 * i + goodrand::getRand(1.0, 2.0) : goodrand::getRand(0.0, 1.0));

    }
  };

  /**
  * Constructor initializing container to r rows and c columns to value v
  * @param r - int value to set the number of rows
  * @param c - int value to set the number of columns
  * @param rand = true - bool to initialize with random numbers
  */
  matrix<T>(int r, int c, T v): container(new T*[r]), col(c), row(r){
    for(int i = 0; i < r; i++){
      container[i] = new T[c];
      for(int j = 0; j < c; j++){
        container[i][j] = v;
      }
    }
  };

  /**
  * Constructor accepting an initializer list
  * @details Allows assignment of matrix like matrix<int> arr = {{1,2,3},{4,5,6}};
  * @param c - an initializer list (i.e {{0,1,2},{3,4,5}})
  */
  matrix<T>(std::initializer_list<std::initializer_list<T>> c){
    container = new T*[c.size()];
    row = c.size();
    col = c.begin()->size();
    int i = 0;
    for(auto c_sub : c){
      if((int)c_sub.size() != col) throw std::runtime_error("All rows must have the same number of elements");
      container[i] = new T[c_sub.size()];
      std::copy(c_sub.begin(), c_sub.end(), container[i]);
      i++;
    }
  }

  /**
  * Overload of matrix index operators. Same as get(i)
  * @param i - index of element. i < mysize
  */
  T* operator[](std::size_t i){ return container[i]; };

  /**
  * Overload of matrix index operators. Same as get(i)
  * @param i - index of element. i < mysize
  */
  T* operator[](std::size_t i) const { return container[i]; };

  bool hasPivoted = false;

  /**
  * Get value at location r,c
  * @param r - row position
  * @param c - column position
  */
  T get(int r, int c){
    return container[r][c];
  }

  /**
  * Set the value at r,c to v
  * @param r - row position
  * @param c - column position
  * @param v - new value
  */
  void set(int r, int c, T v){
    container[r][c] = v;
  }

  /**
  * Check if matrix is symmetric
  */
  bool is_symmetric(){
    for(int i = 0; i < row; i++){
      for(int j = 0; j < col; j++){
        if(container[i][j] != container[j][i])
          return false;
      }
    }

    return true;
  }

  /**
  * Find a pivot using the parital pivoting strategy
  * @param k - the current pass for GE or LU Decomposition
  */
  int find_pivot(int k){
    // Find the best pivot (best is max)
    T qmax = std::abs(container[k][k]);
    int kpiv = k;
    for(int i = k + 1; i < row; i++){
      T qtemp = std::abs(container[i][k]);
      if(qtemp > qmax){
        kpiv = i;
        qmax = qtemp;
        hasPivoted = true;
      }
    }

    return kpiv;
  }

  /**
  * Find a pivot using the scaled parital pivoting strategy
  * @param k - the current pass for GE or LU Decomposition
  */
  int find_scaled_pivot(int k){
    // Find the scale
    array<T> s(row);
    for(int i = 0; i < row; i++){
      s[i] = std::abs(container[i][0]);
      for(int j = 0; j < col; j++){
        if(std::abs(container[i][j]) > s[i])
          s[i] = std::abs(container[i][j]);
      }
    }

    // Find the best pivot (best is max)
    T qmax = std::abs(container[k][k]) / s[k];
    int kpiv = k;
    for(int i = k + 1; k < row; k++){
      T qtmp = std::abs(container[i][k]) / s[i];
      if(qtmp > qmax){
        kpiv = i;
        qmax = qtmp;
        hasPivoted = true;
      }
    }

    return kpiv;
  }

  /**
  * Swap two rows of the matrix
  * @param r1 - row one position
  * @param r2 - row two position
  */
  void swap(int r1, int r2){
    std::swap(container[r1], container[r2]);
  }

  /**
  * Get the number columns in the matrix
  */
  int cols(){ return col; };

  /**
  * Get the number of rows in the matrix
  */
  int rows(){ return row; };

  /**
  * Calculate the one-norm of the matrix
  */
  T one_norm(){
    T max = 0;
    for(int j = 0; j < col; j++){
      T sum = 0;
      for(int i = 0; i < row; i++)
        sum += container[i][j];

      if(sum > max)
        max = sum;
    }

    return max;
  }

  /**
  * Calculate the infinity-norm of the matrix
  */
  T infinity_norm(){
    T max = 0;
    for(int i = 0; i < row; i++){
      T sum = 0;
      for(int j = 0; j < col; j++)
        sum += container[i][j];
      if(sum > max)
        max = sum;
    }

    return max;
  }

  /**
  * Returns a string representation of the matrix
  */
  std::string to_string(){
    std::stringstream ss;
    for(int i = 0; i < row; i++){
      for(int j = 0; j < col; j++)
        ss << std::setw(10) << std::left << container[i][j] << " ";
      ss << std::endl;
    }

    return ss.str();
  }
};

};

#endif
