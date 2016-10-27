#ifndef MATRIX_HPP
#define MATRIX_HPP

#include "goodrand.hpp"

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
          container[i][j] = (i == j ? 10 * i + goodrand::getRand(0.0, 1.0) : goodrand::getRand(0.0, 1.0));

    }
  };

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
};

};

#endif
