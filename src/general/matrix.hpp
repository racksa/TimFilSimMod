// matrix.hpp

// =============================================================================
// Include guard
#ifndef MY_MATRIX_HEADER_INCLUDED
#define MY_MATRIX_HEADER_INCLUDED

// =============================================================================
// Forward declared dependencies

// =============================================================================
// Included dependencies
#include <iostream>
#include <vector>

/*

This matrix class is intended to facilitate replacing our reliance on external
libraries (like Armadillo, Eigen etc.) with direct use of BLAS and LAPACK routines.
Whilst this gives us greater flexibility in terms of the machines we can run on,
it involves re-inventing the wheel a bit.

I'm also not going to include bound-checking etc., so it's on us to make sure we
don't do silly things like try to add two matrices of different sizes.

*/

class matrix {

public:

  int num_rows;
  int num_cols;

  // Storage is guaranteed to be contiguous in std::vector, so this is fine.
  // It has the advantage over "new + delete[]" that the standard library handles all memory for us.
  std::vector<double> data;

  ~matrix();
  matrix();
  matrix(const int num_rows_input, const int num_cols_input);
  matrix(const int num_rows_input, const int num_cols_input, const double *const memptr);
  matrix(const matrix& M);

  //
  // OPERATOR OVERLOADING
  //

  // Element access
  double& operator()(const int row, const int col);
  double operator()(const int row, const int col) const;
  double& operator()(const int elem);
  double operator()(const int elem) const;

  // Assignment
  matrix& operator =(const matrix& X); // M = X

  // In-place scalar multiplication
  matrix& operator *=(const double s); // M *= s

  // In-place scalar division
  matrix& operator /=(const double s); // M /= s

  // Sign swap
  matrix operator -() const; // For calls of the form A = -B;

  // N.B. No in-place matrix multiplication because the size could change!

  // In-place matrix addition
  matrix& operator +=(const matrix& rhs); // M += rhs

  // In-place matrix subtraction
  matrix& operator -=(const matrix& rhs); // M -= rhs

  // Conversion of 1x1 matrix to a double. Could be useful for dot products.
  operator double() const;

  //
  // OTHER METHODS
  //

  void zero();
  void identity();
  void invert();
  void swap(matrix& M);
  bool is_finite() const;
  double trace() const;

  // 2-index block methods
  matrix get_block(const int start_row, const int start_col, const int rows_in_block, const int cols_in_block) const;
  void get_block(const int start_row, const int start_col, const int rows_in_block, const int cols_in_block, matrix& block) const;
  void set_block(const int start_row, const int start_col, const int rows_in_block, const int cols_in_block, const matrix& block);
  void set_block(const int start_row, const int start_col, const int rows_in_block, const int cols_in_block, const double block_constant);
  void add_to_block(const int start_row, const int start_col, const int rows_in_block, const int cols_in_block, const matrix& block);
  void subtract_from_block(const int start_row, const int start_col, const int rows_in_block, const int cols_in_block, const matrix& block);
  void multiply_block(const int start_row, const int start_col, const int rows_in_block, const int cols_in_block, const double block_constant);
  void divide_block(const int start_row, const int start_col, const int rows_in_block, const int cols_in_block, const double block_constant);

  // Single-index block methods. Useful for vectors.
  matrix get_block(const int start_elem, const int elems_in_block) const;
  void get_block(const int start_elem, const int elems_in_block, matrix& block) const;
  void set_block(const int start_elem, const int elems_in_block, const matrix& block);
  void set_block(const int start_elem, const int elems_in_block, const double block_constant);
  void add_to_block(const int start_elem, const int elems_in_block, const matrix& block);
  void subtract_from_block(const int start_elem, const int elems_in_block, const matrix& block);
  void multiply_block(const int start_elem, const int elems_in_block, const double block_constant);
  void divide_block(const int start_elem, const int elems_in_block, const double block_constant);

  // Special cases of the above block methods for when the block is a given row or column
  matrix get_col(const int col) const;
  void get_col(const int col, matrix& block) const;
  void set_col(const int col, const matrix& block);
  void set_col(const int col, const double block_constant);
  void add_to_col(const int col, const matrix& block);
  void subtract_from_col(const int col, const matrix& block);
  void multiply_col(const int col, const double block_constant);
  void divide_col(const int col, const double block_constant);
  matrix get_row(const int row) const;
  void get_row(const int row, matrix& block) const;
  void set_row(const int row, const matrix& block);
  void set_row(const int row, const double block_constant);
  void add_to_row(const int row, const matrix& block);
  void subtract_from_row(const int row, const matrix& block);
  void multiply_row(const int row, const double block_constant);
  void divide_row(const int row, const double block_constant);

}; // End of matrix class definition

//
// BINARY OPERATOR OVERLOADING
//

matrix operator /(matrix M, const double s);
matrix operator *(matrix M, const double s);
matrix operator *(const double s, matrix M);
matrix operator *(const matrix& lhs, const matrix& rhs);
matrix operator +(matrix lhs, const matrix& rhs);
matrix operator -(matrix lhs, const matrix& rhs);
std::ostream& operator <<(std::ostream& stream, const matrix& M);

//
// OTHER FUNCTIONS ASSOCIATED WITH MATRICES
//

void rcross(matrix& M, const double *const v);
matrix rcross(const double *const v);
void rcross(matrix& M, const matrix& v);
matrix rcross(const matrix& v);
matrix inverse(matrix M);
double dot(const matrix& A, const matrix& B);
double norm(const matrix& M);
matrix identity(const int dim);
matrix zero(const int dim);
matrix zero(const int num_rows_input, const int num_cols_input);
matrix cross(const double *const v, const double *const w);
matrix cross(const matrix& v, const double *const w);
matrix cross(const double *const v, const matrix& w);
matrix cross(const matrix& v, const matrix& w);
void transpose(matrix& At, const matrix& A);
matrix transpose(const matrix& A);
double trace(const matrix& A);

#endif // MY_MATRIX_HEADER_INCLUDED
