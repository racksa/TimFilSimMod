// matrix.cpp

#include <cmath>
#include "matrix.hpp"

matrix::~matrix(){

}

matrix::matrix(){

  num_rows = 1;
  num_cols = 1;

}

matrix::matrix(const int num_rows_input, const int num_cols_input){

  num_rows = num_rows_input;
  num_cols = num_cols_input;

  data = std::vector<double>(num_rows*num_cols);
  data.shrink_to_fit();

}

matrix::matrix(const int num_rows_input, const int num_cols_input, const double *const memptr){

  num_rows = num_rows_input;
  num_cols = num_cols_input;

  data = std::vector<double>(num_rows*num_cols);
  data.shrink_to_fit();

  for (int i = 0; i < num_rows*num_cols; i++){

    data[i] = memptr[i];

  }

}

matrix::matrix(const matrix& M){

  if (&M != this){ // Make sure we're not trying to make a copy of ourselves...

    num_rows = M.num_rows;
    num_cols = M.num_cols;

    data = M.data;
    data.shrink_to_fit();

  }

}

//
// OPERATOR OVERLOADING
//

// Element access
double& matrix::operator()(const int row, const int col){

  // Everything must be column-major to be consistent with BLAS, LAPACK etc.
  return data[row + col*num_rows];

}

double matrix::operator()(const int row, const int col) const {

  // Everything must be column-major to be consistent with BLAS, LAPACK etc.
  return data[row + col*num_rows];

}

double& matrix::operator()(const int elem){

  return data[elem];

}

double matrix::operator()(const int elem) const {

  return data[elem];

}

// Assignment
matrix& matrix::operator =(const matrix& X){

  if (&X != this){ // Make sure we're not trying to make a copy of ourselves...

    num_rows = X.num_rows;
    num_cols = X.num_cols;

    data = X.data;
    data.shrink_to_fit();

  }

  return *this;

}

// Apparently there are no BLAS routines for these basic arithmetic operations,
// so the best thing to do is write them naively and let the compiler optimise everything.
// Of course, we can help by ensuring that memory access is contiguous.

// In-place scalar multiplication
matrix& matrix::operator *=(const double s){

  for (int i = 0; i < num_rows*num_cols; i++){

    data[i] *= s;

  }

  return *this;

}

// In-place scalar division
matrix& matrix::operator /=(const double s){

  for (int i = 0; i < num_rows*num_cols; i++){

    data[i] /= s;

  }

  return *this;

}

// Sign swap
matrix matrix::operator -() const {

  matrix A = *this;
  A *= -1.0;
  return A;

}

// In-place matrix addition
matrix& matrix::operator +=(const matrix& rhs){

  for (int i = 0; i < num_rows*num_cols; i++){

    data[i] += rhs.data[i];

  }

  return *this;

}

// In-place matrix subtraction
matrix& matrix::operator -=(const matrix& rhs){

  for (int i = 0; i < num_rows*num_cols; i++){

    data[i] -= rhs.data[i];

  }

  return *this;

}

// Conversion of 1x1 matrix to a double
matrix::operator double() const{

  if ((num_rows != 1) || (num_cols != 1)){

    std::cout << std::endl << "ERROR: Cannot convert a " << num_rows << "-by-" << num_cols << " matrix to a double." << std::endl << std::endl;
    exit(-1);

  }

  return data[0];

}

//
// OTHER METHODS
//

void matrix::zero(){

  for (int n = 0; n < num_rows*num_cols; n++){

    data[n] = 0.0;

  }

}

void matrix::identity(){

  for (int n = 0; n < num_cols; n++){ // Outer loop should be over columns for contiguous access to column-major ordered data.
    for (int m = 0; m < num_rows; m++){

      (*this)(m,n) = (n == m) ? 1.0 : 0.0;

    }

  }

}

void matrix::swap(matrix& M){

  const int Mrows = M.num_rows;
  const int Mcols = M.num_cols;

  M.num_rows = num_rows;
  M.num_cols = num_cols;

  num_rows = Mrows;
  num_cols = Mcols;

  data.swap(M.data);

}

bool matrix::is_finite() const {

  bool out = true;

  int n = 0;

  while (out && (n < num_rows*num_cols)){

    out = out && std::isfinite(data[n++]);


  }

  return out;

}

double matrix::trace() const {

  double tr = 0.0;

  for (int n = 0; n < num_rows; n++){ // As expected, this method is meaningless if num_rows != num_cols. Indeed, it will even segfault if num_rows > num_cols.

    tr += (*this)(n,n);

  }

  return tr;

}

// 2-index block methods
// Some of these could be implemented using the single-index versions, but at the very least
// it won't work for the matrix-returning get_block, so I haven't bothered.
matrix matrix::get_block(const int start_row, const int start_col, const int rows_in_block, const int cols_in_block) const {

  matrix block(rows_in_block, cols_in_block);

  get_block(start_row, start_col, rows_in_block, cols_in_block, block);

  return block;

}

void matrix::get_block(const int start_row, const int start_col, const int rows_in_block, const int cols_in_block, matrix& block) const {

  for (int col = 0; col < cols_in_block; col++){
    for (int row = 0; row < rows_in_block; row++){

      block(row, col) = (*this)(row + start_row, col + start_col);

    }
  }

}

void matrix::set_block(const int start_row, const int start_col, const int rows_in_block, const int cols_in_block, const matrix& block){

  for (int col = 0; col < cols_in_block; col++){
    for (int row = 0; row < rows_in_block; row++){

      (*this)(row + start_row, col + start_col) = block(row, col);

    }
  }

}

void matrix::set_block(const int start_row, const int start_col, const int rows_in_block, const int cols_in_block, const double block_constant){

  for (int col = 0; col < cols_in_block; col++){
    for (int row = 0; row < rows_in_block; row++){

      (*this)(row + start_row, col + start_col) = block_constant;

    }
  }

}

void matrix::add_to_block(const int start_row, const int start_col, const int rows_in_block, const int cols_in_block, const matrix& block){

  for (int col = 0; col < cols_in_block; col++){
    for (int row = 0; row < rows_in_block; row++){

      (*this)(row + start_row, col + start_col) += block(row, col);

    }
  }

}

void matrix::subtract_from_block(const int start_row, const int start_col, const int rows_in_block, const int cols_in_block, const matrix& block){

  for (int col = 0; col < cols_in_block; col++){
    for (int row = 0; row < rows_in_block; row++){

      (*this)(row + start_row, col + start_col) -= block(row, col);

    }
  }

}

void matrix::multiply_block(const int start_row, const int start_col, const int rows_in_block, const int cols_in_block, const double block_constant){

  for (int col = 0; col < cols_in_block; col++){
    for (int row = 0; row < rows_in_block; row++){

      (*this)(row + start_row, col + start_col) *= block_constant;

    }
  }

}

void matrix::divide_block(const int start_row, const int start_col, const int rows_in_block, const int cols_in_block, const double block_constant){

  this->multiply_block(start_row, start_col, rows_in_block, cols_in_block, 1.0/block_constant);

}

// Single-index block methods. Useful for vectors.
matrix matrix::get_block(const int start_elem, const int elems_in_block) const {

  // Returning a vector is the only guaranteed-safe option.
  // We're only likely to call this on vectors anyway, but if we pass a general
  // matrix, it returns a (probably useless) column vector.
  const int output_rows = (num_rows != 1) ? elems_in_block : 1;
  const int output_cols = (output_rows == 1) ? elems_in_block : 1;

  matrix block(output_rows, output_cols);

  get_block(start_elem, elems_in_block, block);

  return block;

}

void matrix::get_block(const int start_elem, const int elems_in_block, matrix& block) const {

  for (int elem = 0; elem < elems_in_block; elem++){

    block(elem) = data[start_elem + elem];

  }

}

void matrix::set_block(const int start_elem, const int elems_in_block, const matrix& block){

  for (int elem = 0; elem < elems_in_block; elem++){

    data[start_elem + elem] = block(elem);

  }

}

void matrix::set_block(const int start_elem, const int elems_in_block, const double block_constant){

  for (int elem = 0; elem < elems_in_block; elem++){

    data[start_elem + elem] = block_constant;

  }

}

void matrix::add_to_block(const int start_elem, const int elems_in_block, const matrix& block){

  for (int elem = 0; elem < elems_in_block; elem++){

    data[start_elem + elem] += block(elem);

  }

}

void matrix::subtract_from_block(const int start_elem, const int elems_in_block, const matrix& block){

  for (int elem = 0; elem < elems_in_block; elem++){

    data[start_elem + elem] -= block(elem);

  }

}

void matrix::multiply_block(const int start_elem, const int elems_in_block, const double block_constant){

  for (int elem = 0; elem < elems_in_block; elem++){

    data[start_elem + elem] *= block_constant;

  }

}

void matrix::divide_block(const int start_elem, const int elems_in_block, const double block_constant){

  this->multiply_block(start_elem, elems_in_block, 1.0/block_constant);

}

// Special cases of the above block methods for when the block is a given row or column
matrix matrix::get_col(const int col) const {

  return this->get_block(0, col, num_rows, 1);

}

void matrix::get_col(const int col, matrix& block) const {

  this->get_block(0, col, num_rows, 1, block);

}

void matrix::set_col(const int col, const matrix& block){

  this->set_block(0, col, num_rows, 1, block);

}

void matrix::set_col(const int col, const double block_constant){

  this->set_block(0, col, num_rows, 1, block_constant);

}

void matrix::add_to_col(const int col, const matrix& block){

  this->add_to_block(0, col, num_rows, 1, block);

}

void matrix::subtract_from_col(const int col, const matrix& block){

  this->subtract_from_block(0, col, num_rows, 1, block);

}

void matrix::multiply_col(const int col, const double block_constant){

  this->multiply_block(0, col, num_rows, 1, block_constant);

}

void matrix::divide_col(const int col, const double block_constant){

  this->divide_block(0, col, num_rows, 1, block_constant);

}

matrix matrix::get_row(const int row) const {

  return this->get_block(row, 0, 1, num_cols);

}

void matrix::get_row(const int row, matrix& block) const {

  this->get_block(row, 0, 1, num_cols, block);

}

void matrix::set_row(const int row, const matrix& block){

  this->set_block(row, 0, 1, num_cols, block);

}

void matrix::set_row(const int row, const double block_constant){

  this->set_block(row, 0, 1, num_cols, block_constant);

}

void matrix::add_to_row(const int row, const matrix& block){

  this->add_to_block(row, 0, 1, num_cols, block);

}

void matrix::subtract_from_row(const int row, const matrix& block){

  this->subtract_from_block(row, 0, 1, num_cols, block);

}

void matrix::multiply_row(const int row, const double block_constant){

  this->multiply_block(row, 0, 1, num_cols, block_constant);

}

void matrix::divide_row(const int row, const double block_constant){

  this->divide_block(row, 0, 1, num_cols, block_constant);

}

extern "C" {

  // LU decomoposition of a general matrix. Overwrites the input with the output.
  void dgetrf_(int*, int*, double*, int*, int*, int*);

  // Generates the inverse of a matrix given its LU decomposition from dgetrf_. Overwrites the input with the output.
  void dgetri_(int*, double*, int*, int*, double*, int*, int*);

}

void matrix::invert(){

  // Allocate everything we need for the LU-decomposition LAPACK call:
  int N = num_rows;
  double *ptr = &data[0];
  int *pivots = new int[N];
  int info;

  // Call the LAPACK routine:
  dgetrf_(&N, &N, ptr, &N, pivots, &info);

  if (info != 0){

    std::cout << std::endl << "ERROR: Matrix is singular." << std::endl << std::endl;
    exit(-1);

  }

  // Allocate additional stuff for the inversion LAPACK call:
  int work_size = N*N;
  double *work = new double[work_size];

  // Call the LAPACK inversion routine:
  dgetri_(&N, ptr, &N, pivots, work, &work_size, &info);

  // Clean up:
  delete[] pivots;
  delete[] work;

}

//
// BINARY OPERATOR OVERLOADING
//

matrix operator /(matrix M, const double s){

  M /= s;
  return M;

};

matrix operator *(matrix M, const double s){

  M *= s;
  return M;

};

matrix operator *(const double s, matrix M){

  M *= s;
  return M;

};

matrix operator +(matrix lhs, const matrix& rhs){

  lhs += rhs;
  return lhs;

};

matrix operator -(matrix lhs, const matrix& rhs){

  lhs -= rhs;
  return lhs;

};

extern "C" {

  // Returns C = alpha*op(A)*op(B) + beta*C -- i.e. overwrites C -- where op(X) = X or X^T and can be different for A and B.
  // The special case op(X) = X and beta = 0 gives us C = A*B.
  void dgemm_(char*, char*, int*, int*, int*, double*, const double*, int*, const double*, int*, double*, double*, int*);

}

matrix operator *(const matrix& lhs, const matrix& rhs){

  // Allocate space for the solution:
  matrix lhs_times_rhs(lhs.num_rows, rhs.num_cols);

  // Allocate everything we need for the LAPACK call:
  char Nchar = 'N';
  int num_rows = lhs_times_rhs.num_rows;
  int num_cols = lhs_times_rhs.num_cols;
  int inner_dim = lhs.num_cols;
  double alpha = 1.0;
  const double *lhs_ptr = &lhs.data[0];
  int ld_lhs = lhs.num_rows;
  const double *rhs_ptr = &rhs.data[0];
  int ld_rhs = rhs.num_rows;
  double beta = 0.0;
  double *soln_ptr = &lhs_times_rhs.data[0];
  int ld_soln = num_rows;

  // Call the LAPACK routine:
  dgemm_(&Nchar, &Nchar, &num_rows, &num_cols, &inner_dim, &alpha, lhs_ptr, &ld_lhs, rhs_ptr, &ld_rhs, &beta, soln_ptr, &ld_soln);

  // Return the solution:
  return lhs_times_rhs;

};

std::ostream& operator <<(std::ostream& stream, const matrix& M){

  for (int row = 0; row < M.num_rows; row++){

    for (int col = 0; col < M.num_cols; col++){

      stream << M(row, col) << " ";

    }

    stream << std::endl;

  }

  return stream;

};

//
// OTHER FUNCTIONS ASSOCIATED WITH MATRICES
//

void rcross(matrix& M, const double *const v){

  M(0,0) = 0.0;
  M(1,0) = -v[2];
  M(2,0) = v[1];

  M(0,1) = v[2];
  M(1,1) = 0.0;
  M(2,1) = -v[0];

  M(0,2) = -v[1];
  M(1,2) = v[0];
  M(2,2) = 0.0;

};

matrix rcross(const double *const v){

  matrix M(3,3);

  rcross(M, v); // This way, if there's a bug I only have to change one set of expressions.

  return M;

};

void rcross(matrix& M, const matrix& v){

  rcross(M, &v.data[0]);

};

matrix rcross(const matrix& v){

  matrix M(3,3);

  rcross(M, v); // This way, if there's a bug I only have to change one set of expressions.

  return M;

};

matrix inverse(matrix M){

  // Copy-construct by taking input by value not reference, and then call in-place method .invert()
  M.invert();
  return M;

};

double dot(const matrix& A, const matrix& B){

  // Returns the Frobenius inner product, which reduces to the vector dot product for singleton dimensions.
  double out = 0.0;

  // Only makes sense if they're the same size.
  for (int n = 0; n < A.num_rows*A.num_cols; n++){

    out += A(n)*B(n);

  }

  return out;

};

double norm(const matrix& M){

  // Returns the Frobenius norm, which reduces to the Euclidean norm for singleton dimensions.
  return std::sqrt(dot(M,M));

};

matrix identity(const int dim){

  matrix I(dim,dim);
  I.identity();
  return I;

};

matrix zero(const int dim){

  matrix Z(dim,dim);
  Z.zero();
  return Z;

};

matrix zero(const int num_rows_input, const int num_cols_input){

  matrix Z(num_rows_input, num_cols_input);
  Z.zero();
  return Z;

};

matrix cross(const double *const v, const double *const w){

  matrix out(3,1);

  out(0) = v[1]*w[2] - v[2]*w[1];
  out(1) = v[2]*w[0] - v[0]*w[2];
  out(2) = v[0]*w[1] - v[1]*w[0];

  return out;

};

matrix cross(const matrix& v, const double *const w){

  return cross(&v.data[0], w);

};

matrix cross(const double *const v, const matrix& w){

  return cross(v, &w.data[0]);

};

matrix cross(const matrix& v, const matrix& w){

  return cross(&v.data[0], &w.data[0]);

};

void transpose(matrix& At, const matrix& A){

  // As always, read from A in column major order for cached access.
  for (int col = 0; col < A.num_cols; col++){
    for (int row = 0; row < A.num_rows; row++){

      At(col,row) = A(row,col);

    }
  }

};

matrix transpose(const matrix& A){

  // This function probably shouldn't be called for one-off calculations of
  // the form transpose(A)*B, A*transpose(B) or transpose(A)*transpose(B) as dgemm_,
  // the LAPACK routine we use for matrix multiplication, can handle these by changing
  // 'N' flags to 'T' flags. This is bound to be more efficient, and if we find
  // ourselves doing these kinds of multiplications frequently we should define
  // routines like Atranspose_times_B(const matrix&A, const matrix& B). For cases
  // where we repeatedly multiply by the same transposed matrix, this function should
  // be fine.
  //
  // Furthermore, because LAPACK assumes (in general) that matrices are stored in
  // column major order, we don't want to cheat and just use a 'transposed' flag inside
  // the matrix class which switches how the element access operators work and which
  // character is passed to dgemm_ etc. Indeed, because LAPACK assumes it, so do I, and
  // various methods in this class would need to be changed to check for such a flag
  // before e.g. adding two matrices.
  //
  // N.B. Don't use the 'take input by value to use the copy to store the output'
  // trick here; the void-return version doesn't use any intermediate storage and
  // so we can't pass it the same matrix for both inputs.
  matrix At(A.num_cols, A.num_rows);
  transpose(At, A);
  return At;

};

double trace(const matrix& A){

  return A.trace();

};
