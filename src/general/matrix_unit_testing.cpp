// matrix_unit_testing.cpp

#include "matrix.hpp"

int main(){

    std::cout << std::endl  << "Running tests for the matrix class..." << std::endl << std::endl;

    //
    // Test 1: Can we access all of the memory we should have allocated?
    //
    int num_rows = 100;
    int num_cols = 200;

    matrix A(num_rows, num_cols);

    for (int elem = 0; elem < num_rows*num_cols; elem++){

        A(elem) = double(elem);

    }

    // If we haven't segfaulted, we say this test has passed.
    // I don't think we can try to catch a segfault because it's lower level than a C++ exception.
    std::cout << "Passed test 1" << std::endl;

    //
    // Test 2: Is the memory organised in column-major order?
    //
    double cumulative_error = 0.0;

    for (int row = 0; row < num_rows; row++){

        cumulative_error += std::abs(A(row,0) - double(row)); // Check the first column...
        cumulative_error += std::abs(A(row, num_cols-1) - double(row + num_rows*(num_cols-1))); // ...and the last column.

    }

    if (cumulative_error > 1e-16){ // It should actually be recognised as 0, but we'll allow some room for round-off error etc.

        std::cout << "Failed test 2: ";
        std::cout << "The total error compared to the expected column-major ordered values is " << cumulative_error << std::endl;
        return -1;

    } else {

        std::cout << "Passed test 2" << std::endl;

    }

    


    std::cout << std::endl << "...all tests passed!" << std::endl;

    return 0;

}