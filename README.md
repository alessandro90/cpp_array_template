# C++ Array and Matrix Template Library

This is a header-only C++ small template library. 
It contains three main classes. An Array class and a Matrix class; both of them inherits from AbstractArray (an abstract class).
To test the library,  you will need a C++14 compiler.

The classes have their own iterators and support the range-based for
loop syntax.

Among the supported operations there are algebraic addition, subtraction and multiplication for array-array, array-matrix and matrix-matrix. It is also possible to stack arrays with matrices, arrays with arrays and matrices with matrices.

#### Dependencies
* [MKL](https://software.intel.com/en-us/mkl) library.