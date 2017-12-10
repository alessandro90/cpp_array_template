This is a header-only c++ small template library. 
It contains two main classes. An Array class and a Matrix class.
Some of the functions make use of the mkl intel library functions.
To test the library, other than the mkl, you will need a C++14
compiler.
I don't have properly tested all the functions yet.
The current version doesn't make use of smart pointers to allocate memory. 

The classes have their own iterators and support the range-based for
loop syntax.