#include <iostream>
#include "Array.h"

using namespace act;

Matrix<float> build_matrix()
{
    return {
        {1.2f, 5.6f},
        {-8.6f, 8.8f},
    };
}


int main()
{
    Array<int> a{};
    Array b = {5, 6, 8};
    a = b;
    std::cout << "Array a:\n" << a << "\n\n";

    Matrix<float> m{};
    m = build_matrix();
    std::cout << "Matrix m:\n" << m << "\n\n";

    m.hstack({{5.5f, 6.9f}, {9.99f, 2.0f}});

    std::cout << "Stacked matrix m:\n" << m << "\n";

    Matrix m2 {
        {5.6, 9.},
        {-99., 15.88},
    };

    std::cout << "Matrix m2:\n" << m2 << "\n";

    std::cout << "m2.trace() = " << m2.trace() << "\n";
}