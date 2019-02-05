#include <plll.hpp>
#include <iostream>

void printShortestVector(const plll::linalg::math_matrix<plll::arithmetic::Integer> & A, 
                         unsigned index, const plll::arithmetic::Integer & sqnorm)
{
    std::cout << "The currently shortest vector is " << A.row(index) 
              << " with squared norm " << sqnorm << " (Integer)\n";
}

void printShortestVector_LI(const plll::linalg::math_matrix<plll::arithmetic::NInt<long> > & A, 
                            unsigned index, const plll::arithmetic::NInt<long> & sqnorm)
{
    std::cout << "The currently shortest vector is " << A.row(index) 
              << " with squared norm " << sqnorm << " (NInt<long>)\n";
}

int main()
{
    plll::linalg::math_matrix<plll::arithmetic::Integer> A;
    std::cin >> A;
    
    plll::LatticeReduction lr;
    lr.setLattice(A);
    lr.setMinCallbackFunction(printShortestVector, printShortestVector_LI);
    lr.bkz(0.99, 40);
}
