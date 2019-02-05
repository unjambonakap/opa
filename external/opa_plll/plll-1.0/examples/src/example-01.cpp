#include <plll.hpp>
#include <iostream>

int main()
{
    plll::linalg::math_matrix<plll::arithmetic::Integer> A;
    std::cin >> A;
    
    plll::LatticeReduction lr;
    lr.setLattice(A);
    lr.lll();
    
    std::cout << lr.getLattice() << "\n";
}
