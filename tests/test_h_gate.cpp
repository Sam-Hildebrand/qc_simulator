#include "../qc.hpp"
#include <iostream>
#include <cassert>
#include <cmath>

void test_h_creates_superposition() {
    qubit q0("q0");
    circuit c({&q0});
    
    c.H(q0);
    
    std::cout << "Test: H creates equal superposition\n";
    std::cout << c.print_circuit();
    std::cout << c.print_amplitudes() << "\n";
    
    double expected = 1.0 / std::sqrt(2.0);
    assert(std::abs(c.amplitudes[0] - expected) < 1e-10);
    assert(std::abs(c.amplitudes[1] - expected) < 1e-10);
    std::cout << "✓ Passed\n\n";
}

void test_double_h() {
    qubit q0("q0");
    circuit c({&q0});
    
    c.H(q0);
    c.H(q0);
    
    std::cout << "Test: HH = I (returns to original state)\n";
    std::cout << c.print_circuit();
    std::cout << c.print_amplitudes() << "\n";
    
    // Should be back to |0>
    assert(std::abs(c.amplitudes[0] - 1.0) < 1e-10);
    assert(std::abs(c.amplitudes[1]) < 1e-10);
    std::cout << "✓ Passed\n\n";
}

void test_h_on_one() {
    qubit q0("q0");
    circuit c({&q0});
    
    c.X(q0);
    c.H(q0);
    
    std::cout << "Test: H on |1> creates superposition with phase\n";
    std::cout << c.print_circuit();
    std::cout << c.print_amplitudes() << "\n";
    
    double sqrt2 = 1.0 / std::sqrt(2.0);
    assert(std::abs(c.amplitudes[0] - sqrt2) < 1e-10);
    assert(std::abs(c.amplitudes[1] + sqrt2) < 1e-10); // negative
    std::cout << "✓ Passed\n\n";
}

int main() {
    test_h_creates_superposition();
    test_double_h();
    test_h_on_one();
    std::cout << "All H gate tests passed!\n";
    return 0;
}