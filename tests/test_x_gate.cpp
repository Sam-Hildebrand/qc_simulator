#include "../qc.hpp"
#include <iostream>
#include <cassert>
#include <cmath>

void test_single_x() {
    qubit q0("q0");
    circuit c({&q0});
    
    c.X(q0);
    
    std::cout << "Test: Single X gate\n";
    std::cout << c.print_circuit();
    std::cout << c.print_amplitudes() << "\n";
    
    // Should be |1> with probability 1
    assert(std::abs(c.amplitudes[0]) < 1e-10); // |0> amplitude ~0
    assert(std::abs(c.amplitudes[1] - 1.0) < 1e-10); // |1> amplitude ~1
    std::cout << "✓ Passed\n\n";
}

void test_double_x() {
    qubit q0("q0");
    circuit c({&q0});
    
    c.X(q0);
    c.X(q0);
    
    std::cout << "Test: Double X gate (should return to |0>)\n";
    std::cout << c.print_circuit();
    std::cout << c.print_amplitudes() << "\n";
    
    // Should be back to |0>
    assert(std::abs(c.amplitudes[0] - 1.0) < 1e-10);
    assert(std::abs(c.amplitudes[1]) < 1e-10);
    std::cout << "✓ Passed\n\n";
}

void test_x_on_two_qubits() {
    qubit q0("q0");
    qubit q1("q1");
    circuit c({&q0, &q1});
    
    c.X(q0);
    c.X(q1);
    
    std::cout << "Test: X on both qubits\n";
    std::cout << c.print_circuit();
    std::cout << c.print_amplitudes() << "\n";
    
    // Should be |11>
    assert(std::abs(c.amplitudes[3] - 1.0) < 1e-10); // |11> = index 3
    std::cout << "✓ Passed\n\n";
}

int main() {
    test_single_x();
    test_double_x();
    test_x_on_two_qubits();
    std::cout << "All X gate tests passed!\n";
    return 0;
}