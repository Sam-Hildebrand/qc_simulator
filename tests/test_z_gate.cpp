#include "../qc.hpp"
#include <iostream>
#include <cassert>
#include <cmath>

void test_z_on_zero() {
    qubit q0("q0");
    circuit c({&q0});
    
    c.Z(q0);
    
    std::cout << "Test: Z gate on |0> (no effect)\n";
    std::cout << c.print_circuit();
    std::cout << c.print_amplitudes() << "\n";
    
    // Should still be |0>
    assert(std::abs(c.amplitudes[0] - 1.0) < 1e-10);
    std::cout << "✓ Passed\n\n";
}

void test_z_on_one() {
    qubit q0("q0");
    circuit c({&q0});
    
    c.X(q0);
    c.Z(q0);
    
    std::cout << "Test: Z gate on |1> (phase flip)\n";
    std::cout << c.print_circuit();
    std::cout << c.print_amplitudes() << "\n";
    
    // Should be -|1>
    assert(std::abs(c.amplitudes[0]) < 1e-10);
    assert(std::abs(c.amplitudes[1] + 1.0) < 1e-10); // negative amplitude
    std::cout << "✓ Passed\n\n";
}

void test_z_on_superposition() {
    qubit q0("q0");
    circuit c({&q0});
    
    c.H(q0);
    c.Z(q0);
    c.H(q0);
    
    std::cout << "Test: HZH = X (phase flip in superposition)\n";
    std::cout << c.print_circuit();
    std::cout << c.print_amplitudes() << "\n";
    
    // HZH = X, so should flip to |1>
    assert(std::abs(c.amplitudes[0]) < 1e-10);
    assert(std::abs(c.amplitudes[1] - 1.0) < 1e-10);
    std::cout << "✓ Passed\n\n";
}

int main() {
    test_z_on_zero();
    test_z_on_one();
    test_z_on_superposition();
    std::cout << "All Z gate tests passed!\n";
    return 0;
}