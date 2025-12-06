#include "../qc.hpp"
#include <iostream>
#include <cassert>
#include <cmath>

void test_cnot_00() {
    qubit q0("control");
    qubit q1("target");
    circuit c({&q0, &q1});
    
    c.CNOT(q0, q1);
    
    std::cout << "Test: CNOT on |00> (no flip)\n";
    std::cout << c.print_circuit();
    std::cout << c.print_amplitudes() << "\n";
    
    // Should remain |00>
    assert(std::abs(c.amplitudes[0] - 1.0) < 1e-10);
    std::cout << "✓ Passed\n\n";
}

void test_cnot_10() {
    qubit q0("control");
    qubit q1("target");
    circuit c({&q0, &q1});
    
    c.X(q0);
    c.CNOT(q0, q1);
    
    std::cout << "Test: CNOT on |10> (flips to |11>)\n";
    std::cout << c.print_circuit();
    std::cout << c.print_amplitudes() << "\n";
    
    // Should be |11>
    assert(std::abs(c.amplitudes[3] - 1.0) < 1e-10);
    std::cout << "✓ Passed\n\n";
}

void test_cnot_01() {
    qubit q0("control");
    qubit q1("target");
    circuit c({&q0, &q1});
    
    c.X(q1);
    c.CNOT(q0, q1);
    
    std::cout << "Test: CNOT on |01> (no flip, control is 0)\n";
    std::cout << c.print_circuit();
    std::cout << c.print_amplitudes() << "\n";
    
    // Should remain |01>
    assert(std::abs(c.amplitudes[1] - 1.0) < 1e-10);
    std::cout << "✓ Passed\n\n";
}

void test_cnot_11() {
    qubit q0("control");
    qubit q1("target");
    circuit c({&q0, &q1});
    
    c.X(q0);
    c.X(q1);
    c.CNOT(q0, q1);
    
    std::cout << "Test: CNOT on |11> (flips to |10>)\n";
    std::cout << c.print_circuit();
    std::cout << c.print_amplitudes() << "\n";
    
    // Should be |10>
    assert(std::abs(c.amplitudes[2] - 1.0) < 1e-10);
    std::cout << "✓ Passed\n\n";
}

int main() {
    test_cnot_00();
    test_cnot_10();
    test_cnot_01();
    test_cnot_11();
    std::cout << "All CNOT gate tests passed!\n";
    return 0;
}