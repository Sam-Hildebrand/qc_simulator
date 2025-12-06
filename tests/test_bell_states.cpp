#include "../qc.hpp"
#include <iostream>
#include <cassert>
#include <cmath>

void test_bell_phi_plus() {
    qubit q0("alice");
    qubit q1("bob");
    circuit c({&q0, &q1});
    
    c.H(q0);
    c.CNOT(q0, q1);
    
    std::cout << "Test: Bell state |Φ+> = (|00> + |11>)/√2\n";
    std::cout << c.print_circuit();
    std::cout << c.print_amplitudes() << "\n";
    
    quantum_computer qc(c);
    std::cout << "Probabilities:\n" << qc.print_prob_dist(qc.calculate_probabilities()) << "\n";
    
    double sqrt2 = 1.0 / std::sqrt(2.0);
    assert(std::abs(c.amplitudes[0] - sqrt2) < 1e-10); // |00>
    assert(std::abs(c.amplitudes[3] - sqrt2) < 1e-10); // |11>
    std::cout << "✓ Passed\n\n";
}

void test_bell_phi_minus() {
    qubit q0("alice");
    qubit q1("bob");
    circuit c({&q0, &q1});
    
    c.X(q1);
    c.H(q0);
    c.CNOT(q0, q1);
    
    std::cout << "Test: Bell state |Φ-> = (|01> + |10>)/√2\n";
    std::cout << c.print_circuit();
    std::cout << c.print_amplitudes() << "\n";
    
    quantum_computer qc(c);
    std::cout << "Probabilities:\n" << qc.print_prob_dist(qc.calculate_probabilities()) << "\n";
    
    double sqrt2 = 1.0 / std::sqrt(2.0);
    assert(std::abs(c.amplitudes[1] - sqrt2) < 1e-10); // |01>
    assert(std::abs(c.amplitudes[2] - sqrt2) < 1e-10); // |10>
    std::cout << "✓ Passed\n\n";
}

void test_bell_psi_plus() {
    qubit q0("alice");
    qubit q1("bob");
    circuit c({&q0, &q1});
    
    c.X(q0);
    c.H(q0);
    c.CNOT(q0, q1);
    
    std::cout << "Test: Bell state |Ψ+> = (|00> - |11>)/√2\n";
    std::cout << c.print_circuit();
    std::cout << c.print_amplitudes() << "\n";
    
    quantum_computer qc(c);
    std::cout << "Probabilities:\n" << qc.print_prob_dist(qc.calculate_probabilities()) << "\n";
    
    double sqrt2 = 1.0 / std::sqrt(2.0);
    assert(std::abs(c.amplitudes[0] - sqrt2) < 1e-10); // |00>
    assert(std::abs(c.amplitudes[3] + sqrt2) < 1e-10); // -|11>
    std::cout << "✓ Passed\n\n";
}

int main() {
    test_bell_phi_plus();
    test_bell_phi_minus();
    test_bell_psi_plus();
    std::cout << "All Bell state tests passed!\n";
    return 0;
}