#include "../qc.hpp"
#include <iostream>
#include <cassert>

void test_three_qubit_entanglement() {
    qubit q0("q0");
    qubit q1("q1");
    qubit q2("q2");
    circuit c({&q0, &q1, &q2});
    
    // Create GHZ state: (|000> + |111>)/√2
    c.H(q0);
    c.CNOT(q0, q1);
    c.CNOT(q1, q2);
    
    std::cout << "Test: Three-qubit GHZ state\n";
    std::cout << c.print_circuit();
    std::cout << c.print_amplitudes() << "\n";
    
    quantum_computer qc(c);
    std::cout << "Probabilities:\n" << qc.print_prob_dist(qc.calculate_probabilities()) << "\n";
    
    // Should have 50% |000> and 50% |111>
    auto probs = qc.calculate_probabilities();
    assert(probs.size() == 2);
    std::cout << "✓ Passed\n\n";
}

void test_quantum_interference_hzh() {
    qubit q0("q0");
    circuit c({&q0});
    
    // HZH sequence creates X gate (bit flip)
    c.H(q0);
    c.Z(q0);
    c.H(q0);
    
    std::cout << "Test: Quantum interference (HZH = X)\n";
    std::cout << c.print_circuit();
    std::cout << c.print_amplitudes() << "\n";
    
    // Should end up in |1> due to interference
    assert(std::abs(c.amplitudes[0]) < 1e-10);
    assert(std::abs(c.amplitudes[1] - 1.0) < 1e-10);
    std::cout << "✓ Passed\n\n";
}

void test_hxh_identity() {
    qubit q0("q0");
    circuit c({&q0});
    
    // HXH returns to |0>
    c.H(q0);
    c.X(q0);
    c.H(q0);
    
    std::cout << "Test: HXH returns to |0>\n";
    std::cout << c.print_circuit();
    std::cout << c.print_amplitudes() << "\n";
    
    // Should be back at |0>
    assert(std::abs(c.amplitudes[0] - 1.0) < 1e-10);
    assert(std::abs(c.amplitudes[1]) < 1e-10);
    std::cout << "✓ Passed\n\n";
}

void test_xhx_phase_flip() {
    qubit q0("q0");
    circuit c({&q0});
    
    // XHX = H (but with global phase)
    c.X(q0);
    c.H(q0);
    c.X(q0);
    
    std::cout << "Test: XHX creates superposition\n";
    std::cout << c.print_circuit();
    std::cout << c.print_amplitudes() << "\n";
    
    // Should create equal superposition (with possible phase)
    double prob0 = c.amplitudes[0] * c.amplitudes[0];
    double prob1 = c.amplitudes[1] * c.amplitudes[1];
    assert(std::abs(prob0 - 0.5) < 1e-10);
    assert(std::abs(prob1 - 0.5) < 1e-10);
    std::cout << "✓ Passed\n\n";
}

int main() {
    test_three_qubit_entanglement();
    test_quantum_interference_hzh();
    test_hxh_identity();
    test_xhx_phase_flip();
    std::cout << "All complex circuit tests passed!\n";
    return 0;
}