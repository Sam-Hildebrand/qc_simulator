#include "../qc.hpp"
#include <iostream>
#include <cassert>
#include <cmath>

void test_basic_apply() {
    std::cout << "\n=== Test 1: Basic Single Gate Application ===" << std::endl;
    
    // Create a simple 1-qubit circuit with H gate
    qubit q1, q2;
    circuit c1({&q1, &q2});
    
    qubit sub_q;
    circuit sub_c({&sub_q}, "H-gate");
    sub_c.H(sub_q);
    
    std::cout << "Sub-circuit:" << std::endl;
    std::cout << sub_c.print_circuit() << std::endl;
    
    // Apply it to q1 in main circuit
    c1.apply_circuit(sub_c, {&q1});
    
    std::cout << "Main circuit after applying sub-circuit to q1:" << std::endl;
    std::cout << c1.print_circuit() << std::endl;
    std::cout << "Amplitudes: " << c1.print_amplitudes() << std::endl;
    
    // Verify: q1 should be in superposition
    double expected = 1.0 / std::sqrt(2.0);
    assert(std::abs(c1.amplitudes[0] - expected) < 1e-10);
    assert(std::abs(c1.amplitudes[2] - expected) < 1e-10); // q1 is MSB
    std::cout << "✓ Test passed" << std::endl;
}

void test_cnot_apply() {
    std::cout << "\n=== Test 2: CNOT Application ===" << std::endl;
    
    // Create a 2-qubit sub-circuit with CNOT
    qubit sq0, sq1;
    circuit sub_c({&sq0, &sq1}, "CNOT");
    sub_c.CNOT(sq0, sq1);
    
    std::cout << "Sub-circuit:" << std::endl;
    std::cout << sub_c.print_circuit() << std::endl;
    
    // Create main circuit, prepare |10> state
    qubit q0, q1, q2;
    circuit c1({&q0, &q1, &q2});
    c1.X(q0);  // q0 = 1
    
    std::cout << "Main circuit before:" << std::endl;
    std::cout << c1.print_circuit() << std::endl;
    std::cout << "Amplitudes: " << c1.print_amplitudes() << std::endl;
    
    // Apply CNOT to q0, q1
    c1.apply_circuit(sub_c, {&q0, &q1});
    
    std::cout << "\nMain circuit after applying CNOT to q0,q1:" << std::endl;
    std::cout << c1.print_circuit() << std::endl;
    std::cout << "Amplitudes: " << c1.print_amplitudes() << std::endl;
    
    // Should now be |110> (q0 controls q1)
    assert(std::abs(c1.amplitudes[0b110] - 1.0) < 1e-10);
    std::cout << "✓ Test passed" << std::endl;
}

void test_multiple_cnots_same_column() {
    std::cout << "\n=== Test 3: Multiple CNOTs in Same Column ===" << std::endl;
    
    // Create a 4-qubit sub-circuit with 2 CNOTs
    qubit sq0, sq1, sq2, sq3;
    circuit sub_c({&sq0, &sq1, &sq2, &sq3}, "Double-CNOT");
    sub_c.CNOT(sq0, sq1);  // These will be in the same column
    sub_c.CNOT(sq2, sq3);  // because they're independent
    
    std::cout << "Sub-circuit (2 CNOTs):" << std::endl;
    std::cout << sub_c.print_circuit() << std::endl;
    
    // Create main circuit, prepare |1010> state
    qubit q0, q1, q2, q3;
    circuit c1({&q0, &q1, &q2, &q3});
    c1.X(q0);  // q0 = 1
    c1.X(q2);  // q2 = 1
    
    std::cout << "Main circuit before:" << std::endl;
    std::cout << c1.print_circuit() << std::endl;
    std::cout << "Amplitudes: " << c1.print_amplitudes() << std::endl;
    
    // Apply both CNOTs
    c1.apply_circuit(sub_c, {&q0, &q1, &q2, &q3});
    
    std::cout << "\nMain circuit after:" << std::endl;
    std::cout << c1.print_circuit() << std::endl;
    std::cout << "Amplitudes: " << c1.print_amplitudes() << std::endl;
    
    // Should now be |1111> (both CNOTs fired)
    assert(std::abs(c1.amplitudes[0b1111] - 1.0) < 1e-10);
    std::cout << "✓ Test passed" << std::endl;
}

void test_mixed_gates() {
    std::cout << "\n=== Test 4: Mixed Gates (H, X, CNOT) ===" << std::endl;
    
    // Create a Bell state sub-circuit
    qubit sq0, sq1;
    circuit bell({&sq0, &sq1}, "Bell");
    bell.H(sq0);
    bell.CNOT(sq0, sq1);
    
    std::cout << "Bell sub-circuit:" << std::endl;
    std::cout << bell.print_circuit() << std::endl;
    std::cout << "Bell amplitudes: " << bell.print_amplitudes() << std::endl;
    
    // Apply to main circuit
    qubit q0, q1, q2;
    circuit c1({&q0, &q1, &q2});
    c1.apply_circuit(bell, {&q1, &q2});
    
    std::cout << "\nMain circuit after applying Bell to q1,q2:" << std::endl;
    std::cout << c1.print_circuit() << std::endl;
    std::cout << "Amplitudes: " << c1.print_amplitudes() << std::endl;
    
    // Should have |000> and |011> with equal amplitude
    double expected = 1.0 / std::sqrt(2.0);
    assert(std::abs(c1.amplitudes[0b000] - expected) < 1e-10);
    assert(std::abs(c1.amplitudes[0b011] - expected) < 1e-10);
    std::cout << "✓ Test passed" << std::endl;
}

void test_sequential_applications() {
    std::cout << "\n=== Test 5: Sequential Applications ===" << std::endl;
    
    // Create two sub-circuits
    qubit sq1;
    circuit h_gate({&sq1}, "H");
    h_gate.H(sq1);
    
    qubit sq2;
    circuit x_gate({&sq2}, "X");
    x_gate.X(sq2);
    
    // Apply them sequentially
    qubit q0, q1;
    circuit c1({&q0, &q1});
    
    c1.apply_circuit(h_gate, {&q0});
    c1.apply_circuit(x_gate, {&q0});
    c1.apply_circuit(h_gate, {&q1});
    
    std::cout << "Circuit after H→X on q0, H on q1:" << std::endl;
    std::cout << c1.print_circuit() << std::endl;
    std::cout << "Amplitudes: " << c1.print_amplitudes() << std::endl;
    
    std::cout << "✓ Test passed" << std::endl;
}

void test_unnamed_circuit() {
    std::cout << "\n=== Test 6: Unnamed Circuit (no box) ===" << std::endl;
    
    // Create unnamed sub-circuit
    qubit sq0;
    circuit unnamed({&sq0}); // No name
    unnamed.H(sq0);
    
    qubit q0, q1;
    circuit c1({&q0, &q1});
    c1.apply_circuit(unnamed, {&q0});
    
    std::cout << "Circuit after applying unnamed sub-circuit:" << std::endl;
    std::cout << c1.print_circuit() << std::endl;
    std::cout << "(No box should appear, just the gates)" << std::endl;
    
    std::cout << "✓ Test passed" << std::endl;
}

void test_complex_circuit() {
    std::cout << "\n=== Test 7: Complex Multi-Application ===" << std::endl;
    
    // Create a swap-like circuit
    qubit sq0, sq1;
    circuit swap_circuit({&sq0, &sq1}, "SWAP");
    swap_circuit.CNOT(sq0, sq1);
    swap_circuit.CNOT(sq1, sq0);
    swap_circuit.CNOT(sq0, sq1);
    
    std::cout << "SWAP circuit:" << std::endl;
    std::cout << swap_circuit.print_circuit() << std::endl;
    
    // Main circuit with initial state
    qubit q0, q1, q2, q3;
    circuit c1({&q0, &q1, &q2, &q3});
    c1.X(q0);  // Start with |1000>
    
    std::cout << "\nBefore:" << std::endl;
    std::cout << c1.print_circuit() << std::endl;
    std::cout << "Amplitudes: " << c1.print_amplitudes() << std::endl;
    
    // Apply swap to q0 and q1
    c1.apply_circuit(swap_circuit, {&q0, &q1});
    
    std::cout << "\nAfter SWAP on q0,q1:" << std::endl;
    std::cout << c1.print_circuit() << std::endl;
    std::cout << "Amplitudes: " << c1.print_amplitudes() << std::endl;
    
    // Should now be |0100> (swapped)
    assert(std::abs(c1.amplitudes[0b0100] - 1.0) < 1e-10);
    std::cout << "✓ Test passed" << std::endl;
}

int main() {
    std::cout << "Testing apply_circuit functionality..." << std::endl;
    
    test_basic_apply();
    test_cnot_apply();
    test_multiple_cnots_same_column();
    test_mixed_gates();
    test_sequential_applications();
    test_unnamed_circuit();
    test_complex_circuit();
    
    std::cout << "\n=== All tests passed! ===" << std::endl;
    return 0;
}