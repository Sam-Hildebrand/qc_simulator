#include <iostream>
#include <stdexcept>
#include <cassert>

#include "../qc.hpp"     // qubit, gate, circuit

int main() {
    std::cout << "[TEST] qubit_in_circuit validation...\n";

    // --- Setup ---
    qubit qa("a");
    qubit qb("b");
    qubit qc("c_outside");

    circuit c({&qa, &qb}, "testCircuit");

    bool caught = false;

    // ================================================================
    // 1. Valid qubit should NOT throw
    // ================================================================
    try {
        c.H(qa);  // should work fine
    } catch (...) {
        std::cerr << "FAIL: H(qa) threw unexpectedly.\n";
        return 1;
    }

    try {
        c.X(qb);  // should also work
    } catch (...) {
        std::cerr << "FAIL: X(qb) threw unexpectedly.\n";
        return 1;
    }

    // ================================================================
    // 2. Invalid qubit should THROW for single-qubit gates
    // ================================================================
    caught = false;
    try {
        c.Z(qc); // qc not in circuit
    } catch (const std::runtime_error&) {
        caught = true;
    }

    if (!caught) {
        std::cerr << "FAIL: Z(qc) did NOT throw for out-of-circuit qubit.\n";
        return 1;
    }

    // ================================================================
    // 3. Invalid qubits should THROW for CNOT
    // ================================================================

    // Case A: control not in circuit
    caught = false;
    try {
        c.CNOT(qc, qa);
    } catch (const std::runtime_error&) {
        caught = true;
    }
    if (!caught) {
        std::cerr << "FAIL: CNOT(qc, qa) did NOT throw.\n";
        return 1;
    }

    // Case B: target not in circuit
    caught = false;
    try {
        c.CNOT(qa, qc);
    } catch (const std::runtime_error&) {
        caught = true;
    }
    if (!caught) {
        std::cerr << "FAIL: CNOT(qa, qc) did NOT throw.\n";
        return 1;
    }

    // Case C: both invalid
    caught = false;
    try {
        c.CNOT(qc, qc);
    } catch (const std::runtime_error&) {
        caught = true;
    }
    if (!caught) {
        std::cerr << "FAIL: CNOT(qc, qc) did NOT throw.\n";
        return 1;
    }

    // ================================================================
    // PASSED
    // ================================================================
    std::cout << "PASS: qubit_in_circuit validation works correctly.\n";
    return 0;
}
