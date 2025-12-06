#include "qc.hpp"
#include <iostream>

int main() {
    int num_runs = 1000;

    qubit q0, q1, q2, q3;

    circuit c_bell({&q0, &q1, &q2, &q3}, "Bell");

    // Entangling Alice and Bob's qubits
    c_bell.H(q0);
    c_bell.CNOT(q0, q1);

    std::cout << "Bell circuit:\n";
    std::cout << c_bell.print_circuit(); 

    qubit q_alice("alice");
    qubit q_bob("bob");
    qubit q4;
    qubit q5;
    qubit q6;
    qubit q7;

    circuit c({&q_alice, &q_bob, &q4, &q5, &q6, &q7}, "Quantum Teleportation");

    c.CNOT(q6, q7);

    c.apply_circuit(c_bell, {&q_alice, &q_bob, &q4, &q5});

    std::cout << "\nCircuit after applying gates:\n";
    std::cout << c.print_circuit(); 

    std::cout << "\nAmplitudes:\n";
    std::cout << c.print_amplitudes(2) << "\n";

    quantum_computer qc(c);

    std::cout << "\nProbabilistic Values:\n";
    std::cout << qc.print_prob_dist(qc.calculate_probabilities(), 2) << "\n";

    std::cout << "\nSimulation Results (" << num_runs << " runs):\n";
    std::cout << qc.print_prob_dist(qc.simulate(num_runs), 2) << "\n";

    return 0;
}