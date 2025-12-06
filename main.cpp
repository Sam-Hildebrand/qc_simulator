#include "qc.hpp"
#include <iostream>

int main() {
    int num_runs = 1000;

    qubit q_alice("alice");
    qubit q_bob("bob");

    circuit c({&q_alice, &q_bob});

    // Entangling Alice and Bob's qubits
    c.H(q_alice);
    c.CNOT(q_alice, q_bob);

    std::cout << "Circuit after applying gates:\n";
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