#include "qc.hpp"
#include <iostream>

int main() {
    int num_runs = 1000;

    qubit qA, qB, qC, qD;

    circuit c_4bell({&qA, &qB, &qC, &qD}, "4 Qubit Bell");

    c_4bell.H(qA);
    c_4bell.CNOT(qA, qB);

    c_4bell.H(qC);
    c_4bell.CNOT(qC, qD);

    c_4bell.CNOT(qA, qC);

    std::cout << "\n4 Qubit Bell circuit:\n";
    std::cout << c_4bell.print_circuit(); 

    qubit q0;
    qubit q1;
    qubit q2;
    qubit q3;

    circuit c({&q0, &q1, &q2, &q3}, "4-Qubit Bell test Circuit");


    c.apply_circuit(c_4bell, {&q0, &q1, &q2, &q3});

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