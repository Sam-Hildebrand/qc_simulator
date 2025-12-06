#include "computer.cpp"
#include <iostream>

int main() {
    int num_runs = 10000;

    qubit q1;
    qubit q2;

    circuit c1({H});
    circuit c2({X, H});

    c1.apply_to(q1);
    std::cout  << "Qubit state after applying " << c1.print_circuit() << ": " << print_amplitudes(q1) << std::endl;
    q1.reset();
    std::cout << "Probability of 1 after "<< num_runs << " runs: " << simulate_probability(10000, c1, q1) * 100.0 << "%" << std::endl;

    c2.apply_to(q2);
    std::cout  << "\nQubit state after applying " << c2.print_circuit() << ": " << print_amplitudes(q2) << std::endl;
    q2.reset();
    std::cout << "Probability of 1 after "<< num_runs << " runs: " << simulate_probability(10000, c2, q2) * 100.0 << "%" << std::endl;

    return 0;
}