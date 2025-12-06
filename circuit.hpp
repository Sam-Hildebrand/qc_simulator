#include <cmath>
#include <stdexcept>
#include <vector>
#include <string>

enum class gate {
    None,
    I,        // identity / placeholder
    H,        // Hadamard
    X,        // Pauli-X  (bit flip)
    Z,        // Pauli-Z  (phase flip)
    Control,  // control dot for CNOT
    Target    // target + for CNOT
};



class qubit {
    private:
        std::vector<gate> gates;
        
    public:
        int measured : 1;
        int value    : 1;

        // Index of the qubit in the circuit
        size_t index = 0;

        // Each qubit tracks how many gates have been applied to it
        size_t gateCount = 0;

        // Qubit name
        std::string name;

        // Constructors
        qubit(const std::string& n = "") : measured(0), value(0), name(n) {}

        void addGate(gate g) {
            gates.push_back(g);
            gateCount++;
        }

        std::string printGates() const {
            std::string seq;
            for (const auto& g : gates) {
                switch (g) {
                    case gate::H: seq += "-H"; break;
                    case gate::X: seq += "-X"; break;
                    case gate::Z: seq += "-Z"; break;
                    case gate::I: seq += "--"; break;
                    case gate::Control: seq += "-•"; break;
                    case gate::Target: seq += "-+"; break;
                    default: break;
                }
            }
            return seq + "-";
        }
};

class circuit {

    private:
        std::vector<qubit*> qubits;

    public:
        std::vector<double> amplitudes;

        circuit(const std::vector<qubit*>& q) : qubits(q){
            amplitudes.resize(1 << qubits.size());
            amplitudes[0] = 1.0;

            for(size_t i = 0; i < qubits.size(); ++i) {
                qubits[i]->index = i;
            }
        }

        size_t size() const {
            return qubits.size();
        }

        qubit& operator[](size_t index) {
            return *qubits[index];   // dereference pointer
        }

        std::string print_circuit() const {
            std::string result;
            // Step 1: determine the maximum length of qubit name + index for alignment
            size_t max_name_len = 0;
            for (const auto& q : qubits) {
                const std::string& name = q->name.empty() ? "q" + std::to_string(q->index) : q->name;
                size_t len = std::to_string(q->index).size() + 2 + name.size();
                if (len > max_name_len) max_name_len = len;
            }
            // Step 2: determine the maximum gate count
            size_t max_gate_count = 0;
            for (const auto& q : qubits) {
                if (q->gateCount > max_gate_count) max_gate_count = q->gateCount;
            }
            // Step 3: print each qubit line, pad name and gates properly
            for (size_t i = 0; i < qubits.size(); ++i) {
                const std::string& name = qubits[i]->name.empty() ? "q" + std::to_string(qubits[i]->index) : qubits[i]->name;
                std::string line = std::to_string(qubits[i]->index) + ": " + name;
                // pad the name so first '-' lines up, plus 1 extra space
                line += std::string(max_name_len - line.size() + 1, ' ');
                std::string gates = qubits[i]->printGates();
                // pad the gates on the right based on gate count difference
                size_t padding_needed = (max_gate_count - qubits[i]->gateCount) * 2; // 2 chars per gate (-X, -H, etc)
                gates += std::string(padding_needed, '-');
                line += gates;
                result += line + "\n";
            }
            return result;
        }

        std::string print_amplitudes(int precision = 2) const {
            std::string result;
            size_t n = qubits.size();
            bool first = true;

            for (size_t i = 0; i < amplitudes.size(); ++i) {
                double amp = amplitudes[i];
                if (std::abs(amp) < 1e-12) continue;

                // Build basis string: leftmost = qubits[0], rightmost = qubits[N-1]
                std::string ket;
                ket.reserve(n);
                for (size_t q = 0; q < n; ++q) {  // <-- LEFTMOST qubit first
                    size_t bit = 1ULL << (n - q - 1); // MSB = qubits[0]
                    ket.push_back((i & bit) ? '1' : '0');
                }

                // Round and remove trailing zeros
                double factor = std::pow(10, precision);
                double rounded = std::round(amp * factor) / factor;
                std::string amp_str = std::to_string(rounded);
                amp_str.erase(amp_str.find_last_not_of('0') + 1, std::string::npos);
                if (amp_str.back() == '.') amp_str.pop_back();

                if (!first) result += " + ";
                result += "(" + amp_str + ")|" + ket + ">";
                first = false;
            }

            return result.empty() ? "0" : result;
        }

        void I(qubit& q) {
            q.addGate(gate::I);
        }

        void H(qubit& q) {
            size_t n = qubits.size();
            size_t bit = 1ULL << (n - q.index - 1); // MSB-first mask
            std::vector<double> new_amplitudes(amplitudes.size());

            const double s = M_SQRT1_2; // 1/sqrt(2)

            for (size_t i = 0; i < amplitudes.size(); ++i) {
                if ((i & bit) == 0) {
                    size_t partner = i | bit; // set the qubit bit to 1
                    new_amplitudes[i]       = s * (amplitudes[i] + amplitudes[partner]);
                    new_amplitudes[partner] = s * (amplitudes[i] - amplitudes[partner]);
                }
            }

            amplitudes = std::move(new_amplitudes);

            q.addGate(gate::H);
        }


        void X(qubit& q) {
            size_t n = qubits.size();
            size_t bit = 1ULL << (n - q.index - 1); // MSB-first mask
            std::vector<double> new_amplitudes = amplitudes;

            for (size_t i = 0; i < amplitudes.size(); ++i) {
                size_t flipped = i ^ bit;
                new_amplitudes[flipped] = amplitudes[i]; // swap amplitude
            }

            amplitudes = std::move(new_amplitudes);

            q.addGate(gate::X);
        }

        void Z(qubit& q) {
            size_t n = qubits.size();
            size_t bit = 1ULL << (n - q.index - 1); // MSB-first mask

            for (size_t i = 0; i < amplitudes.size(); ++i) {
                if (i & bit) // qubit = 1
                    amplitudes[i] = -amplitudes[i];
            }

            q.addGate(gate::Z);
        }

        void CNOT(qubit& control, qubit& target) {
            size_t n = qubits.size();
            size_t control_bit = 1ULL << (n - control.index - 1);
            size_t target_bit  = 1ULL << (n - target.index - 1);

            std::vector<double> new_amplitudes = amplitudes;

            // Flip target qubit amplitudes when control is 1
            for (size_t i = 0; i < amplitudes.size(); ++i) {
                if (i & control_bit) {
                    size_t flipped = i ^ target_bit;
                    new_amplitudes[flipped] = amplitudes[i];
                } else {
                    new_amplitudes[i] = amplitudes[i];
                }
            }

            amplitudes = std::move(new_amplitudes);

            // Pad all qubits with I so they line up before adding CNOT symbols
            size_t max_len = 0;
            for (auto* q : qubits) max_len = std::max(max_len, q->gateCount + 1); // +1 for new gate

            for (auto* q : qubits) {
                while (q->gateCount < max_len - 1) // pad until right before the new column
                    I(*q);  // this increments gateCount
            }

            // Add CNOT gate symbols (automatically increments gateCount)
            control.addGate(gate::Control);
            target.addGate(gate::Target);
        }

};
