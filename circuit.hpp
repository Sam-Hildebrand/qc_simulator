#include <vector>
#include "gates.hpp"

class circuit{
    
    private:
        std::vector<gate> c;

    public:
        circuit(const std::vector<gate>& circuit_gates) : c(circuit_gates) {}

        std::string print_circuit() const {
            std::string result = "";
            for (size_t i = 0; i < c.size(); ++i) {
                result += "-";
                result += c[i].symbol();
            }
            return result + "-";
        }

        void apply_to(qubit& q) const {
            for (const auto& gate : c) {
                gate(q);
            }
        }
};
