#include "qubit.hpp"
#include <string>
#include <functional>

class gate {
public:
    using GateFunction = std::function<void(qubit&)>;

private:
    std::string symbol_;   // e.g. "H"
    std::string name_;     // e.g. "Hadamard"
    GateFunction func_;    // function(H)

public:
    gate(const std::string& symbol,
         const std::string& name,
         GateFunction func)
        : symbol_(symbol), name_(name), func_(func) {}

    void operator()(qubit& q) const {
        func_(q);
    }

    const std::string& symbol() const { return symbol_; }
    const std::string& name() const { return name_; }
};

gate H("H", "Hadamard", [](qubit& q) {
    // Compute new theta/phi based on current theta/phi
    double theta_old = q.theta();
    double phi_old = q.phi();

    // Conversion to Cartesian coordinates on Bloch sphere
    double x = std::sin(theta_old) * std::cos(phi_old);
    double y = std::sin(theta_old) * std::sin(phi_old);
    double z = std::cos(theta_old);

    // Apply Hadamard rotation on Bloch sphere
    double x_new = z;
    double y_new = -y;
    double z_new = x;

    // Convert back to theta/phi
    q.adjust_angles(std::acos(z_new), std::atan2(y_new, x_new));
});

gate X("X", "Pauli-X", [](qubit& q) {
    double theta_old = q.theta();
    double phi_old = q.phi();

    // Convert to Cartesian coordinates
    double x = std::sin(theta_old) * std::cos(phi_old);
    double y = std::sin(theta_old) * std::sin(phi_old);
    double z = std::cos(theta_old);

    // Pauli-X rotation on Bloch sphere
    double x_new = x;
    double y_new = -y;
    double z_new = -z;

    // Convert back to theta/phi
    q.adjust_angles(std::acos(z_new), std::atan2(y_new, x_new));
});