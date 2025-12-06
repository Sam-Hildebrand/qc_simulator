#include <random>
#include <cmath>

class qubit {
    private:
        struct QubitState {
            double theta;
            double phi;
            unsigned int measured : 1;
            unsigned int value : 1;
        } state;

        // Random number generator
        static double random_zero_to_one() {
            static std::random_device rd;
            static std::mt19937 gen(rd());
            static std::uniform_real_distribution<double> dist(0.0, 1.0);
            return dist(gen);
        }

    public:
        qubit() {
            state.theta = 0.0;
            state.phi = 0.0;
            state.measured = 0;
            state.value = 0;
        }

        double theta() const {
            return state.theta;
        }

        double phi() const {
            return state.phi;
        }

        void adjust_angles(double theta, double phi) {
            state.theta = theta;
            state.phi = phi;
        }

        int measure() {
            if (!state.measured) {
                double prob_zero = std::pow(std::cos(state.theta / 2.0), 2.0);
                state.value = (random_zero_to_one() < prob_zero) ? 0 : 1;
                state.measured = 1;
            }
            return state.value;
        }

        // Reset qubit to pre-measurement state
        void reset() {
            state.theta = 0.0;
            state.phi = 0.0;
            state.measured = 0;
            state.value = 0;
        }
};
