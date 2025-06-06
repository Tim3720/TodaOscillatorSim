#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>

void rk4Step(double* state, double* helper, double t, double A);
double driver(double t, double A);

// external constants:
const double L = 64e-3;
const double C0 = 100e-12;
const double R = 10000;
const double US = 0.5;

const double a = R * sqrt(C0 / L) * (2 * M_PI);
const double phi0 = 0;

const double f0 = 55.4e3;
const double omega = f0 * 2 * M_PI * sqrt(L * C0) * 2 * M_PI;

const int nPeriods = 2000;
const int stepsPerPeriod = 1001;
const double dt = 1 / (stepsPerPeriod * omega / (2 * M_PI));
const size_t steps = nPeriods * stepsPerPeriod;

#define BIFURKATIONMODE

int main()
{
    double dA = 0.0001;
    double A0 = 3.05;
    double A1 = 3.15;
    size_t ASteps = (A1 - A0) / dA;

#pragma omp parallel for
    for (size_t i = 0; i < ASteps; i++) {
        const double A = A0 + i * dA;
#ifdef BIFURKATIONMODE
        std::ofstream file("BifurkationFull/A_" + std::to_string(A) + ".dat",
            std::ios::trunc);
#else
        std::ofstream file("Simulations/A_" + std::to_string(A) + ".dat",
            std::ios::trunc);
#endif

        std::cout << A << std::endl;

        double state[] = {-0.1, -1};  // abritrary starting state
        double helper[8];
        double t = 0;

        for (size_t step = 0; step < steps; step++) {
            t += dt;
            rk4Step(state, helper, t, A);

#ifdef BIFURKATIONMODE
            if (step > 3 * steps / 4 && step % stepsPerPeriod == 0) {
                // file << std::setprecision(15) << t << "\t" << driver(t, A) << "\t"
                //      << state[0] << "\t" << state[1] << std::endl;
                file << std::setprecision(15) << state[0] << std::endl;
            }
#else
            if (step % 10 == 0) {
                file << std::setprecision(15) << t << "\t" << driver(t, A) << "\t"
                     << state[0] << "\t" << state[1] << std::endl;
            }
#endif
        }
    }
};

double driver(double t, double A)
{
    return A / US * cos(omega * t + phi0) * 4 * M_PI * M_PI;
}

void diff(const double* state, double* dst, double t, double A)
{
    double I = state[0];
    double Q = state[1];

    dst[0] = -a * I - 4 * M_PI * M_PI * (exp(Q) - 1) - driver(t, A);
    dst[1] = I;
}

void rk4Step(double* state, double* helper, double t, double A)
{
    double copy[] = {state[0], state[1]};

    double* k1 = helper;
    diff(state, k1, t, A);
    state[0] = copy[0] + 0.5 * dt * k1[0];
    state[1] = copy[1] + 0.5 * dt * k1[1];

    double* k2 = helper + 2;
    diff(state, k2, t + 0.5 * dt, A);
    state[0] = copy[0] + 0.5 * dt * k2[0];
    state[1] = copy[1] + 0.5 * dt * k2[1];

    double* k3 = helper + 4;
    diff(state, k3, t + 0.5 * dt, A);
    state[0] = copy[0] + dt * k3[0];
    state[1] = copy[1] + dt * k3[1];

    double* k4 = helper + 6;
    diff(state, k4, t + dt, A);
    state[0] = copy[0] + dt / 6.0 * (k1[0] + 2 * k2[0] + 2 * k3[0] + k4[0]);
    state[1] = copy[1] + dt / 6.0 * (k1[1] + 2 * k2[1] + 2 * k3[1] + k4[1]);
}
