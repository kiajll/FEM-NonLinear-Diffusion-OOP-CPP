#include <iostream>
#include <vector>
#include <cmath>
#include <Eigen/Dense>

using namespace std;
using namespace Eigen;

// Base class for a general FEM Solver
class FEMSolver {
protected:
    int nx;        // Number of nodes
    double L;      // Length of domain
    double dx;     // Spatial step size
    double dt;     // Time step size
    int nt;        // Number of time steps
    VectorXd u;    // Solution vector
    MatrixXd M;    // Mass matrix

public:
    FEMSolver(int nx_, double L_, double dt_, int nt_)
        : nx(nx_), L(L_), dt(dt_), nt(nt_) {
        dx = L / (nx - 1);
        u = VectorXd::Ones(nx);  // Initialize solution vector to ones
        M = assemble_mass_matrix();  // Assemble mass matrix
    }

    // Virtual functions for derived classes to implement specific behavior
    virtual MatrixXd assemble_stiffness_matrix() = 0;
    virtual void apply_boundary_conditions(VectorXd& u_new) = 0;

    // Encapsulated mass matrix assembly (same for all solvers)
    MatrixXd assemble_mass_matrix() {
        MatrixXd M = MatrixXd::Zero(nx, nx);
        for (int i = 1; i < nx - 1; ++i) {
            M(i, i) = 2.0 / 3.0 * dx;
            M(i, i - 1) = 1.0 / 6.0 * dx;
            M(i - 1, i) = 1.0 / 6.0 * dx;
        }
        M(0, 0) = 1.0;  // Dirichlet condition at x = 0
        M(nx - 1, nx - 1) = 1.0;  // Dirichlet condition at x = L
        return M;
    }

    // Function to run the simulation (solving the system over time)
    void solve() {
        for (int n = 0; n < nt; ++n) {
            MatrixXd K = assemble_stiffness_matrix();
            VectorXd rhs = M * u;
            VectorXd u_new = M.colPivHouseholderQr().solve(rhs - dt * K * u);
            apply_boundary_conditions(u_new);
            u = u_new;
        }
    }

    // Display the solution
    void display_solution() {
        for (int i = 0; i < nx; ++i) {
            cout << "x[" << i << "] = " << i * dx << ", u[" << i << "] = " << u[i] << endl;
        }
    }
};

// Derived class for nonlinear diffusion (inherits from FEMSolver)
class NonlinearDiffusionSolver : public FEMSolver {
public:
    NonlinearDiffusionSolver(int nx_, double L_, double dt_, int nt_)
        : FEMSolver(nx_, L_, dt_, nt_) {}

    // Nonlinear diffusion coefficient as a function of u
    double D(double u) {
        return 1.0 + 0.5 * u;  // Example: linear dependence on u
    }

    // Override the stiffness matrix assembly specific to nonlinear diffusion
    MatrixXd assemble_stiffness_matrix() override {
        MatrixXd K = MatrixXd::Zero(nx, nx);
        for (int i = 1; i < nx - 1; ++i) {
            double diffusion_coeff = D(u[i]);  // Nonlinear diffusion coefficient
            K(i, i) = 2 * diffusion_coeff / dx;
            K(i, i - 1) = -diffusion_coeff / dx;
            K(i - 1, i) = -diffusion_coeff / dx;
        }
        K(0, 0) = 1.0;  // Dirichlet condition at x = 0
        K(nx - 1, nx - 1) = 1.0;  // Dirichlet condition at x = L
        return K;
    }

    // Override boundary condition handling (Dirichlet BCs in this case)
    void apply_boundary_conditions(VectorXd& u_new) override {
        u_new(0) = 1.0;     // Dirichlet condition at x = 0
        u_new(nx - 1) = 1.0;  // Dirichlet condition at x = L
    }
};

int main() {
    // Create an instance of NonlinearDiffusionSolver with nx=20, L=2, dt=0.001, nt=100
    NonlinearDiffusionSolver solver(20, 2.0, 0.001, 100);
    
    // Solve the nonlinear diffusion equation
    solver.solve();
    
    // Display the results
    solver.display_solution();
    
    return 0;
}
