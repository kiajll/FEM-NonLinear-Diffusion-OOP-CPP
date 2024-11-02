/************************************************************************
* Governing equation: 1D Nonlinear Diffusion 
* Spatial Numerical Discretization: Finite Element Method
* Temporal integrating Technique: Finite Difference Forward Euler Method
* Structure paradigm: Object Oriented Programming
* Author: Kiarash Jalali
* Licensed under the Creative Commons Attribution-NoDerivatives 4.0, 
* International License.
************************************************************************/
#include <iostream>
#include <vector>
#include <cmath>
#include <Eigen/Dense>

using namespace std;
using namespace Eigen;

// AbstractFemSolver: Pure abstract base class for an FEM Solver
class AbstractFemSolver {
public:
    // Pure virtual methods for abstraction
    virtual MatrixXd assemble_mass_matrix() = 0;
    virtual MatrixXd assemble_stiffness_matrix() = 0;
    virtual void apply_boundary_conditions(VectorXd& u_new) = 0;
    virtual void solve() = 0;
    virtual void display_solution() = 0;

    // Virtual destructor
    virtual ~AbstractFemSolver() = default;
};

// FEMSolver class inheriting from AbstractFemSolver, implementing general FEM functionality
class FEMSolver: public AbstractFemSolver {
protected:
    int nx;        // Number of nodes
    double L;      // Length of domain
    double dx;     // Spatial step size
    double dt;     // Time step size
    int nt;        // Number of time steps
    VectorXd u;    // Solution vector
    MatrixXd M;    // Mass matrix

public:
    // Constructor initializing FEMSolver with parameters and setting up the initial mass matrix
    FEMSolver(int nx_, double L_, double dt_, int nt_)
        : nx(nx_), L(L_), dt(dt_), nt(nt_) {
        dx = L / (nx - 1);
        u = VectorXd::Ones(nx);  // Initialize solution vector to ones
        M = assemble_mass_matrix();  // Assemble mass matrix
    }

    // Encapsulated mass matrix assembly (same for all solvers)
    MatrixXd assemble_mass_matrix() override {
        MatrixXd M = MatrixXd::Zero(nx, nx);
        for (int i = 1; i < nx - 1; ++i) {
            M(i, i) = 2.0 / 3.0 * dx;
            M(i, i - 1) = 1.0 / 6.0 * dx;
            M(i - 1, i) = 1.0 / 6.0 * dx;
        }
        M(0, 0) = 1.0;  // Dirichlet boundary condition at x = 0
        M(nx - 1, nx - 1) = 1.0;  // Dirichlet boundary condition at x = L
        return M;
    }

    // Virtual methods for derived classes
    virtual MatrixXd assemble_stiffness_matrix() = 0;
    virtual void apply_boundary_conditions(VectorXd& u_new) = 0;

    // Function to run the simulation (solving the system over time)
    void solve() override {
        for (int n = 0; n < nt; ++n) {
            MatrixXd K = assemble_stiffness_matrix();
            VectorXd rhs = M * u;
            VectorXd u_new = M.colPivHouseholderQr().solve(rhs - dt * K * u);
            apply_boundary_conditions(u_new);
            u = u_new;
        }
    }

    // Display the solution
    void display_solution() override {
        for (int i = 0; i < nx; ++i) {
            cout << "x[" << i << "] = " << i * dx << ", u[" << i << "] = " << u[i] << endl;
        }
    }
};

// NonlinearDiffusionSolver class inherits from FEMSolver, handling nonlinear-specific logic
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
            K(i, i) = 2 * diffusion_coeff / dx;  // Diagonal elements
            K(i, i - 1) = -diffusion_coeff / dx;  // Off-diagonal elements
            K(i - 1, i) = -diffusion_coeff / dx;  // Off-diagonal elements
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
    AbstractFemSolver* solver = new NonlinearDiffusionSolver(20, 2.0, 0.001, 100);
    
    // Solve the nonlinear diffusion equation
    solver->solve();
    
    // Display the results
    solver->display_solution();
    
    delete solver;
    return 0;
}
