# 1D Nonlinear Diffusion Solver Using Finite Element Method (FEM)
**1. Introduction**
<br/> This project implements a 1D nonlinear diffusion equation solver using the Finite Element Method (FEM). Diffusion phenomena are commonly observed in fluid flow, heat transfer, and other physical systems where quantities such as temperature, concentration, or velocity diffuse over time due to molecular interactions.
In this particular implementation, we model diffusion in the presence of a nonlinear diffusion coefficient $D(u)$, which depends on the solution field $u(x,t)$. The numerical solution is computed using the FEM approach, where the governing equation is discretized in space and evolved in time using a forward Euler scheme.

**2. Mathematical Definition**
<br/> The general form of the 1D nonlinear diffusion equation is given by:

$$\frac{\partial u}{\partial t}= \frac {\partial}{\partial x}( D(u) \frac{\partial u}{\partial x} )$$

Where:
* $u(x,t)$ is the scalar field representing the physical quantity being diffused (e.g., temperature or concentration).
* $D(u)$ is the nonlinear diffusion coefficient that depends on the solution $u$.
* The term $\frac{\partial u}{\partial t}$ represents the rate of change of the field over time.
* The term $\frac {\partial}{\partial x}( D(u) \frac{\partial u}{\partial x} )$ describes the spatial diffusion, with the diffusion coefficient varying based on the solution.

*2.1. Boundary and Initial Conditions*
* Dirichlet Boundary Conditions: Fixed values for $u$ are applied at both boundaries of the 1D domain, i.e., $u(0,t)=u_0$ and $u(L,t)=u_L$.
* Initial Condition: The solution is initialized with a uniform value across the domain, $u(x,0)=1.0$, or other specified initial distributions.

**3. Methodology**
<br/> The solution to the nonlinear diffusion equation is computed using the Finite Element Method (FEM). This method involves discretizing the domain into small elements and approximating the solution over these elements using basis functions.

*3.1. Weak Formulation*
<br/> To derive the weak form, we multiply the original equation by a test function $v(x)$ and integrate by parts:

$$\int_{\Omega}{\left(\frac{\partial u}{\partial t} N_j+c \frac{\partial u}{\partial x}N_j\right)dx}=0$$
