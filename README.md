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
* Initial Condition: The solution is initialized with a uniform value across the domain, $u(x,0)=1.0$.

**3. Methodology**
<br/> The solution to the nonlinear diffusion equation is computed using the Finite Element Method (FEM). This method involves discretizing the domain into small elements and approximating the solution over these elements using basis functions.

*3.1. Weak Formulation*
<br/> To derive the weak form, we multiply the original equation by a test function $v(x)$ and integrate by parts:

$$\int_{\Omega}{\upsilon(x)\frac{\partial u}{\partial t}}dx = \int_{\Omega}{\upsilon(x)\frac{\partial}{\partial x}(D(u)\frac{\partial u}{\partial x})}dx  $$

Using Galerkin FEM, the solution $u(x,t)$ is approximated by a linear combination of basis functions $N_i(x)$:

$$ u\left(x\right) \approx \sum_{i=1}^{n}{u_i N_i(x)} $$

*3.2. Explanation of Basis Functions*
<br/> The domain is discretized into $N$ nodes, and between each pair of adjacent nodes, we define a linear basis function. These basis functions $N_i(x)$ have the following properties:
* Local support: Each basis function is associated with a specific node, and its value is non-zero only on the two elements that share this node.
* Piecewise linear: The function is linear between the nodes, meaning it varies linearly within the interval between two adjacent nodes.

<br/> For node $i$, the basis function $N_i(x)$ is defined as:

This leads to a system of equations of the form:

$$M^{T} \frac{du}{dt}=-K^{T}u$$

Where:
* $M$ is the mass matrix.
* $𝐾$ is the stiffness matrix that depends on $u$ and the nonlinear diffusion coefficient $D(u)$.

<br/> The linear basis function is as follows

$$ N_i(x)=\begin{cases}
    \frac{x-x_{i-1}}{x_i-x_{i-1}} & \text{$x_{i-1}~ \leq  x \leq x_{i}$}\\
    \frac{x_{i+1}-x}{x_{i+1}-x_{i}} & \text{$x_{i}~ \leq  x \leq x_{i+1}$}\\
    0, & \text{otherwise}
       \end{cases} $$

Where:

* $N_i(x)$ takes the value of 1 at node $i$ (at $x=x_i$).
* $N_i(x)$ is 0 at nodes $i−1$ and $i+1$.
* $N_i(x)$ is linear between the nodes, i.e., it rises from 0 to 1 between $x_{i-1}$ and $x_i$, and then falls from 1 to 0 between $x_i$ and $x_{i+1}$.
To form the mass matrix, we integrate over the domain using the basis functions:

𝑀
𝑖
𝑗
=
∫
Ω
𝑁
𝑖
(
𝑥
)
𝑁
𝑗
(
𝑥
)
 
𝑑
𝑥
M 
ij
​
 =∫ 
Ω
​
 N 
i
​
 (x)N 
j
​
 (x)dx
Where:

𝑀
𝑖
𝑗
M 
ij
​
  represents the entry in the mass matrix corresponding to basis functions 
𝑁
𝑖
(
𝑥
)
N 
i
​
 (x) and 
𝑁
𝑗
(
𝑥
)
N 
j
​
 (x).
𝑁
𝑖
(
𝑥
)
N 
i
​
 (x) and 
𝑁
𝑗
(
𝑥
)
N 
j
​
 (x) are the linear basis functions (hat functions) defined over the elements.
Formula for 1D Linear Basis Functions:
For linear basis functions in 1D, the mass matrix entries are calculated using Gaussian quadrature (or analytically) over each element. For an element with length 
ℎ
h (which in 1D corresponds to 
ℎ
=
Δ
𝑥
h=Δx), the mass matrix contributions are:

𝑀
𝑖
𝑖
=
2
ℎ
6
,
𝑀
𝑖
𝑗
=
ℎ
6
for 
𝑖
≠
𝑗
M 
ii
​
 = 
6
2h
​
 ,M 
ij
​
 = 
6
h
​
 for i

=j
Where 
ℎ
h is the element length (or 
Δ
𝑥
Δx).

Mass Matrix Implementation:
In the code, the mass matrix is assembled using these contributions. Here's the formula used:

Diagonal terms 
𝑀
𝑖
𝑖
=
2
3
⋅
ℎ
M 
ii
​
 = 
3
2
​
 ⋅h.
Off-diagonal terms 
𝑀
𝑖
,
𝑖
−
1
=
𝑀
𝑖
−
1
,
𝑖
=
1
6
⋅
ℎ
M 
i,i−1
​
 =M 
i−1,i
​
 = 
6
1
​
 ⋅h.
*3.3. Properties of the Linear Basis Functions*
<br/> 3.3.1. Partition of Unity: The sum of all basis functions at any point $x$ is equal to 1, i.e., $\sum_{i}N_i(x)=1$.
<br/> 3.3.2. Locality: Each basis function is non-zero only in the elements adjacent to the node $i$. This results in a sparse system of equations, which is computationally efficient.
<br/> 3.3.2. Continuity: The piecewise linear basis functions are continuous, meaning they ensure that the solution $u(x)$ is continuous across the elements, but the derivative may not be.

*3.4. Time Discretization*
<br/> For time integration, we use the forward Euler method:

$$\frac{u^{n+1}-u^n}{\Delta t}=-(M^T)^{-1}K^Tu^n$$

This results in the update formula:

$$u^{n+1}=u^n-\Delta t(M^T)^{-1}K^Tu^n$$
 
The solution is iteratively updated over the specified time steps.

*3.5. Spatial Discretization*
<br/> The domain is discretized into $N$ elements of equal size $Δx$, and linear basis functions are used to approximate the solution over each element. The mass matrix $M$ and stiffness matrix $K$ are assembled based on the basis functions and the nonlinear diffusion coefficient.
* Mass Matrix: The mass matrix is constant and set initially. It represents how quantities are distributed over the elements.
* Stiffness Matrix: The stiffness matrix is assembled at each time step based on the current value of the solution $u$, making it dependent on the nonlinear diffusion coefficient $D(u)$.

