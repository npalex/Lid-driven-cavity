# **Lid-driven cavity flow of an incompressible, Newtonian fluid**

&emsp; This program solves the mass conservation and incompressible Navier Stokes equations in 2D, given by

$$ \frac{\partial u}{\partial x} + \frac{\partial v}{\partial y} = 0, $$

$$ Re\left( \frac{\partial u}{\partial t} + \frac{\partial u^2}{\partial x} + \frac{\partial (u v)}{\partial y} + \frac{\partial p}{\partial x}\right) = 
    \left(\frac{\partial^2 u}{\partial x^2} + \frac{\partial^2 u}{\partial y^2} \right), $$

and

$$ Re\left( \frac{\partial v}{\partial t} + \frac{\partial (u v)}{\partial x} + \frac{\partial v^2}{\partial y} + \frac{\partial p}{\partial y}\right) = 
    \left(\frac{\partial^2 v}{\partial x^2} + \frac{\partial^2 v}{\partial y^2} \right). $$

&emsp; Here, the Reynolds number is defined according to 

$$Re = \frac{\rho U L}{\eta},$$

where $\rho$, $\eta$, $U$, and $L$ are the fluid density, fluid viscosity, lid speed, and lid length, respectively, and $u$ and $v$ are the non-dimensional components of the fluid velocity in the $x$ and $y$ directions. The time scale $\bar t$ for this problem is defined by the characteristic shear rate according to $\bar t = \frac{L}{U}$ and the inertial pressure scale was chosen, equal to $\bar p = \rho U^2$. The boundary conditions are no-slip and no-flow at the cavity walls and the fluid is initially at rest. 

&emsp; Following Lee and Leveque (2003), a *fractional step* approach is used to solve the system of equations above, which decomposes the problem into the following steps:

**Step 1. Solve the convection equation:**

&emsp;A Roe solver $\textemdash$ i.e., a locally linear, approximate Riemann solver based on the Godunov method $\textemdash$ is employed to evaulate $q = (u,v)$ satisfying:

$$ q_t + \hat A \cdot q_x + \hat B \cdot q_y = 0 . $$

Dimensional splitting via the donor cell upwind method (DCU) is used to advanced the cell-centered velocities $Q=(U,V)$ forward in time via sweeps in the x-direction

$$ Q_{i,j}^{*}= Q_{i,j}^n - \Delta t (\hat A \cdot Q_x^{n})_{i,j}$$

followed by sweeps in the y-direction

$$ Q_{i,j}^{**}= Q_{i,j}^{*} - \Delta t ( \hat B \cdot Q^{*}_{i,j} )_{i,j}$$

$$ Q_{i,j}^{**}= Q_{i,j}^{*} - \Delta t (\hat B \cdot Q_y^{*})_{i,j}$$

where the Roe matrices $\hat A$ and $\hat B$ are apporoximate Jacobian matrices given by

$$ \hat A =        \begin{bmatrix} 
                                2\hat u & 0 \\
                                \hat v & \hat u \end{bmatrix}, $$

and

$$ \hat B =        \begin{bmatrix} 
                                \hat v & \hat u \\
                                0 & 2 \hat v \end{bmatrix}. $$

Here, $\hat u$ and $\hat v$ are Roe averages (in this case, they are linear interpolations of cell-centered velocities) defined at the edge of each grid cell. For example, Roe averages for the x-sweeps used to evaluate the matrix $A$ are given by

$$ \hat u_{i-\frac{1}{2},j} = \frac{U_{i,j} + U_{i-1,j}}{2}$$

and

$$ \hat v_{i-\frac{1}{2},j} = \frac{V_{i,j} + V_{i-1,j}}{2}$$

In addition, monotenzied central flux limiters are used to achieve second order accuracy for this step where the solution is smooth. Note, a transverse Riemann solver has also been defined for this problem so that the corner-transport upwind method can be used instead of the DCU method if desired.

**Step 2. Solve the diffusion equation:** 

&emsp; An alternating direction implicit (ADI) method is employed to update the numerical solution for diffusion via

$$ q_t = \frac{1}{Re}(q_{xx} + q_{yy}). $$

Two difference equations are used to advance successive time steps of $\frac{\Delta t}{2}$. The first equation is implict in the x-direction

$$ Q_{i,j}^{***} = Q_{i,j}^{**} + \frac{\alpha}{2}\left(\delta^2_x Q^{***} + \delta^2_y Q^{**}\right) $$

and the second equation is implicit in the y-direction

$$ \widetilde{Q}_{i,j} = Q_{i,j}^{***} + \frac{\alpha}{2}\left(\delta^2_x Q^{***} + \delta^2_y \widetilde{Q}\right). $$

The parameter $\alpha$ is defined by

$$ \alpha = \frac{\Delta t}{(\Delta x)^2Re}.$$

and $\delta^2_x$ denotes the central difference of the 2nd partial derivative with respect to $x$.

**Step 3. Update the edge velocities for diffusion:**

&emsp; The velocities at the edges of each grid cell are updated via linear interpolation:

$$ \widetilde{q}_{i-\frac{1}{2},j} = \frac{\widetilde{Q}_{i-1,j} + \widetilde{Q}_{i,j}}{2} $$

**Step 4. Compute the pressure distribution:**

&emsp;The divergence of the equation above provides a Laplacian equation for the pressure,

$$ \nabla^2 p^{n+1} = \frac{1}{\Delta t} \nabla \cdot \widetilde{q}_{i,j} ,$$
which is then discretized with Nuemann boundary conditions to produce a system of linear equations, $Ax = b$. However, the matrix $A$ is singular because the equation set has an inifinite number of solutions within an arbitrary reference pressure. Hence, a ficticious source term $C_0 p^{n+1}$ has been added with proportionality constant $C_0$, which is defined on the order of 1e-9 to render the influence of the source negligble, so that $A$ is non-singular.  

&emsp;So far, the solution $\widetilde{q}$ is not divergence free. In order to satisfy continuity, $\widetilde{q}$ is projected into a divergence-free space by correcting the result for pressure-driven flow via

$$ q^{n+1}_{i,j} = \widetilde{q}^{}_{i,j} -\Delta t\nabla p^{n+1}$$

**Step 5. Update the edge velocities for pressure-driven flow:**

&emsp; The cell-centered velocities are determined by using central differences for the pressure gradient via

$$ u^{n+1}_{i-\frac{1}{2},j} = \widetilde{u}_{i-\frac{1}{2},j} - \Delta t \left( \frac{p^{n+1}_{i,j} - p^{n+1}_{i-1,j}}{\Delta x}\right) $$

and

$$ v^{n+1}_{i,j-\frac{1}{2}} = \widetilde{v}_{i,j-\frac{1}{2}} - \Delta t \left( \frac{p^{n+1}_{i,j} - p^{n+1}_{i,j-1}}{\Delta y}\right) $$

**Step 6. Update the cell-centered velocities for pressure-driven flow:**

&emsp;Finally, the cell-centered velocities are determined by using central differences for the pressure gradient according to

$$ U^{n+1}_{i,j} = \widetilde{U}_{i,j} - \Delta t \left( \frac{p^{n+1}_{i+1,j} - p^{n+1}_{i-1,j}}{2 \Delta x}\right) $$

and

$$ V^{n+1}_{i,j} = \widetilde{V}_{i,j} - \Delta t \left( \frac{p^{n+1}_{i,j+1} - p^{n+1}_{i,j-1})}{2 \Delta y}\right), $$

which produces checkerboard oscillations in pressure. However, these oscillations do not affect the accuracy of the velocity field.

## **Results**:

**References**:

Clawpack Development Team (2023), Clawpack Version 5.9.2,
    http://www.clawpack.org, doi: 10.5281/zenodo.10076317

R. J. LeVeque, 1997. Wave propagation algorithms for multi-dimensional 
    hyperbolic systems. J. Comput. Phys. 131, 327–353.

R. J. LeVeque. Finite Volume Methods for Hyperbolic Problems. Cambridge 
    University Press, Cambridge, UK, 2002.