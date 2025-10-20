# rk2expint
Rational Krylov Runge--Kutta exponential integrators

This repository implements (RK)$^2$EXPINT, a technique leveraging rational Krylov subspace approximations to the matrix exponential acting on vectors for the efficient and scalable solution of large stiff systems of ordinary differential equations (ODEs) with exponential Runge--Kutta integrators. Furthermore, it contains Matlab codes reproducing the numerical experiments from the

**Preprint:**
[1] K. Bergermann and M. Stoll. [Adaptive rational Krylov methods for exponential Runge–Kutta integrators](https://doi.org/10.1137/23M1559439), SIAM Journal on Matrix Analysis and Applications, 45(1), p.744-770, 2024.

**Requirements:**
All codes have been tested with Matlab version R2020b on Ubuntu 20.04.6 LTS. For everything to work, we depend on several (free) external software packages:

**External packages:**
- AGMG (free academic license required, we used v.3.3.5): http://agmg.eu/
- kiops: https://gitlab.com/stephane.gaudreault/kiops
- phipm: http://www1.maths.leeds.ac.uk/~jitse/software.html
- rktoolbox: http://guettel.com/rktoolbox/
- pardiso (free academic license required, we used v6, for installation, you may find this helpful: https://github.com/blechta/pardiso-matlab-recipes): https://www.pardiso-project.org/

Please obtain these packages as indicated on the respective website and copy the obtained codes into the empty placeholder directories (and get pardiso running. You can avoid the dependence on pardiso by using the linear_system_solver='lu_Matlab' option. Our numerical experiments showed a much better performance of pardiso though).

**Data:**

[1, Subsections 6.3 and 6.4] present numerical experiments on networks for which network data is required. The directory 'network_data' contains some of those files, but in order to keep the size of the repository reasonable, we have excluded the larger network files of the 'loc-Brightkite', 'ny2010', and 'roadNet-PA' networks. If you want to use them, please download them in .mat format using the links below and move them into the 'network_data' directory.
- https://sparse.tamu.edu/SNAP/loc-Brightkite
- https://sparse.tamu.edu/DIMACS10/ny2010
- https://sparse.tamu.edu/SNAP/roadNet-PA

**Additional directories:**
- pole_files: contains .mat files of the poles for the rational Krylov methods discussed in [1, Subsection 4.1]
- subroutines: contains some subroutines listed below
- system_solve_fhs: contains functions for solving the shifted linear systems by the strategies discussed in [1, Subsection 4.2]
- time_integrators: implements three example explicit exponential Runge--Kutta integrators, cf. [1, Section 3], using very similar syntax to Matlab's ode solvers

**Main scripts:**
- main.m: reproduces the numerical experiments from [1, Section 6]. Specifies poles and the time discretization and calls 'run_problem.m'
- compute_rkfit_poles.m: Computes optimal poles for the function $f(x)=e^{-x}$ using the RKFIT algorithm [2]

**Functions:**
- subroutines/compute_factorizations.m: precomputes LU or Cholesky decompositions of the matrices required in the direct solution of block linear systems
- subroutines/max_connected_component.m: detects largest connected component of a graph
- subroutines/rk2expint.m: Rational Krylov Runge--Kutta exponential integrator (RK)^2EXPINT routine as introduced in [1] to adaptively compute linear combinations of $\varphi$-functions by means of approximating the action of the matrix exponential on a certain vector [1,3,4,5]. This code is an adaptation of the KIOPS routine (https://gitlab.com/stephane.gaudreault/kiops) [3], which, in turn, is based on the PHIPM and EXPMVP codes (http://www1.maths.leeds.ac.uk/~jitse/software.html) [4,5].
- subroutines/run_problem.m: efficiently solves a stiff system of ODEs by exponential integrators, cf. [1,3,4,5,6]
- subroutines/setup_problem.m: sets up stiff systems of ODEs considered in [1, Section 6]. We consider the Allen--Cahn (AC) equation as well as the Gierer--Meinhardt (GM) equations on finite difference discretizations of two-dimensional unit squares as well as graph Laplacians of inherently discrete network domains
- system_solve_fhs/solve_agmg_mult.m: solves the block linear system [1, Subsection 4.2] by backsubstitution and a flexible conjugate gradient preconditioned with aggregation-based multigrid (AGMG)
- system_solve_fhs/solve_lu_Matlab.m: solves the block linear system [1, Subsection 4.2] by backsubstitution and an LU or Cholesky decomposition precomputed by Matlab's amd and lu functions
- system_solve_fhs/solve_lu_pardiso.m: solves the block linear system [1, Subsection 4.2] by backsubstitution and an LU or Cholesky decomposition precomputed by pardiso
- time_integrators/ETD3RK.m: ETD3RK exponential integrator of stage 3 and stiff order 3 [7]
- time_integrators/Krogstad4.m: Krogstad exponential integrator of stage 4 and stiff order 4 [8]
- time_integrators/SW2.m: Strehmel & Weiner exponential integrator of stage 2 and stiff order 2 [9]

**License:**
 - LICENSE: GNU General Public License v2.0

**References:**
- [1] K. Bergermann and M. Stoll. Adaptive rational Krylov methods for exponential Runge–Kutta integrators, SIAM Journal on Matrix Analysis and Applications, 45(1), p.744-770, 2024.
- [2] M. Berljafa and S. Güttel, The RKFIT algorithm for nonlinear rational approximation, SIAM J. Sci. Comput., 39 (2017), pp. A2049–A2071.
- [3] Gaudreault, S., Rainwater, G. and Tokman, M., 2018. KIOPS: A fast adaptive Krylov subspace solver for exponential integrators. Journal of Computational Physics.
- [4] Niesen, J. and Wright, W.M., 2011. A Krylov subspace method for option pricing. SSRN 1799124
- [5] Niesen, J. and Wright, W.M., 2012. Algorithm 919: A Krylov subspace algorithm for evaluating the \varphi-functions appearing in exponential integrators. ACM Transactions on Mathematical Software (TOMS), 38(3), p.22
- [6] A. H. Al-Mohy and N. J. Higham, Computing the action of the matrix exponential, with an application to exponential integrators, SIAM J. Sci. Comput., 33 (2011), pp. 488–511.
- [7] S. M. Cox and P. C. Matthews, Exponential time differencing for stiff systems, J. Comput. Phys., 176 (2002), pp. 430–455.
- [8] S. Krogstad, Generalized integrating factor methods for stiff PDEs, J. Comput. Phys., 203 (2005), pp. 72–88.
- [9] R. Weiner, Linear-implizite Runge-Kutta-Methoden und ihre Anwendung, vol. 127, Springer-Verlag, 2013.

**Contact:**

Kai Bergermann ([kai.bergermann@math.tu-chemnitz.de](mailto:kai.bergermann@math.tu-chemnitz.de))

