
### Solid-Phase

We will now provide the complete theoretical background used within the
tidal models. Admittedly, this section will get quite complicated, but I
will try to formulate things in a clear way! We will start with the
original solid tides model, namely "Lovepy" or as we refer to it
"solid1d". But first let us establish some core concepts.

In gravito-elastic tidal theory the response of a medium is represented
by so called $y$-functions. These functions form a state vector at every
radius $\pmb{y}(r)$. These functions also depend on the $(n,m)$-pair;
$\pmb{y_{n,m}}(r)$. The solid tide modules incorporated into Obliqua
make use of the spheroidal 6--vector solution and the spheroidal
8--vector solution, where the latter is used specifically for the
poro-visco elastic models. Generally the first six $y$-functions are
labeled as follows

$$\pmb{y}_{n,m} = (U_{n,m} , V_{n,m} , X_{n,m} , Y_{n,m}, \Phi_{n,m} , \Psi_{n,m} )^T$$

where the first and second components are the radial and tangential
displacements, the third and fourth components the radial and tangential
stresses, the fifth component the potential and the sixth component the
so called 'potential stress'. By extension, the seventh and eighth
$y$-functions are

$$\pmb{y}_{n,m} = (U_{n,m} , V_{n,m} , X_{n,m} , Y_{n,m}, \Phi_{n,m} , \Psi_{n,m}, P_{n,m}, R_{n,m})^T$$

with $P_{n,m}$ the pore pressure, and $R_{n,m}$ the relative tangential
displacement. Notably, these $y$-functions come in pairs, although this
may not be generally true. We can define two sets, both sets are bounded
from below (e.g. core or internal source), but only one set is bounded
from above (e.g. surface or internal sink). The sets are

$$\text{lower} = U_{n,m} , V_{n,m}, \Phi_{n,m}, P_{n,m}$$ 

and

$$\text{upper} = X_{n,m} , Y_{n,m}, \Psi_{n,m}, R_{n,m}.$$ 

It is worth
noting that there is a direct relation between the lower bound
$y$-functions at the surface and the Lovenumbers 

$$\begin{aligned}
    h_n &= y_1(R) \\
    l_n &= y_2(R) \\
    k_n &= y_5(R) - 1
\end{aligned}$$ 

assuming the $y$-functions are dimensionless (note the
index $i$ in $y_i$ corresponds to the $i$th $y$-function in the state
vector $\pmb{y}_{n,m}$ given above). The relations are identical for the
load Lovenumbers. 

$$\begin{aligned}
    h'_n &= y_1(R) \\
    l'_n &= y_2(R) \\
    k'_n &= y_5(R) - 1
\end{aligned}$$ 

However, in the case of a pressure load (atmosphere),
the relations change slightly 

$$\begin{aligned}
    h^P_n &= y_1(R) \\
    l^P_n &= y_2(R) \\
    k^P_n &= y_5(R)
\end{aligned}$$ 

At the same time the upper $y$-functions at the surface
are subject to the cause specific boundary conditions. For the general
case, if we denote $P$ as an external pressure, $\zeta$ as a surface
mass load, $U$ as an external potential, and $\tau$ as a tangential
component of traction, we can write the general boundary conditions for
a given $n$th degree as 

$$\begin{aligned}
\begin{array}{c|cccc}
 & U & \zeta_n & \tau & P \\ \hline
y_{3}(a) & 0 & -g_e \zeta_n & 0 & -P_n \\
y_{4}(a) & 0 & 0 & \tau_n & 0 \\
\dfrac{n+1}{a} y_{5}(a) + y_{6}(a)
           & \dfrac{2n+1}{a} U_n & 4\pi G \zeta_n & 0 & 0
\end{array}
\end{aligned}$$ 

where $g_e$ is the norm of surface gravity. The surface
mass load can also be written as an external potential $U'$ such that
$\zeta_n = [(2n + 1)/4 \pi G a] U'_n$. 

$$\begin{aligned}
y_{3}(R)   &= - \frac{(2n + 1)g_e}{4 \pi G R} U'_n - P_n \\
y_{4}(R)   &= \tau_n \\
\dfrac{n+1}{R} y_{5}(R) + y_{6}(R) 
            &= \frac{2n+1}{R} (U_n + U'_n)
\end{aligned}$$ 

For now we only use the following: the tidal Love number
corresponding to the case of an external potential perturbation $U$ and
the load Love number computed for a surface mass load perturbation
$\zeta$ or equivalently a potential $U'$. They can be directly
determined setting, respectively $(U, U', \tau, P)$ to $(1,0,0,0)$ for
TLN and $(0,1,0,0)$ for LLN in system. Interestingly, the tidal
potential imposed here $U = 1$ is identical to the tidal potential
$U_{n,m,k}$ listed earlier (up to some normalization factor). We simply
set $U = 1$ here to nondimensionalize the Lovenumbers, and scale the
solution later on with $U_{n,m,k}$.

It should be stated that all of these terms ($U, \zeta_n, \tau, P$) are
non-dimensional here (i.e. 1 if present or 0 if not present). They each
give rise to their own set of Lovenumbers ($h_n , l_n , k_n$), which can
be combined into a distortion potential of the planet through

$$U_n^D = k_n^T + (1+k_n^L)U_n^L + K_n^P U_n^P,$$

where $U_n^L$ and
$U_n^P$ are, respectively, the mass-load and pressure potentials:

$$U_n^L = \frac{3}{2n+1}\frac{g \mathcal{S}_n}{\rho_b}, \qquad U_n^P = \frac{3}{2n+1}\frac{\mathcal{Q_n}}{\rho_b}$$

with surface gravity $g$, bulk density $\rho_b$, surface density of the
mass load $\mathcal{S}_n$, and surface pressure $\mathcal{Q_n}$. In
principle, we can account for the stacking of different layers with
there own Lovenumbers by manipulating this expression. For example,
Farhat et al. (2025) provides the forms for $U_n^L$ for a magma ocean
and $\mathcal{Q_n}$ for a uniform lithosphere/crust for the Andrade
rheology. M. Beuthe (2016) provides its explicit form for a crust
characterized by a Maxwell rheology, allowing for radial variations in
the physical properties. I don't know the exact form we should use for
atmospheric pressure, but it can probably be found. Nevertheless, we
should note that each of these terms is effectively just a correction to
the global Lovenumber (and as shown by Farhat et al., 2025 for the
crust, often has only a small effect). As such, we may consider
implementing these corrections in the future, but they should not alter
the results drastically. Given this, we currently combine the
Lovenumbers of the solid and fluid using the simple form

$$k_n = k^T_n + (1 + k^L_n) k_n^{T, \text{(fluid)}}$$ 

to get the global
Lovenumber $k_n$. For this, the determination of tidal and load
lovenumber $k_n$ has been implemented in Obliqua. The other terms are
currently set to zero, but can easily also be determined, allowing us to
further expand on this in the future. Finally, note that the loading
Lovenumbers (both gravity and pressure) are only defined for the solid
part of the mantle.

Evidently, these $y$-functions are functions of radius, hence we wrote
$\pmb{y_{n,m}}(r)$ before. They are effectively the solutions at each
radii to the following first order differential equation

$$\frac{d\pmb{y}_{n,m}(r)}{dr} = \pmb{A}_n (r) \pmb{y}_{n,m}(r) - \pmb{f}_{n,m}(r)$$

where $\pmb{A}_n$ represents the motion matrix with unit $\sim1/r$, it
is a differential operator. The correction term $\pmb{f}_{n,m}(r)$
represents sources and sinks that may be present within the interior
(e.g. earthquakes, actually we will use it for porous interfaces!).

In principle $\pmb{A}_n$ can be anything you want, as long as it
represents the medium that you are considering. There are many different
motion matrices in the literature, but for our purposes one is of
specific interest 

$$\small
\pmb{A}_n(r) = \\
\begin{pmatrix}
-\dfrac{2\lambda}{r\beta} & \dfrac{\ell(\ell+1)\lambda}{r\beta} & \dfrac{1}{\beta} & 0 & 0 & 0 \\[1.2em]

-\dfrac{1}{r} & \dfrac{1}{r} & 0 & \frac{1}{\mu} & 0 & 0 \\[1.2em]

\dfrac{4}{r}\!\left( \dfrac{3\kappa\mu}{r\beta} - \rho_{0}g \right) { - \rho_0 \omega^2} &
\dfrac{\ell(\ell+1)}{r}\!\left(\rho_{0}g - \dfrac{6\kappa\mu}{r\beta}\right) &
-\dfrac{4\mu}{r\beta} & \dfrac{\ell(\ell+1)}{r} & - \dfrac{\rho_{0}(\ell+1)}{r} &
\rho_{0} \\[1.2em]
\dfrac{1}{r}\!\left(\rho_{0}g - \dfrac{6\mu\kappa}{r\beta}\right) &
\dfrac{2\mu}{r^{2}}\!\left[\ell(\ell+1)\!\left(1+\dfrac{\lambda}{\beta}\right)-1\right] { - \rho_0 \omega^2} &
-\dfrac{\lambda}{r\beta} & - \dfrac{3}{r} &
\dfrac{\rho_{0}}{r} & 0 \\[1.2em] 
-4\pi G \rho_{0} & 0 & 0 & 0 & -\dfrac{\ell+1}{r} & 1 \\[1.2em]
-\dfrac{4\pi G \rho_{0} (\ell+1)}{r} & \dfrac{4\pi G \rho_{0} \ell(\ell+1)}{r} & 0 & 0 & 0 & \dfrac{\ell-1}{r} 
\end{pmatrix}$$ 

with ($\ell = n$, same thing, different notation) and

$$\beta = \lambda + 2 \mu$$

and $\lambda$ is the second Lamé parameter
that is expressed in terms of the shear modulus $\mu$ (also known as
first Lamé parameter) and the bulk modulus $\kappa$

$$\lambda = \kappa - \frac{2}{3}\mu.$$

As you can see, this matrix has size $6\times 6$, implying that it
couples six $y$-functions. Specifically these are
$U_{n,m} , V_{n,m} , X_{n,m} , Y_{n,m}, \Phi_{n,m} , \Psi_{n,m}$. We may
extend the motion matrix to also include corrections from porosity, the
corresponding a matrix will have size $8\times 8$ 

$$\tiny
A =
\begin{pmatrix}
\cdot & \cdot & \cdot & \cdot & \cdot & \cdot & \alpha \beta^{-1} & 0 \\
\cdot & \cdot & \cdot & \cdot & \cdot & \cdot & 0 & 0 \\
1i k \rho_\ell^2 g^2 n(n+1)\! \dfrac{r^{-2}}{\omega \eta_\ell} & \cdot & \cdot & \cdot & 
-\! (n+1) r^{-1} 1i k \rho_\ell^2 g n \dfrac{r^{-1}}{\omega \eta_\ell} &
\cdot & 1i k \rho_\ell g n(n+1)\!\dfrac{r^{-2}}{\omega \eta_\ell} - 4\mu \alpha \beta^{-1} r^{-1} & 1i k \rho_\ell^2 g^2 n(n+1)\!\dfrac{r^{-2}}{\omega \eta_\ell} - 4 \phi \rho_\ell g r^{-1} \\
0 & 0 & 0 & 0 & 0 & 0 & 2\alpha\mu r^{-1}\beta^{-1} & \phi \rho_\ell g r^{-1} \\
0 & 0 & 0 & 0 & 0 & 0 & 0 & 4\pi G \rho_\ell \phi \\
-1i 4\pi G n(n+1) r^{-1} k\rho_\ell^2 g\dfrac{r^{-1}}{\omega\eta_\ell} & 0 & 0 & 0 &
1i 4\pi n(n+1)G\rho_\ell^2 k\dfrac{r^{-2}}{\omega\eta_\ell} & 0 &
-1i 4\pi n(n+1)G\rho_\ell k\dfrac{r^{-2}}{\omega\eta_\ell} &
4\pi G (n+1)r^{-1}\!\left(\phi\rho_\ell - 1i n k\rho_\ell^2 g\dfrac{r^{-1}}{\omega\eta_\ell}\right) \\
\rho_\ell g r^{-1}\!\left(4 - 1i k\rho_\ell g n(n+1)\dfrac{r^{-1}}{\omega \phi \eta_\ell}\right) &
-\rho_\ell n(n+1) g r^{-1} &
0 & 0 &
-\rho_\ell (n+1)r^{-1}\!\left(1 - 1i k\rho_\ell g n\dfrac{r^{-1}}{\omega \phi \eta_\ell}\right) &
\rho_\ell &
-1i k\rho_\ell g n(n+1)\!\dfrac{r^{-2}}{\omega \phi \eta_\ell} &
-1i \omega \phi \eta_\ell / k - 4\pi G(\rho - \phi \rho_\ell)\rho_\ell + \rho_\ell g r^{-1}\!\left(4 - 1i k\rho_\ell g n(n+1)\dfrac{r^{-1}}{\omega \phi \eta_\ell}\right) \\
r^{-1}\!\left(1i k\rho_\ell g n(n+1)\dfrac{r^{-1}}{\omega \phi \eta_\ell} - \dfrac{\alpha}{\phi} 4\mu\beta^{-1}\right) &
\dfrac{\alpha}{\phi} 2n(n+1)\mu \beta^{-1} r^{-1} &
-\dfrac{\alpha}{\phi}\beta^{-1} &
0 &
-1i k \rho_\ell n(n+1)\!\dfrac{r^{-2}}{\omega \phi \eta_\ell} &
0 &
1i k n(n+1)\!\dfrac{r^{-2}}{\omega \phi \eta_\ell} - \dfrac{1}{\phi}(S + \alpha^2 \beta^{-1}) &
1i k \rho_\ell g n(n+1)\!\dfrac{r^{-2}}{\omega \phi \eta_\ell} - 2r^{-1}
\end{pmatrix}$$

where 

$$\beta^{-1} = \frac{1}{2\mu + \lambda}, \quad 
S = \frac{\phi}{K_\ell} + \frac{\alpha - \phi}{K}, \quad 
\lambda = K_d - \frac{2\mu}{3}.$$

The terms involving $k$, $\eta_\ell$,
and $\omega$ represent viscous coupling between the solid and liquid
phases, while the terms with $G$ describe self-gravitational feedback.
When $\phi = 0$, the system reduces to the single-phase viscoelastic
formulation.

It is important to note that the ordering of the elements in
$\pmb{A}_n(r)$ and $\pmb{y}_{n,m}(r)$ must be consistent. Different
papers will order the $y$-functions differently (specifically $y_2$ and
$y_3$ are often flipped). We will also perform a reorganization, where
we order $\pmb{y}_{n,m}(r)$ as $[\text{lower}, \text{upper}]^T$.
However, for now we can omit this as there is currently no benefit to
doing this yet.

This puts us in a bit of an awkward spot, since we have a coupled system
of ODEs, but we can only constrain half of the $y$-functions at the
surface. The missing piece is the core solution vector! Interestingly,
we can determine also tides acting within the core if we really want to,
since it would merely require another $\pmb{A}_n(r)$ for that regime. We
can then simply impose a constraint on the lower $y$-functions at the
center of the planet leaving us with six unknowns and six knows, solving
the system. There is a catch, as this particular approach has a
singularity at $r=0$. One can potentially omit this by starting the
solution at some $\delta r > 0$, but this still leaves us with a core
that has no radially resolved structure. Instead we can make use of the
analytically derived solutions for a homogeneous core! The solutions
place constraints on all six $y$-functions at the Core-Mantle Boundary,
hence over constraining the system. To combat this, we now have three
solution vectors (instead of just 1) at the CMB. These solution vectors
are so-called elementary solutions ($\mathbf{y}^{(i)}, i = 1,2,3$) of an
isotropic homogeneous elastic-gravitational sphere. The true solution is
naturally some weighted average of the three solutions, hence
introducing three new unknowns: we are saved the system is no longer
overestimated!

Let us define $\pmb{C}$ the vector of constants of integration (weights)
$$\pmb{C} = (a_1, a_2, a_3),$$ these constants of integration must be
determined using the boundary condition at the Earth's surface for the
stress components of the spheroidal vector solution. The (true) modal
solution vector can be written as a superposition of the elementary
solutions $y^{(i)}$ which are regular at the centre of the body,

$$\pmb{y}_n(r_C^+) = a_1 \mathbf{y}^{(1)} + a_2 \mathbf{y}^{(2)} + a_3 \mathbf{y}^{(3)},$$

implying that as our true solution $\pmb{y}_n(r_C^+)$ approaches the
core from above ($r_C^+$) it must smoothly (for mathematicians:
smoothness class $C^0$; continuous) transition into the core solution
vector (which is then the weighted superposition of the elementary
solutions). So just to be clear, the core also only has one solution,
but it can be expressed as three elementary solutions. The moment we fix
$\pmb{C}$ we implicitly fix also the core solution!

We can now introduce the exact elementary solutions at the CMB.
Equations (95) for toroidal modes and (98) and (100) for spheroidal
modes of Takeuchi & Saito (1972) are adopted for the inner boundary
conditions. We express the spheroidal vector solution at the bottom of
the solid mantle as 

$$\pmb{y}_n(r_C^+) = \pmb{I}_C \pmb{C}$$

where
$\pmb{I}_C$ is the core-mantle boundary (CMB) matrix for a fluid core
(quantities correspond to core density/radius/surface gravity!) 

$$\small
\pmb{I}_C = \\
\begin{pmatrix}
-r^n/g & 0 & 1 \\[1.2em]
0 & 1 & 0 \\[1.2em]
0 & 0 & g/\rho \\[1.2em]
0 & 0 & 0 \\[1.2em]
r^n & 0 & 0 \\[1.2em]
2(n-1)r^{(n-1)} & 0 & 4 \pi G \rho
\end{pmatrix}$$ 

The three columns correspond to the elementary
solutions. Similar to how we were able to define different
$\pmb{A}_l(r)$, we may also impose different effects on our core,
yielding different elementary solutions. Specifically, we may impose
that the core is a incompressible solid, or that it is a rotating fluid
core, taking into account inertia effects at high forcing frequencies.
For now, this simple and elegant solution will suffice, though the
inclusion of inertia effects will be essential for the Earth-Moon
system, this does not change the method so we can proceed like this for
now.

Given the above you may skip this subsection, here I will briefly
introduce $\pmb{I}_C$ for the compressible fluid core with inertial
effects (assuming shear is zero.) We start with the aforementioned
elementary solutions from Equations from Takeuchi & Saito (1972)

$$\begin{aligned}
    y_1 (r) &= n r^{n-1} \\
    y_2 (r) &= r^{n-1} \\
    y_3 (r) &= 0 \\
    y_4 (r) &= 0 \\
    y_5 (r) &= -(n\gamma - \omega^2 ) r^n \\
    y_6 (r) &= -[2(n-1)n\gamma - (2n+1)\omega^2] r^{n-1} 
\end{aligned}$$ 

where $\gamma = 4\pi G\rho /3$. The second elementary
solution is more troublesome 

$$\begin{aligned}
    y_1 (r) &= - \frac{r^{l+1}}{2n+3} \left[ \frac12 n h \psi_n(x) + f \phi_{n+1} (x) \right] \\
    y_2 (r) &= - \frac{r^{l+1}}{2n+3} \left[ \frac12 h \psi_n(x) - \phi_{n+1} (x) \right] \\
    y_3 (r) &= - \lambda r^n f \phi_n (x) \\
    y_4 (r) &= 0 \\
    y_5 (r) &= - r^{n+2} \left[ \frac{\alpha^2 f}{r^2} - \frac{3\gamma f}{2(2n+3)} \psi_n (x) \right] \\
    y_6 (r) &= - r^{n+1} \left[ \frac{(2n+1) \alpha^2 f}{r^2} - \frac{3\gamma [ (2n+1)f - nh]}{2(2n+3)}\psi_n(x)\right]
\end{aligned}$$ 

where 

$$\begin{aligned}
    x &= kr \\
    k^2 &= \frac{1}{\alpha^2} \left[ \omega^2 + 4\gamma - \frac{n(n+1)\gamma^2}{\omega^2} \right] \\
    f &= -\frac{\omega^2}{\gamma} \\
    h &= f - (n+1) \\
    \phi_n (x) &= \frac{(2n+1)!!}{x^n} j_n (x) \\
    \psi_n (x) &= \frac{2(2n+3)}{x^2} [1-\phi_n (x)].
\end{aligned}$$ 

Here $\alpha$ is the compressional wave speed and
$j_n (x)$ is the spherical Bessel function of the first kind. The issue
is that as we decrease the forcing frequency $\omega$, $k^2$ becomes
negative, and $x$ becomes large and imaginary. The spherical Bessel
function of the first kind grows exponentially, and the elementary
solution rapidly grows beyond numerical precision. One potential fix
would be to use non-dimensional scaling, however this is currently not
working properly. Another way to resolve this is to note that the
elementary solution may be scaled entirely by some arbitrary constant
($C_2$), as noted earlier in the text. As such, we may simply divide by
$j_n(x)$ to remain around unity. Note that without this the solution
diverges and the system cannot be resolved. With this in mind we need
the recursion relations for the spherical Bessel function of the first
kind, define 

$$z_{n} (x) = x j_{n+1} (x) / j_n (x)$$ 

such that

$$z_{n-1} (x) = \frac{x^2}{(2n+1) - z_n(x)}$$ 

This recursion calculation
should be carried out in order of decreasing $n$ starting from a
sufficiently large wavenumber with an initial value
$z_n(x) = x^2/(2n+3)$. We can now do some algebraic manipulation of the
elementary solution to obtain a form with only terms $\propto 1/j_n(x)$
or unity. Since, as mentioned, $j_n(x)$ will diverge for large imaginary
$x$, the newly obtained solution will be well-behaved around unity. Let
us define 

$$\begin{aligned}
    \bar{\phi}_n (x) &\equiv \frac{\phi_n (x)}{j_n (x)} = \frac{(2n+1)!!}{x^n} \\
    \bar{\phi}_{n+1} (x) &\equiv \frac{\phi_{n+1} (x)}{j_n (x)} = \frac{(2n+3)!!}{x^{n+1}} \frac{j_{n+1}(x)}{j_n (x)} = \frac{(2n+3)!!}{x^{n+2}}z_n(x) \\
    \bar{\psi}_n (x) &= \frac{\psi_n(x)}{j_n(x)} = \frac{2(2n+3)}{x^2} \left[ \frac{1}{j_n (x)} - \bar{\phi}_n (x)\right]
\end{aligned}$$ 

The elementary solution is then written:

$$\begin{aligned}
    \bar{y}_1 (r) &= - \frac{r^{n+1}}{2n+3} \left[ \frac{1}{2} n h \bar{\psi}_n(x) + f \frac{(2n+3)!!}{x^{n+2}} z_n(x) \right] \\
    \bar{y}_2 (r) &= - \frac{r^{n+1}}{2n+3} \left[ \frac{1}{2} h \bar{\psi}_n(x) - \frac{(2n+3)!!}{x^{n+2}} z_n(x) \right] \\
    \bar{y}_3 (r) &= - \lambda r^n f \bar{\phi}_n (x) \\
    \bar{y}_4 (r) &= 0 \\
    \bar{y}_5 (r) &= - r^{n+2} \left[ \frac{\alpha^2 f}{r^2} \frac{1}{j_n(x)} - \frac{3\gamma f}{2(2n+3)} \bar{\psi}_n (x) \right] \\
    \bar{y}_6 (r) &= - r^{n+1} \left[ \frac{(2n+1) \alpha^2 f}{r^2} \frac{1}{j_n(x)} - \frac{3\gamma [ (2n+1)f - nh]}{2(2n+3)}\bar{\psi}_n(x)\right]
\end{aligned}$$ 

Now we note that we can simplify the system of equations
when $k^2 \ll -1$, which implies 

$$\begin{aligned}
    \omega^2 + 4\gamma &\ll \frac{n(n+1)\gamma^2}{\omega^2} \\
    \omega^4 + 4\gamma \omega^2 &\ll n(n+1)\gamma^2 
\end{aligned}$$

Noting that the LHS is dominated by $4\gamma \omega^2$ as $\omega \to 0$, we have

$$\begin{aligned}
    4\gamma \omega^2 &\ll n(n+1)\gamma^2 \\
    4\omega^2 &\ll n(n+1)\gamma \qquad \wedge \qquad \gamma \neq 0
\end{aligned}$$ 

We can now find a typical forcing frequency for the
Earth's iron core, below which we may further simplify the system. Let's
assume for the outer liquid core
$\rho \approx 10000 \text{ to } 12000 \text{ kg/m}^3$. Then, plugging in
to obtain $\gamma$, we have

$$\gamma \approx \frac{4}{3} \pi (6.67 \times 10^{-11}) (11000) \approx 3 \times 10^{-6} \text{ s}^{-2}$$

For a degree $n=2$, that becomes

$$\frac{n(n+1)}{4} \gamma = \frac{6}{4} (3 \times 10^{-6}) = 4.5 \times 10^{-6} \text{ s}^{-1}$$

So that would imply that $\omega \ll 10^{-3} \text{ s}^{-2}$. In this
regime we should use 

$$\begin{aligned}
    \bar{y}_1 (r) &= - \frac{r^{n+1}}{2n+3} \left[ \frac{n h (2n+3)}{x^2} (-\bar{\phi}_n) + f \frac{(2n+3)!!}{x^{n+2}} z_n(x) \right] \\
    \bar{y}_2 (r) &= - \frac{r^{n+1}}{2n+3} \left[ \frac{h (2n+3)}{x^2} (-\bar{\phi}_n) - \frac{(2n+3)!!}{x^{n+2}} z_n(x) \right] \\
    \bar{y}_3 (r) &= - \lambda r^n f \bar{\phi}_n (x) \\
    \bar{y}_4 (r) &= 0 \\
    \bar{y}_5 (r) &= - r^{n+2} \left[ - \frac{3\gamma f}{2(2n+3)} \bar{\psi}_n (x) \right] \\
    \bar{y}_6 (r) &= - r^{n+1} \left[ - \frac{3\gamma [ (2n+1)f - nh]}{2(2n+3)}\bar{\psi}_n(x)\right]
\end{aligned}$$ 

which further simplify to 

$$\begin{aligned}
    \bar{y}_1 (r) &= - \frac{r^{n+1}(2n+1)!!}{x^{n+2}} \left[ f z_n(x) - nh \right] \\
    \bar{y}_2 (r) &= - \frac{r^{n+1}(2n+1)!!}{x^{n+2}} \left[ -z_n(x) -h \right] \\
    \bar{y}_3 (r) &= - \frac{\lambda f r^n (2n+1)!!}{x^n} \\
    \bar{y}_4 (r) &= 0 \\
    \bar{y}_5 (r) &= - \frac{3\gamma f r^{n+2} (2n+1)!!}{x^{n+2}} \\
    \bar{y}_6 (r) &= - \frac{3\gamma [(2n+1)f - nh] r^{n+1} (2n+1)!!}{x^{n+2}}
\end{aligned}$$ 

Every term now shares the common factor
$\mathcal{K} = \frac{(2n+1)!!}{x^{n+2}}$. If you wish, you can even
divide the entire solution by $\mathcal{K}$ (since the solution is
scalable by an arbitrary constant $C_2$), which leaves us with a
remarkably clean system: 

$$\begin{aligned}
    \bar{y}_1 (r) &= -r^{n+1} [f z_n(x) - nh] \\
    \bar{y}_2 (r) &= -r^{n+1} [-z_n(x) - h] \\
    \bar{y}_3 (r) &= -\lambda f r^n x^2 \\
    \bar{y}_4 (r) &= 0 \\
    \bar{y}_5 (r) &= -3\gamma f r^{n+2} \\
    \bar{y}_6 (r) &= -3\gamma [(2n+1)f - nh] r^{n+1}
\end{aligned}$$

Finally, we should note that the porosity related $y$-functions need
also be bounded from below and partially from above. Moreover, since a
porous layer is likely to form at some point midway through the mantle,
and not just at the CMB, we will need to account for these boundary
conditions in a slightly different way. Think of the lower boundary as
an extension to the already existing elementary solutions, instead of
three there are now four such elementary solutions. As long as you know
the current form (at some radius) of the three elementary solutions, you
can simply add one more, then forward propagate the four solutions to
the surface (or where ever the porous layer ends) and apply the boundary
conditions directly to obtain the weights
$\pmb{C} = (a_1, a_2, a_3, a_4)$. This is basically the approach taken
in the so-called "shooting method". However, as we will see, this
becomes more difficult once we move to solution/transfer space where we
do not have direct access to the current form of the three elementary
solutions, but rather only know how the solution vector $\pmb{y}_n(r)$
(this can be any solution vector, not just the true solution) transforms
between layers. This is a central aspect of the "relaxation method" (or
Henyey style solver) that we will employ as it is substantially more
stable then the shooting method.

For now we will finish this section by stating the elementary boundary
conditions to be imposed on the poro-viscoelastic solution vector. At
the lower boundary of the porous layer the pore pressure is nonzero
$b_{(7,4)} = 1$ and we have zero radial Darcy flux $b_{(8,4)} = 1$. At
the upper part of the porous layer we impose again zero radial Darcy
flux $b_{(8)} = 1$, without any constraint on the pore pressure.

We will now discuss the different solver schemes.
