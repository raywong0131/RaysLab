# Supplementary Theory



> [!abstract] Draft status and theoretical spine
> This note is a working framework for the complete Supplementary Theory.
> Sections S1--S4 contain the derivation spine and the equations to retain.
> Section S5 develops the finite-cavity wave-vector decomposition and its connection to radiation.
> Editorial callouts and Obsidian comments will be removed before submission.

The theory separates the closed-system band problem from the open-system radiation problem. The Hermitian Maxwell operator determines the band ordering and gap. A distinct radiation operator determines the outgoing amplitudes. The target condition is

$$
\boxed{
D(\Gamma)|p_y\rangle=\mathbf 0,
\qquad
\Delta_y(\Gamma)\neq 0,
}
$$

where $|p_y\rangle$ denotes the full electromagnetic mode with dominant $y$-polarized radiation. The quantity $\Delta_y$ is its symmetry-resolved band separation from the coupled $d$-like state.

The first component of $D(\Gamma)|p_y\rangle$ vanishes by $C_{2v}$ symmetry. The second component is symmetry allowed. Geometry-dependent interference can tune that component to zero without closing $\Delta_y$.

This distinction provides the working definition used here:

> [!important] Working definition of the dipolar singularity
> A dipolar singularity is a radiation zero of a dipole-like $\Gamma$-point Bloch state in a $C_{2v}$-anisotropic photonic Dirac lattice. One electric-field component is forbidden by symmetry. The remaining allowed component vanishes through geometric interference. The corresponding Hermitian band gap remains open.

## Drafting map

| Section | Central question | Main source | Current status |
|---|---|---|---|
| S1 | Why may the photonic crystal be treated with Bloch, tight-binding, and $k\cdot p$ methods? | Maxwell equations and existing Dirac notes | Derivation spine written |
| S2 | How do $\eta$ and $\zeta$ change couplings, eigenvectors, and symmetry-resolved gaps? | [[Dirac_C2_dipolar_singularity_theory_SI]] | Framework written; microscopic label map needs validation |
| S3 | How do the normal-radiation channels evolve from $C_{6v}$ to $C_{2v}$, and why can an allowed $p_y$ channel cross zero? | [[C2v_py_mode_Ex_Ey_symmetry_supplementary_note]] | Channel definition, selection rules, and zero-existence argument written |
| S4 | How does geometry tune the allowed radiation amplitude to zero? | [[SI_theory_WuHu_radiation_integral]] | Exact reduced-model identity and response framework written |
| S5 | How do nearby Bloch wave vectors form the finite-cavity far field? | Finite-lattice Fourier and Bloch-decomposition notes | Decomposition framework written; boundary mixing requires validation |

> [!warning] Basis convention to lock before final typesetting
> The imported notes use two different $p_x/p_y$ conventions.
> One convention labels the scalar $H_z$ pattern.
> The other labels the full mode by its dominant electric-field polarization.
> This manuscript uses the second convention.
>
> For sites numbered $n=1,\ldots,6$ at $\theta_n=(n-1)\pi/3$, the target $p_y$ mode has $H_z\propto\cos\theta_n$ under the coordinates used below.
> Its full-vector representation is $B_2$ and its dominant allowed radiation is $E_y$.
> Any microscopic formula imported from an $H_z$-orbital basis must first be translated into this convention.

%%
SOURCE SYNTHESIS
- Strict result: C2v symmetry forces c_x(0,0)=0 for the full-vector p_y mode at Gamma.
- Strict result: c_y is allowed and may vanish only through an additional interference condition.
- Strict result at $\Gamma$ in the 2D TE model: the area integral reduces to an interface line integral. Finite $\mathbf k$ adds an explicit unit-cell area term.
- Model result: local boundary moments may be linearized in eta and zeta.
- Model result: common local susceptibilities can predict opposite px/py shifts, but this requires full-wave validation.
- Do not identify a microscopic target mass formula until the H_z-orbital labels are translated.
%%

## S1. Maxwell eigenproblem and the photonic-crystal master equation

We consider a linear, nonmagnetic, and source-free dielectric. The relative permittivity $\varepsilon(\mathbf r)$ is real and periodic in the lossless reference problem. We use the time dependence $e^{-i\omega t}$. Maxwell's equations are

$$
\nabla\times\mathbf E
=i\omega\mu_0\mathbf H,
\qquad
\nabla\times\mathbf H
=-i\omega\varepsilon_0\varepsilon(\mathbf r)\mathbf E.
$$

Eliminating $\mathbf E$ gives the magnetic-field master equation

$$
\boxed{
\nabla\times
\frac{1}{\varepsilon(\mathbf r)}
\nabla\times\mathbf H(\mathbf r)
=
\left(\frac{\omega}{c}\right)^2
\mathbf H(\mathbf r).
}
\tag{S1}
$$

For real positive $\varepsilon$, the operator on the left-hand side of Eq. (S1) is Hermitian on the transverse-field subspace with periodic boundary conditions. Its eigenvalues are $(\omega/c)^2$.

%% CITATION NEEDED: standard photonic-crystal master equation and Hermiticity. %%

The dielectric function obeys

$$
\varepsilon(\mathbf r+\mathbf R)
=\varepsilon(\mathbf r)
$$

for every Bravais vector $\mathbf R$. The eigenfields can therefore be written as

$$
\mathbf H_{n\mathbf k}(\mathbf r)
=e^{i\mathbf k\cdot\mathbf r}
\mathbf u_{n\mathbf k}(\mathbf r),
$$

where $\mathbf u_{n\mathbf k}$ is lattice periodic. Substitution into Eq. (S1) yields

$$
\boxed{
\hat\Theta_{\mathbf k}\mathbf u_{n\mathbf k}
=
\left(\frac{\omega_{n\mathbf k}}{c}\right)^2
\mathbf u_{n\mathbf k},
}
\tag{S2}
$$

with

$$
\hat\Theta_{\mathbf k}
=
(\nabla+i\mathbf k)\times
\frac{1}{\varepsilon(\mathbf r)}
(\nabla+i\mathbf k)\times.
$$

Equation (S2) permits Bloch-band analysis, localized-basis projection, and a $k\cdot p$ expansion near $\Gamma$.

For a TE mode, the nonzero field components are

$$
\mathbf H=H_z(x,y)\hat{\mathbf z},
\qquad
\mathbf E=E_x\hat{\mathbf x}+E_y\hat{\mathbf y}.
$$

The corresponding scalar master equation is

$$
\boxed{
-\nabla_{\parallel}\cdot
\left[
\frac{1}{\varepsilon(x,y)}
\nabla_{\parallel}H_z
\right]
=
\left(\frac{\omega}{c}\right)^2H_z.
}
\tag{S3}
$$

For all subsequent derivations, we set $c=1$ and measure frequency relative to the Dirac-point frequency, so that $\omega_D=0$.

## S2. Six-site tight-binding and $k\cdot p$ theory near $\Gamma$

The radial modulation $\eta$ controls the common Dirac mass of the $C_{6v}$ lattice. The $C_{2v}$ modulation $\zeta$ resolves this mass into two anisotropic masses and shifts the two standing-wave sectors in opposite directions.

### S2.1 Six-site Hamiltonian

We use the site basis

$$
\mathcal B_{\mathrm{site}}
=\{|1\rangle,|2\rangle,|3\rangle,
|4\rangle,|5\rangle,|6\rangle\}.
$$

The Bravais vectors are

$$
\mathbf a_1=a(1,0),
\qquad
\mathbf a_2=\frac{a}{2}(1,\sqrt3),
\qquad
\mathbf a_3=\mathbf a_2-\mathbf a_1.
$$

The parameters $t_0(\eta)$ and $t_1(\eta)$ denote the intracell and intercell nearest-neighbour couplings. We choose negative off-diagonal hopping amplitudes. This convention places the $p$ doublet below the $d$ doublet for $t_0>t_1$.

In the $C_{6v}$ structure, all six sites and all symmetry-related bonds are equivalent. With $\mathcal B_{\mathrm{site}}=\{|1\rangle,|2\rangle,|3\rangle,|4\rangle,|5\rangle,|6\rangle\}$ as the ordered basis, the Bloch Hamiltonian is

$$
\boxed{
H_{\mathrm{site}}^{C_{6v}}(\mathbf k)
=
\begin{pmatrix}
0&-t_0&0&-t_1e^{i\mathbf k\cdot\mathbf a_1}&0&-t_0\\
-t_0&0&-t_0&0&-t_1e^{i\mathbf k\cdot\mathbf a_2}&0\\
0&-t_0&0&-t_0&0&-t_1e^{i\mathbf k\cdot\mathbf a_3}\\
-t_1e^{-i\mathbf k\cdot\mathbf a_1}&0&-t_0&0&-t_0&0\\
0&-t_1e^{-i\mathbf k\cdot\mathbf a_2}&0&-t_0&0&-t_0\\
-t_0&0&-t_1e^{-i\mathbf k\cdot\mathbf a_3}&0&-t_0&0
\end{pmatrix}.
}
\tag{S4}
$$

To introduce the $C_{2v}$ modulation, the hole sizes at sites 1 and 4 are changed while the remaining four hole sizes compensate the total area:

$$

b_{1,4}=\zeta b,
\qquad
b_{2,3,5,6}=\sqrt{\frac{3-\zeta^2}{2}}b,
\qquad

\tag{S5}
$$

The geometric parameter $\zeta$ controls three bond-coupling ratios, $\alpha(\zeta)$, $\beta(\zeta)$, and $\gamma(\zeta)$, defined relative to their $C_{6v}$ values. These ratios vary monotonically over the modulation range considered. As $\zeta$ increases, $\alpha$ and $\beta$ increase, whereas $\gamma$ decreases. Their local variation rates satisfy $\partial_\zeta\alpha>0$, $\partial_\zeta\beta>0$, $\partial_\zeta\gamma<0$, and $|\partial_\zeta\alpha|>|\partial_\zeta\gamma|>|\partial_\zeta\beta|$. The onsite parameter $\mu(\zeta)$ represents a comparatively independent onsite modulation and obeys $\mu(\zeta)\geq0$ under the adopted convention. These inequalities specify qualitative trends rather than fixed functional forms.

Equation (S5) preserves the total triangular-hole area. The implemented hole orientation retains two in-plane mirrors, so the point group changes from $C_{6v}$ to $C_{2v}$. The onsite response is $2\mu$ on sites 1 and 4 and $-\mu$ on the other four sites.

With $\mathcal B_{\mathrm{site}}=\{|1\rangle,|2\rangle,|3\rangle,|4\rangle,|5\rangle,|6\rangle\}$ as the ordered basis, the complete $C_{2v}$ Bloch Hamiltonian is

$$
\boxed{
H_{\mathrm{site}}^{C_{2v}}(\mathbf k)
=
\begin{pmatrix}
2\mu &-t_0\beta&0&-t_1\alpha e^{i\mathbf k\cdot\mathbf a_1}&0&-t_0\beta\\
-t_0\beta&-\mu&-t_0\gamma&0&-t_1\gamma e^{i\mathbf k\cdot\mathbf a_2}&0\\
0&-t_0\gamma&-\mu&-t_0\beta&0&-t_1\gamma e^{i\mathbf k\cdot\mathbf a_3}\\
-t_1\alpha e^{-i\mathbf k\cdot\mathbf a_1}&0&-t_0\beta&2\mu&-t_0\beta&0\\
0&-t_1\gamma e^{-i\mathbf k\cdot\mathbf a_2}&0&-t_0\beta&-\mu&-t_0\gamma\\
-t_0\beta&0&-t_1\gamma e^{-i\mathbf k\cdot\mathbf a_3}&0&-t_0\gamma&-\mu
\end{pmatrix}.
}
\tag{S6}
$$

When the structure returns to $C_{6v}$, $\alpha=\beta=\gamma=1$ and $\mu=0$, so Eq. (S6) reduces exactly to Eq. (S4). The additional condition $\eta=1$ gives $t_0=t_1$.

### S2.2 Six-site standing-wave basis

At $\Gamma$, the $C_{6v}$ Hamiltonian is diagonal in six real standing-wave vectors. The coefficients of these vectors describe the scalar $H_z$ distribution on the six sites:

$$
\boxed{
\begin{aligned}
|s\rangle&=\frac{1}{\sqrt6}(1,1,1,1,1,1)^{\mathsf T},\\
|p_x\rangle&=\frac{1}{2}(0,1,1,0,-1,-1)^{\mathsf T},\\
|p_y\rangle&=\frac{1}{\sqrt{12}}(-2,-1,1,2,1,-1)^{\mathsf T},\\
|d_{x^2-y^2}\rangle&=\frac{1}{\sqrt{12}}(-2,1,1,-2,1,1)^{\mathsf T},\\
|d_{xy}\rangle&=\frac{1}{2}(0,1,-1,0,1,-1)^{\mathsf T},\\
|f\rangle&=\frac{1}{\sqrt6}(1,-1,1,-1,1,-1)^{\mathsf T}.
\end{aligned}
}
\tag{S7}
$$

The $p_x$ and $p_y$ labels follow the dominant electric-field polarization. Consequently, $H_z^{(p_x)}\propto\sin\theta_n$ and $H_z^{(p_y)}\propto-\cos\theta_n$. This convention agrees with the component analysis in Section S3.

Let the ordered six-state basis and its transformation matrix be

$$
\mathcal B_6=
\{|s\rangle,|p_x\rangle,|p_y\rangle,
|d_{x^2-y^2}\rangle,|d_{xy}\rangle,|f\rangle\},
\qquad
U_6=
\big(
|s\rangle,|p_x\rangle,|p_y\rangle,
|d_{x^2-y^2}\rangle,|d_{xy}\rangle,|f\rangle
\big).
$$

For the $C_{6v}$ structure, the Hamiltonian in the basis $\mathcal B_6$ is

$$
\boxed{
U_6^\dagger H_{\mathrm{site}}^{C_{6v}}(\Gamma)U_6
=
\operatorname{diag}
\left(
-2t_0-t_1,-m,-m,m,m,2t_0+t_1
\right),
\qquad
m=t_0-t_1.
}
\tag{S8}
$$

The $p$ and $d$ doublets become fourfold degenerate when $m=0$. This occurs at $t_0=t_1$ in the nearest-neighbour model.

### S2.3 Evolution of the four-band Hamiltonian from $C_{6v}$ to $C_{2v}$

The six sites are numbered counter-clockwise from the positive $x$ axis, with

$$
\mathbf R_n=R(\cos\theta_n,\sin\theta_n),
\qquad
\theta_n=\frac{(n-1)\pi}{3},
\qquad
n=1,\ldots,6.
$$

The radial modulation is $\eta=3R/a$. It preserves $C_{6v}$ while changing the common mass through $t_0(\eta)-t_1(\eta)$. The four-band Hamiltonian of this familiar $C_{6v}$ structure is stated first. Appendix A gives the Löwdin reduction from the six-site model and the frequency-based parameter extraction. The rotating states are defined from the standing waves in Eq. (S7) by

$$
\boxed{
\begin{aligned}
|d_\pm\rangle
&=\frac{|d_{x^2-y^2}\rangle\pm i|d_{xy}\rangle}{\sqrt2},\\
|p_\pm\rangle
&=-\frac{|p_x\rangle\mp i|p_y\rangle}{\sqrt2}.
\end{aligned}
}
\tag{S9}
$$

The phases in Eq. (S9) follow the electric-polarization labeling adopted in Eq. (S7). Alternative overall phases change some off-diagonal signs but leave every eigenfrequency and Dirac mass unchanged. The ordered rotating basis is

$$
\mathcal B_{\mathrm{spin}}
=\{|d_+\rangle,|p_+\rangle,
|d_-\rangle,|p_-\rangle\},
$$

whereas the ordered standing-wave basis remains

$$
\mathcal B_{\mathrm{sw}}
=\{|d_{xy}\rangle,|d_{x^2-y^2}\rangle,
|p_y\rangle,|p_x\rangle\}.
$$

With these explicit orderings, the coordinate transformation and the corresponding operator transformation are

$$
U_{\mathrm{sw}\leftarrow\mathrm{spin}}
=\frac{1}{\sqrt2}
\begin{pmatrix}
i&0&-i&0\\
1&0&1&0\\
0&i&0&-i\\
0&-1&0&-1
\end{pmatrix},
\qquad
H^{\mathrm{sw}}
=U_{\mathrm{sw}\leftarrow\mathrm{spin}}
H^{\mathrm{spin}}
U_{\mathrm{sw}\leftarrow\mathrm{spin}}^{\dagger}.
\tag{S10}
$$

In the $C_{6v}$ structure, the two rotating sectors are independent. With $m=t_0-t_1$, $v=at_1/2$, and $k_\pm=k_x\pm ik_y$, the Hamiltonian with $\mathcal B_{\mathrm{spin}}=\{|d_+\rangle,|p_+\rangle,|d_-\rangle,|p_-\rangle\}$ as the ordered basis is

$$
\boxed{
H_{C_{6v}}^{\mathrm{spin}}(\mathbf k)
=
\begin{pmatrix}
m&-vk_+&0&0\\
-vk_-&-m&0&0\\
0&0&m&vk_-\\
0&0&vk_+&-m
\end{pmatrix}.
}
\tag{S11}
$$

The same operator, with $\mathcal B_{\mathrm{sw}}=\{|d_{xy}\rangle,|d_{x^2-y^2}\rangle,|p_y\rangle,|p_x\rangle\}$ as the ordered basis, is

$$
\boxed{
H_{C_{6v}}^{\mathrm{sw}}(\mathbf k)
=
\begin{pmatrix}
m&0&-ivk_y&ivk_x\\
0&m&ivk_x&ivk_y\\
ivk_y&-ivk_x&-m&0\\
-ivk_x&-ivk_y&0&-m
\end{pmatrix}.
}
\tag{S12}
$$

Equations (S11) and (S12) describe the same $C_{6v}$ operator. The rotating basis exposes two time-reversed Dirac sectors, while the standing-wave basis exposes the real $p/d$ field patterns. At $\Gamma$, the $p$ and $d$ doublets have eigenfrequencies $-m$ and $+m$, respectively. A single Dirac mass $m$ therefore controls both polarization-resolved $p/d$ separations.

The projected Maxwell eigenproblem is written in frequency units, and all coupling parameters below use the convention fixed at the end of S1.

The $C_{2v}$ modulation removes the equivalence of the two standing-wave sectors. The two symmetry-resolved Dirac masses are

$$
\boxed{
m_x=\gamma(t_0-t_1),
\qquad
m_y=\frac{(4\beta-\gamma)t_0-(2\alpha+\gamma)t_1}{3}.
}
\tag{S13}
$$

It is useful to define

$$
\boxed{
M=\frac{m_x+m_y}{2},
\qquad
\Delta=\frac{m_y-m_x}{2},
\qquad
v_+=\frac{\chi+\gamma}{2}v,
\qquad
v_-=\frac{\chi-\gamma}{2}v,
\qquad
\chi=\frac{4\alpha-\gamma}{3}.
}
\tag{S14}
$$

When the structure returns to $C_{6v}$, these parameters reduce to

$$
\boxed{
\alpha=\beta=\gamma=1,
\quad \mu=0,
\quad m_x=m_y=m,
\quad M=m,
\quad \Delta=0,
\quad v_+=v,
\quad v_-=0.
}
\tag{S15}
$$

With $\mathcal B_{\mathrm{spin}}=\{|d_+\rangle,|p_+\rangle,|d_-\rangle,|p_-\rangle\}$ as the ordered basis, the first-order $C_{2v}$ Hamiltonian is

$$
\boxed{
H_{C_{2v}}^{\mathrm{spin}}(\mathbf k)
=
\begin{pmatrix}
M&-v_+k_x-i\gamma vk_y&\Delta+\mu&v_-k_x\\
-v_+k_x+i\gamma vk_y&-M&-v_-k_x&\Delta-\mu\\
\Delta+\mu&-v_-k_x&M&v_+k_x-i\gamma vk_y\\
v_-k_x&\Delta-\mu&v_+k_x+i\gamma vk_y&-M
\end{pmatrix}
+\mathcal O(k^2).
}
\tag{S16}
$$

The $k$-independent entries $\Delta+\mu$ and $\Delta-\mu$ couple the two rotating $d$ states and the two rotating $p$ states, respectively. Hence, the pseudospin block structure visible in Eq. (S11) is no longer preserved. The additional coefficient $v_-$ and the unequal $k_x$ and $k_y$ couplings record the velocity anisotropy allowed by $C_{2v}$.

With $\mathcal B_{\mathrm{sw}}=\{|d_{xy}\rangle,|d_{x^2-y^2}\rangle,|p_y\rangle,|p_x\rangle\}$ as the ordered basis, transforming Eq. (S16) gives

$$
\boxed{
H_{C_{2v}}^{\mathrm{sw}}(\mathbf k)
=
\begin{pmatrix}
m_x-\mu&0&-i\gamma vk_y&i\gamma vk_x\\
0&m_y+\mu&i\chi vk_x&i\gamma vk_y\\
i\gamma vk_y&-i\chi vk_x&-m_y+\mu&0\\
-i\gamma vk_x&-i\gamma vk_y&0&-m_x-\mu
\end{pmatrix}
+\mathcal O(k^2).
}
\tag{S17}
$$

The real standing waves diagonalize Eq. (S17) at $\Gamma$, but not at a generic two-dimensional wave vector. Along $k_x$, the coupled pairs are $(d_{xy},p_x)$ and $(d_{x^2-y^2},p_y)$. Along $k_y$, the coupled pairs are $(d_{xy},p_y)$ and $(d_{x^2-y^2},p_x)$.

Symmetry reduction therefore changes a single mass into two independent masses rather than merely perturbing one isotropic gap. The corresponding $\Gamma$-point separations are

$$
\boxed{
\Delta_x
=\left|\omega_{d_{xy}}-\omega_{p_x}\right|
=2|m_x|,
\qquad
\Delta_y
=\left|\omega_{d_{x^2-y^2}}-\omega_{p_y}\right|
=2|m_y|.
}
\tag{S18}
$$

The target electric-polarized $p_y$ state belongs to the $m_y$ sector. Section S2 establishes how the $C_{2v}$ modulation changes its Hermitian band separation. Radiation is a separate open-system question addressed in Sections S3 and S4.

## S3. Normal-radiation selection rules and the tunable $\Gamma$-point zero

### S3.1 Open radiation channels and the $C_{6v}$ structure

Below the first diffraction threshold, a periodic slab at $\Gamma$ has only the zeroth-order plane-wave radiation channels. Let $|\mathcal R_\alpha\rangle$, with $\alpha\in\{x,y\}$, denote the normalized outgoing plane wave with in-plane electric polarization $\hat{\mathbf e}_\alpha$ and the required vertical parity. Its complex coupling to a periodic Bloch mode $|\psi\rangle$ is

$$
\boxed{
c_\alpha(\Gamma)
=\langle\mathcal R_\alpha|\psi\rangle
\propto
\frac{1}{A_{\mathrm{cell}}}
\int_{\mathrm{cell}}E_\alpha(\boldsymbol\rho,z_0)\,d^2\rho.
}
\tag{S19}
$$

Equation (S19) defines the physical meaning of $c_x$ and $c_y$: each is a complex zeroth-order outgoing plane-wave amplitude on an exterior plane. An in-plane symmetry permits coupling only when the mode and radiation-channel representations can form the identity representation,

$$
\boxed{
\Gamma_{\mathrm{mode}}\otimes
\Gamma_{\mathrm{rad}}
\supset A_1.
}
\tag{S20}
$$

Compatibility with the vertical mirror parity is an additional condition and is assumed for the TE-like modes considered here. At $\Gamma$, the zeroth-order plane wave is spatially uniform within the unit cell. Its scalar spatial part therefore belongs to $A_1$, while its two in-plane polarization vectors together transform as $E_1$. The complete normal-radiation channel consequently belongs to $E_1$.

In the $C_{6v}$ structure, the dipolar doublet belongs to $E_1$ and the quadrupolar doublet belongs to $E_2$. Substituting these representations into Eq. (S20) gives

$$
\boxed{
\begin{aligned}
&\Gamma_{\mathrm{rad}}^{C_{6v}}=E_1,
\qquad
\Gamma(p_x,p_y)=E_1,
\qquad
\Gamma(d_{x^2-y^2},d_{xy})=E_2,\\[2pt]
&E_1\otimes E_1
=A_1\oplus A_2\oplus E_2
\supset A_1,\\
&E_2\otimes E_1
=B_1\oplus B_2\oplus E_1
\not\supset A_1.
\end{aligned}
}
\tag{S21}
$$

The first product contains $A_1$, so the $p$ doublet can couple to the normal free-space channel. The second product does not contain $A_1$, so the $d$ doublet cannot couple to that channel at $\Gamma$. More generally, for the present vertical parity and diffraction order, a $C_{6v}$ mode can radiate normally only through its $E_1$ component. The familiar $C_{6v}$ case therefore shows that a $p$-like state is not generically dark.

### S3.2 Symmetry descent from $C_{6v}$ to $C_{2v}$

The surviving operations are defined by

$$
C_2(z):(x,y,z)\rightarrow(-x,-y,z),
$$

$$
\sigma_x\equiv\sigma_{yz}:(x,y,z)\rightarrow(-x,y,z),
$$

$$
\sigma_y\equiv\sigma_{xz}:(x,y,z)\rightarrow(x,-y,z).
$$

The corresponding character table is

| Irrep | $E$ | $C_2(z)$ | $\sigma_x$ | $\sigma_y$ | Representative basis |
| ----- | --: | -------: | ---------: | ---------: | -------------------- |
| $A_1$ |   1 |        1 |          1 |          1 | $1,z,x^2,y^2$        |
| $A_2$ |   1 |        1 |         -1 |         -1 | $R_z,xy$             |
| $B_1$ |   1 |       -1 |         -1 |          1 | $x,R_y$              |
| $B_2$ |   1 |       -1 |          1 |         -1 | $y,R_x$              |

Reducing the $C_{6v}$ representations to $C_{2v}$ gives

$$
\boxed{
E_1\downarrow C_{2v}=B_1\oplus B_2,
\qquad
E_2\downarrow C_{2v}=A_1\oplus A_2.
}
\tag{S22}
$$

Consequently, $p_x$ and the $x$-polarized outgoing channel transform as $B_1$, while $p_y$ and the $y$-polarized outgoing channel transform as $B_2$. The $C_{2v}$ modulation splits the $E_1$ doublet but does not forbid its matched radiation channels. In particular, the $p_y$ mode remains symmetry allowed to emit through $E_y$ at $\Gamma$.

Indeed, the matched and crossed products are

$$
B_1\otimes B_1=A_1,
\qquad
B_2\otimes B_2=A_1,
\qquad
B_1\otimes B_2=A_2.
$$

Thus, $p_x$ can couple to the $x$-polarized channel and $p_y$ can couple to the $y$-polarized channel, whereas the crossed polarization couplings are forbidden at $\Gamma$.

### S3.3 Component selection and existence of the allowed-channel zero

The full target $p_y$ mode belongs to $B_2$. Because $H_z$ is an axial component with intrinsic representation $A_2$, its scalar spatial distribution transforms as

$$
\boxed{
\Gamma(H_z)=B_2\otimes A_2=B_1.
}
\tag{S23}
$$

The dielectric function transforms as $A_1$, and $\partial_x$ and $\partial_y$ transform as $B_1$ and $B_2$, respectively. Within the TE-like reduction,

$$
\boxed{
E_x\propto\frac{1}{\varepsilon}\partial_yH_z,
\qquad
E_y\propto-\frac{1}{\varepsilon}\partial_xH_z.
}
\tag{S24}
$$

The scalar electric-field distributions therefore satisfy

$$
\boxed{
\Gamma(E_x)=A_2,
\qquad
\Gamma(E_y)=A_1.
}
\tag{S25}
$$

Only an $A_1$ scalar distribution can have a nonzero unit-cell average. Combining Eqs. (S19) and (S25) gives

$$
\boxed{
c_x(\Gamma;\lambda)=0
\quad\text{for every mirror-preserving }C_{2v}\text{ geometry},
\qquad
c_y(\Gamma;\lambda)=c_{y0}(\lambda)
\quad\text{is symmetry allowed}.
}
\tag{S26}
$$

Here, $\lambda$ denotes a continuous $C_{2v}$-preserving geometric modulation, and $c_{y0}(\lambda)$ is the physical $y$-polarized outgoing amplitude defined by Eq. (S19). Equation (S26) distinguishes two different statements. The zero of $c_x$ is enforced by the mirrors. Symmetry places no zero constraint on $c_{y0}$.

An allowed matrix element can nevertheless vanish because different parts of the periodic field contribute with opposite signs. For a lossless reciprocal slab, a symmetry-compatible phase convention makes $c_{y0}(\lambda)$ a continuous signed coefficient. Once two $C_{2v}$ geometries bracket a sign reversal, continuity proves the existence of an intermediate zero:

$$
\boxed{
c_{y0}(\lambda_-)c_{y0}(\lambda_+)<0
\quad\Longrightarrow\quad
\exists\,\lambda_*\in(\lambda_-,\lambda_+)
\text{ such that }c_{y0}(\lambda_*)=0.
}
\tag{S27}
$$

Locally, if $\partial_\lambda c_{y0}|_{\lambda_0}\neq0$, its position is

$$
\lambda_*
\simeq\lambda_0-
\frac{c_{y0}(\lambda_0)}
{\left.\partial_\lambda c_{y0}\right|_{\lambda_0}}.
$$

The geometry sweep and the interface formulation in S4 supply the required sign and slope. The apparently counterintuitive result is therefore precise: $E_y$ is an allowed radiation channel for the periodic $p_y$ state, but its overlap can be tuned continuously through zero. The zero is an interference zero inside an allowed channel, not a second symmetry prohibition and not a band-degeneracy condition.

For comparison, the two dipolar modes obey

| Full-vector mode | Scalar $H_z$ | Scalar $E_x$ | Scalar $E_y$ | Symmetry-allowed normal component |
|---|---|---|---|---|
| $p_x$ ($B_1$) | $B_2$ | $A_1$ | $A_2$ | $x$ |
| $p_y$ ($B_2$) | $B_1$ | $A_2$ | $A_1$ | $y$ |

### S3.4 Symmetry-constrained expansion of the radiation matrix elements

The $\Gamma$-point selection rule extends to a symmetry-constrained $k\cdot p$ expansion of the outgoing matrix elements. Let $|\mathcal R_\alpha(\mathbf k)\rangle$ denote the zero-order outgoing channel with in-plane polarization $\alpha\in\{x,y\}$. The radiation coefficients are

$$
c_\alpha(\mathbf k)
=\langle\mathcal R_\alpha(\mathbf k)|\psi_{p_y}(\mathbf k)\rangle.
$$

For an isolated target band in an analytic, symmetry-adapted Bloch gauge,

$$
\mathbf c(g\mathbf k)
=\chi_{B_2}(g)R_g\mathbf c(\mathbf k),
\qquad g\in C_{2v},
\tag{S27a}
$$

where $R_g$ acts on the in-plane polar vector. Since the outgoing $x$ and $y$ channels transform as $B_1$ and $B_2$, respectively, the method of invariants gives

$$
\boxed{
\Gamma(c_x)=B_1\otimes B_2=A_2,
\qquad
\Gamma(c_y)=B_2\otimes B_2=A_1.
}
\tag{S27b}
$$

The lowest-order wave-vector polynomials in the present coordinate convention are

| Irrep | Symmetry-adapted polynomials near $\Gamma$ |
|---|---|
| $A_1$ | $1$, $k_x^2$, $k_y^2$, $\ldots$ |
| $A_2$ | $k_xk_y$, $k_x^3k_y$, $k_xk_y^3$, $\ldots$ |
| $B_1$ | $k_x$, $k_x^3$, $k_xk_y^2$, $\ldots$ |
| $B_2$ | $k_y$, $k_y^3$, $k_x^2k_y$, $\ldots$ |

Equivalently, Eq. (S27a) gives the mirror parities

$$
\begin{aligned}
c_x(-k_x,k_y)&=-c_x(k_x,k_y),
&c_x(k_x,-k_y)&=-c_x(k_x,k_y),\\
c_y(-k_x,k_y)&=+c_y(k_x,k_y),
&c_y(k_x,-k_y)&=+c_y(k_x,k_y).
\end{aligned}
\tag{S27c}
$$

Analyticity around $\Gamma$ and the representations in Eq. (S27b) therefore fix the leading Cartesian expansion:

$$
\boxed{
\begin{aligned}
c_x(k_x,k_y)
&=a_{xy}k_xk_y+\mathcal O(|\mathbf k|^4),\\
c_y(k_x,k_y)
&=c_{y0}+a_{xx}k_x^2+a_{yy}k_y^2
+\mathcal O(|\mathbf k|^4),
\end{aligned}
}
\tag{S27d}
$$

where $c_{y0}=c_y(0,0;\lambda)$ is the symmetry-allowed coefficient in Eq. (S26). The coefficients can be defined without a model-dependent fit as

$$
a_{xy}
=\left.\frac{\partial^2c_x}{\partial k_x\partial k_y}\right|_{\mathbf k=0},
\qquad
a_{xx}
=\left.\frac12\frac{\partial^2c_y}{\partial k_x^2}\right|_{\mathbf k=0},
\qquad
a_{yy}
=\left.\frac12\frac{\partial^2c_y}{\partial k_y^2}\right|_{\mathbf k=0}.
$$

At the tuned periodic-structure zero $\lambda=\lambda_*$ established by Eq. (S27), $c_{y0}=0$. Equation (S27d) then contains no constant or linear radiation term. This result follows from the representation of the radiation matrix elements and does not rely on a topological classification. Section S4 independently recovers the same expansion from the Maxwell interface and unit-cell integrals.

### S3.5 Why $C_2$ alone is insufficient

Under $C_2$ alone, both in-plane polarizations are odd. The two mirror parities are then unavailable. A strict statement such as $c_x(0,0)=0$ for the $p_y$ mode requires the surviving $C_{2v}$ mirrors.

The final geometry and numerical boundary conditions must preserve these mirrors. Small fabrication asymmetries can weakly reopen the nominally forbidden component.

## S4. Normal out-of-plane radiation from interface line integrals

The central result is an exact area-to-interface reduction for the two-dimensional TE model at $\Gamma$. The zero-order electric-field average has no independent contribution from a homogeneous hole interior. It is determined by the magnetic field on the dielectric interfaces and by their oriented normals.

This result turns normal radiation into a boundary-design problem. Translating or deforming a hole perturbs the integration contour, the interface normal, and the self-consistent boundary field. These responses determine the effective radiative dipole, emitted intensity, and radiation zero.

### S4.1 Outgoing coefficient in a photonic-crystal slab

Let

$$
\mathbf E_{\mathbf k}(\boldsymbol\rho,z)
=e^{i\mathbf k\cdot\boldsymbol\rho}
\mathbf u_{E,\mathbf k}(\boldsymbol\rho,z).
$$

On an exterior plane $z=z_0$, define the zero-order tangential coefficient

$$
\boxed{
\mathbf c(\mathbf k;z_0)
=\frac{1}{A_{\mathrm{cell}}}
\int_{\mathrm{cell}}
\mathbf u_{E,\mathbf k,\parallel}(\boldsymbol\rho,z_0)
\,d^2\rho.
}
\tag{S28}
$$

Below the first diffraction threshold, $\mathbf c$ determines the open plane-wave amplitude with in-plane momentum $\mathbf k$. Equation (S28) is the primary full-wave definition used for numerical verification.

### S4.2 Detailed integration by parts at $\Gamma$

Maxwell's equation gives

$$
E_x
=\frac{i}{\omega\varepsilon_0\varepsilon}
\partial_yH_z,
\qquad
E_y
=-\frac{i}{\omega\varepsilon_0\varepsilon}
\partial_xH_z.
\tag{S29}
$$

At $\Gamma$, define

$$
\langle E_\alpha\rangle
=\frac{1}{A_{\mathrm{cell}}}
\int_{\mathrm{cell}}E_\alpha\,d^2\rho.
$$

The following cancellation applies only at $\mathbf k=\mathbf0$. At this point, $H_z$ and $1/\varepsilon$ are periodic over the photonic-crystal unit cell. We show the integration by parts explicitly for $E_y$:

$$
\begin{aligned}
\langle E_y\rangle
&=-\frac{i}{\omega\varepsilon_0A_{\mathrm{cell}}}
\int_{\mathrm{cell}}
\frac{1}{\varepsilon}\partial_xH_z\,d^2\rho\\
&=-\frac{i}{\omega\varepsilon_0A_{\mathrm{cell}}}
\left[
\int_{\mathrm{cell}}
\partial_x\!\left(\frac{H_z}{\varepsilon}\right)d^2\rho
-\int_{\mathrm{cell}}
H_z\partial_x\!\left(\frac{1}{\varepsilon}\right)d^2\rho
\right]\\
&=-\frac{i}{\omega\varepsilon_0A_{\mathrm{cell}}}
\left[
\oint_{\partial\mathrm{cell}}
\frac{H_z}{\varepsilon}N_x\,dl
-\int_{\mathrm{cell}}
H_z\partial_x\!\left(\frac{1}{\varepsilon}\right)d^2\rho
\right]\\
&=\frac{i}{\omega\varepsilon_0A_{\mathrm{cell}}}
\int_{\mathrm{cell}}
H_z\partial_x\!\left(\frac{1}{\varepsilon}\right)d^2\rho.
\end{aligned}
\tag{S30}
$$

Here, $\hat{\mathbf N}$ is the outward normal of the unit-cell boundary. Let $C_+$ and $C_-$ be opposite edges related by a lattice translation $\mathbf R$. The combined contribution is

$$
\int_{C_+}\frac{H_z(\mathbf r)}{\varepsilon(\mathbf r)}N_x\,dl
+\int_{C_-}\frac{H_z(\mathbf r)}{\varepsilon(\mathbf r)}N_x\,dl
=0.
$$

At $\Gamma$, $H_z(\mathbf r+\mathbf R)=H_z(\mathbf r)$ and $\varepsilon(\mathbf r+\mathbf R)=\varepsilon(\mathbf r)$. The paired outward normals have opposite signs. Summing all opposite-edge pairs therefore gives a zero unit-cell boundary term.

The same calculation for $E_x$ gives

$$
\langle E_x\rangle
=-\frac{i}{\omega\varepsilon_0A_{\mathrm{cell}}}
\int_{\mathrm{cell}}
H_z\partial_y\!\left(\frac{1}{\varepsilon}\right)d^2\rho.
\tag{S31}
$$

Equations (S30) and (S31) move the derivative from $H_z$ to $1/\varepsilon$. Its gradient vanishes inside every homogeneous region. Its distributional support lies only on the air--dielectric interfaces.

### S4.3 Exact reduction to an interface line integral

Choose each interface normal $\hat{\mathbf n}$ from dielectric to air. Introduce a signed normal coordinate $u$ and a tangential coordinate $s$. The inverse permittivity changes only along $u$. For a stepwise dielectric profile,

$$
\partial_x\!\left(\frac{1}{\varepsilon}\right)
=\Delta\!\left(\frac{1}{\varepsilon}\right)n_x\delta_S,
\qquad
\partial_y\!\left(\frac{1}{\varepsilon}\right)
=\Delta\!\left(\frac{1}{\varepsilon}\right)n_y\delta_S,
$$

where

$$
\Delta\!\left(\frac{1}{\varepsilon}\right)
=\frac{1}{\varepsilon_{\mathrm{air}}}
-\frac{1}{\varepsilon_{\mathrm{diel}}}.
$$

Equations (S30) and (S31) become

$$
\boxed{
\langle E_x\rangle
=-\frac{i\Delta(1/\varepsilon)}{\omega\varepsilon_0A_{\mathrm{cell}}}
\oint_S H_z n_y\,dl,
}
\tag{S32}
$$

$$
\boxed{
\langle E_y\rangle
=\frac{i\Delta(1/\varepsilon)}{\omega\varepsilon_0A_{\mathrm{cell}}}
\oint_S H_z n_x\,dl.
}
\tag{S33}
$$

The vector identity is

$$
\boxed{
\langle\mathbf E\rangle
=\frac{i\Delta(1/\varepsilon)}{\omega\varepsilon_0A_{\mathrm{cell}}}
\oint_S
H_z(\hat{\mathbf z}\times\hat{\mathbf n})\,dl.
}
\tag{S34}
$$

Here, $S$ is the union of all air--dielectric interfaces in one unit cell. The factors $n_x$ and $n_y$ are Cartesian projections of the interface-normal delta function, rather than phenomenological bond weights.

Equation (S34) is the compact area-to-line conversion. At $\Gamma$, normal radiation is controlled by the interface line integral of $H_z$. Hole interiors influence the result through the interface geometry and the self-consistent boundary field.

### S4.4 Radiation at $\mathbf k=\delta\mathbf k$

The interface-only expression in Eq. (S34) is specific to $\mathbf k=\mathbf0$. At finite wave vector, the complete Bloch field is not periodic over one unit cell. We therefore write

$$
H_{z,\mathbf k}(\boldsymbol\rho)
=e^{i\mathbf k\cdot\boldsymbol\rho}
u_{H,\mathbf k}(\boldsymbol\rho),
$$

where $u_{H,\mathbf k}$ is periodic. The periodic electric-field components satisfy

$$
\begin{aligned}
(u_E)_x
&=\frac{i}{\omega_{\mathbf k}\varepsilon_0\varepsilon}
(\partial_y+ik_y)u_{H,\mathbf k},\\
(u_E)_y
&=-\frac{i}{\omega_{\mathbf k}\varepsilon_0\varepsilon}
(\partial_x+ik_x)u_{H,\mathbf k}.
\end{aligned}
$$

The unit-cell boundary term cancels for the periodic function $u_{H,\mathbf k}/\varepsilon$. Differentiating the Bloch phase nevertheless produces an explicit wave-vector term. The zero-order radiation coefficient becomes

$$
\boxed{
\begin{aligned}
\mathbf c(\mathbf k)
={}&\frac{i\Delta(1/\varepsilon)}
{\omega_{\mathbf k}\varepsilon_0A_{\mathrm{cell}}}
\oint_Su_{H,\mathbf k}
(\hat{\mathbf z}\times\hat{\mathbf n})\,dl\\
&+\frac{\hat{\mathbf z}\times\mathbf k}
{\omega_{\mathbf k}\varepsilon_0A_{\mathrm{cell}}}
\int_{\mathrm{cell}}
\frac{u_{H,\mathbf k}}{\varepsilon}\,d^2\rho.
\end{aligned}
}
\tag{S35}
$$

In Cartesian components,

$$
\begin{aligned}
c_x(\mathbf k)
&=-\frac{i\Delta(1/\varepsilon)}
{\omega_{\mathbf k}\varepsilon_0A_{\mathrm{cell}}}
\oint_Su_{H,\mathbf k}n_y\,dl
-\frac{k_y}{\omega_{\mathbf k}\varepsilon_0A_{\mathrm{cell}}}
\int_{\mathrm{cell}}\frac{u_{H,\mathbf k}}{\varepsilon}\,d^2\rho,\\
c_y(\mathbf k)
&=\frac{i\Delta(1/\varepsilon)}
{\omega_{\mathbf k}\varepsilon_0A_{\mathrm{cell}}}
\oint_Su_{H,\mathbf k}n_x\,dl
+\frac{k_x}{\omega_{\mathbf k}\varepsilon_0A_{\mathrm{cell}}}
\int_{\mathrm{cell}}\frac{u_{H,\mathbf k}}{\varepsilon}\,d^2\rho.
\end{aligned}
$$

The first contribution in each component is an interface term. The second is a unit-cell area integral generated by the Bloch phase. Both contributions are required away from exact $\Gamma$.

For $\mathbf k=\delta\mathbf k$, expand

$$
u_{H,\delta\mathbf k}
=u_{H,0}+\delta u_H+\mathcal O(|\delta\mathbf k|^2),
\qquad
\delta u_H
=\sum_{j=x,y}\delta k_j
\left.\partial_{k_j}u_{H,\mathbf k}\right|_{\mathbf k=0}.
$$

This derivative assumes a continuous Bloch-mode gauge. For a degenerate manifold, the projected $k\cdot p$ Hamiltonian must first define analytic eigenstates. For an isolated reciprocal state at $\Gamma$, $\omega_{\delta\mathbf k}=\omega_0+\mathcal O(|\delta\mathbf k|^2)$. Equation (S35) then gives

$$
\boxed{
\begin{aligned}
\mathbf c(\delta\mathbf k)
={}&\mathbf c(\mathbf0)
+\frac{i\Delta(1/\varepsilon)}
{\omega_0\varepsilon_0A_{\mathrm{cell}}}
\oint_S\delta u_H
(\hat{\mathbf z}\times\hat{\mathbf n})\,dl\\
&+\frac{\hat{\mathbf z}\times\delta\mathbf k}
{\omega_0\varepsilon_0A_{\mathrm{cell}}}
\int_{\mathrm{cell}}
\frac{u_{H,0}}{\varepsilon}\,d^2\rho
+\mathcal O(|\delta\mathbf k|^2).
\end{aligned}
}
\tag{S36}
$$

Equation (S36) is the general first-order result before the band-specific constraints derived in S3 are imposed. At the tuned radiation zero, $\mathbf c(\mathbf0)=\mathbf0$. Section S3 predicts that every linear coefficient vanishes for the target $p_y$ band. The following second-order expansion of the complete Maxwell expression shows how the interface and unit-cell area terms reproduce the symmetry-constrained quadratic coefficients.

### S4.5 Maxwell recovery of the symmetry-constrained expansion

To expand Eq. (S35) without hiding either contribution, define the smooth prefactor, the interface-integral vector, and the unit-cell area integral as

$$
\begin{aligned}
\mathcal K(\mathbf k)
&=\frac{1}{\omega_{\mathbf k}\varepsilon_0A_{\mathrm{cell}}},\\
\mathbf L(\mathbf k)
&=\Delta(1/\varepsilon)
\oint_Su_{H,\mathbf k}
(\hat{\mathbf z}\times\hat{\mathbf n})\,dl,\\
\mathcal A(\mathbf k)
&=\int_{\mathrm{cell}}
\frac{u_{H,\mathbf k}}{\varepsilon}\,d^2\rho.
\end{aligned}
\tag{S36a}
$$

Equation (S35) is then

$$
\mathbf c(\mathbf k)
=\mathcal K(\mathbf k)
\left[
i\mathbf L(\mathbf k)
+(\hat{\mathbf z}\times\mathbf k)\mathcal A(\mathbf k)
\right].
$$

For $i,j\in\{x,y\}$, introduce the derivatives

$$
\mathbf L_i
=\left.\partial_{k_i}\mathbf L\right|_{\mathbf k=0},
\qquad
\mathbf L_{ij}
=\left.\partial_{k_i}\partial_{k_j}\mathbf L\right|_{\mathbf k=0},
\qquad
\mathcal A_i
=\left.\partial_{k_i}\mathcal A\right|_{\mathbf k=0}.
$$

For an isolated reciprocal state, $\mathcal K(\mathbf k)=\mathcal K_0+\mathcal O(|\mathbf k|^2)$. At $\Gamma$, the selection rule gives $L_x(\mathbf0)=0$, while the allowed coefficient in Eq. (S27d) is $c_{y0}=i\mathcal K_0L_y(\mathbf0)$. Tuning $c_{y0}=0$ therefore gives $\mathbf L(\mathbf0)=\mathbf0$. In addition, Eq. (S19) gives $u_{H,0}/\varepsilon\sim B_1$, so its unit-cell average vanishes:

$$
\mathcal A_0
=\int_{\mathrm{cell}}\frac{u_{H,0}}{\varepsilon}\,d^2\rho
=0.
$$

Expanding Eq. (S35) through second order gives

$$
\boxed{
\begin{aligned}
\mathbf c(\mathbf k)
={}&\mathcal K_0
\sum_i k_i
\left[
i\mathbf L_i
+(\hat{\mathbf z}\times\hat{\mathbf e}_i)\mathcal A_0
\right]\\
&+\mathcal K_0
\left[
\frac{i}{2}\sum_{ij}k_ik_j\mathbf L_{ij}
+(\hat{\mathbf z}\times\mathbf k)
\sum_i k_i\mathcal A_i
\right]
+\mathcal O(|\mathbf k|^3).
\end{aligned}
}
\tag{S36b}
$$

The absence of linear terms in Eq. (S27d) requires the first line of Eq. (S36b) to satisfy

$$
\boxed{
i\mathbf L_i
+(\hat{\mathbf z}\times\hat{\mathbf e}_i)\mathcal A_0
=\mathbf0,
\qquad i=x,y.
}
\tag{S36c}
$$

For the target $p_y$ band, $\mathcal A_0=0$, so Eq. (S36c) further gives $\mathbf L_x=\mathbf L_y=\mathbf0$. The complete linear order therefore vanishes term by term. The interface response and the Bloch-phase area contribution first combine nontrivially at second order.

Write $L_{\alpha,ij}=\hat{\mathbf e}_\alpha\cdot\mathbf L_{ij}$. The forbidden quadratic coefficients provide the additional consistency relations

$$
L_{x,xx}=0,
\qquad
\frac{i}{2}L_{x,yy}-\mathcal A_y=0,
\qquad
iL_{y,xy}+\mathcal A_y=0.
$$

The remaining second-order terms give the following result. The mirror parities in Eq. (S27c) also forbid all cubic terms, so the remainder begins at fourth order:

$$
\boxed{
\begin{aligned}
c_x(k_x,k_y)
&=a_{xy}k_xk_y+\mathcal O(|\mathbf k|^4),\\
c_y(k_x,k_y)
&=a_{xx}k_x^2+a_{yy}k_y^2
+\mathcal O(|\mathbf k|^4),
\end{aligned}
}
\tag{S36d}
$$

with the microscopic coefficient identities

$$
\boxed{
\begin{aligned}
a_{xy}
&=\mathcal K_0\left(iL_{x,xy}-\mathcal A_x\right),\\
a_{xx}
&=\mathcal K_0\left(\frac{i}{2}L_{y,xx}+\mathcal A_x\right),\\
a_{yy}
&=\mathcal K_0\frac{i}{2}L_{y,yy}.
\end{aligned}
}
$$

Equation (S36d) is identical to Eq. (S27d) after setting $c_{y0}=0$. The S3 result follows from representations and analyticity; the S4 result follows from a second-order expansion of the full Maxwell radiation integral. The coefficient identities show that the same $a_{xy}$, $a_{xx}$, and $a_{yy}$ contain both interface-field response and finite-wave-vector unit-cell contributions.

Let the smooth outgoing-channel normalization be $\mathcal N_I(\mathbf k)=\mathcal N_I(0)+\mathcal O(|\mathbf k|^2)$, with $\mathcal N_I(0)>0$. The mutually verified radiation vector gives the anisotropic Cartesian intensity

$$
\boxed{
\begin{aligned}
I(k_x,k_y)
={}&\mathcal N_I(0)
\left[
|a_{xy}|^2k_x^2k_y^2
+|a_{xx}k_x^2+a_{yy}k_y^2|^2
\right]
+\mathcal O(|\mathbf k|^6)\\
={}&\mathcal N_I(0)
\left[
|a_{xx}|^2k_x^4
+|a_{yy}|^2k_y^4
+\left(|a_{xy}|^2+2\operatorname{Re}[a_{xx}a_{yy}^*]\right)k_x^2k_y^2
\right]
+\mathcal O(|\mathbf k|^6).
\end{aligned}
}
\tag{S36e}
$$

In particular,

$$
I(k_x,0)=\mathcal N_I(0)|a_{xx}|^2k_x^4+\mathcal O(k_x^6),
\qquad
I(0,k_y)=\mathcal N_I(0)|a_{yy}|^2k_y^4+\mathcal O(k_y^6).
$$

The quartic order is fixed after tuning, while the unequal Cartesian coefficients retain the $C_{2v}$ anisotropy.

### S4.6 Boundary first moments and the effective radiative dipole

For the $n$th triangular hole, define

$$
\boxed{
\mathbf q_n
=\oint_{\partial h_n}H_z(\mathbf r)\hat{\mathbf n}(\mathbf r)\,dl.
}
\tag{S37}
$$

Equation (S34) becomes

$$
\langle\mathbf E\rangle
=
\frac{i\Delta(1/\varepsilon)}{\omega\varepsilon_0A_{\mathrm{cell}}}
\hat{\mathbf z}\times
\sum_{n=1}^{6}\mathbf q_n.
\tag{S38}
$$

For a triangular hole with edges $e_{n,m}$, this moment must be evaluated edge by edge:

$$
\mathbf q_n
=\sum_{m=1}^{3}
\hat{\mathbf n}_{n,m}
\int_{e_{n,m}}H_z\,dl.
$$

Thus, the boundary moment retains both the signed field on each edge and the orientation of that edge. Replacing a triangle by its perimeter or its center value discards the information required for radiation.

If $H_z$ is constant on the complete boundary of one hole, then

$$
\mathbf q_n
=H_{z,n}\oint_{\partial h_n}\hat{\mathbf n}\,dl
=\mathbf 0.
\tag{S39}
$$

The six site amplitudes alone cannot determine a nonzero radiation amplitude. The edge-to-edge variation of $H_z$ supplies the required local boundary form factor.

Below the first diffraction threshold, $\langle\mathbf E\rangle$ is proportional to the effective in-plane radiative electric-dipole amplitude per cell. The proportionality contains the vertical overlap and the outgoing-channel normalization. We may therefore write

$$
\mathbf p_{\parallel}^{\mathrm{rad}}
=\mathcal N_p\langle\mathbf E\rangle
=\mathcal C_p\hat{\mathbf z}\times
\sum_{n=1}^{6}\mathbf q_n,
$$

and

$$
I_{\perp}
=\mathcal N_I\left(|c_x|^2+|c_y|^2\right),
\qquad
\mathcal N_I>0.
$$

The constants $\mathcal N_p$, $\mathcal C_p$, and $\mathcal N_I$ depend on channel normalization. The position of a radiation zero does not depend on these nonzero factors.

For the target mode, define $r_y$ as a shorthand for the physical $\Gamma$-point coefficient $c_y(\Gamma)=c_{y0}$ introduced in S3. Then

$$
r_y\equiv c_y(\Gamma)=c_{y0}
\propto
\langle E_y\rangle
\propto
\sum_{n=1}^{6}
\oint_{\partial h_n}H_z n_x\,dl.
\tag{S40}
$$

The $C_{2v}$ selection rule gives $c_x(0,0)=0$ for the target $p_y$ state. The dipolar singularity therefore occurs when the complex, signed edge contributions in Eq. (S40) cancel exactly.

### S4.7 Connection to the tight-binding eigenvector

Expand the magnetic field as

$$
H_z(\mathbf r)
\simeq
\sum_{n=1}^{6}a_n\phi_n(\mathbf r;\eta,\zeta).
$$

Equation (S40) can then be written as

$$
\boxed{
r_y(\eta,\zeta)
=\mathbf g_y^\dagger(\eta,\zeta)
\mathbf a_{p_y}(\eta,\zeta).
}
\tag{S41}
$$

The vector $\mathbf a_{p_y}$ is determined by the Hermitian tight-binding problem. The vector $\mathbf g_y$ contains the interface moments and the outgoing-mode normalization. Its site-resolved entries have the form

$$
(g_y)_n
\propto
\sum_{j=1}^{6}
\oint_{\partial h_j}
\phi_n(\mathbf r;\eta,\zeta)n_x\,dl
\simeq
\oint_{\partial h_n}
\phi_n n_x\,dl.
$$

The final approximation uses localized orbitals and neglects cross-site boundary overlaps. Thus, $\mathbf a_{p_y}$ specifies the six-site standing wave. The vector $\mathbf g_y$ specifies how the boundaries convert that standing wave into normal radiation.

For either geometry parameter $\lambda\in\{\eta,\zeta\}$,

$$
\boxed{
\frac{\partial r_y}{\partial\lambda}
=
\left(\frac{\partial\mathbf g_y^\dagger}{\partial\lambda}\right)
\mathbf a_{p_y}
+
\mathbf g_y^\dagger
\frac{\partial\mathbf a_{p_y}}{\partial\lambda}.
}
\tag{S42}
$$

The first term describes direct changes of the boundary radiation operator. The second describes redistribution of the Hermitian eigenvector. This identity explains why the band Hamiltonian alone cannot predict the radiation zero.

### S4.8 Boundary perturbation and construction of the radiation zero

Changing $\eta$ translates the six interfaces and changes the boundary field. Changing $\zeta$ deforms the interfaces while preserving the total hole area. Both perturbations modify $\mathbf g_y$ and $\mathbf a_{p_y}$ through Eq. (S42).

For either parameter $\lambda$, the first variation of one triangular boundary moment is

$$
\delta_\lambda\mathbf q_n
=\delta\lambda\sum_{m=1}^{3}
\left[
(\partial_\lambda\hat{\mathbf n}_{n,m})
\int_{e_{n,m}}H_z\,dl
+\hat{\mathbf n}_{n,m}
\partial_\lambda
\left(\int_{e_{n,m}}H_z\,dl\right)
\right].
$$

The second derivative acts on the moving edge, its length, and the self-consistent field. The present translations and size changes preserve each edge orientation, so $\partial_\lambda\hat{\mathbf n}_{n,m}=0$. The first term is retained for general boundary deformations. This expression makes the boundary character of the geometric perturbation explicit.

The area-preserving size response has the first-order weights

$$
\gamma_n=
\begin{cases}
1,&n=1,4,\\
-\dfrac12,&n=2,3,5,6.
\end{cases}
$$

The factor $-1/2$ follows from $\left.\partial_\zeta\sqrt{(3-\zeta^2)/2}\right|_{\zeta=1}=-1/2$ for the compensating holes in Eq. (S5). At the undeformed folded-lattice reference, the baseline moments satisfy $\sum_n\mathbf q_n^{(0)}=\mathbf0$. Individual baseline moments need not vanish. A symmetry-reduced boundary perturbation is

$$
\boxed{
\mathbf q_n
=\mathbf q_n^{(0)}
+a_n\left[
\chi_\eta(\eta-1)
+\chi_\zeta\gamma_n(\zeta-1)
\right]\hat{\mathbf r}_n
+\mathcal O(\delta^2).
}
\tag{S43}
$$

Here, $\hat{\mathbf r}_n=\mathbf R_n/R$, and $\delta$ collectively denotes the two geometric perturbations. The amplitudes $a_n$ are evaluated at the reference geometry to this order. The coefficients $\chi_\eta$ and $\chi_\zeta$ are derivatives of the complete three-edge boundary operator. Each coefficient includes contour motion, normal rotation, and changes of the local basis-field profile. Neither coefficient is determined by the hole perimeter alone. Eigenvector redistribution remains the second term in Eq. (S42).

For the $p_y$ vector in Eq. (S7), direct summation gives

$$
\sum_{n=1}^{6}(p_y)_n\hat{\mathbf r}_n
=-\sqrt3\,\hat{\mathbf x},
\qquad
\sum_{n=1}^{6}(p_y)_n\gamma_n\hat{\mathbf r}_n
=-\frac{\sqrt3}{2}\hat{\mathbf x}.
$$

Within this fixed-eigenvector common-response model, Eqs. (S38) and (S43) yield the allowed scalar radiation amplitude

$$
r_y
=\mathcal C_y
\left[
\chi_\eta(\eta-1)
+\frac{\chi_\zeta}{2}(\zeta-1)
\right]
+\mathcal O(\delta^2),
\qquad
c_x(0,0)=0,
$$

where $\mathcal C_y$ is a nonzero channel-normalization factor. This expression shows how the two boundary perturbations can cancel in the allowed $y$-polarized channel.

Around an arbitrary reference geometry $(\eta_0,\zeta_0)$, the general first-order response is

$$
r_y(\eta,\zeta)
=r_{y0}
+B_\eta(\eta-\eta_0)
+B_\zeta(\zeta-\zeta_0)
+\mathcal O(\delta^2),
$$

with $B_\lambda=\partial_\lambda r_y$ from Eq. (S42). At fixed $\eta$, the predicted zero is

$$
\zeta^*(\eta)
\simeq
\zeta_0-
\frac{r_y(\eta,\zeta_0)}
{\left.\partial_\zeta r_y\right|_{(\eta,\zeta_0)}}.
\tag{S44}
$$

Along a mirror-preserving lossless tuning path, a continuous eigenmode gauge gives a signed scalar $r_y$. The complex phase must remain collinear before Eq. (S44) is applied to full-wave data. The intensity and radiative decay then satisfy

$$
I_\perp=\mathcal N_I|r_y|^2,
\qquad
\gamma_{\mathrm{rad}}=\mathcal N_\gamma|r_y|^2,
\qquad
Q_{\mathrm{rad}}^{-1}=\mathcal N_Q|r_y|^2,
$$

with positive normalization factors. Near a simple zero, $I_\perp\propto|\partial_\zeta r_y|^2(\zeta-\zeta^*)^2$. The singularity must therefore be located from the complex signed amplitude, rather than from an intensity minimum alone.

In practice, the construction starts from the eighteen edge contributions in Eq. (S40). Finite differences of basis-resolved edge integrals determine $\partial_\lambda\mathbf g_y$. Projections of the full-wave modes determine $\partial_\lambda\mathbf a_{p_y}$. Without this decomposition, finite differences of the total edge sum determine $B_\lambda$ directly. Equation (S44) then locates the signed zero crossing.

An ideal common-response model may predict opposite $\zeta$ shifts for the $p_x$ and $p_y$ zeros. This sign relation is a model prediction. It is not a group-theory identity.

%%
NUMERICAL INSERTS REQUIRED:
- complex c_x and c_y at Gamma versus zeta for fixed eta;
- continuous phase convention and signed zero crossing;
- comparison of the exterior-plane Fourier coefficient with the interface integral;
- edge-resolved contributions for the three sides of all six holes;
- Gamma-point band gap at the fitted radiation zero;
- uncertainty or mesh-convergence estimate for zeta*.
%%

### S4.9 Open-system summary

The closed and open descriptions can be combined as

$$
H_{\mathrm{eff}}
=H_{\mathrm{Herm}}
+\Delta_{\mathrm{rad}}
-\frac{i}{2}D^\dagger D.
\tag{S45}
$$

At the dipolar singularity,

$$
D(\Gamma)|p_y\rangle=\mathbf0,
$$

while

$$
m_y\neq0.
$$

The anti-Hermitian loss of the target state vanishes at exact $\Gamma$. The Hermitian band separation remains finite. These two statements define different operator conditions.

### S4.10 Scope of the line-integral model

Equations (S30)--(S35) are exact within the two-dimensional TE model with a stepwise dielectric profile. Equation (S34) is the interface-only identity at $\Gamma$. Equation (S35) is its finite-$\mathbf k$ extension and contains the additional unit-cell area term. Equation (S36) is the corresponding first-order expansion for an isolated reciprocal state. Equations (S36a)--(S36e) expand the complete Maxwell expression through second order, establish the absence of all linear terms for the target state, and recover the symmetry-constrained radiation vector and its anisotropic intensity.

Equations (S37)--(S40) define the exact $\Gamma$-point boundary moments within the TE model. Equation (S41) is exact for a complete basis and approximate after the six-orbital truncation. Equation (S42) is an exact derivative identity within the reduced model. Equation (S43) is a symmetry-reduced first-order boundary response.

The actual slab is three dimensional and TE-like. Its outgoing coefficient must be evaluated from Eq. (S28) or an equivalent scattering projection.

The line-integral model provides the exact-$\Gamma$ in-plane interference mechanism. At finite wave vector, the explicit area term in Eq. (S35) must also be retained. The proportionality to the slab radiation amplitude includes a vertical-overlap factor. Full-wave agreement is therefore required before quantitative use.

## S5. Brillouin-zone wave-vector weighting of the finite cavity

The finite cavity does not possess exact in-plane translation symmetry, so its eigenmode is not labeled by a single Bloch wave vector. A Bloch basis nevertheless provides a complete representation of the field on a selected collection of complete unit cells. The decomposition proceeds in two distinct stages. A finite-lattice Fourier transform first resolves the cell-index dependence into discrete wave vectors. The periodic cell content at each wave vector is then projected onto the Bloch eigenfields of an infinite reference crystal. The resulting complex coefficients give the wave-vector- and band-resolved composition of the finite mode.

### S5.1 Finite-mode data on a lattice of complete unit cells

Let $\boldsymbol\Psi_f(\mathbf r)$ denote a selected complex field of the finite eigenmode, such as the magnetic field used in the master equation. Choose a bulk-like analysis region containing $N$ complete unit cells with Bravais origins

$$
S=\{\mathbf R_1,\ldots,\mathbf R_N\}\subset\Lambda.
$$

Writing $\mathbf r=(\mathbf R+\boldsymbol\rho,z)$, define the field attached to each cell as

$$
\boldsymbol\Psi_{\mathbf R}(\boldsymbol\rho,z)
=\boldsymbol\Psi_f(\mathbf R+\boldsymbol\rho,z),
\qquad \mathbf R\in S.
$$

The Fourier transform below acts only on the discrete cell label $\mathbf R$. It retains the complete intra-cell coordinate $\boldsymbol\rho$, vertical coordinate $z$, vector components, and complex phase. A positive cell metric $B$ defines

$$
\langle\mathbf f,\mathbf g\rangle_B
=\int_{\Omega_{\mathrm{cell}}}
\mathbf f^\dagger B\mathbf g\,d^2\rho\,dz,
\qquad
\|\mathbf f\|_B^2=\langle\mathbf f,\mathbf f\rangle_B.
$$

The metric must match the field representation used for the reference Bloch modes. Examples include the ordinary magnetic-field inner product, an $\varepsilon$-weighted electric-field inner product, or the corresponding finite-element mass matrix.

### S5.2 Exact finite-lattice Fourier basis

An exact orthogonal plane-wave basis exists when $S$ is a complete representative set of a finite lattice quotient $\Lambda/\Lambda_{\mathrm{sup}}$, where $\Lambda_{\mathrm{sup}}$ is a finite-index superlattice. The dual finite wave-vector set is

$$
K=\Lambda_{\mathrm{sup}}^*/\Lambda^*,
\qquad |K|=|S|=N.
$$

The characters of this finite Abelian group obey

$$
\boxed{
\frac{1}{N}
\sum_{\mathbf R\in S}
e^{i(\mathbf k-\mathbf k')\cdot\mathbf R}
=\delta_{\mathbf k,\mathbf k'},
\qquad \mathbf k,\mathbf k'\in K.
}
\tag{S46}
$$

For a block $\mathbf R_{mn}=m\mathbf a_1+n\mathbf a_2$ with $0\le m<M_1$ and $0\le n<M_2$, the corresponding wave vectors are

$$
\mathbf k_{\ell_1\ell_2}
=\frac{\ell_1}{M_1}\mathbf b_1
+\frac{\ell_2}{M_2}\mathbf b_2.
$$

The construction is not restricted to parallelograms. A centered hexagonal region on a triangular lattice,

$$
S_L=\{m\mathbf a_1+n\mathbf a_2:
|m|\le L,\ |n|\le L,\ |m+n|\le L\},
$$

contains $N=3L(L+1)+1$ cells and forms an exact quotient representative set for the superlattice generated by

$$
\mathbf T_1=(L+1)\mathbf a_1+L\mathbf a_2,
\qquad
\mathbf T_2=-L\mathbf a_1+(2L+1)\mathbf a_2.
$$

The associated $N$ wave vectors satisfy $e^{i\mathbf k\cdot\mathbf T_1}=e^{i\mathbf k\cdot\mathbf T_2}=1$, modulo the reciprocal lattice. The exactness follows from the quotient property, not from the visual shape of the region.

If the selected cells form an irregular mask, contain missing cells, or do not represent a finite quotient, Eq. (S46) does not hold. The plane-wave dictionary can still be used, but the coefficients must then be obtained from its Gram matrix or from a least-squares fit. A conventional Fourier transform of an arbitrary mask is a windowed spectrum and contains wave-vector leakage.

### S5.3 Extraction of the wave-vector-resolved periodic cell content

This manuscript uses the Bloch convention of S1,

$$
\boldsymbol\Psi_{b\mathbf k}(\mathbf r)
=e^{i\mathbf k\cdot\mathbf r}
\mathbf u_{b\mathbf k}(\boldsymbol\rho,z).
$$

The unitary cell-index Fourier pair consistent with this convention is

$$
\boxed{
\begin{aligned}
\mathbf F(\mathbf k,\boldsymbol\rho,z)
&=\frac{1}{\sqrt N}
\sum_{\mathbf R\in S}
e^{-i\mathbf k\cdot\mathbf R}
\boldsymbol\Psi_{\mathbf R}(\boldsymbol\rho,z),\\
\boldsymbol\Psi_{\mathbf R}(\boldsymbol\rho,z)
&=\frac{1}{\sqrt N}
\sum_{\mathbf k\in K}
e^{i\mathbf k\cdot\mathbf R}
\mathbf F(\mathbf k,\boldsymbol\rho,z).
\end{aligned}
}
\tag{S47}
$$

The function $\mathbf F$ still contains the intra-cell Bloch phase. Removing it defines the periodic cell content

$$
\mathbf G(\mathbf k,\boldsymbol\rho,z)
=e^{-i\mathbf k\cdot\boldsymbol\rho}
\mathbf F(\mathbf k,\boldsymbol\rho,z).
$$

The finite mode is then reconstructed exactly on the selected cells as

$$
\boxed{
\boldsymbol\Psi_f(\mathbf R+\boldsymbol\rho,z)
=\frac{1}{\sqrt N}
\sum_{\mathbf k\in K}
e^{i\mathbf k\cdot(\mathbf R+\boldsymbol\rho)}
\mathbf G(\mathbf k,\boldsymbol\rho,z).
}
\tag{S48}
$$

Equation (S48) is already an exact Bloch-component decomposition of the finite field on $S$. It does not assume a slowly varying envelope and does not require prior knowledge of any band.

The relation to an ordinary spatial Fourier spectrum follows by expanding the periodic function in reciprocal vectors:

$$
\mathbf G(\mathbf k,\boldsymbol\rho,z)
=\sum_{\mathbf G_r}
\mathbf G_{\mathbf G_r}(\mathbf k,z)
e^{i\mathbf G_r\cdot\boldsymbol\rho}.
$$

Substitution into Eq. (S48) gives spatial harmonics at

$$
\boxed{
\mathbf q=\mathbf k+\mathbf G_r.
}
\tag{S49}
$$

Consequently, Fourier peaks separated by reciprocal lattice vectors are replicas of the same reduced Bloch wave vector. Folding the global Fourier spectrum into the first Brillouin zone identifies $\mathbf k$, while $\mathbf G(\mathbf k,\boldsymbol\rho,z)$ retains the full intra-cell field profile that a simple Fourier-intensity map discards.

Parseval's identity gives

$$
\boxed{
\sum_{\mathbf R\in S}
\|\boldsymbol\Psi_{\mathbf R}\|_B^2
=\sum_{\mathbf k\in K}
\|\mathbf G(\mathbf k)\|_B^2.
}
\tag{S50}
$$

The total wave-vector weight is therefore

$$
W(\mathbf k)=\|\mathbf G(\mathbf k)\|_B^2,
\qquad
\eta(\mathbf k)=
\frac{W(\mathbf k)}{\sum_{\mathbf k'\in K}W(\mathbf k')}.
$$

The map $W(\mathbf k)$ identifies where the finite mode resides in the Brillouin zone. It does not by itself identify which Bloch band contributes at a given wave vector.

If fields are exported using the COMSOL convention $e^{-i\mathbf k\cdot\mathbf r}$, the preprocessing transform uses the opposite sign and the reported wave vector is mapped by $\mathbf k\rightarrow-\mathbf k$ relative to Eqs. (S47)--(S49). The same convention must be used for the finite field and the reference Bloch modes.

### S5.4 Projection onto infinite-crystal Bloch bands

Let $\mathbf u_{b\mathbf k}^{\infty}(\boldsymbol\rho,z)$ be the periodic profile of band $b$ in an infinite reference crystal at the same reduced wave vector $\mathbf k$. For a complete $B$-orthonormal cell eigenbasis,

$$
\boxed{
\mathbf G(\mathbf k,\boldsymbol\rho,z)
=\sum_b
A_b(\mathbf k)
\mathbf u_{b\mathbf k}^{\infty}(\boldsymbol\rho,z),
\qquad
A_b(\mathbf k)
=\left\langle
\mathbf u_{b\mathbf k}^{\infty},
\mathbf G(\mathbf k)
\right\rangle_B.
}
\tag{S51}
$$

Combining Eqs. (S48) and (S51) yields the finite-mode Bloch-band expansion

$$
\boxed{
\boldsymbol\Psi_f(\mathbf r)
=\frac{1}{\sqrt N}
\sum_{\mathbf k\in K}
\sum_b
A_b(\mathbf k)
e^{i\mathbf k\cdot\mathbf r}
\mathbf u_{b\mathbf k}^{\infty}(\boldsymbol\rho,z).
}
\tag{S52}
$$

With all bands retained, Eq. (S52) is complete on the selected cells. It is a representation statement: the finite, terminated, or defective structure is not being identified with an infinite-crystal eigenproblem. The coefficients only specify how the finite field decomposes in the chosen complete basis.

The band-resolved weights are

$$
W_b(\mathbf k)=|A_b(\mathbf k)|^2,
\qquad
\eta_b(\mathbf k)
=\frac{|A_b(\mathbf k)|^2}
{\sum_{\mathbf k',b'}|A_{b'}(\mathbf k')|^2}.
$$

If the retained Bloch profiles are not exactly $B$-orthonormal, define the overlap matrix $S_{bb'}(\mathbf k)=\langle\mathbf u_{b\mathbf k}^{\infty},\mathbf u_{b'\mathbf k}^{\infty}\rangle_B$ and the data vector $d_b(\mathbf k)=\langle\mathbf u_{b\mathbf k}^{\infty},\mathbf G(\mathbf k)\rangle_B$. The coefficients follow from $\mathbf S(\mathbf k)\mathbf A(\mathbf k)=\mathbf d(\mathbf k)$ rather than from independent scalar projections.

### S5.5 Complex phases, gauge, and reconstruction tests

The coefficient $A_b(\mathbf k)$ is complex. Its phase is essential for field reconstruction and for any later coherent radiation calculation. A Bloch-gauge transformation acts as

$$
\boxed{
\mathbf u_{b\mathbf k}^{\infty}
\rightarrow e^{i\chi_b(\mathbf k)}
\mathbf u_{b\mathbf k}^{\infty},
\qquad
A_b(\mathbf k)
\rightarrow e^{-i\chi_b(\mathbf k)}A_b(\mathbf k).
}
\tag{S53}
$$

Only the product $A_b\mathbf u_{b\mathbf k}^{\infty}$ is gauge invariant. A continuous gauge may be fixed by parallel transport from $\Gamma$ or by requiring the overlap with a reference profile to be real and positive. Intensity-only data cannot determine these complex coefficients.

For a retained set $\mathcal K_{\mathrm{keep}}$ and band window $\mathcal B_{\mathrm{keep}}$, reconstruct

$$
\boldsymbol\Psi_{\mathrm{rec}}(\mathbf r)
=\frac{1}{\sqrt N}
\sum_{\mathbf k\in\mathcal K_{\mathrm{keep}}}
\sum_{b\in\mathcal B_{\mathrm{keep}}}
A_b(\mathbf k)e^{i\mathbf k\cdot\mathbf r}
\mathbf u_{b\mathbf k}^{\infty},
$$

and report

$$
\epsilon_{\mathrm{rec}}
=\frac{\|\boldsymbol\Psi_f-\boldsymbol\Psi_{\mathrm{rec}}\|_B}
{\|\boldsymbol\Psi_f\|_B}.
$$

A large residual can indicate an insufficient band window, inconsistent normalization, boundary-cell contamination, or a mismatch between the finite structure and the selected reference crystal.

### S5.6 Incomplete Bloch dictionaries and data-driven cell profiles

If the infinite-crystal band dictionary is unavailable, the exact object supplied by the finite Fourier transform remains $\mathbf G(\mathbf k,\boldsymbol\rho,z)$. A normalized profile

$$
\mathbf p_{\mathbf k}
=\frac{\mathbf G(\mathbf k)}{\sqrt{W(\mathbf k)}}
$$

can be compared across wave vectors through the gauge-invariant similarity

$$
\mathcal S(\mathbf k,\mathbf k')
=\left|
\langle\mathbf p_{\mathbf k},\mathbf p_{\mathbf k'}\rangle_B
\right|^2.
$$

Wave vectors with similar periodic profiles may be grouped and represented by a weighted singular-value decomposition,

$$
\mathbf G(\mathbf k,\boldsymbol\rho,z)
\simeq
\sum_{s=1}^{r}
A_s(\mathbf k)\mathbf u_s(\boldsymbol\rho,z).
$$

A rank-one group identifies one common unit-cell profile with a wave-vector-dependent complex coefficient. A higher rank indicates a varying mixture of cell profiles. Without infinite-crystal reference modes, these data-driven profiles should not be assigned unique band labels. Reference Bloch modes, symmetry labels, or several finite eigenmodes analyzed jointly are required to resolve multiple bands contributing at the same $\mathbf k$.

### S5.7 Envelope limit and connection to the radiation vector

The familiar carrier-envelope form is a controlled compression of Eq. (S52), rather than its starting assumption. If one target band dominates, its periodic profile varies weakly over a narrow set of wave vectors around $\mathbf k_0$, and a continuous gauge has been fixed, then

$$
\mathbf G(\mathbf k,\boldsymbol\rho,z)
\simeq A(\mathbf k)\mathbf u_0(\boldsymbol\rho,z).
$$

Equation (S48) reduces to

$$
\boxed{
\boldsymbol\Psi_f(\mathbf r)
\simeq
e^{i\mathbf k_0\cdot\mathbf r}
\mathbf u_0(\boldsymbol\rho,z)
\psi(\mathbf r),
\qquad
\psi(\mathbf r)
=\frac{1}{\sqrt N}
\sum_{\mathbf k}
A(\mathbf k)e^{i(\mathbf k-\mathbf k_0)\cdot\mathbf r}.
}
\tag{S54}
$$

The envelope is slowly varying only when the relative spectrum $\mathbf k-\mathbf k_0$ is narrow compared with the reciprocal-lattice scale. A spatially localized mode therefore occupies a finite region of the Brillouin zone, with a characteristic width that decreases as the cavity size increases.

Let $\mathbf c_b(\mathbf k)=(c_{b,x},c_{b,y})$ denote the outgoing radiation vector of the infinite-crystal Bloch state, as derived in S3 and S4. In a translation-invariant reference channel, the radiation carried by one reduced wave vector is

$$
\boxed{
\mathbf d_f(\mathbf k)
=\sum_b A_b(\mathbf k)\mathbf c_b(\mathbf k).
}
\tag{S55}
$$

The gauge phases in Eq. (S53) cancel between $A_b$ and $\mathbf c_b$, so $\mathbf d_f$ is gauge invariant. Equation (S55) also shows why the finite-mode composition must be extracted as complex amplitudes rather than as $W(\mathbf k)$ alone.

Different reduced wave vectors radiate into different in-plane far-field momenta in the translation-invariant reference problem. A coherent sum over distinct $\mathbf k$ values at one observed momentum requires an additional finite-boundary or aperture mixing kernel. The general form is

$$
\mathbf E_{\mathrm{far}}(\mathbf K)
=\sum_{\mathbf k,b}
\mathcal M(\mathbf K,\mathbf k)
A_b(\mathbf k)\mathbf c_b(\mathbf k),
$$

where $\mathcal M$ must be obtained from the actual cavity truncation, interface scattering, or a direct Fourier transform of the complete finite field. It must not be replaced by an assumed convolution without validation. The decomposition in Eqs. (S47)--(S53) supplies the finite mode's complex Bloch content required for that subsequent radiation calculation.

## Appendix A. Four-band reduction and $\Gamma$-point parameter extraction

### A.1 Löwdin reduction from the six-site model

The $s$- and $f$-like states remain spectrally separated from the four $p/d$ states near $\Gamma$. Define

$$
U_{pd}
=\big(
|d_{xy}\rangle,
|d_{x^2-y^2}\rangle,
|p_y\rangle,
|p_x\rangle
\big),
\qquad
P=U_{pd}U_{pd}^\dagger,
\qquad
Q=I_6-P.
$$

The exact energy-dependent Löwdin downfolding is

$$
\boxed{
H_{\mathrm{eff}}(E)
=U_{pd}^\dagger H_{\mathrm{site}}U_{pd}
+U_{pd}^\dagger H_{\mathrm{site}}Q
\left(E-QH_{\mathrm{site}}Q\right)^{-1}
QH_{\mathrm{site}}U_{pd}.
}
\tag{A1}
$$

The inverse in Eq. (A1) is restricted to the $Q$ subspace. For the $C_{6v}$ structure at $\Gamma$, symmetry separates the $p/d$ and $s/f$ subspaces, so direct projection is exact. After the $C_{2v}$ modulation, $s$ and $d_{x^2-y^2}$ may share one irreducible representation, as may $f$ and $p_y$. The leading off-block elements of the six-site Hamiltonian are

$$
\begin{aligned}
\kappa_s
&=\langle s|H_{\mathrm{site}}|d_{x^2-y^2}\rangle
=\frac{\sqrt2}{3}
\left(\alpha t_1+\beta t_0-\gamma t_0-\gamma t_1-3\mu\right),\\
\kappa_f
&=\langle f|H_{\mathrm{site}}|p_y\rangle
=\frac{\sqrt2}{3}
\left(-\alpha t_1-\beta t_0+\gamma t_0+\gamma t_1-3\mu\right).
\end{aligned}
$$

Both matrix elements vanish when the structure returns to $C_{6v}$ but are generally nonzero under $C_{2v}$. Their influence enters through the self-energy in Eq. (A1). If the remote-state separation $\Delta_{sf}$ exceeds these couplings, the static corrections scale as $|\kappa_{s,f}|^2/\Delta_{sf}$.

Let $\delta_{2v}$ denote the off-block modulation scale. Retaining first order in $\mathbf k$ and first order in $\delta_{2v}$ gives

$$
H_{pd}(\mathbf k)
=U_{pd}^\dagger H_{\mathrm{site}}(\mathbf k)U_{pd}
+\mathcal O\!\left(
\frac{(vk)^2}{\Delta_{sf}},
\frac{\delta_{2v}vk}{\Delta_{sf}},
\frac{\delta_{2v}^2}{\Delta_{sf}}
\right).
$$

Spectral isolation of the $p/d$ manifold is the basis of the $6\times6$ to $4\times4$ reduction. Equation (A1) must be retained when the remote-state correction is not negligible.

### A.2 $\Gamma$-point eigenfrequencies and effective parameters

At $\Gamma$, Eq. (S17) gives

$$
\boxed{
\begin{aligned}
\omega_{d_{x^2-y^2}}&=m_y+\mu,\\
\omega_{p_y}&=-m_y+\mu,\\
\omega_{d_{xy}}&=m_x-\mu,\\
\omega_{p_x}&=-m_x-\mu.
\end{aligned}
}
\tag{A2}
$$

The onsite response shifts the centers of the two sectors by $\pm\mu$ and cancels from each $p/d$ separation. The hopping anisotropy controls $m_x$ and $m_y$. The parameters reduce as follows:

$$
\boxed{
\begin{aligned}
\zeta=1:\quad
&\alpha=\beta=\gamma=1,
\quad \mu=0,
\quad m_x=m_y=t_0-t_1,\\
(\eta,\zeta)=(1,1):\quad
&t_0=t_1,
\quad m_x=m_y=0.
\end{aligned}
}
\tag{A3}
$$

The three independent effective parameters can be extracted directly from full-wave $\Gamma$-point frequencies:

$$
\boxed{
\begin{aligned}
m_y&=\frac{\omega_{d_{x^2-y^2}}-\omega_{p_y}}{2},\\
m_x&=\frac{\omega_{d_{xy}}-\omega_{p_x}}{2},\\
\mu&=\frac{\omega_{d_{x^2-y^2}}+\omega_{p_y}
-\omega_{d_{xy}}-\omega_{p_x}}{4}.
\end{aligned}
}
\tag{A4}
$$

When fitted to full-wave frequencies, Eq. (A4) defines renormalized parameters that include remote-state shifts. The factors $\alpha$, $\beta$, and $\gamma$ require additional near-$\Gamma$ velocity information or microscopic overlap calculations. The factorized form in Eq. (S6) remains a nearest-neighbour approximation.

For the target sector,

$$
\boxed{
\Delta_y
=\left|\omega_{d_{x^2-y^2}}-\omega_{p_y}\right|
=2|m_y|.
}
\tag{A5}
$$

Equation (A5) is the frequency-domain verification of the mass definition in Eqs. (S13) and (S18).

## Appendix B. Localized-basis projection and the $k\cdot p$ expansion

Let $\{\phi_n\}_{n=1}^{6}$ be six localized Wannier-like functions in one enlarged unit cell, and let

$$
\mathcal B_{\mathrm{loc}}
=\{\phi_{1\mathbf k},\phi_{2\mathbf k},\ldots,\phi_{6\mathbf k}\}.
$$

The lattice-periodic field can be expanded as

$$
u_{H,\mathbf k}(\mathbf r)
\simeq
\sum_{n=1}^{6}a_n(\mathbf k)\phi_{n\mathbf k}(\mathbf r).
$$

With $\mathcal B_{\mathrm{loc}}$ as the ordered basis, projection of Eq. (S2) gives the generalized matrix eigenproblem

$$
\boxed{
H_{\mathrm{TB}}(\mathbf k)\mathbf a
=\lambda(\mathbf k)S(\mathbf k)\mathbf a.
}
\tag{B1}
$$

Here, $H_{\mathrm{TB}}$ and $S$ are the projected operator and overlap matrices in $\mathcal B_{\mathrm{loc}}$. The orthonormal tight-binding approximation sets $S\simeq I$. Expanding the projected operator around $\Gamma$ gives

$$
\boxed{
H_{\mathrm{TB}}(\mathbf k)
=H_{\mathrm{TB}}(\mathbf 0)
+k_i\left.\partial_{k_i}H_{\mathrm{TB}}\right|_{\Gamma}
+\mathcal O(k^2).
}
\tag{B2}
$$

Equations (B1) and (B2) provide the localized-basis origin of the effective tight-binding and $k\cdot p$ Hamiltonians used in S2.

%%
FINAL COMPLETION CHECKLIST
- Replace editorial callouts and comments with submission prose.
- Lock the geometry labels and surviving point group.
- Translate every H_z-orbital label into the electric-polarization convention.
- Add the explicit six-site matrix and fitted coupling values.
- Validate the microscopic target mass formula.
- Insert full-wave radiation coefficients and interface contributions.
- Specify the analyzed unit-cell set and verify its finite-quotient condition.
- Insert the complex $A_b(\mathbf k)$, Parseval error, and field-reconstruction error.
- Validate the far-field momentum-mixing kernel $\mathcal M(\mathbf K,\mathbf k)$ against the complete finite field.
- Add citations after the derivation and claims are stable.
- Keep the topological classification outside Sections S1--S4.
%%
