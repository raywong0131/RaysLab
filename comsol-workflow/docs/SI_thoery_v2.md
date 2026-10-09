---
title: "SI Theory"
aliases:
  - "SI_thoery"
  - "Supplementary Theory"
type: "supporting-information"
status: "framework"
created: 2026-08-06
tags:
  - supporting-information
  - manuscript
  - theory
related:
  - "[[SI_thoery_cn]]"
  - "[[00_Supporting Information Index]]"
  - "[[Supporting Information Draft]]"
  - "[[C2v_py_mode_Ex_Ey_symmetry_supplementary_note]]"
  - "[[Dirac_C2_dipolar_singularity_theory_SI]]"
  - "[[SI_theory_WuHu_radiation_integral]]"
  - "[[I. Tight Binding C6v to C2v Dipolar Singularity]]"
  - "[[II. Group Theory Symmetry Radiation Dipolar Singularity]]"
  - "[[Anisotropic Dirac Cavity Derivation]]"
---

# Supplementary Theory

Chinese companion: [[SI_thoery_cn|SI Theory 中文版]]

> [!abstract] Draft status and theoretical spine
> This note is a working framework for the complete Supplementary Theory.
> Sections S1--S4 contain the derivation spine and the equations to retain.
> Section S5 develops the finite-cavity wave-vector decomposition and its connection to radiation.
> Editorial callouts and Obsidian comments will be removed before submission.

The theory separates the closed-system band problem from the open-system radiation problem. The Hermitian Maxwell operator determines the band ordering and gap. A distinct radiation operator determines the outgoing amplitudes. The target condition is

$$
D(\Gamma)|p_y\rangle=\mathbf 0,
\qquad
\Delta_y(\Gamma)\neq 0,

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
\nabla\times
\frac{1}{\varepsilon(\mathbf r)}
\nabla\times\mathbf H(\mathbf r)
=
\left(\frac{\omega}{c}\right)^2
\mathbf H(\mathbf r).
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
\hat\Theta_{\mathbf k}\mathbf u_{n\mathbf k}
=
\left(\frac{\omega_{n\mathbf k}}{c}\right)^2
\mathbf u_{n\mathbf k},
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
-\nabla_{\parallel}\cdot
\left[
\frac{1}{\varepsilon(x,y)}
\nabla_{\parallel}H_z
\right]
=
\left(\frac{\omega}{c}\right)^2H_z.
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
\tag{S6}
$$

When the structure returns to $C_{6v}$, $\alpha=\beta=\gamma=1$ and $\mu=0$, so Eq. (S6) reduces exactly to Eq. (S4). The additional condition $\eta=1$ gives $t_0=t_1$.

### S2.2 Six-site standing-wave basis

At $\Gamma$, the $C_{6v}$ Hamiltonian is diagonal in six real standing-wave vectors. The coefficients of these vectors describe the scalar $H_z$ distribution on the six sites:

$$
\begin{aligned}
|s\rangle&=\frac{1}{\sqrt6}(1,1,1,1,1,1)^{\mathsf T},\\
|p_x\rangle&=\frac{1}{2}(0,1,1,0,-1,-1)^{\mathsf T},\\
|p_y\rangle&=\frac{1}{\sqrt{12}}(-2,-1,1,2,1,-1)^{\mathsf T},\\
|d_{x^2-y^2}\rangle&=\frac{1}{\sqrt{12}}(-2,1,1,-2,1,1)^{\mathsf T},\\
|d_{xy}\rangle&=\frac{1}{2}(0,1,-1,0,1,-1)^{\mathsf T},\\
|f\rangle&=\frac{1}{\sqrt6}(1,-1,1,-1,1,-1)^{\mathsf T}.
\end{aligned}
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
U_6^\dagger H_{\mathrm{site}}^{C_{6v}}(\Gamma)U_6
=
\operatorname{diag}
\left(
-2t_0-t_1,-m,-m,m,m,2t_0+t_1
\right),
\qquad
m=t_0-t_1.
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
\begin{aligned}
|d_\pm\rangle
&=\frac{|d_{x^2-y^2}\rangle\pm i|d_{xy}\rangle}{\sqrt2},\\
|p_\pm\rangle
&=-\frac{|p_x\rangle\mp i|p_y\rangle}{\sqrt2}.
\end{aligned}
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
H_{C_{6v}}^{\mathrm{spin}}(\mathbf k)
=
\begin{pmatrix}
m&-vk_+&0&0\\
-vk_-&-m&0&0\\
0&0&m&vk_-\\
0&0&vk_+&-m
\end{pmatrix}.
\tag{S11}
$$

The same operator, with $\mathcal B_{\mathrm{sw}}=\{|d_{xy}\rangle,|d_{x^2-y^2}\rangle,|p_y\rangle,|p_x\rangle\}$ as the ordered basis, is

$$
H_{C_{6v}}^{\mathrm{sw}}(\mathbf k)
=
\begin{pmatrix}
m&0&-ivk_y&ivk_x\\
0&m&ivk_x&ivk_y\\
ivk_y&-ivk_x&-m&0\\
-ivk_x&-ivk_y&0&-m
\end{pmatrix}.
\tag{S12}
$$

Equations (S11) and (S12) give the same $C_{6v}$ operator through first order in $\mathbf k$. The rotating basis exposes two time-reversed Dirac sectors, while the standing-wave basis exposes the real $p/d$ field patterns. At $\Gamma$, the $p$ and $d$ doublets have eigenfrequencies $-m$ and $+m$, respectively. A single Dirac mass $m$ therefore controls both polarization-resolved $p/d$ separations.

The projected Maxwell eigenproblem is written in frequency units, and all coupling parameters below use the convention fixed at the end of S1.

The $C_{2v}$ modulation removes the equivalence of the two standing-wave sectors. The two symmetry-resolved Dirac masses are

$$
m_x=\gamma(t_0-t_1),
\qquad
m_y=\frac{(4\beta-\gamma)t_0-(2\alpha+\gamma)t_1}{3}.
\tag{S13}
$$

It is useful to define

$$
M=\frac{m_x+m_y}{2},
\qquad
\Delta=\frac{m_y-m_x}{2},
\qquad
v_+=\frac{\chi+\gamma}{2}v,
\qquad
v_-=\frac{\chi-\gamma}{2}v,
\qquad
\chi=\frac{4\alpha-\gamma}{3}.
\tag{S14}
$$

When the structure returns to $C_{6v}$, these parameters reduce to

$$
\alpha=\beta=\gamma=1,
\quad \mu=0,
\quad m_x=m_y=m,
\quad M=m,
\quad \Delta=0,
\quad v_+=v,
\quad v_-=0.
\tag{S15}
$$

With $\mathcal B_{\mathrm{spin}}=\{|d_+\rangle,|p_+\rangle,|d_-\rangle,|p_-\rangle\}$ as the ordered basis, the first-order $C_{2v}$ Hamiltonian is

$$
H_{C_{2v}}^{\mathrm{spin}}(\mathbf k)
=
\begin{pmatrix}
M&-v_+k_x-i\gamma vk_y&\Delta+\mu&v_-k_x\\
-v_+k_x+i\gamma vk_y&-M&-v_-k_x&\Delta-\mu\\
\Delta+\mu&-v_-k_x&M&v_+k_x-i\gamma vk_y\\
v_-k_x&\Delta-\mu&v_+k_x+i\gamma vk_y&-M
\end{pmatrix}
+\mathcal O(k^2).
\tag{S16}
$$

The $k$-independent entries $\Delta+\mu$ and $\Delta-\mu$ couple the two rotating $d$ states and the two rotating $p$ states, respectively. Hence, the pseudospin block structure visible in Eq. (S11) is no longer preserved. The coefficient $v_-$ and the unequal $k_x$ and $k_y$ couplings record the first-order velocity anisotropy.

With $\mathcal B_{\mathrm{sw}}=\{|d_{xy}\rangle,|d_{x^2-y^2}\rangle,|p_y\rangle,|p_x\rangle\}$ as the ordered basis, transforming Eq. (S16) gives

$$
H_{C_{2v}}^{\mathrm{sw}}(\mathbf k)
=
\begin{pmatrix}
m_x-\mu&0&-i\gamma vk_y&i\gamma vk_x\\
0&m_y+\mu&i\chi vk_x&i\gamma vk_y\\
i\gamma vk_y&-i\chi vk_x&-m_y+\mu&0\\
-i\gamma vk_x&-i\gamma vk_y&0&-m_x-\mu
\end{pmatrix}
+\mathcal O(k^2).
\tag{S17}
$$

The real standing waves diagonalize Eq. (S17) at $\Gamma$, but not at a generic two-dimensional wave vector. Along $k_x$, the first-order coupled pairs are $(d_{xy},p_x)$ and $(d_{x^2-y^2},p_y)$. Along $k_y$, they are $(d_{xy},p_y)$ and $(d_{x^2-y^2},p_x)$.

Symmetry reduction therefore changes a single mass into two independent masses rather than merely perturbing one isotropic gap. The corresponding $\Gamma$-point separations are

$$
\Delta_x
=\left|\omega_{d_{xy}}-\omega_{p_x}\right|
=2|m_x|,
\qquad
\Delta_y
=\left|\omega_{d_{x^2-y^2}}-\omega_{p_y}\right|
=2|m_y|.
\tag{S18}
$$

The target electric-polarized $p_y$ state belongs to the $m_y$ sector. Section S2 establishes how the $C_{2v}$ modulation changes its Hermitian band separation. Radiation is a separate open-system question addressed in Sections S3 and S4.

## S3. Group theory analysis of Dipolar singularity

### S3.1 Conditions for coupling to free space and mode symmetry

Below the first diffraction threshold, a periodic slab at $\Gamma$ couples only to the zeroth-order outgoing plane waves. For polarization $\alpha\in\{x,y\}$, the complex radiation amplitude is

$$
c_\alpha(\Gamma)
=\langle\mathcal R_\alpha|\psi\rangle
\propto
\frac{1}{A_{\mathrm{cell}}}
\int_{\mathrm{cell}}E_\alpha(\boldsymbol\rho,z_0)\,d^2\rho.
\tag{S19}
$$

The vertical parities of the mode and outgoing wave must first be compatible. For the TE-like modes considered here, this condition is assumed. The remaining in-plane selection rule requires the direct product to contain the identity representation:

$$
\Gamma_{\mathrm{mode}}\otimes
\Gamma_{\mathrm{rad}}
\supset A_1.
\tag{S20}
$$

At $\Gamma$, the $C_{6v}$ dipolar and quadrupolar doublets and the normal-radiation channel transform as

$$
\Gamma_{\mathrm{rad}}^{C_{6v}}=E_1,
\qquad
\Gamma(p_x,p_y)=E_1,
\qquad
\Gamma(d_{x^2-y^2},d_{xy})=E_2.
\tag{S21}
$$

Under the $C_{2v}$ modulation, these representations reduce as

$$
E_1\downarrow C_{2v}=B_1\oplus B_2,
\qquad
E_2\downarrow C_{2v}=A_1\oplus A_2.
\tag{S22}
$$

where $x$ and $p_x$ transform as $B_1$, while $y$ and $p_y$ transform as $B_2$. Because $H_z$ is an axial component with intrinsic representation $A_2$, the scalar magnetic-field pattern of the target full-vector $p_y$ mode satisfies

$$
\Gamma(H_z)=B_2\otimes A_2=B_1.
\tag{S23}
$$

The component analysis starts from Maxwell's equations. For any in-plane point-group symmetry, a TE mode obeys

$$
E_x=\frac{i}{\omega\varepsilon_0\varepsilon}\partial_yH_z,
\qquad
E_y=-\frac{i}{\omega\varepsilon_0\varepsilon}\partial_xH_z.
\tag{S24}
$$

Equation (S24) is symmetry independent. The point group enters only when the transformation properties of $H_z$, $1/\varepsilon$, and the derivatives are applied to this Maxwell relation.

### S3.2 Symmetry descent from $C_{6v}$ to $C_{2v}$

In the $C_{6v}$ structure, $E_1\otimes E_1$ contains $A_1$, whereas $E_2\otimes E_1$ does not. The $p$ doublet can therefore radiate normally, while the $d$ doublet is forbidden from the same zeroth-order channel. Thus, a dipolar state is not generically dark in the unmodulated structure.

The $C_{2v}$ modulation retains $C_2(z)$ and the two mirrors $\sigma_x\equiv\sigma_{yz}$ and $\sigma_y\equiv\sigma_{xz}$. Equation (S22) shows that the two-dimensional $E_1$ representation separates into the nondegenerate $B_1$ and $B_2$ sectors. The free-space channel separates in the same way. Consequently, $p_x$ remains allowed to couple to the $x$-polarized channel, and $p_y$ remains allowed to couple to the $y$-polarized channel. The crossed couplings are forbidden because $B_1\otimes B_2=A_2$ does not contain $A_1$.

The essential change is therefore polarization resolution, rather than complete radiation suppression. Lowering the symmetry lifts the dipolar degeneracy and allows the two matched radiation amplitudes to evolve independently. For the target $p_y$ mode, $C_{2v}$ forbids the $x$-polarized component but still permits the $y$-polarized component. Any zero of the latter must arise from interference within an allowed channel.

### S3.3 $\Gamma$-point selection rule and radiation expansion

For the target $p_y$ mode, $1/\varepsilon$ belongs to $A_1$, while $\partial_x$ and $\partial_y$ belong to $B_1$ and $B_2$, respectively. Applying Eq. (S24) to Eq. (S23) gives

$$
\Gamma(E_x)=A_2,
\qquad
\Gamma(E_y)=A_1.
\tag{S25}
$$

Only the $A_1$ scalar distribution can have a nonzero unit-cell average. Equations (S19) and (S25) therefore give the $\Gamma$-point selection rule

$$
c_x(\Gamma;\lambda)=0,
\qquad
c_y(\Gamma;\lambda)=c_{y0}(\lambda),
\tag{S26}
$$

for every mirror-preserving $C_{2v}$ geometry. The first zero is enforced by symmetry, whereas $c_{y0}$ is allowed. Here, $\lambda$ is a continuous $C_{2v}$-preserving geometric parameter. For a lossless reciprocal slab, a symmetry-compatible phase convention makes $c_{y0}(\lambda)$ a continuous signed amplitude. Oppositely signed field contributions can therefore tune the allowed amplitude through zero:

$$
c_{y0}(\lambda_-)c_{y0}(\lambda_+)<0
\quad\Longrightarrow\quad
\exists\,\lambda_*\in(\lambda_-,\lambda_+)
\text{ such that }c_{y0}(\lambda_*)=0.
\tag{S27}
$$

This is an interference zero within the allowed $y$-polarized channel. It is neither a second symmetry prohibition nor a band-degeneracy condition. Section S4 relates this sign change to the triangular-hole boundaries.

The same selection rule constrains the radiation amplitudes near $\Gamma$. For an isolated target band in an analytic, symmetry-adapted Bloch gauge,

$$
\mathbf c(g\mathbf k)
=\chi_{B_2}(g)R_g\mathbf c(\mathbf k),
\qquad g\in C_{2v},
\tag{S27a}
$$

where $R_g$ acts on the in-plane polar vector. Hence,

$$
\Gamma(c_x)=B_1\otimes B_2=A_2,
\qquad
\Gamma(c_y)=B_2\otimes B_2=A_1.
\tag{S27b}
$$

The two mirror operations impose

$$
\begin{aligned}
c_x(-k_x,k_y)&=-c_x(k_x,k_y),
&c_x(k_x,-k_y)&=-c_x(k_x,k_y),\\
c_y(-k_x,k_y)&=+c_y(k_x,k_y),
&c_y(k_x,-k_y)&=+c_y(k_x,k_y).
\end{aligned}
\tag{S27c}
$$

Analyticity then fixes the lowest-order Cartesian expansion:

$$
\begin{aligned}
c_x(k_x,k_y)
&=a_{xy}k_xk_y+\mathcal O(|\mathbf k|^4),\\
c_y(k_x,k_y)
&=c_{y0}+a_{xx}k_x^2+a_{yy}k_y^2
+\mathcal O(|\mathbf k|^4).
\end{aligned}
\tag{S27d}
$$

At the tuned periodic-structure zero, $c_{y0}=0$. Equation (S27d) then contains neither a constant nor a linear radiation term. This result follows from group theory and does not require a topological classification.

### S3.4 Why $C_2$ is insufficient

Under $C_2$ alone, the $x$- and $y$-polarized outgoing waves belong to the same odd representation. The target dipolar mode can then couple to both channels. Without the two mirror parities, $C_2$ cannot distinguish the crossed and matched polarizations and cannot enforce $c_x(\Gamma)=0$.

The component-resolved selection rule therefore requires the full $C_{2v}$ symmetry of the geometry and numerical boundary conditions. Mirror-breaking perturbations can reopen the nominally forbidden component, even when $C_2$ remains intact.

## S4. Normal out-of-plane radiation from triangular-hole boundaries

This section relates the normal outgoing amplitude to the triangular-hole boundaries. The derivation is exact for the two-dimensional TE model at $\Gamma$; its use for the three-dimensional TE-like slab requires the corresponding vertical-overlap factor.

### S4.1 From a unit-cell area integral to interface line integrals

Below the first diffraction threshold, the normal outgoing coefficient is proportional to the zero-order electric-field component,

$$
\mathbf c(\mathbf0)\propto\langle\mathbf E\rangle
=\frac{1}{A_{\mathrm{cell}}}
\int_{\mathrm{cell}}\mathbf E(\boldsymbol\rho)\,d^2\rho.
\tag{S28}
$$

For a TE mode, Maxwell's equations give

$$
E_x=\frac{i}{\omega\varepsilon_0\varepsilon}\partial_yH_z,
\qquad
E_y=-\frac{i}{\omega\varepsilon_0\varepsilon}\partial_xH_z.
\tag{S29}
$$

The area-to-line conversion follows by integration by parts. For $E_y$,

$$
\begin{aligned}
\langle E_y\rangle
&=-\frac{i}{\omega\varepsilon_0A_{\mathrm{cell}}}
\int_{\mathrm{cell}}\frac{1}{\varepsilon}\partial_xH_z\,d^2\rho\\
&=-\frac{i}{\omega\varepsilon_0A_{\mathrm{cell}}}
\left[
\int_{\mathrm{cell}}\partial_x\!\left(\frac{H_z}{\varepsilon}\right)d^2\rho
-\int_{\mathrm{cell}}H_z\partial_x\!\left(\frac{1}{\varepsilon}\right)d^2\rho
\right]\\
&=-\frac{i}{\omega\varepsilon_0A_{\mathrm{cell}}}
\left[
\oint_{\partial\mathrm{cell}}\frac{H_z}{\varepsilon}N_x\,dl
-\int_{\mathrm{cell}}H_z\partial_x\!\left(\frac{1}{\varepsilon}\right)d^2\rho
\right]\\
&=\frac{i}{\omega\varepsilon_0A_{\mathrm{cell}}}
\int_{\mathrm{cell}}H_z\partial_x\!\left(\frac{1}{\varepsilon}\right)d^2\rho.
\end{aligned}
\tag{S30}
$$

The last equality holds only at $\mathbf k=\mathbf0$: opposite unit-cell edges have opposite outward normals, while $H_z/\varepsilon$ has identical values because it is periodic at $\Gamma$. For a stepwise dielectric profile, $\nabla(1/\varepsilon)$ is supported only on the air--dielectric interfaces. Choosing $\hat{\mathbf n}$ from dielectric to air and defining $\Delta(1/\varepsilon)=1/\varepsilon_{\mathrm{air}}-1/\varepsilon_{\mathrm{diel}}$ gives

$$
\begin{aligned}
\langle E_x\rangle
&=-\mathcal K_\Gamma\oint_S H_z n_y\,dl,
&\langle E_y\rangle
&=\mathcal K_\Gamma\oint_S H_z n_x\,dl,\\
\langle\mathbf E\rangle
&=\mathcal K_\Gamma\oint_S
H_z(\hat{\mathbf z}\times\hat{\mathbf n})\,dl,
&\mathcal K_\Gamma
&=\frac{i\Delta(1/\varepsilon)}
{\omega\varepsilon_0A_{\mathrm{cell}}}.
\end{aligned}
\tag{S31}
$$

Thus, at $\Gamma$, a homogeneous hole interior supplies no independent contribution to normal radiation. The hole affects the outgoing amplitude through the position and length of its interfaces and through the self-consistent value of $H_z$ sampled on them.

### S4.2 Boundary-field extraction under $C_{6v}$ and $C_{2v}$ modulation

For the $n$th triangular hole, write its boundary contribution edge by edge as

$$
\mathbf q_n
=\oint_{\partial h_n}H_z\hat{\mathbf n}\,dl
=\sum_{m=1}^{3}\ell_{n,m}\overline H_{n,m}\hat{\mathbf n}_{n,m},
\qquad
\overline H_{n,m}=\frac{1}{\ell_{n,m}}
\int_{e_{n,m}}H_z\,dl.
\tag{S32}
$$

Equation (S32) shows the extraction action of a triangular boundary: each side selects a signed segment of the magnetic-field distribution and projects it through its normal. Increasing a side length changes both the amount of field sampled and the position at which it is sampled. A constant $H_z$ around a closed triangle gives $\mathbf q_n=H_z\oint\hat{\mathbf n}dl=\mathbf0$; radiation therefore depends on the edge-to-edge variation of the field, not on the perimeter alone.

We first consider the $C_{6v}$ structure, in which the six identical holes lie at $\mathbf R_n=(\eta a/3)\hat{\mathbf r}_n$. At $\eta=1$, the six-site cell returns to the undeformed honeycomb structure. The primitive-translation structure factor cancels the zero-order boundary sum at the folded $\Gamma$ point, producing the Dirac singularity with $\mathbf c(\mathbf0)=\mathbf0$. A radial contraction or expansion, $\eta\neq1$, preserves $C_{6v}$ but removes this translation-induced cancellation. For small $|\eta-1|$, the symmetry-reduced boundary response is

$$
\mathbf q_n(\eta)
=\mathbf q_n^{(0)}
+\chi_\eta(\eta-1)a_n\hat{\mathbf r}_n
+\mathcal O[(\eta-1)^2],
\qquad
\sum_n\mathbf q_n^{(0)}=\mathbf0,
\tag{S33}
$$

where $a_n$ is the six-site amplitude and $\chi_\eta$ is obtained from the complete three-edge response. Using the standing-wave vectors in Eq. (S7),

$$
\begin{aligned}
\sum_n(p_x)_n\hat{\mathbf r}_n&=\sqrt3\,\hat{\mathbf y},
&\langle\mathbf E\rangle_{p_x}
&=C_\eta(\eta-1)\hat{\mathbf x},\\
\sum_n(p_y)_n\hat{\mathbf r}_n&=-\sqrt3\,\hat{\mathbf x},
&\langle\mathbf E\rangle_{p_y}
&=C_\eta(\eta-1)\hat{\mathbf y}.
\end{aligned}
\tag{S34}
$$

With the normal and time-harmonic conventions used above, $C_\eta=-\sqrt3\mathcal K_\Gamma\chi_\eta$. Hence both $p$ modes are dark at $\eta=1$, become normally radiative when $\eta\neq1$, and satisfy

$$
I_{\Gamma,p}\propto|\eta-1|^2
\qquad (C_{6v}).
\tag{S35}
$$

We next introduce the $C_{2v}$ size modulation of Eq. (S5). Holes 1 and 4 form class $A$ with side length $b_A=\zeta b$, while the other four form class $B$ with $b_B=s(\zeta)b$ and $s(\zeta)=\sqrt{(3-\zeta^2)/2}$. For the target $p_y$ mode, $H_z(-x,y)=-H_z(x,y)$ and $H_z(x,-y)=H_z(x,y)$. The two mirrors therefore cancel the $x$-polarized amplitude, whereas the allowed $y$-polarized amplitude reduces exactly to

$$
c_x(\Gamma)\propto\langle E_x\rangle=0,
\qquad
r_y\equiv c_y(\Gamma)
\propto\langle E_y\rangle
=\mathcal K_\Gamma\left[2Q_x^{(A)}+4Q_x^{(B)}\right],
\qquad
Q_x^{(\nu)}=\oint_{\partial h_\nu}H_z n_x\,dl.
\tag{S36}
$$

The $E_y$ channel is symmetry allowed; its zero is produced when the two signed hole classes satisfy $Q_x^{(A)}+2Q_x^{(B)}=0$. If the $H_z$ distribution and the edge averages change weakly while the triangle sizes are varied, a common-phase convention makes these averages signed scalars and Eq. (S32) gives the frozen-field estimate

$$
\begin{aligned}
r_y^{\mathrm{fr}}(\zeta)
&\propto\mathcal K_\Gamma b
\left[2\zeta F_A+4s(\zeta)F_B\right],\\
F_\nu&=\sum_{m=1}^{3}n_{x,\nu m}\overline H_{\nu m},
&\rho&=-\frac{2F_B}{F_A}>0,
&\zeta_{\mathrm{fr}}^*&=\sqrt{\frac{3\rho^2}{2+\rho^2}}.
\end{aligned}
\tag{S37}
$$

The analytic root requires $F_AF_B<0$ and should be evaluated using edge integrals from the simulated $p_y$ field. Self-consistent field redistribution changes the coefficients but not the cancellation mechanism. Locally, the complete response may be written as

$$
r_y(\eta,\zeta)
=C_{\eta y}(\eta-1)+C_{\zeta y}(\zeta-1)+\mathcal O(\delta^2),
\qquad
\zeta_y^*(\eta)
=1-\frac{C_{\eta y}}{C_{\zeta y}}(\eta-1).
\tag{S38}
$$

Here, $\delta$ collectively denotes $\eta-1$ and $\zeta-1$. Equation (S38) summarizes the two steps: the $C_{6v}$ radial displacement first opens the $p_y$ radiation amplitude, and the $C_{2v}$ size modulation subsequently tunes the signed boundary extraction until that allowed amplitude vanishes.

### S4.3 Radiation away from $\Gamma$

The interface-only identity in Eq. (S31) is specific to $\mathbf k=\mathbf0$. At finite wave vector, write $H_{z,\mathbf k}=e^{i\mathbf k\cdot\boldsymbol\rho}u_{H,\mathbf k}$ with periodic $u_{H,\mathbf k}$. Differentiating the Bloch phase produces an additional unit-cell area term:

$$
\begin{aligned}
\mathbf c(\mathbf k)
={}&\frac{i\Delta(1/\varepsilon)}
{\omega_{\mathbf k}\varepsilon_0A_{\mathrm{cell}}}
\oint_Su_{H,\mathbf k}
(\hat{\mathbf z}\times\hat{\mathbf n})\,dl\\
&+\frac{\hat{\mathbf z}\times\mathbf k}
{\omega_{\mathbf k}\varepsilon_0A_{\mathrm{cell}}}
\int_{\mathrm{cell}}\frac{u_{H,\mathbf k}}{\varepsilon}\,d^2\rho.
\end{aligned}
\tag{S39}
$$

The first term remains the interface contribution, whereas the explicit wave-vector prefactor makes the second term vanish at exact $\Gamma$. Expanding both terms and imposing the $C_{2v}$ mirror parities obtained in S3 gives

$$
\begin{aligned}
c_x(k_x,k_y)
&=a_{xy}k_xk_y+\mathcal O(|\mathbf k|^4),\\
c_y(k_x,k_y)
&=c_{y0}+a_{xx}k_x^2+a_{yy}k_y^2
+\mathcal O(|\mathbf k|^4),\\
I(k_x,k_y)\big|_{c_{y0}=0}
&=\mathcal N_I(0)
\left[|a_{xy}|^2k_x^2k_y^2
+|a_{xx}k_x^2+a_{yy}k_y^2|^2\right]
+\mathcal O(|\mathbf k|^6).
\end{aligned}
\tag{S40}
$$

Thus, at the tuned $p_y$ zero, neither radiation component contains a constant or linear term. The quadratic, anisotropic Cartesian dependence in Eq. (S40) is identical to the group-theoretical result in Eq. (S27d), while Eq. (S39) identifies its Maxwell origin in the combined interface and finite-$\mathbf k$ area contributions.

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
\frac{1}{N}
\sum_{\mathbf R\in S}
e^{i(\mathbf k-\mathbf k')\cdot\mathbf R}
=\delta_{\mathbf k,\mathbf k'},
\qquad \mathbf k,\mathbf k'\in K.
\tag{S41}
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

If the selected cells form an irregular mask, contain missing cells, or do not represent a finite quotient, Eq. (S41) does not hold. The plane-wave dictionary can still be used, but the coefficients must then be obtained from its Gram matrix or from a least-squares fit. A conventional Fourier transform of an arbitrary mask is a windowed spectrum and contains wave-vector leakage.

### S5.3 Extraction of the wave-vector-resolved periodic cell content

This manuscript uses the Bloch convention of S1,

$$
\boldsymbol\Psi_{b\mathbf k}(\mathbf r)
=e^{i\mathbf k\cdot\mathbf r}
\mathbf u_{b\mathbf k}(\boldsymbol\rho,z).
$$

The unitary cell-index Fourier pair consistent with this convention is

$$
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
\tag{S42}
$$

The function $\mathbf F$ still contains the intra-cell Bloch phase. Removing it defines the periodic cell content

$$
\mathbf G(\mathbf k,\boldsymbol\rho,z)
=e^{-i\mathbf k\cdot\boldsymbol\rho}
\mathbf F(\mathbf k,\boldsymbol\rho,z).
$$

The finite mode is then reconstructed exactly on the selected cells as

$$
\boldsymbol\Psi_f(\mathbf R+\boldsymbol\rho,z)
=\frac{1}{\sqrt N}
\sum_{\mathbf k\in K}
e^{i\mathbf k\cdot(\mathbf R+\boldsymbol\rho)}
\mathbf G(\mathbf k,\boldsymbol\rho,z).
\tag{S43}
$$

Equation (S43) is already an exact Bloch-component decomposition of the finite field on $S$. It does not assume a slowly varying envelope and does not require prior knowledge of any band.

The relation to an ordinary spatial Fourier spectrum follows by expanding the periodic function in reciprocal vectors:

$$
\mathbf G(\mathbf k,\boldsymbol\rho,z)
=\sum_{\mathbf G_r}
\mathbf G_{\mathbf G_r}(\mathbf k,z)
e^{i\mathbf G_r\cdot\boldsymbol\rho}.
$$

Substitution into Eq. (S43) gives spatial harmonics at

$$
\mathbf q=\mathbf k+\mathbf G_r.
\tag{S44}
$$

Consequently, Fourier peaks separated by reciprocal lattice vectors are replicas of the same reduced Bloch wave vector. Folding the global Fourier spectrum into the first Brillouin zone identifies $\mathbf k$, while $\mathbf G(\mathbf k,\boldsymbol\rho,z)$ retains the full intra-cell field profile that a simple Fourier-intensity map discards.

Parseval's identity gives

$$
\sum_{\mathbf R\in S}
\|\boldsymbol\Psi_{\mathbf R}\|_B^2
=\sum_{\mathbf k\in K}
\|\mathbf G(\mathbf k)\|_B^2.
\tag{S45}
$$

The total wave-vector weight is therefore

$$
W(\mathbf k)=\|\mathbf G(\mathbf k)\|_B^2,
\qquad
\eta(\mathbf k)=
\frac{W(\mathbf k)}{\sum_{\mathbf k'\in K}W(\mathbf k')}.
$$

The map $W(\mathbf k)$ identifies where the finite mode resides in the Brillouin zone. It does not by itself identify which Bloch band contributes at a given wave vector.

If fields are exported using the COMSOL convention $e^{-i\mathbf k\cdot\mathbf r}$, the preprocessing transform uses the opposite sign and the reported wave vector is mapped by $\mathbf k\rightarrow-\mathbf k$ relative to Eqs. (S42)--(S44). The same convention must be used for the finite field and the reference Bloch modes.

### S5.4 Projection onto infinite-crystal Bloch bands

Let $\mathbf u_{b\mathbf k}^{\infty}(\boldsymbol\rho,z)$ be the periodic profile of band $b$ in an infinite reference crystal at the same reduced wave vector $\mathbf k$. For a complete $B$-orthonormal cell eigenbasis,

$$
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
\tag{S46}
$$

Combining Eqs. (S43) and (S46) yields the finite-mode Bloch-band expansion

$$
\boldsymbol\Psi_f(\mathbf r)
=\frac{1}{\sqrt N}
\sum_{\mathbf k\in K}
\sum_b
A_b(\mathbf k)
e^{i\mathbf k\cdot\mathbf r}
\mathbf u_{b\mathbf k}^{\infty}(\boldsymbol\rho,z).
\tag{S47}
$$

With all bands retained, Eq. (S47) is complete on the selected cells. It is a representation statement: the finite, terminated, or defective structure is not being identified with an infinite-crystal eigenproblem. The coefficients only specify how the finite field decomposes in the chosen complete basis.

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
\mathbf u_{b\mathbf k}^{\infty}
\rightarrow e^{i\chi_b(\mathbf k)}
\mathbf u_{b\mathbf k}^{\infty},
\qquad
A_b(\mathbf k)
\rightarrow e^{-i\chi_b(\mathbf k)}A_b(\mathbf k).
\tag{S48}
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

The familiar carrier-envelope form is a controlled compression of Eq. (S47), rather than its starting assumption. If one target band dominates, its periodic profile varies weakly over a narrow set of wave vectors around $\mathbf k_0$, and a continuous gauge has been fixed, then

$$
\mathbf G(\mathbf k,\boldsymbol\rho,z)
\simeq A(\mathbf k)\mathbf u_0(\boldsymbol\rho,z).
$$

Equation (S43) reduces to

$$
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
\tag{S49}
$$

The envelope is slowly varying only when the relative spectrum $\mathbf k-\mathbf k_0$ is narrow compared with the reciprocal-lattice scale. A spatially localized mode therefore occupies a finite region of the Brillouin zone, with a characteristic width that decreases as the cavity size increases.

Let $\mathbf c_b(\mathbf k)=(c_{b,x},c_{b,y})$ denote the outgoing radiation vector of the infinite-crystal Bloch state, as derived in S3 and S4. In a translation-invariant reference channel, the radiation carried by one reduced wave vector is

$$
\mathbf d_f(\mathbf k)
=\sum_b A_b(\mathbf k)\mathbf c_b(\mathbf k).
\tag{S50}
$$

The gauge phases in Eq. (S48) cancel between $A_b$ and $\mathbf c_b$, so $\mathbf d_f$ is gauge invariant. Equation (S50) also shows why the finite-mode composition must be extracted as complex amplitudes rather than as $W(\mathbf k)$ alone.

Different reduced wave vectors radiate into different in-plane far-field momenta in the translation-invariant reference problem. A coherent sum over distinct $\mathbf k$ values at one observed momentum requires an additional finite-boundary or aperture mixing kernel. The general form is

$$
\mathbf E_{\mathrm{far}}(\mathbf K)
=\sum_{\mathbf k,b}
\mathcal M(\mathbf K,\mathbf k)
A_b(\mathbf k)\mathbf c_b(\mathbf k),
$$

where $\mathcal M$ must be obtained from the actual cavity truncation, interface scattering, or a direct Fourier transform of the complete finite field. It must not be replaced by an assumed convolution without validation. The decomposition in Eqs. (S42)--(S48) supplies the finite mode's complex Bloch content required for that subsequent radiation calculation.

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
H_{\mathrm{eff}}(E)
=U_{pd}^\dagger H_{\mathrm{site}}U_{pd}
+U_{pd}^\dagger H_{\mathrm{site}}Q
\left(E-QH_{\mathrm{site}}Q\right)^{-1}
QH_{\mathrm{site}}U_{pd}.
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

Let $\delta_{2v}$ denote the off-block modulation scale. Equations (S16) and (S17) retain the direct four-band projection through first order in $\mathbf k$. Spectral isolation of the $p/d$ manifold is the basis of the $6\times6$ to $4\times4$ reduction. The leading omitted self-energy scales are $(vk)^2/\Delta_{sf}$, $\delta_{2v}vk/\Delta_{sf}$, and $\delta_{2v}^2/\Delta_{sf}$; Eq. (A1) must be evaluated when these corrections are not negligible. The direct second-order momentum expansion is given separately in Appendix C.

### A.2 $\Gamma$-point eigenfrequencies and effective parameters

At $\Gamma$, Eq. (S17) gives

$$
\begin{aligned}
\omega_{d_{x^2-y^2}}&=m_y+\mu,\\
\omega_{p_y}&=-m_y+\mu,\\
\omega_{d_{xy}}&=m_x-\mu,\\
\omega_{p_x}&=-m_x-\mu.
\end{aligned}
\tag{A2}
$$

The onsite response shifts the centers of the two sectors by $\pm\mu$ and cancels from each $p/d$ separation. The hopping anisotropy controls $m_x$ and $m_y$. The parameters reduce as follows:

$$
\begin{aligned}
\zeta=1:\quad
&\alpha=\beta=\gamma=1,
\quad \mu=0,
\quad m_x=m_y=t_0-t_1,\\
(\eta,\zeta)=(1,1):\quad
&t_0=t_1,
\quad m_x=m_y=0.
\end{aligned}
\tag{A3}
$$

The three independent effective parameters can be extracted directly from full-wave $\Gamma$-point frequencies:

$$
\begin{aligned}
m_y&=\frac{\omega_{d_{x^2-y^2}}-\omega_{p_y}}{2},\\
m_x&=\frac{\omega_{d_{xy}}-\omega_{p_x}}{2},\\
\mu&=\frac{\omega_{d_{x^2-y^2}}+\omega_{p_y}
-\omega_{d_{xy}}-\omega_{p_x}}{4}.
\end{aligned}
\tag{A4}
$$

When fitted to full-wave frequencies, Eq. (A4) defines renormalized parameters that include remote-state shifts. The factors $\alpha$, $\beta$, and $\gamma$ require additional near-$\Gamma$ velocity information or microscopic overlap calculations. The factorized form in Eq. (S6) remains a nearest-neighbour approximation.

For the target sector,

$$
\Delta_y
=\left|\omega_{d_{x^2-y^2}}-\omega_{p_y}\right|
=2|m_y|.
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
H_{\mathrm{TB}}(\mathbf k)\mathbf a
=\lambda(\mathbf k)S(\mathbf k)\mathbf a.
\tag{B1}
$$

Here, $H_{\mathrm{TB}}$ and $S$ are the projected operator and overlap matrices in $\mathcal B_{\mathrm{loc}}$. The orthonormal tight-binding approximation sets $S\simeq I$. Expanding the projected operator around $\Gamma$ gives

$$
H_{\mathrm{TB}}(\mathbf k)
=H_{\mathrm{TB}}(\mathbf 0)
+\sum_{i=x,y}k_i
\left.\partial_{k_i}H_{\mathrm{TB}}\right|_{\Gamma}
+\frac12\sum_{i,j=x,y}k_i k_j
\left.\partial_{k_i}\partial_{k_j}H_{\mathrm{TB}}\right|_{\Gamma}
+\mathcal O(k^3).
\tag{B2}
$$

Equations (B1) and (B2) provide the localized-basis origin of the first- and second-order effective tight-binding and $k\cdot p$ Hamiltonians used in S2 and Appendix C.

## Appendix C. Second-order four-band expansion under $C_{2v}$ modulation

The main derivation in S2 retains only terms through first order in $\mathbf k$. This appendix records the direct quadratic correction obtained from the same six-site Hamiltonian, without adding the separate Löwdin self-energy generated by the remote $s/f$ subspace.

### C.1 Quadratic expansion of the intercell Bloch phases

Each intercell phase in Eq. (S6) is expanded as

$$
e^{i\mathbf k\cdot\mathbf a_j}
=1+i\mathbf k\cdot\mathbf a_j
-\frac12(\mathbf k\cdot\mathbf a_j)^2
+\mathcal O(k^3).
\tag{C1}
$$

Direct projection onto the four $p/d$ states gives the real standing-wave coefficients

$$
\begin{aligned}
\Lambda_x(\mathbf k)
&=\frac{a^2t_1\gamma}{8}(k_x^2+3k_y^2),\\
\Lambda_y(\mathbf k)
&=\frac{a^2t_1}{24}\left[(8\alpha+\gamma)k_x^2+3\gamma k_y^2\right],\\
\Lambda_{xy}(\mathbf k)
&=\frac{a^2t_1\gamma}{4}k_xk_y,\\
\Lambda_0(\mathbf k)
&=\frac{\Lambda_x+\Lambda_y}{2}
=\frac{a^2t_1}{12}
\left[(2\alpha+\gamma)k_x^2+3\gamma k_y^2\right],\\
\Lambda_2(\mathbf k)
&=\frac{\Lambda_y-\Lambda_x}{2}-i\Lambda_{xy}\\
&=\frac{a^2t_1}{24}
\left[(4\alpha-\gamma)k_x^2-3\gamma k_y^2
-6i\gamma k_xk_y\right].
\end{aligned}
\tag{C2}
$$

The parameters $\beta$ and $\mu$ do not enter Eq. (C2), because they modify only intracell or onsite matrix elements and therefore carry no Bloch phase. When the structure returns to $C_{6v}$, $\alpha=\gamma=1$, and

$$
\begin{aligned}
\Lambda_x&=\frac{a^2t_1}{8}(k_x^2+3k_y^2),
&\Lambda_y&=\frac{a^2t_1}{8}(3k_x^2+k_y^2),\\
\Lambda_{xy}&=\frac{a^2t_1}{4}k_xk_y,
&\Lambda_0&=\frac{a^2t_1}{4}(k_x^2+k_y^2),\\
\Lambda_2&=\frac{a^2t_1}{8}(k_x-ik_y)^2.
\end{aligned}
\tag{C3}
$$

### C.2 Second-order matrices in the rotating and standing-wave bases

With $\mathcal B_{\mathrm{spin}}=\{|d_+\rangle,|p_+\rangle,|d_-\rangle,|p_-\rangle\}$ as the ordered basis, the direct second-order projection is

$$
H_{C_{2v}}^{\mathrm{spin},(2)}(\mathbf k)
=
\begin{pmatrix}
M+\Lambda_0&-v_+k_x-i\gamma vk_y&\Delta+\mu+\Lambda_2&v_-k_x\\
-v_+k_x+i\gamma vk_y&-M-\Lambda_0&-v_-k_x&\Delta-\mu+\Lambda_2^*\\
\Delta+\mu+\Lambda_2^*&-v_-k_x&M+\Lambda_0&v_+k_x-i\gamma vk_y\\
v_-k_x&\Delta-\mu+\Lambda_2&v_+k_x+i\gamma vk_y&-M-\Lambda_0
\end{pmatrix}
+\mathcal O(k^3).
\tag{C4}
$$

The scalar quadratic term $\Lambda_0$ shifts the two $d$ states upward and the two $p$ states downward. The complex term $\Lambda_2$ is the symmetry-allowed momentum-dependent coupling between opposite rotating sectors.

With $\mathcal B_{\mathrm{sw}}=\{|d_{xy}\rangle,|d_{x^2-y^2}\rangle,|p_y\rangle,|p_x\rangle\}$ as the ordered basis, the same operator is

$$
H_{C_{2v}}^{\mathrm{sw},(2)}(\mathbf k)
=
\begin{pmatrix}
m_x-\mu+\Lambda_x&\Lambda_{xy}&-i\gamma vk_y&i\gamma vk_x\\
\Lambda_{xy}&m_y+\mu+\Lambda_y&i\chi vk_x&i\gamma vk_y\\
i\gamma vk_y&-i\chi vk_x&-m_y+\mu-\Lambda_y&\Lambda_{xy}\\
-i\gamma vk_x&-i\gamma vk_y&\Lambda_{xy}&-m_x-\mu-\Lambda_x
\end{pmatrix}
+\mathcal O(k^3).
\tag{C5}
$$

The term $\Lambda_{xy}\propto k_xk_y$ mixes the two $d$ standing waves and, with the same sign, the two $p$ standing waves. It vanishes along the two principal axes and is nonzero for a generic off-axis wave vector.

### C.3 Scope of the direct second-order projection

Using $U_{pd}$ defined in Appendix A, let $\Delta_{sf}$ denote the spectral separation from the remote $s/f$ subspace and let $\delta_{2v}$ denote the corresponding off-block modulation scale. Define the direct four-band projection through second order by

$$
\begin{aligned}
H_{pd}^{(2)}(\mathbf k)
&:=\left[
U_{pd}^\dagger H_{\mathrm{site}}(\mathbf k)U_{pd}
\right]_{\text{through }k^2},\\
H_{\mathrm{eff}}(\mathbf k)
&=H_{pd}^{(2)}(\mathbf k)
+\mathcal O(k^3)
+\mathcal O\!\left(
\frac{(vk)^2}{\Delta_{sf}},
\frac{\delta_{2v}vk}{\Delta_{sf}},
\frac{\delta_{2v}^2}{\Delta_{sf}}
\right).
\end{aligned}
\tag{C6}
$$

Equations (C4) and (C5) contain every direct-projection term through $k^2$, but not the additional Löwdin self-energy from virtual coupling to the $s/f$ subspace. In particular, $(vk)^2/\Delta_{sf}$ is itself quadratic in momentum and must be included for a quantitatively complete second-order model whenever the remote-state separation is not sufficiently large.

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
