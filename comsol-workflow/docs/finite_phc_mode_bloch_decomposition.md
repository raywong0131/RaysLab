# Bloch-Component Decomposition of a Finite 2D Photonic-Crystal Eigenmode

This document describes a reusable post-processing algorithm for decomposing a
finite 2D photonic-crystal eigenmode into a discrete set of Bloch-like
components.

The ideal complete-basis picture is simple. If the selected finite cell region
supports an exact finite Fourier basis, and if all infinite-crystal Bloch
profiles are known at the corresponding discrete wave vectors, then any field
on that finite region has an exact expansion

$$
E_f(\mathbf R+\boldsymbol\rho)
=
\frac{1}{\sqrt N}
\sum_{\mathbf k\in K}
\sum_b
c_b(\mathbf k)
e^{-i\mathbf k\cdot(\mathbf R+\boldsymbol\rho)}
u_b(\mathbf k,\boldsymbol\rho),
$$

on the selected cells. This is a basis expansion, not a slow-envelope ansatz.

The practical problem is that the complete Bloch dictionary is often unknown,
or only some specific $(b,\mathbf k)$ entries are known. The algorithm therefore
first extracts, without band knowledge, the per-$\mathbf k$ periodic cell
content

$$
G(\mathbf k,\boldsymbol\rho).
$$

Known bands can then be used by projection. If no band dictionary is available,
the extracted functions $G(\mathbf k,\boldsymbol\rho)$ are analyzed directly by
weights, profile similarity, clustering, and SVD. The slow-envelope picture is
recovered only as a later compression, when many nearby $\mathbf k$ components
share one nearly identical periodic unit-cell profile.

## COMSOL-compatible convention

COMSOL frequency-domain fields use the time convention

$$
\operatorname{Re}\{\widetilde E(\mathbf r)e^{+i\omega t}\}.
$$

With this convention, a Bloch wave with COMSOL wave vector $\mathbf k$ is
written as

$$
E_{b\mathbf k}(\mathbf r)
=
e^{-i\mathbf k\cdot\mathbf r}
u_b(\mathbf k,\boldsymbol\rho).
$$

Here

$$
\mathbf r=\mathbf R+\boldsymbol\rho,
$$

where $\mathbf R$ is a lattice vector labeling the unit cell and
$\boldsymbol\rho$ is the coordinate inside one unit cell.

The periodic part satisfies

$$
u_b(\mathbf k,\boldsymbol\rho+\mathbf R)=u_b(\mathbf k,\boldsymbol\rho).
$$

With this sign convention, a lattice Fourier transform with phase
$e^{+i\mathbf k\cdot\mathbf R}$ detects the COMSOL wave vector $\mathbf k$
directly.

## Data model

Choose a bulk-like analysis region consisting of complete unit cells. Let

$$
\mathbf R_{mn}=m\mathbf a_1+n\mathbf a_2,
$$

with reciprocal vectors satisfying

$$
\mathbf a_i\cdot\mathbf b_j=2\pi\delta_{ij}.
$$

Define the field restricted to cell $\mathbf R$:

$$
E_{\mathbf R}(\boldsymbol\rho)
=
E_f(\mathbf R+\boldsymbol\rho).
$$

Use a unit-cell inner product

$$
\langle f,g\rangle_B
=
\int_{\Omega_{\rm cell}}
f^*(\boldsymbol\rho)
B(\boldsymbol\rho)
g(\boldsymbol\rho)
d\boldsymbol\rho.
$$

Common choices are $B=\varepsilon$, a finite-element mass matrix, or $B=1$ for
a scalar diagnostic. The induced norm is

$$
\|f\|_B^2=\langle f,f\rangle_B.
$$

## Core finite-region assumption

The exact version of this algorithm assumes that the selected analysis cells

$$
S=\{\mathbf R_1,\ldots,\mathbf R_N\}\subset\Lambda
$$

are a complete representative set of a finite lattice quotient

$$
\Lambda/\Lambda_{\rm super},
$$

where $\Lambda_{\rm super}\subset\Lambda$ is a finite-index superlattice.
Equivalently, every lattice point can be written uniquely as

$$
\mathbf R=\mathbf R_S+\mathbf T,
\qquad
\mathbf R_S\in S,
\qquad
\mathbf T\in\Lambda_{\rm super}.
$$

Then there is a finite set of wave vectors

$$
K=\Lambda_{\rm super}^*/\Lambda^*
$$

such that

$$
\frac{1}{N}
\sum_{\mathbf R\in S}
e^{+i(\mathbf q-\mathbf k)\cdot\mathbf R}
=
\delta_{\mathbf q,\mathbf k},
\qquad
\mathbf q,\mathbf k\in K.
$$

This makes the cell-index Fourier transform unitary. The finite-lattice basis
and its completeness are summarized in
[`finite_lattice_fourier_bases.md`](finite_lattice_fourier_bases.md).

A regular block in lattice-index coordinates is the simplest example:

$$
\mathbf R_{mn}=m\mathbf a_1+n\mathbf a_2,
\qquad
m=0,\ldots,M_1-1,
\qquad
n=0,\ldots,M_2-1,
$$

with $N=M_1M_2$ cells. The corresponding lattice Fourier grid is

$$
\mathbf q_{\ell_1\ell_2}
=
\frac{\ell_1}{M_1}\mathbf b_1
+
\frac{\ell_2}{M_2}\mathbf b_2,
\qquad
\ell_1=0,\ldots,M_1-1,
\qquad
\ell_2=0,\ldots,M_2-1.
$$

This is a rectangle in the integer lattice indices $(m,n)$, not necessarily a
Cartesian rectangle in real space. The primitive vectors $\mathbf a_1$ and
$\mathbf a_2$ may be nonorthogonal, so the same statement applies to a
hexagonal or triangular Bravais lattice if the selected cells form such an
$M_1\times M_2$ block.

Other quotient-compatible regions are also possible. For example, a standard
centered hexagonal shell patch on a triangular lattice is not a parallelogram,
but it can still be a complete representative set of a suitable superlattice
quotient.

If the selected core cut is not a quotient representative set, for example
because cells are missing or the boundary is an arbitrary mask, the exact
orthogonality is lost. Then replace the exact discrete Fourier transform below
by a least-squares or Gram-matrix version.

### Strip with one Bloch-periodic direction

A strip simulation with one finite direction and one Bloch-periodic direction
is an exact special case of the same assumption. Let
$\mathbf a_\parallel$ denote the periodic lattice vector and let
$\mathbf a_\perp$ label rows across the strip. For a simulation in a fixed
Bloch sector $k_\parallel$,

$$
E_f(\mathbf r+\mathbf a_\parallel)
=
e^{-i k_\parallel\cdot\mathbf a_\parallel}E_f(\mathbf r).
$$

Thus the periodic direction is already quotiented out by the boundary
condition. Equivalently, the gauge-transformed field

$$
\widetilde E_f(\mathbf r)=e^{+i k_\parallel\cdot\mathbf r}E_f(\mathbf r)
$$

is periodic along $\mathbf a_\parallel$. The remaining finite cell index is the
row coordinate $n$ in

$$
\mathbf R_n=n\mathbf a_\perp
\quad
\text{mod } \langle\mathbf a_\parallel\rangle.
$$

If the selected bulk rows are a complete interval

$$
S_N=\{n_0,n_0+1,\ldots,n_0+N-1\},
$$

then they are a complete representative set of $\mathbb Z/N\mathbb Z$. The
finite Fourier wave vectors in this fixed Bloch sector are

$$
\mathbf k_m
=
\mathbf k_\parallel+\frac{m}{N}\mathbf b_\perp,
\qquad
m=0,\ldots,N-1,
$$

where $\mathbf b_\perp\cdot\mathbf a_\perp=2\pi$ and
$\mathbf b_\perp\cdot\mathbf a_\parallel=0$. The row transform is

$$
F_m(\boldsymbol\rho)
=
\frac{1}{\sqrt N}
\sum_{n\in S_N}
\exp\left(+\frac{2\pi i mn}{N}\right)
E_n(\boldsymbol\rho).
$$

The same intra-cell phase removal used below then becomes

$$
G_m(\boldsymbol\rho)
=
e^{+i\mathbf k_m\cdot\boldsymbol\rho}
F_m(\boldsymbol\rho).
$$

When a row cell is split by the periodic seam, the pieces must be assembled
before the row transform. A piece shifted by $s\mathbf a_\parallel$ carries the
COMSOL-convention phase

$$
e^{-i k_\parallel\cdot s\mathbf a_\parallel}
$$

for the physical field $E_f$. In the periodic gauge $\widetilde E_f$, this
phase is removed. The $k_\parallel=0$ case is the special case where seam pieces
are glued without any phase.

For the current strip geometry with bulk radius $L=3$, the bulk rows are

$$
n=-3,-2,-1,0,1,2,3,
$$

so the exact finite transform has $N=7$ row wave vectors. This is the correct
specialization of the finite-region theory for the strip. The two-dimensional
hexagonal quotient with $N=3L(L+1)+1=37$ applies to the later full finite
hexagonal region, not to this strip quotient.

## Complete-band picture and the unknown-band problem

Under the finite-region assumption above, every field on the selected cells has
the exact cell-index Fourier expansion

$$
E_{\mathbf R}(\boldsymbol\rho)
=
\frac{1}{\sqrt N}
\sum_{\mathbf k\in K}
e^{-i\mathbf k\cdot\mathbf R}
F(\mathbf k,\boldsymbol\rho).
$$

The forward transform is

$$
F(\mathbf k,\boldsymbol\rho)
=
\frac{1}{\sqrt N}
\sum_{\mathbf R\in S}
e^{+i\mathbf k\cdot\mathbf R}
E_{\mathbf R}(\boldsymbol\rho).
$$

Remove the trivial intra-cell phase by defining

$$
G(\mathbf k,\boldsymbol\rho)
=
e^{+i\mathbf k\cdot\boldsymbol\rho}
F(\mathbf k,\boldsymbol\rho).
$$

Then the exact reconstruction on the selected cells is

$$
E_{\mathbf R}(\boldsymbol\rho)
=
\frac{1}{\sqrt N}
\sum_{\mathbf k\in K}
e^{-i\mathbf k\cdot(\mathbf R+\boldsymbol\rho)}
G(\mathbf k,\boldsymbol\rho).
$$

If the complete infinite-crystal Bloch profiles
$u_b(\mathbf k,\boldsymbol\rho)$ are known for every $\mathbf k\in K$, then the
periodic cell content has the complete band expansion

$$
G(\mathbf k,\boldsymbol\rho)
=
\sum_b
c_b(\mathbf k)
u_b(\mathbf k,\boldsymbol\rho).
$$

Substituting gives the complete Bloch-band representation

$$
E_{\mathbf R}(\boldsymbol\rho)
=
\frac{1}{\sqrt N}
\sum_{\mathbf k\in K}
\sum_b
c_b(\mathbf k)
e^{-i\mathbf k\cdot(\mathbf R+\boldsymbol\rho)}
u_b(\mathbf k,\boldsymbol\rho).
$$

This is the clean complete-dictionary case. The rest of this document explains
what to do when the band profiles are unknown or only partially known.

When no band dictionary is available, do not try to assign band indices first.
Use the extracted functions

$$
G(\mathbf k,\boldsymbol\rho)
$$

themselves. For each lattice wave vector $\mathbf k$,
$G(\mathbf k,\boldsymbol\rho)$ is the periodic unit-cell content of the finite
mode at that $\mathbf k$. It is the key extracted object in the unknown-band
setting.

If only one band/profile contributes appreciably at that $\mathbf k$, then

$$
G(\mathbf k,\boldsymbol\rho)
\approx
c(\mathbf k)u(\mathbf k,\boldsymbol\rho).
$$

In that case $G(\mathbf k,\boldsymbol\rho)$ directly gives the periodic Bloch
profile at that $\mathbf k$, up to an arbitrary complex phase and normalization.

## Exact reconstruction from the extracted components

The inverse cell-index transform is

$$
E_{\mathbf R}(\boldsymbol\rho)
=
\frac{1}{\sqrt N}
\sum_{\mathbf k\in K}
e^{-i\mathbf k\cdot\mathbf R}
F(\mathbf k,\boldsymbol\rho).
$$

Using $F(\mathbf k,\boldsymbol\rho)
=e^{-i\mathbf k\cdot\boldsymbol\rho}G(\mathbf k,\boldsymbol\rho)$ gives

$$
E_{\mathbf R}(\boldsymbol\rho)
=
\frac{1}{\sqrt N}
\sum_{\mathbf k\in K}
e^{-i\mathbf k\cdot(\mathbf R+\boldsymbol\rho)}
G(\mathbf k,\boldsymbol\rho).
$$

Equivalently, on the selected cells,

$$
E_f(\mathbf r)
=
\frac{1}{\sqrt N}
\sum_{\mathbf k\in K}
e^{-i\mathbf k\cdot\mathbf r}
G(\mathbf k,\boldsymbol\rho).
$$

This is already a Bloch-component decomposition. No envelope has been assumed
and no band dictionary has been used.

If only a subset $K_{\rm keep}$ of significant wave vectors is retained, define

$$
E_{\rm rec}(\mathbf r)
=
\frac{1}{\sqrt N}
\sum_{\mathbf k\in K_{\rm keep}}
e^{-i\mathbf k\cdot\mathbf r}
G(\mathbf k,\boldsymbol\rho).
$$

The truncation error is

$$
\epsilon_{\rm direct}
=
\frac{\|E_f-E_{\rm rec}\|_B}{\|E_f\|_B}.
$$

## Per-k weight and normalized profile

Define the Bloch-component weight

$$
W(\mathbf k)
=
\|G(\mathbf k)\|_B^2.
$$

For $W(\mathbf k)>0$, define a normalized per-$\mathbf k$ profile

$$
p_{\mathbf k}(\boldsymbol\rho)
=
\frac{G(\mathbf k,\boldsymbol\rho)}{\sqrt{W(\mathbf k)}}.
$$

Then

$$
\|p_{\mathbf k}\|_B=1.
$$

If $G(\mathbf k,\boldsymbol\rho)=c(\mathbf k)u(\mathbf k,\boldsymbol\rho)$ and
$\|u(\mathbf k)\|_B=1$, then

$$
W(\mathbf k)=|c(\mathbf k)|^2,
$$

and

$$
p_{\mathbf k}(\boldsymbol\rho)
=
e^{i\phi(\mathbf k)}u(\mathbf k,\boldsymbol\rho),
$$

where $e^{i\phi(\mathbf k)}$ is an arbitrary gauge phase. The phase split
between $c(\mathbf k)$ and $u(\mathbf k,\boldsymbol\rho)$ is not unique.

This gauge ambiguity is harmless for profile clustering if the similarity uses
the absolute value of the inner product.

## Gauge fixing and complex coefficients

The factorization

$$
G(\mathbf k,\boldsymbol\rho)
\approx
c(\mathbf k)u(\mathbf k,\boldsymbol\rho)
$$

has a phase gauge freedom:

$$
c(\mathbf k)u(\mathbf k,\boldsymbol\rho)
=
\left[c(\mathbf k)e^{-i\chi(\mathbf k)}\right]
\left[e^{+i\chi(\mathbf k)}u(\mathbf k,\boldsymbol\rho)\right].
$$

Therefore the complex coefficient $c(\mathbf k)$ is not uniquely defined until a
gauge is chosen for $u(\mathbf k,\boldsymbol\rho)$.

One simple per-$\mathbf k$ convention is

$$
u(\mathbf k,\boldsymbol\rho)=p_{\mathbf k}(\boldsymbol\rho),
\qquad
c(\mathbf k)=\sqrt{W(\mathbf k)}.
$$

This reconstructs $G(\mathbf k,\boldsymbol\rho)$ exactly, but it stores all
phase information in $u(\mathbf k,\boldsymbol\rho)$ rather than in
$c(\mathbf k)$.

For a cluster of similar profiles, a more useful convention is to choose a
reference profile $u_{\rm ref}$ and align each profile so that

$$
\langle u_{\rm ref},u(\mathbf k)\rangle_B
$$

is real and positive. After this alignment, the complex coefficient is

$$
c(\mathbf k)
=
\langle u(\mathbf k),G(\mathbf k)\rangle_B,
\qquad
\|u(\mathbf k)\|_B=1.
$$

If the cluster is rank 1, one can instead use the cluster SVD to choose a common
profile $u_\alpha$ and define

$$
A_\alpha(\mathbf k)
=
\langle u_\alpha,G(\mathbf k)\rangle_B.
$$

These $A_\alpha(\mathbf k)$ are the complex coefficients used for reconstruction
and, optionally, for forming an envelope after choosing a representative carrier.

## Profile similarity and automatic clustering

Compare two retained wave vectors using

$$
S(\mathbf k,\mathbf k')
=
\left|
\langle p_{\mathbf k},p_{\mathbf k'}\rangle_B
\right|^2.
$$

Because both profiles are normalized,

$$
0\le S(\mathbf k,\mathbf k')\le1.
$$

If

$$
S(\mathbf k,\mathbf k')\approx1,
$$

then the two wave vectors have nearly the same periodic cell profile, up to a
global complex phase. They likely belong to the same Bloch-profile family or
band-like branch.

If

$$
S(\mathbf k,\mathbf k')\ll1,
$$

then the corresponding periodic profiles are different, and the two components
should not be compressed into one common cell profile.

A practical clustering procedure is:

1. Keep only wave vectors with $W(\mathbf k)$ above a noise threshold.
2. Build the similarity matrix $S(\mathbf k,\mathbf k')$.
3. Build a graph with an edge when $S(\mathbf k,\mathbf k')>S_{\rm min}$.
4. Use connected components or hierarchical clustering to obtain clusters
   $\mathcal C_\alpha$.

This clustering is based on unit-cell profile similarity, not on distance in
$\mathbf k$ space. Nearby $\mathbf k$ points may belong to different profile
families, and far-apart $\mathbf k$ points may still have the same periodic
profile symmetry or band character.

## Cluster-level decomposition

For a cluster $\mathcal C_\alpha$, the partial field is

$$
E_\alpha(\mathbf r)
=
\frac{1}{\sqrt N}
\sum_{\mathbf k\in\mathcal C_\alpha}
e^{-i\mathbf k\cdot\mathbf r}
G(\mathbf k,\boldsymbol\rho).
$$

The cluster weight is

$$
W_\alpha
=
\sum_{\mathbf k\in\mathcal C_\alpha}
W(\mathbf k).
$$

The weight fraction is

$$
\eta_\alpha
=
\frac{W_\alpha}{\sum_{\mathbf k}W(\mathbf k)}.
$$

Within each cluster, perform a $B$-weighted SVD of the functions
$G(\mathbf k,\boldsymbol\rho)$:

$$
G(\mathbf k,\boldsymbol\rho)
\approx
\sum_s
A_{\alpha s}(\mathbf k)
u_{\alpha s}(\boldsymbol\rho),
\qquad
\mathbf k\in\mathcal C_\alpha.
$$

The rank-$r$ retained weight is

$$
\eta_{\alpha,r}
=
\frac{
\sum_{s=1}^{r}\sigma_{\alpha s}^2
}{
\sum_s\sigma_{\alpha s}^2
}.
$$

If $\eta_{\alpha,1}\approx1$, then all significant $\mathbf k$ components in
the cluster share one common periodic profile. If multiple singular values are
important, the cluster represents a low-dimensional unit-cell subspace rather
than a single profile.

The rank-$r$ cluster reconstruction is

$$
E_{\alpha,r}(\mathbf r)
=
\frac{1}{\sqrt N}
\sum_{\mathbf k\in\mathcal C_\alpha}
\sum_{s=1}^{r}
e^{-i\mathbf k\cdot\mathbf r}
A_{\alpha s}(\mathbf k)
u_{\alpha s}(\boldsymbol\rho).
$$

The full finite-mode reconstruction is the sum over clusters:

$$
E_{\rm rec}(\mathbf r)
=
\sum_\alpha E_{\alpha,r_\alpha}(\mathbf r).
$$

## Discussion: mixed profiles and output interpretation

The extraction does not require every $\mathbf k$ component to be a pure
single-band Bloch mode. The primary extracted object is always
$G(\mathbf k,\boldsymbol\rho)$. Clustering and SVD determine the appropriate
level of description.

### Fixed mixture as one composite profile

Suppose two profiles contribute with a fixed relative amplitude and phase:

$$
G(\mathbf k,\boldsymbol\rho)
=
A(\mathbf k)
\left[
\alpha u_1(\boldsymbol\rho)
+
\beta u_2(\boldsymbol\rho)
\right].
$$

Then all retained $\mathbf k$ components share the same composite profile

$$
u_{\rm mix}(\boldsymbol\rho)
=
\alpha u_1(\boldsymbol\rho)+\beta u_2(\boldsymbol\rho).
$$

The similarity clustering will likely group these $\mathbf k$ values together,
and the cluster SVD will be rank 1:

$$
\eta_{\alpha,1}\approx1.
$$

This is a valid output. The algorithm has found a stable composite unit-cell
profile. Later infinite-crystal validation may reveal that this composite
profile overlaps with several true infinite-crystal modes, but the finite-mode
decomposition itself is still meaningful.

### Varying mixture as a rank-r subspace

If the relative amplitude or phase varies with $\mathbf k$,

$$
G(\mathbf k,\boldsymbol\rho)
=
c_1(\mathbf k)u_1(\boldsymbol\rho)
+
c_2(\mathbf k)u_2(\boldsymbol\rho),
$$

then the normalized apparent profile

$$
p_{\mathbf k}(\boldsymbol\rho)
=
\frac{G(\mathbf k,\boldsymbol\rho)}{\|G(\mathbf k)\|_B}
$$

can rotate inside the subspace spanned by $u_1$ and $u_2$. If the clustering
keeps these $\mathbf k$ values together, the cluster SVD will show rank 2:

$$
\eta_{\alpha,1}<1,
\qquad
\eta_{\alpha,2}\approx1.
$$

In that case, the correct output is not one profile but a two-dimensional
profile subspace:

$$
G(\mathbf k,\boldsymbol\rho)
\approx
A_{\alpha1}(\mathbf k)u_{\alpha1}(\boldsymbol\rho)
+
A_{\alpha2}(\mathbf k)u_{\alpha2}(\boldsymbol\rho).
$$

More generally, significant singular values indicate the dimension of the
profile subspace needed to represent the data. The SVD basis
$u_{\alpha s}$ is a data-driven orthonormal basis for this subspace. It need not
be identical to the true infinite-crystal band basis; it may differ by a unitary
rotation within the same subspace. Infinite-crystal modes, operator residuals,
or symmetry labels can be used afterward to rotate or identify the physical
band basis.

### Strongly varying mixtures as per-k profiles

If the apparent profiles vary too strongly, similarity clustering may split the
data into small clusters or even one cluster per $\mathbf k$. This is also a
valid output. It means the finite mode should be described at the per-$\mathbf k$
level:

$$
E_f(\mathbf r)
\approx
\frac{1}{\sqrt N}
\sum_{\mathbf k}
e^{-i\mathbf k\cdot\mathbf r}
G(\mathbf k,\boldsymbol\rho),
$$

without compressing several $\mathbf k$ values into one shared profile family.

For a singleton cluster, SVD cannot separate multiple hidden Bloch modes at that
same $\mathbf k$; it only returns the extracted vector
$G(\mathbf k,\boldsymbol\rho)$ itself. If further band separation is needed,
project $G(\mathbf k,\boldsymbol\rho)$ onto infinite-crystal modes or analyze
multiple finite eigenmodes jointly.

Thus a fragmented clustering result is not necessarily a failure. It says that
the apparent unit-cell profile is strongly $\mathbf k$ dependent, so the most
faithful description is the uncompressed per-$\mathbf k$ decomposition.

## Relation to the envelope picture

The envelope picture is optional. It should be used only after the
Bloch-component decomposition shows that a cluster is close to rank 1.

Suppose cluster $\mathcal C_\alpha$ has

$$
G(\mathbf k,\boldsymbol\rho)
\approx
A_\alpha(\mathbf k)u_\alpha(\boldsymbol\rho),
\qquad
\mathbf k\in\mathcal C_\alpha.
$$

Then

$$
E_\alpha(\mathbf r)
\approx
\frac{1}{\sqrt N}
\sum_{\mathbf k\in\mathcal C_\alpha}
A_\alpha(\mathbf k)
e^{-i\mathbf k\cdot\mathbf r}
u_\alpha(\boldsymbol\rho).
$$

Choose a representative carrier $\mathbf k_\alpha$. Define

$$
\psi_\alpha(\mathbf r)
=
\frac{1}{\sqrt N}
\sum_{\mathbf k\in\mathcal C_\alpha}
A_\alpha(\mathbf k)
e^{-i(\mathbf k-\mathbf k_\alpha)\cdot\mathbf r}.
$$

Then

$$
E_\alpha(\mathbf r)
\approx
e^{-i\mathbf k_\alpha\cdot\mathbf r}
u_\alpha(\boldsymbol\rho)
\psi_\alpha(\mathbf r).
$$

Thus a carrier-envelope form is a compression of the Bloch-component
decomposition. It is not an input assumption.

The envelope is slow only if the relative spectrum
$\mathbf q=\mathbf k-\mathbf k_\alpha$ is narrow. A useful diagnostic is

$$
\delta_{\alpha,j}^2
=
\frac{
\sum_{\mathbf k\in\mathcal C_\alpha}
|A_\alpha(\mathbf k)|^2
|e^{-i(\mathbf k-\mathbf k_\alpha)\cdot\mathbf a_j}-1|^2
}{
\sum_{\mathbf k\in\mathcal C_\alpha}
|A_\alpha(\mathbf k)|^2
}.
$$

If $\delta_{\alpha,j}\ll1$, the envelope varies slowly along lattice vector
$\mathbf a_j$. If not, the Bloch-component decomposition may still be valid,
but the slow-envelope interpretation is not.

## What this algorithm can output

For each finite eigenmode, the algorithm can report:

1. finite eigenfrequency $\omega_f$;
2. selected bulk-core cell region $S$;
3. quotient/spectral-basis check for $S$;
4. lattice Fourier wave-vector set $K$ and Parseval error;
5. Bloch-component weight map $W(\mathbf k)$;
6. normalized extracted profiles $p_{\mathbf k}(\boldsymbol\rho)$;
7. profile similarity matrix $S(\mathbf k,\mathbf k')$;
8. automatically detected profile clusters $\mathcal C_\alpha$;
9. cluster weights $W_\alpha$ and fractions $\eta_\alpha$;
10. cluster SVD spectra $\sigma_{\alpha s}$;
11. retained rank $r_\alpha$ for each cluster;
12. cluster profiles $u_{\alpha s}(\boldsymbol\rho)$;
13. complex coefficients $A_{\alpha s}(\mathbf k)$;
14. partial reconstructions $E_{\alpha,r_\alpha}$;
15. total reconstruction error $\epsilon_{\rm direct}$;
16. optional known-band weights and coefficients $c_j(\mathbf k)$;
17. optional known, tested-residual, and untested weight fractions;
18. optional envelope $\psi_\alpha$ for rank-1 clusters;
19. optional infinite-crystal validation overlaps.

## Validation with infinite-crystal modes

After extraction, solve a small number of infinite-crystal modes near each
important $\mathbf k$ and frequency $\omega_f$.

For a normalized extracted profile $p_{\mathbf k}$ and an infinite-crystal
periodic profile $u_j^\infty(\mathbf k)$, compute

$$
\mathcal O_j(\mathbf k)
=
\frac{
|\langle p_{\mathbf k},u_j^\infty(\mathbf k)\rangle_B|^2
}{
\|p_{\mathbf k}\|_B^2
\|u_j^\infty(\mathbf k)\|_B^2
}.
$$

If $\mathcal O_j(\mathbf k)\approx1$, then the extracted per-$\mathbf k$
profile agrees with a true infinite-crystal Bloch mode.

For a cluster subspace $U_{\rm ext}$ and an infinite-crystal subspace
$U_{\rm inf}$, first $B$-orthonormalize both sets of columns. Then compute

$$
C=U_{\rm ext}^\dagger \mathcal B U_{\rm inf}.
$$

The singular values of $C$ are the cosines of the principal angles between the
two subspaces.

## Known infinite-crystal modes

The blind decomposition above does not require any known band structure. If
some infinite-crystal Bloch profiles are known, use them as an optional
dictionary after extracting

$$
G(\mathbf k,\boldsymbol\rho).
$$

Let $K$ be the finite quotient wave-vector set used for the extraction. Suppose
the known data are a set of specific mode entries

$$
\mathcal D_{\rm known}
=
\{(j,\mathbf k): u_j^\infty(\mathbf k,\boldsymbol\rho)
\text{ is known}\}.
$$

No Cartesian-product assumption is made. Knowing $(j_1,\mathbf k_1)$ and
$(j_2,\mathbf k_2)$ does not imply that $(j_1,\mathbf k_2)$ or
$(j_2,\mathbf k_1)$ is known.

For each finite-grid wave vector, define the available known band indices

$$
\mathcal J_{\rm known}(\mathbf k)
=
\{j:(j,\mathbf k)\in\mathcal D_{\rm known}\}.
$$

### Complete known-band setup

In the complete setup, the known dictionary contains the full infinite-crystal
cell eigenbasis for every $\mathbf k\in K$ in the chosen field space. Then the
band expansion of $G(\mathbf k,\boldsymbol\rho)$ is exact.

If the known profiles are $B$-orthonormal, compute

$$
c_j(\mathbf k)
=
\langle u_j^\infty(\mathbf k),G(\mathbf k)\rangle_B,
$$

and

$$
G_{\rm known}(\mathbf k,\boldsymbol\rho)
=
\sum_j
c_j(\mathbf k)
u_j^\infty(\mathbf k,\boldsymbol\rho).
$$

The band-resolved weight is

$$
W_j(\mathbf k)=|c_j(\mathbf k)|^2.
$$

The finite mode can then be reconstructed as

$$
E_{\rm rec}(\mathbf r)
=
\frac{1}{\sqrt N}
\sum_{\mathbf k\in K}
\sum_j
c_j(\mathbf k)
e^{-i\mathbf k\cdot\mathbf r}
u_j^\infty(\mathbf k,\boldsymbol\rho).
$$

With the full basis this reconstruction equals the finite-mode data on the
selected cells. In practice one often keeps only a finite band window. Then the
residual

$$
G_{\rm res}(\mathbf k)
=
G(\mathbf k)-G_{\rm known}(\mathbf k)
$$

measures the part not represented by the retained bands. A large residual means
that the band window is too small, the mode normalization or inner product is
inconsistent, the field spaces do not match, or the selected finite region
contains non-bulk content.

### Sparse known-pair setup

In the sparse setup, only some specific $(j,\mathbf k)$ entries are known. For
each $\mathbf k$, use only the modes that actually exist in
$\mathcal J_{\rm known}(\mathbf k)$.

If $\mathcal J_{\rm known}(\mathbf k)$ is nonempty, form the known subspace

$$
U_{\rm known}(\mathbf k)
=
\operatorname{span}
\{u_j^\infty(\mathbf k):j\in\mathcal J_{\rm known}(\mathbf k)\}.
$$

Let $P_{\rm known}(\mathbf k)$ be the $B$-orthogonal projector onto this
subspace. Then

$$
G_{\rm known}(\mathbf k)
=
P_{\rm known}(\mathbf k)G(\mathbf k),
$$

and

$$
G_{\rm res}(\mathbf k)
=
\left[I-P_{\rm known}(\mathbf k)\right]G(\mathbf k).
$$

Here $G_{\rm res}$ means "tested at this $\mathbf k$, but not explained by the
known modes at this $\mathbf k$."

If $\mathcal J_{\rm known}(\mathbf k)$ is empty, do not call the result a
residual. There was no known mode to test against. Instead define

$$
G_{\rm untested}(\mathbf k)=G(\mathbf k).
$$

A simple weight summary is

$$
W_{\rm known}
=
\sum_{\mathbf k:\mathcal J_{\rm known}(\mathbf k)\ne\emptyset}
\|G_{\rm known}(\mathbf k)\|_B^2,
$$

$$
W_{\rm res}
=
\sum_{\mathbf k:\mathcal J_{\rm known}(\mathbf k)\ne\emptyset}
\|G_{\rm res}(\mathbf k)\|_B^2,
$$

and

$$
W_{\rm untested}
=
\sum_{\mathbf k:\mathcal J_{\rm known}(\mathbf k)=\emptyset}
\|G(\mathbf k)\|_B^2.
$$

These three numbers should be reported separately. $W_{\rm known}$ is the part
explained by the known $(j,\mathbf k)$ entries. $W_{\rm res}$ is the part at
tested $\mathbf k$ values that those entries did not explain.
$W_{\rm untested}$ is the part at $\mathbf k$ values where no known mode was
available, so no band assignment should be made.

Normalize these by

$$
W_{\rm total}
=
\sum_{\mathbf k\in K}\|G(\mathbf k)\|_B^2
$$

to obtain $\eta_{\rm known}$, $\eta_{\rm res}$, and
$\eta_{\rm untested}$. In the complete known-band setup,
$\eta_{\rm untested}=0$. If all bands in the chosen field space are retained,
$\eta_{\rm res}$ is zero up to numerical error. With a finite band window,
$\eta_{\rm res}$ is the truncated-away content.

If the known modes at a given $\mathbf k$ are not exactly $B$-orthonormal, use
the subspace projector instead of summing $|c_j(\mathbf k)|^2$ directly.

## Important limitations

This method extracts the Bloch-like content present in one finite eigenmode. It
does not always uniquely identify every underlying infinite-crystal band.

The exact Fourier version assumes the selected cell set is a complete
representative set of a finite lattice quotient. If the selected cells are an
arbitrary mask, have missing cells, or otherwise do not support an exact finite
Fourier basis, use a Gram-matrix or least-squares version of the cell-index
decomposition.

The complete Bloch-band expansion is a representation statement. It does not
mean that a mode of an open, defective, or boundary-terminated finite system is
itself an eigenmode of the infinite periodic system. It only means that the
finite data on the selected cells can be expanded in that basis.

The main failure mode occurs when multiple Bloch profiles contribute at the same
wave vector:

$$
G(\mathbf k,\boldsymbol\rho)
=
c_1(\mathbf k)u_1(\mathbf k,\boldsymbol\rho)
+
c_2(\mathbf k)u_2(\mathbf k,\boldsymbol\rho)
+\cdots.
$$

From a single finite eigenmode at that same $\mathbf k$, this mixture is not
uniquely separable. The extracted object is the linear combination
$G(\mathbf k,\boldsymbol\rho)$.

To separate multiple true bands at the same $\mathbf k$, one may need:

1. several finite eigenmodes analyzed jointly;
2. an operator residual using the unit-cell Maxwell operator;
3. direct infinite-crystal validation modes;
4. symmetry labels or additional physical constraints.

Other practical limitations:

1. If only a finite band window or sparse $(j,\mathbf k)$ dictionary is known,
   the known-mode projection is not complete. Report residual and untested
   weights separately.
2. If $W(\mathbf k)$ is very small, $p_{\mathbf k}$ is noise-sensitive.
3. If comparing $\mathbf k$ points across Brillouin-zone boundaries, reciprocal
   lattice gauge conventions must be aligned.
4. If boundary or defect cells contaminate the selected region, the extracted
   profiles may not be bulk-like.

## Minimal pseudo-code

```text
given finite mode E_f(r)

choose bulk-like complete cells S
verify S is a finite lattice quotient representative set
compute the corresponding discrete wave-vector set K

sample each cell field E_R(rho) = E_f(R + rho)

for each k in K:
    F_k(rho) = (1/sqrt(N)) * sum_R exp(+i k dot R) * E_R(rho)
    G_k(rho) = exp(+i k dot rho) * F_k(rho)
    W_k = ||G_k||_B^2

if full infinite-crystal band dictionary is available:
    for each k in K:
        project G_k onto all band profiles u_j(k, rho)
        c_j(k) = <u_j(k), G_k>_B
    reconstruct from exp(-i k dot r) * c_j(k) * u_j(k, rho)
    report band-resolved weights
    if using only a finite band window, also report truncation residual

else if sparse known pairs (j, k) are available:
    for each k in K:
        if known modes exist at k:
            project G_k onto their B-orthogonal subspace
            store known contribution and tested residual
        else:
            store G_k as untested content
    report known, tested-residual, and untested weights
    optionally run blind analysis on residual or untested content

else:
    keep K_keep = {k : W_k > threshold}

    for k in K_keep:
        p_k(rho) = G_k(rho) / sqrt(W_k)

    build similarity matrix:
        S(k, k') = |<p_k, p_k'>_B|^2

    cluster k values using S(k, k')

    for each cluster C_alpha:
        compute weight W_alpha = sum_{k in C_alpha} W_k
        perform B-weighted SVD of {G_k(rho) : k in C_alpha}
        choose rank r_alpha from singular values
        reconstruct cluster contribution:
            E_alpha(r) =
                (1/sqrt(N)) *
                sum_{k in C_alpha}
                sum_{s <= r_alpha}
                exp(-i k dot r) * A_alpha_s(k) * u_alpha_s(rho)

    sum cluster reconstructions:
        E_rec(r) = sum_alpha E_alpha(r)

    compute direct reconstruction error:
        epsilon = ||E_f - E_rec||_B / ||E_f||_B

    if a cluster is rank 1:
        optionally choose representative carrier k_alpha
        optionally form envelope from A_alpha(k)
```

## Conceptual summary

The exact version starts from a finite cell set $S$ that is a complete
representative set of a finite lattice quotient. This gives a unitary finite
Fourier basis with discrete wave vectors $K$.

If the full infinite-crystal band dictionary is known at those $\mathbf k$
values, the finite-region field has a complete Bloch-band expansion. If the
dictionary is not known, the fundamental extracted object is

$$
G(\mathbf k,\boldsymbol\rho),
$$

the periodic unit-cell content of the finite mode at lattice wave vector
$\mathbf k$.

On the selected cells, the finite mode is reconstructed exactly as

$$
E_f(\mathbf r)
=
\frac{1}{\sqrt N}
\sum_{\mathbf k\in K}
e^{-i\mathbf k\cdot\mathbf r}
G(\mathbf k,\boldsymbol\rho).
$$

Known band profiles can be projected onto $G(\mathbf k,\boldsymbol\rho)$. With
no band dictionary, profile similarity and SVD discover which $\mathbf k$
components share a unit-cell profile family. Only after that, if a cluster is
rank 1 and its relative $\mathbf k$ spectrum is narrow, does the usual
carrier-envelope picture become a useful compressed description.
