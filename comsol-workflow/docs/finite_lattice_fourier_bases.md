# Orthogonal Fourier Bases on Finite Lattice Regions

This note explains when a finite set of lattice points admits an exact
orthogonal basis made of plane waves.

No physics background is required. The only objects are:

1. an infinite lattice;
2. a finite subset of lattice points;
3. a set of wave vectors used to analyze data on that subset.

## Basic problem

A $d$-dimensional lattice is the set

$$
\Lambda
=
\left\{
\sum_{i=1}^{d}m_i\mathbf a_i:
m_i\in\mathbb Z
\right\},
$$

where $\mathbf a_i$ are primitive lattice vectors.

The reciprocal vectors $\mathbf b_i$ satisfy

$$
\mathbf a_i\cdot\mathbf b_j=2\pi\delta_{ij}.
$$

## Sign convention

This note uses the COMSOL-compatible Bloch sign convention. A component with
wave vector $\mathbf k$ has spatial phase

$$
e^{-i\mathbf k\cdot\mathbf R}.
$$

The forward analysis transform therefore uses the opposite phase,

$$
e^{+i\mathbf k\cdot\mathbf R},
$$

so that a signal proportional to $e^{-i\mathbf k_0\cdot\mathbf R}$ is detected
at $\mathbf k=\mathbf k_0$.

With the common frequency-domain time convention
$\operatorname{Re}\{\widetilde f e^{+i\omega t}\}$, this is the same sign used
for COMSOL Bloch wave vectors. If a different convention writes Bloch waves as
$e^{+i\mathbf k\cdot\mathbf R}$, replace $\mathbf k$ by $-\mathbf k$ throughout.

Choose a finite set of lattice points

$$
S=\{\mathbf R_1,\ldots,\mathbf R_N\}\subset\Lambda.
$$

For a wave vector $\mathbf k$, define the plane wave on $S$ by

$$
\phi_{\mathbf k}(\mathbf R)
=
\frac{1}{\sqrt N}e^{-i\mathbf k\cdot\mathbf R},
\qquad
\mathbf R\in S.
$$

We want a set $K$ of $N$ wave vectors such that

$$
\frac{1}{N}
\sum_{\mathbf R\in S}
e^{i(\mathbf k-\mathbf k')\cdot\mathbf R}
=
\delta_{\mathbf k,\mathbf k'},
\qquad
\mathbf k,\mathbf k'\in K.
$$

If this holds, then the plane waves form an orthonormal basis for data sampled
on $S$.

For scalar data $f(\mathbf R)$, the exact finite Fourier transform is

$$
\widehat f(\mathbf k)
=
\frac{1}{\sqrt N}
\sum_{\mathbf R\in S}
e^{+i\mathbf k\cdot\mathbf R}
f(\mathbf R),
$$

and the inverse is

$$
f(\mathbf R)
=
\frac{1}{\sqrt N}
\sum_{\mathbf k\in K}
e^{-i\mathbf k\cdot\mathbf R}
\widehat f(\mathbf k).
$$

The same transform applies componentwise if $f(\mathbf R)$ is a vector or any
other object attached to each lattice point.

## Spectral pairs

Define the structure factor of the finite set $S$:

$$
D_S(\mathbf q)
=
\sum_{\mathbf R\in S}
e^{i\mathbf q\cdot\mathbf R}.
$$

The orthogonality condition is equivalent to

$$
D_S(\mathbf k-\mathbf k')=0
\quad
\text{for all distinct }
\mathbf k,\mathbf k'\in K.
$$

A pair $(S,K)$ satisfying this condition is called a spectral pair.

Not every finite set $S$ has an obvious plane-wave spectrum $K$. An arbitrary
finite-dimensional data space always has some orthonormal basis, but it may not
be made of lattice plane waves.

The cleanest guaranteed construction comes from finite lattice quotients.

## Meaning of the quotient language

The infinite lattice $\Lambda$ is the full set of allowed grid points. The
finite set $S$ is the region where data are sampled.

A sublattice is a coarser lattice inside $\Lambda$. For example,
$\langle 4\mathbf a_1,3\mathbf a_2\rangle$ is a sublattice of
$\langle\mathbf a_1,\mathbf a_2\rangle$.

A quotient means that points separated by a sublattice vector are identified as
the same point. If

$$
\mathbf R-\mathbf R'\in\Lambda_{\rm super},
$$

then $\mathbf R$ and $\mathbf R'$ represent the same quotient element. This is
the lattice version of reducing an integer modulo $N$.

A representative set chooses one point from each quotient class. In one
dimension, for $\mathbb Z/N\mathbb Z$, the usual representatives are

$$
\{0,1,\ldots,N-1\}.
$$

In two dimensions, a representative set may look like a parallelogram, a
hexagon, or another shape. The shape is secondary. The important property is
that exactly one point is chosen from each quotient class.

The quotient is a finite Abelian group because points can be added and then
reduced modulo the sublattice. "Abelian" only means that the order of addition
does not matter.

A character of this finite group is a function $\chi$ satisfying

$$
\chi(\mathbf R+\mathbf R')=\chi(\mathbf R)\chi(\mathbf R').
$$

Plane waves are characters, since

$$
e^{-i\mathbf k\cdot(\mathbf R+\mathbf R')}
=
e^{-i\mathbf k\cdot\mathbf R}
e^{-i\mathbf k\cdot\mathbf R'}.
$$

The useful fact is that distinct characters of a finite Abelian group are
orthogonal. This is the reason quotient representative sets have exact Fourier
bases.

## Finite quotient construction

Choose a finite-index sublattice

$$
\Lambda_{\rm super}\subset\Lambda.
$$

The quotient

$$
G=\Lambda/\Lambda_{\rm super}
$$

is a finite Abelian group. Its size is the index

$$
N=[\Lambda:\Lambda_{\rm super}].
$$

Now choose $S$ to be one complete set of representatives of this quotient. That
means every lattice point can be written uniquely as

$$
\mathbf R=\mathbf R_S+\mathbf T,
\qquad
\mathbf R_S\in S,
\qquad
\mathbf T\in\Lambda_{\rm super}.
$$

Equivalently, translates of $S$ by $\Lambda_{\rm super}$ tile the original
lattice without overlap.

The dual wave vectors are those that are periodic on the superlattice:

$$
e^{i\mathbf k\cdot\mathbf T}=1
\qquad
\text{for every }
\mathbf T\in\Lambda_{\rm super}.
$$

Wave vectors that differ by an original reciprocal lattice vector are
equivalent, because they give the same values on $\Lambda$. Therefore the
finite spectrum is

$$
K=\Lambda_{\rm super}^*/\Lambda^*,
$$

where

$$
\Lambda_{\rm super}^*
=
\{\mathbf k:\mathbf k\cdot\mathbf T\in2\pi\mathbb Z
\text{ for all }\mathbf T\in\Lambda_{\rm super}\}.
$$

This set has exactly $N$ inequivalent wave vectors.

The orthogonality follows from the standard character identity on a finite
group:

$$
\sum_{g\in G}\chi(g)=0
$$

for every nontrivial character $\chi$. The plane waves restricted to $S$ are
exactly the characters of $G$.

Thus every complete representative set of a finite lattice quotient has an
exact orthogonal plane-wave basis.

## Coordinate formula

In lattice coordinates, write a point as an integer vector

$$
\mathbf m=(m_1,\ldots,m_d)\in\mathbb Z^d,
$$

meaning

$$
\mathbf R_{\mathbf m}
=
\sum_i m_i\mathbf a_i.
$$

Let the superlattice be generated by an integer full-rank matrix $M$:

$$
\mathbf T_j
=
\sum_i M_{ij}\mathbf a_i.
$$

Then

$$
N=|\det M|.
$$

The reciprocal coordinates of allowed wave vectors are

$$
\boldsymbol\xi
=
M^{-T}\mathbf h
\pmod{\mathbb Z^d},
\qquad
\mathbf h\in\mathbb Z^d,
$$

and

$$
\mathbf k_{\mathbf h}
=
\sum_i \xi_i\mathbf b_i.
$$

There are $|\det M|$ distinct choices modulo the original reciprocal lattice.

## What the wave vectors look like

The allowed wave vectors are discrete points on the reciprocal torus. They are
defined modulo the original reciprocal lattice $\Lambda^*$, because
$\mathbf k$ and $\mathbf k+\mathbf G$ give the same plane wave on $\Lambda$
when $\mathbf G\in\Lambda^*$.

One may choose representatives of these discrete wave vectors inside the first
Brillouin zone. This is only a convention. Points on the boundary may have
several equivalent representatives.

For a rectangular product block, the points form a rectangular grid in
reciprocal-lattice coordinates. For a skew supercell, they form a skew finite
grid. For a non-parallelogram representative set, such as the hexagonal example
below, the same $N$ wave vectors still exist, but after folding into the first
Brillouin zone they may not look like a simple rectangular mesh.

## Completeness on the finite region

Assume in this section that $S$ is a complete representative set of a finite
lattice quotient. Then the corresponding plane waves are not only orthogonal.
They are complete for scalar data on $S$.

Let

$$
U_{\mathbf R,\mathbf k}
=
\frac{1}{\sqrt N}e^{-i\mathbf k\cdot\mathbf R},
\qquad
\mathbf R\in S,
\qquad
\mathbf k\in K.
$$

Orthogonality says

$$
U^\dagger U=I.
$$

Since $U$ is an $N\times N$ matrix, this also implies

$$
UU^\dagger=I.
$$

Thus $U$ is unitary. Therefore every scalar function $f(\mathbf R)$ on $S$ has
an exact expansion

$$
f(\mathbf R)
=
\frac{1}{\sqrt N}
\sum_{\mathbf k\in K}
e^{-i\mathbf k\cdot\mathbf R}
\widehat f(\mathbf k).
$$

Now suppose each lattice point carries an internal object, such as a vector or a
function of an intra-cell coordinate $\boldsymbol\rho$:

$$
f_{\mathbf R}(\boldsymbol\rho).
$$

The same finite Fourier transform gives

$$
f_{\mathbf R}(\boldsymbol\rho)
=
\frac{1}{\sqrt N}
\sum_{\mathbf k\in K}
e^{-i\mathbf k\cdot\mathbf R}
F(\mathbf k,\boldsymbol\rho).
$$

This is still exact. The Fourier basis handles the finite cell index
$\mathbf R$, while $F(\mathbf k,\boldsymbol\rho)$ keeps all internal degrees of
freedom.

Equivalently, the finite-region space factors as

$$
\ell^2(S)\otimes\mathcal H_{\rm cell}
\cong
\bigoplus_{\mathbf k\in K}\mathcal H_{\rm cell}.
$$

Here $\mathcal H_{\rm cell}$ is the space of internal data attached to one
lattice point.

Now consider a standard self-adjoint periodic eigenproblem on one cell at fixed
$\mathbf k$. The spectral theorem gives a complete orthonormal eigenbasis for
the cell space:

$$
\{\varphi_b(\mathbf k,\boldsymbol\rho)\}_b.
$$

In a finite-dimensional discretization, this is the familiar statement that the
full set of eigenvectors of a Hermitian generalized eigenproblem

$$
A(\mathbf k)u=\lambda B(\mathbf k)u,
\qquad
B(\mathbf k)>0,
$$

forms a complete $B$-orthonormal basis.

Therefore every internal function has the exact expansion

$$
F(\mathbf k,\boldsymbol\rho)
=
\sum_b c_b(\mathbf k)\varphi_b(\mathbf k,\boldsymbol\rho),
$$

and hence

$$
f_{\mathbf R}(\boldsymbol\rho)
=
\frac{1}{\sqrt N}
\sum_{\mathbf k\in K}
\sum_b
c_b(\mathbf k)
e^{-i\mathbf k\cdot\mathbf R}
\varphi_b(\mathbf k,\boldsymbol\rho).
$$

In a Bloch-wave convention, the internal basis function is often split into an
intra-cell plane-wave phase and a periodic profile:

$$
\varphi_b(\mathbf k,\boldsymbol\rho)
=
e^{-i\mathbf k\cdot\boldsymbol\rho}
u_b(\mathbf k,\boldsymbol\rho).
$$

Substituting this definition gives

$$
f(\mathbf R+\boldsymbol\rho)
=
\frac{1}{\sqrt N}
\sum_{\mathbf k\in K}
\sum_b
c_b(\mathbf k)
e^{-i\mathbf k\cdot(\mathbf R+\boldsymbol\rho)}
u_b(\mathbf k,\boldsymbol\rho).
$$

The extra factor $e^{-i\mathbf k\cdot\boldsymbol\rho}$ is therefore not an
additional assumption. It is part of the definition of the internal basis when
one chooses to write it in terms of a periodic Bloch profile $u_b$.

Thus the discrete $\mathbf k$ values associated with the finite quotient,
together with all bands $b$ at those $\mathbf k$ values, are sufficient to
represent arbitrary data on the finite region.

No extra completeness assumption is being added here. The cell-index
completeness comes from the unitary finite Fourier transform, and the band
completeness comes from using the full eigenbasis of the cell eigenproblem.
Keeping only a finite band window is a truncation, and then the representation
is approximate.

This is a statement about representation completeness. It does not mean that a
mode of an open, defective, or boundary-terminated finite system is itself an
eigenmode of the infinite periodic system. It only means that the finite data
can be expanded in this basis, and the coefficients describe how much of each
basis component is present.

## Efficient implementation

The quotient construction also explains when a fast transform exists.

For a general superlattice matrix $M$, compute a Smith normal form

$$
U M V
=
\operatorname{diag}(d_1,\ldots,d_d),
$$

where $U$ and $V$ are integer matrices with determinant $\pm1$. Then

$$
\Lambda/\Lambda_{\rm super}
\cong
\mathbb Z_{d_1}\times\cdots\times\mathbb Z_{d_d}.
$$

After reindexing the representative points by this product of cyclic groups,
the exact Fourier transform is an ordinary multidimensional FFT with sizes
$d_1,\ldots,d_d$. Factors with $d_i=1$ can be ignored.

Special cases:

1. A rectangular $M_1\times M_2$ block gives a standard 2D FFT.
2. A cyclic quotient $\mathbb Z_N$ gives a standard 1D FFT.
3. If the selected set is not a quotient representative set, there is no exact
   quotient FFT. Use a dense transform, Gram matrix, or least-squares fit.

## Example 1: parallelogram block

In two dimensions, choose

$$
\Lambda_{\rm super}
=
\langle M_1\mathbf a_1,\ M_2\mathbf a_2\rangle.
$$

A natural representative set is

$$
S
=
\{
m\mathbf a_1+n\mathbf a_2:
0\le m<M_1,\ 0\le n<M_2
\}.
$$

This is a rectangle in integer lattice coordinates. In real space it is a
parallelogram unless $\mathbf a_1$ and $\mathbf a_2$ are orthogonal.

The allowed wave vectors are

$$
\mathbf k_{\ell_1\ell_2}
=
\frac{\ell_1}{M_1}\mathbf b_1
+
\frac{\ell_2}{M_2}\mathbf b_2,
\qquad
\ell_1=0,\ldots,M_1-1,
\quad
\ell_2=0,\ldots,M_2-1.
$$

The orthogonality sum factorizes into two one-dimensional geometric sums. This
is the ordinary finite discrete Fourier transform.

## Example 2: skew supercell

More generally, choose two integer supercell vectors

$$
\mathbf T_1=p\mathbf a_1+q\mathbf a_2,
\qquad
\mathbf T_2=r\mathbf a_1+s\mathbf a_2.
$$

The index is

$$
N=|ps-qr|.
$$

Any complete set of $N$ representatives of
$\Lambda/\langle\mathbf T_1,\mathbf T_2\rangle$ gives an exact orthogonal
plane-wave basis. The representative set may be a skew parallelogram, but it
does not have to be.

The allowed wave vectors are the $N$ inequivalent solutions of

$$
e^{i\mathbf k\cdot\mathbf T_1}=1,
\qquad
e^{i\mathbf k\cdot\mathbf T_2}=1,
$$

modulo the original reciprocal lattice.

## Example 3: hexagonal region on a triangular lattice

Now specialize to a triangular lattice,

$$
|\mathbf a_1|=|\mathbf a_2|,
\qquad
\angle(\mathbf a_1,\mathbf a_2)=60^\circ.
$$

A centered hexagonal region with $L$ shells is

$$
S_L
=
\{
m\mathbf a_1+n\mathbf a_2:
|m|\le L,\ |n|\le L,\ |m+n|\le L
\}.
$$

It contains

$$
N=1+6+12+\cdots+6L
=
3L(L+1)+1
$$

lattice points.

This region is not a parallelogram. Nevertheless, it can be a complete
representative set for a high-symmetry superlattice quotient. One convenient
choice of superlattice generators is

$$
\mathbf T_1=(L+1)\mathbf a_1+L\mathbf a_2,
$$

$$
\mathbf T_2=-L\mathbf a_1+(2L+1)\mathbf a_2.
$$

Their index is

$$
\left|
\det
\begin{pmatrix}
L+1 & -L \\
L & 2L+1
\end{pmatrix}
\right|
=
3L^2+3L+1
=
N.
$$

The corresponding allowed wave vectors are the $N$ inequivalent solutions of

$$
e^{i\mathbf k\cdot\mathbf T_1}=1,
\qquad
e^{i\mathbf k\cdot\mathbf T_2}=1.
$$

With those wave vectors,

$$
\frac{1}{N}
\sum_{\mathbf R\in S_L}
e^{i(\mathbf k-\mathbf k')\cdot\mathbf R}
=
\delta_{\mathbf k,\mathbf k'}.
$$

The hexagonal shape works because it is a complete quotient representative set,
not because every visually hexagonal finite region has this property.

If a different triangular-lattice coordinate convention is used, the condition
$|m+n|\le L$ may appear as $|m-n|\le L$. This is only a sign convention.

### Efficient transform for the hexagonal region

For the hexagonal region above, the quotient is cyclic:

$$
\Lambda/\Lambda_{\rm super}\cong\mathbb Z_N,
\qquad
N=3L(L+1)+1.
$$

A convenient one-dimensional index for a lattice point

$$
\mathbf R_{mn}=m\mathbf a_1+n\mathbf a_2
$$

is

$$
x(m,n)
=
(2L+1)m+Ln
\pmod N.
$$

For $(m,n)\in S_L$, this map hits each residue

$$
x=0,\ldots,N-1
$$

exactly once. Therefore data on the hexagonal region can be reordered into a
length-$N$ array.

The characters are

$$
\chi_t(m,n)
=
\exp\left[
-\frac{2\pi i t}{N}x(m,n)
\right],
\qquad
t=0,\ldots,N-1.
$$

The corresponding wave vectors are

$$
\mathbf k_t
=
\frac{t}{N}
\left[
(2L+1)\mathbf b_1+L\mathbf b_2
\right]
\pmod{\Lambda^*},
\qquad
t=0,\ldots,N-1.
$$

They satisfy

$$
e^{i\mathbf k_t\cdot\mathbf T_1}=1,
\qquad
e^{i\mathbf k_t\cdot\mathbf T_2}=1.
$$

With the COMSOL-compatible forward-transform convention used in this note,

$$
\widehat f_t
=
\frac{1}{\sqrt N}
\sum_{(m,n)\in S_L}
\exp\left[
+\frac{2\pi i t}{N}x(m,n)
\right]
f(m,n).
$$

Thus the transform can be implemented by reordering the hexagonal data into
$f_x$ and applying a length-$N$ inverse FFT:

```python
F_t = np.fft.ifft(f_x, axis=0, norm="ortho")
```

The inverse transform is

```python
f_x = np.fft.fft(F_t, axis=0, norm="ortho")
```

If each lattice point carries a vector, field profile, or multiple components,
keep those as extra array dimensions and apply the FFT only along the cell
index axis.

## Example 4: 1D strip with one Bloch-periodic direction

A strip simulation can also satisfy the exact finite-Fourier assumptions, even
when the simulated geometry is rectangular and one direction is open. Suppose
the underlying two-dimensional lattice is generated by

$$
\Lambda=\langle \mathbf a_\parallel,\mathbf a_\perp\rangle,
$$

where $\mathbf a_\parallel$ is the direction made periodic in the simulation
and $\mathbf a_\perp$ labels the finite row direction. A strip eigenmode belongs
to one fixed Bloch sector in the periodic direction:

$$
E(\mathbf r+\mathbf a_\parallel)
=
e^{-i k_\parallel\cdot\mathbf a_\parallel}E(\mathbf r).
$$

Equivalently, after the gauge transformation

$$
\widetilde E(\mathbf r)=e^{+i k_\parallel\cdot\mathbf r}E(\mathbf r),
$$

the field is periodic along $\mathbf a_\parallel$. Therefore the periodic
direction may be quotiented out. The remaining lattice is the row lattice

$$
\Lambda/\langle\mathbf a_\parallel\rangle\cong\mathbb Z,
$$

generated by the finite row coordinate. If the selected analysis rows are

$$
S_N=\{n_0,n_0+1,\ldots,n_0+N-1\},
$$

then they form a complete representative set of the quotient
$\mathbb Z/N\mathbb Z$. The exact characters are

$$
\chi_m(n)
=
\exp\left[-\frac{2\pi i mn}{N}\right],
\qquad
m=0,\ldots,N-1.
$$

With the COMSOL-compatible forward-transform convention,

$$
\widehat f_m
=
\frac{1}{\sqrt N}
\sum_{n\in S_N}
\exp\left[+\frac{2\pi i mn}{N}\right]
f_n.
$$

The full Bloch wave vectors represented by this transform are

$$
\mathbf k_m
=
\mathbf k_\parallel+\frac{m}{N}\mathbf b_\perp,
\qquad
m=0,\ldots,N-1,
$$

where $\mathbf b_\perp$ is the reciprocal vector satisfying
$\mathbf b_\perp\cdot\mathbf a_\perp=2\pi$ and
$\mathbf b_\perp\cdot\mathbf a_\parallel=0$.

If a row cell is cut by the simulation's periodic seam, its pieces are not
separate cells. They must first be glued by the Bloch condition. A point mapped
across the seam by $s\mathbf a_\parallel$ contributes the phase

$$
e^{-i k_\parallel\cdot s\mathbf a_\parallel}
$$

for the physical field $E$, or no phase after using the periodic gauge
$\widetilde E$. For $k_\parallel=0$, this phase is one, so seam pieces are
glued directly.

In the current strip construction with bulk radius $L=3$, the bulk analysis
rows are

$$
n=-3,-2,-1,0,1,2,3,
$$

so $N=7$. This is an exact $\mathbb Z_7$ Fourier problem. It is a valid
special case of the same finite-lattice theory, but it is not the same quotient
as the two-dimensional hexagonal $N=3L(L+1)+1=37$ construction above.

## What can fail

The exact construction can fail if the selected set is not a complete set of
representatives of any useful quotient. Common examples are:

1. missing cells;
2. irregular boundary cuts;
3. a partial shell of a hexagonal region;
4. an arbitrary mask chosen by geometry rather than by a lattice quotient.

In such cases the plane waves are generally not orthogonal on the selected set.
The transform becomes a windowed Fourier transform, and different wave vectors
leak into each other.

A practical replacement is to use a Gram matrix or least-squares fit. For a
chosen list of wave vectors, define

$$
G_{\alpha\beta}
=
\sum_{\mathbf R\in S}
e^{i(\mathbf k_\alpha-\mathbf k_\beta)\cdot\mathbf R}.
$$

If $G$ is not diagonal, the plane waves are not orthogonal on $S$. One can still
fit data using this nonorthogonal dictionary, but the coefficients must be
interpreted with the Gram matrix.

## Practical checklist

To get an exact orthogonal finite Fourier basis:

1. Choose a finite lattice set $S$ with $N$ points.
2. Find a finite-index superlattice $\Lambda_{\rm super}$ with index $N$.
3. Check that $S$ contains exactly one representative of each coset of
   $\Lambda/\Lambda_{\rm super}$.
4. Generate the $N$ wave vectors satisfying superlattice periodicity.
5. Verify the matrix

$$
U_{\mathbf R,\mathbf k}
=
\frac{1}{\sqrt N}e^{-i\mathbf k\cdot\mathbf R}
$$

is unitary on the selected $S$ and $K$.

If this test passes, the finite Fourier transform is exact and Parseval's
identity holds. If it fails, use a Gram-matrix or least-squares formulation
instead.
