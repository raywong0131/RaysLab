# Recovering Real Fourier-Subspace Energies from Cyclic-Shift Inner Products

This document is a self-contained reference for the following problem:

> **Given** a vector $x \in \mathbb{R}^N$ and access only to inner products of the form  
> $$
> r_m := \langle x,\;R^m x\rangle,\qquad m=0,1,\dots,N-1,
> $$
> where $R$ is the cyclic shift operator,  
> **recover** the **energy in each real Fourier subspace** (DC, each $(c_k,s_k)$ plane, and Nyquist when $N$ is even).

It also proves why recovering the **individual Fourier coefficients** (with sign/phase) is impossible under these constraints.

---

## 1) Problem definition

### 1.1 Cyclic shift operator

Work in $\mathbb{R}^N$. Let vectors be indexed by $0,1,\dots,N-1$.

Define the **right cyclic shift** $R:\mathbb{R}^N\to\mathbb{R}^N$ by

$$
R(x_0,x_1,\dots,x_{N-1}) = (x_{N-1},x_0,x_1,\dots,x_{N-2}).
$$

Then $R^m$ shifts by $m$ steps, and $R^N = I$.

---

### 1.2 The real orthonormal Fourier basis

Define angles

$$
\theta_n = \frac{2\pi n}{N},\qquad n=0,1,\dots,N-1.
$$

The **real orthonormal Fourier basis** of $\mathbb{R}^N$ consists of:

#### (A) DC (constant) vector
$$
v_0 = \frac{1}{\sqrt{N}}(1,1,\dots,1)^\top.
$$

#### (B) Cosine/sine pairs for each frequency $k$
For
$$
k = 1,2,\dots,\left\lfloor\frac{N-1}{2}\right\rfloor,
$$
define
$$
c_k = \sqrt{\frac{2}{N}}\big(\cos(k\theta_n)\big)_{n=0}^{N-1},\qquad
s_k = \sqrt{\frac{2}{N}}\big(\sin(k\theta_n)\big)_{n=0}^{N-1}.
$$

#### (C) Nyquist vector (only when $N$ is even)
If $N$ is even, include
$$
v_{N/2} = \frac{1}{\sqrt{N}}(1,-1,1,-1,\dots)^\top,
$$
i.e. entry $n$ equals $(-1)^n$.

This is an orthonormal basis; every $x\in\mathbb{R}^N$ decomposes uniquely as:

- If $N$ is odd:
  $$
  x
  = \alpha_0\,v_0
  + \sum_{k=1}^{(N-1)/2} \big(\alpha_{c,k}\,c_k + \alpha_{s,k}\,s_k\big).
  $$
- If $N$ is even:
  $$
  x
  = \alpha_0\,v_0
  + \sum_{k=1}^{N/2-1} \big(\alpha_{c,k}\,c_k + \alpha_{s,k}\,s_k\big)
  + \alpha_{N/2}\,v_{N/2}.
  $$

---

### 1.3 What we are allowed to measure

You are only allowed the scalar measurements

$$
r_m := \langle x,\;R^m x\rangle,\qquad m=0,1,\dots,N-1,
$$

i.e. inner products between $x$ and its cyclic shifts (no other probe vectors).

These $\{r_m\}$ are the **cyclic autocorrelation** of $x$.

---

### 1.4 Goal: recover **energies** of Fourier subspaces

Define the **subspace energies**:

- DC energy:
  $$
  E_0 := \alpha_0^2.
  $$
- For each $k\in\{1,\dots,\lfloor (N-1)/2\rfloor\}$, the energy in the 2D plane $\mathrm{span}\{c_k,s_k\}$:
  $$
  E_k := \alpha_{c,k}^2 + \alpha_{s,k}^2.
  $$
- If $N$ even, Nyquist energy:
  $$
  E_{N/2} := \alpha_{N/2}^2.
  $$

> **Problem:** Recover all $E_k$ using only $\{r_m\}$.

---

## 2) Why individual coefficients cannot be recovered (counterexamples)

Under the constraint “only $\langle x,R^m x\rangle$”, the data $\{r_m\}$ do **not** uniquely determine
$\alpha_0, \alpha_{c,k}, \alpha_{s,k}, \alpha_{N/2}$. Two simple counterexamples:

### 2.1 Global sign ambiguity ($x$ vs $-x$)

Let $x'=-x$. Then for every $m$,

$$
r'_m
= \langle x',R^m x'\rangle
= \langle -x,\;R^m(-x)\rangle
= \langle x,\;R^m x\rangle
= r_m.
$$

So $\{r_m\}$ are identical, but all Fourier coefficients flip sign:
$$
\alpha'_0=-\alpha_0,\quad \alpha'_{c,k}=-\alpha_{c,k},\quad \alpha'_{s,k}=-\alpha_{s,k},\quad \alpha'_{N/2}=-\alpha_{N/2}.
$$
Therefore the coefficients are not uniquely recoverable.

### 2.2 Shift ambiguity ($x$ vs $R^t x$)

Let $x' = R^t x$ for some integer $t$. Then

$$
r'_m
= \langle R^t x,\;R^m R^t x\rangle
= \langle R^t x,\;R^t R^m x\rangle
= \langle x,\;R^m x\rangle
= r_m,
$$
because $R^t$ is an orthogonal matrix (a permutation), hence preserves inner products.

But shifting $x$ typically **rotates** each coefficient pair $(\alpha_{c,k},\alpha_{s,k})$ by a $k$-dependent angle.
So even the pair $(\alpha_{c,k},\alpha_{s,k})$ is not determined—only its radius (energy) is.

> Conclusion: From $\{r_m\}$, you can at best recover **invariants** under sign and shift—namely the subspace energies $E_k$.

---

## 3) How to compute subspace energies from $\{r_m\}$

### 3.1 Define the DFT of the autocorrelation

Let
$$
\omega = e^{-i\frac{2\pi}{N}}.
$$

Define the (non-orthonormal) **Discrete Fourier Transform (DFT)** of the sequence $r_0,\dots,r_{N-1}$ as

$$
R_k := \sum_{m=0}^{N-1} r_m\,\omega^{km},\qquad k=0,1,\dots,N-1.
$$

Facts:
- Because $r_m\in\mathbb{R}$ and $r_m=r_{N-m}$, each $R_k$ is real and $R_k\ge 0$.
- This DFT uses complex numbers only as intermediate values; final outputs $R_k$ are real.

---

### 3.2 The core identity (why DFT of autocorrelation gives power)

Define the **complex Fourier coefficient** (intermediate object) for $x$:

$$
X_k := \sum_{n=0}^{N-1} x_n\,\omega^{kn},\qquad k=0,\dots,N-1.
$$

> We will show:
$$
R_k = |X_k|^2.
$$

**Proof:**

Start from the definition of autocorrelation under the chosen shift:
$$
r_m = \langle x,R^m x\rangle = \sum_{n=0}^{N-1} x_n\,(R^m x)_n.
$$
With right shift, $(R^m x)_n = x_{n-m\ (\mathrm{mod}\ N)}$. Therefore
$$
r_m = \sum_{n=0}^{N-1} x_n\,x_{n-m}.
$$

Now compute the DFT:
$$
R_k = \sum_{m=0}^{N-1} r_m\,\omega^{km}
= \sum_{m=0}^{N-1}\left(\sum_{n=0}^{N-1} x_n x_{n-m}\right)\omega^{km}.
$$
Swap the sums:
$$
R_k = \sum_{n=0}^{N-1} x_n \sum_{m=0}^{N-1} x_{n-m}\,\omega^{km}.
$$
Change variable $p=n-m$ (mod $N$); as $m$ runs over $0,\dots,N-1$, so does $p$:
$$
\sum_{m=0}^{N-1} x_{n-m}\,\omega^{km}
= \sum_{p=0}^{N-1} x_p\,\omega^{k(n-p)}
= \omega^{kn}\sum_{p=0}^{N-1} x_p\,\omega^{-kp}.
$$
So
$$
R_k = \sum_{n=0}^{N-1} x_n \left(\omega^{kn}\sum_{p=0}^{N-1} x_p\,\omega^{-kp}\right)
= \left(\sum_{n=0}^{N-1} x_n\,\omega^{kn}\right)\left(\sum_{p=0}^{N-1} x_p\,\omega^{-kp}\right).
$$
Recognize the two factors:
$$
R_k = X_k \cdot \overline{X_k} = |X_k|^2.
$$
That proves the identity. $\square$

Thus the DFT of the autocorrelation yields the **power spectrum** $|X_k|^2$.

---

### 3.3 Mapping $|X_k|^2$ to real Fourier-subspace energies $E_k$

The real basis groups complex frequencies into real subspaces:

- $k=0$ corresponds to DC (1D).
- For $k=1,\dots,\lfloor (N-1)/2\rfloor$, the real subspace $\mathrm{span}\{c_k,s_k\}$ corresponds to the complex pair $(k,N-k)$.
- If $N$ is even, $k=N/2$ is Nyquist (1D).

With the normalization used in Section 1, the subspace energies are obtained from $R_k=|X_k|^2$ by:

- DC:
  $$
  E_0 = \alpha_0^2 = \frac{R_0}{N}.
  $$
- Nyquist (only if $N$ even):
  $$
  E_{N/2} = \alpha_{N/2}^2 = \frac{R_{N/2}}{N}.
  $$
- For $k=1,\dots,\left\lfloor\frac{N-1}{2}\right\rfloor$:
  $$
  E_k = \alpha_{c,k}^2+\alpha_{s,k}^2 = \frac{2R_k}{N}.
  $$

**Why these scaling factors are correct (proof sketch):**

- Under our DFT convention, DC amplitude is $X_0=\sum_n x_n$.
  The DC coefficient in the real orthonormal basis is
  $$
  \alpha_0 = \langle v_0,x\rangle = \frac{1}{\sqrt{N}}\sum_{n=0}^{N-1}x_n = \frac{X_0}{\sqrt{N}}.
  $$
  Hence
  $$
  E_0=\alpha_0^2=\frac{|X_0|^2}{N}=\frac{R_0}{N}.
  $$

- For $1\le k\le \lfloor (N-1)/2\rfloor$, the complex coefficient $X_k$ encodes both cosine and sine components.
  For real $x$, $X_{N-k}=\overline{X_k}$, so the energy in that conjugate pair is $2|X_k|^2$.
  That conjugate-pair energy matches the real-plane energy $E_k=\alpha_{c,k}^2+\alpha_{s,k}^2$ up to the factor $1/N$,
  giving
  $$
  E_k=\frac{2|X_k|^2}{N}=\frac{2R_k}{N}.
  $$

- If $N$ is even, the Nyquist component is purely real with $X_{N/2}=\sum_n x_n(-1)^n$, and
  $$
  \alpha_{N/2} = \frac{X_{N/2}}{\sqrt{N}}\;\Rightarrow\;E_{N/2}=\frac{|X_{N/2}|^2}{N}=\frac{R_{N/2}}{N}.
  $$

This completes the mapping from shift-inner-products to real-subspace energies.

---

## 4) Full algorithm (uses only allowed measurements + arithmetic)

### Inputs
- $x\in\mathbb{R}^N$
- cyclic shift $R$ (defined above)
- $\omega = e^{-i2\pi/N}$

### Step 1 — measure autocorrelation via allowed inner products
For $m=0,1,\dots,N-1$:
$$
r_m := \langle x,\;R^m x\rangle.
$$

### Step 2 — compute DFT of $\{r_m\}$
For $k=0,1,\dots,N-1$:
$$
R_k := \sum_{m=0}^{N-1} r_m\,\omega^{km}.
$$
(Complex arithmetic allowed internally; output $R_k$ is real $\ge 0$.)

### Step 3 — output real Fourier-subspace energies
- DC:
  $$
  E_0 = \frac{R_0}{N}.
  $$
- For $k=1,\dots,\left\lfloor\frac{N-1}{2}\right\rfloor$:
  $$
  E_k = \frac{2R_k}{N}.
  $$
- If $N$ even:
  $$
  E_{N/2} = \frac{R_{N/2}}{N}.
  $$

These $E_k$ are exactly:
- $E_0=\alpha_0^2$
- $E_k=\alpha_{c,k}^2+\alpha_{s,k}^2$
- $E_{N/2}=\alpha_{N/2}^2$ (if applicable)

---

## 5) Example: $N=6$

Let
$$
x=(1,2,3,4,5,6)^\top,\qquad N=6.
$$
Use the right shift $R(x_0,\dots,x_5)=(x_5,x_0,\dots,x_4)$.

### 5.1 Autocorrelation measurements
Compute:
$$
r_m=\langle x,R^m x\rangle,\quad m=0,\dots,5.
$$

You get:
- $r_0 = 91$
- $r_1 = 76$
- $r_2 = 67$
- $r_3 = 64$
- $r_4 = 67$
- $r_5 = 76$

So:
$$
(r_0,r_1,r_2,r_3,r_4,r_5)=(91,76,67,64,67,76).
$$

### 5.2 DFT of $r$
Let $\omega=e^{-i2\pi/6}$. Compute
$$
R_k=\sum_{m=0}^{5} r_m\,\omega^{km},\quad k=0,\dots,5.
$$

Result:
$$
(R_0,R_1,R_2,R_3,R_4,R_5)=(441,\ 36,\ 12,\ 9,\ 12,\ 36).
$$
(As expected, $R_{6-k}=R_k$ and all are nonnegative reals.)

### 5.3 Convert to real Fourier-subspace energies
For $N=6$, the real subspaces are:
- $E_0$ (DC)
- $E_1$ for the plane $\mathrm{span}\{c_1,s_1\}$
- $E_2$ for the plane $\mathrm{span}\{c_2,s_2\}$
- $E_3$ (Nyquist)

Use the formulas:
$$
E_0=\frac{R_0}{6},\qquad
E_1=\frac{2R_1}{6},\qquad
E_2=\frac{2R_2}{6},\qquad
E_3=\frac{R_3}{6}.
$$

Compute:
- DC:
  $$
  E_0=\frac{441}{6}=73.5
  \quad\Rightarrow\quad |\alpha_0|=\sqrt{73.5}\approx 8.573.
  $$
- $k=1$ plane:
  $$
  E_1=\frac{2\cdot 36}{6}=12
  \quad\Rightarrow\quad \sqrt{\alpha_{c,1}^2+\alpha_{s,1}^2}=\sqrt{12}\approx 3.464.
  $$
- $k=2$ plane:
  $$
  E_2=\frac{2\cdot 12}{6}=4
  \quad\Rightarrow\quad \sqrt{\alpha_{c,2}^2+\alpha_{s,2}^2}=\sqrt{4}=2.
  $$
- Nyquist:
  $$
  E_3=\frac{9}{6}=1.5
  \quad\Rightarrow\quad |\alpha_{3}|=\sqrt{1.5}\approx 1.225.
  $$

These are exactly the energies in the Fourier subspaces, computed using only $\langle x,R^m x\rangle$ measurements.

---

## Summary

- The allowed measurements $\{r_m=\langle x,R^m x\rangle\}$ do **not** uniquely determine the signed Fourier coefficients.
- They **do** uniquely determine the **energy** in each Fourier subspace.
- The computation is:
  1) measure autocorrelation $r_m$,
  2) take DFT to get $R_k$,
  3) map $R_k$ to energies $E_k$ via simple scaling.