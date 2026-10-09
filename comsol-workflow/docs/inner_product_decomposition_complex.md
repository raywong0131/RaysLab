# Recovering Real Fourier-Subspace Energies from Cyclic-Shift Inner Products (Complex $x$)

This document is a self-contained reference for the following problem:

> **Given** a vector $x \in \mathbb{C}^N$ and access only to inner products of the form  
> $$
> r_m := \langle x,\;R^m x\rangle,\qquad m=0,1,\dots,N-1,
> $$
> where $R$ is the cyclic shift operator and $\langle\cdot,\cdot\rangle$ is the **Hermitian** inner product,  
> **recover** the **energy in each real Fourier subspace** (DC, each $(c_k,s_k)$ plane, and Nyquist when $N$ is even).

It also proves why recovering the **individual Fourier coefficients** (including phases) is impossible under these constraints.

---

## 1) Problem definition

### 1.1 Cyclic shift operator

Work in $\mathbb{C}^N$. Let vectors be indexed by $0,1,\dots,N-1$.

Define the **right cyclic shift** $R:\mathbb{C}^N\to\mathbb{C}^N$ by

$$
R(x_0,x_1,\dots,x_{N-1}) = (x_{N-1},x_0,x_1,\dots,x_{N-2}).
$$

Equivalently, in component form (indices mod $N$):

$$
(Rx)_n = x_{n-1},\qquad (R^m x)_n = x_{n-m}.
$$

Then $R^N = I$.

---

### 1.2 Inner product (Hermitian)

Use the standard **Hermitian** inner product on $\mathbb{C}^N$:

$$
\langle a,b\rangle := \sum_{n=0}^{N-1}\overline{a_n}\,b_n.
$$

With this inner product, every cyclic shift $R^m$ is unitary (it preserves inner products and norms).

---

### 1.3 The real orthonormal Fourier basis (used as a real subspace decomposition)

Define angles

$$
\theta_n = \frac{2\pi n}{N},\qquad n=0,1,\dots,N-1.
$$

The **real orthonormal Fourier basis vectors** (they have real entries) are:

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

These vectors are orthonormal in $\mathbb{C}^N$ as well (because they are orthonormal in $\mathbb{R}^N$ and the Hermitian inner product restricts to the real dot product on real vectors).

Therefore every $x\in\mathbb{C}^N$ decomposes uniquely as:

- If $N$ is odd:
  $$
  x
  = \alpha_0\,v_0
  + \sum_{k=1}^{(N-1)/2} \big(\alpha_{c,k}\,c_k + \alpha_{s,k}\,s_k\big),
  $$
- If $N$ is even:
  $$
  x
  = \alpha_0\,v_0
  + \sum_{k=1}^{N/2-1} \big(\alpha_{c,k}\,c_k + \alpha_{s,k}\,s_k\big)
  + \alpha_{N/2}\,v_{N/2}.
  $$

Here the coefficients $\alpha_\bullet$ are generally **complex**.

---

### 1.4 What we are allowed to measure

You are only allowed the scalar measurements

$$
r_m := \langle x,\;R^m x\rangle,\qquad m=0,1,\dots,N-1,
$$

i.e. inner products between $x$ and its cyclic shifts (no other probe vectors).

These $\{r_m\}$ are the **cyclic autocorrelation** of $x$.

Useful identities that always hold:

- $r_0=\langle x,x\rangle=\|x\|^2\in\mathbb{R}_{\ge 0}$,
- $r_{N-m}=\overline{r_m}$ (Hermitian symmetry).

---

### 1.5 Goal: recover **energies** of real Fourier subspaces

Define the **subspace energies** (coefficients are complex, so energy uses modulus):

- DC energy:
  $$
  E_0 := |\alpha_0|^2.
  $$
- For each $k\in\{1,\dots,\lfloor (N-1)/2\rfloor\}$, the energy in the 2D plane $\mathrm{span}\{c_k,s_k\}$:
  $$
  E_k := |\alpha_{c,k}|^2 + |\alpha_{s,k}|^2.
  $$
- If $N$ even, Nyquist energy:
  $$
  E_{N/2} := |\alpha_{N/2}|^2.
  $$

> **Problem:** Recover all $E_k$ using only $\{r_m\}$.

---

## 2) Why individual coefficients cannot be recovered (counterexamples)

Under the constraint “only $\langle x,R^m x\rangle$”, the data $\{r_m\}$ do **not** uniquely determine
$\alpha_0, \alpha_{c,k}, \alpha_{s,k}, \alpha_{N/2}$ (including their phases). Two counterexamples:

### 2.1 Global phase ambiguity ($x$ vs $e^{i\phi}x$)

Let $x' = e^{i\phi}x$ for any real $\phi$. Then for every $m$,

$$
r'_m
= \langle x',R^m x'\rangle
= \langle e^{i\phi}x,\;R^m(e^{i\phi}x)\rangle
= \overline{e^{i\phi}}\,e^{i\phi}\,\langle x,R^m x\rangle
= r_m.
$$

So $\{r_m\}$ are identical, but every Fourier coefficient (in any basis) gets multiplied by $e^{i\phi}$.  
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
because $R^t$ is unitary and preserves inner products.

But shifting $x$ generally changes the coefficients $\alpha_{c,k},\alpha_{s,k}$ (and their phases).  
So even the pair $(\alpha_{c,k},\alpha_{s,k})$ is not determined—only its total energy is.

> Conclusion: From $\{r_m\}$, you can at best recover **invariants** under global phase and shift—namely the subspace energies $E_k$.

---

## 3) How to compute subspace energies from $\{r_m\}$

### 3.1 Define the DFT of the autocorrelation

Let
$$
\omega = e^{-i\frac{2\pi}{N}}.
$$

Define the **Discrete Fourier Transform (DFT)** of the sequence $r_0,\dots,r_{N-1}$ as

$$
R_k := \sum_{m=0}^{N-1} r_m\,\omega^{km},\qquad k=0,1,\dots,N-1.
$$

(Complex arithmetic is allowed internally. We will prove below that each $R_k$ is real and nonnegative.)

---

### 3.2 The core identity (why DFT of autocorrelation gives power)

Define the (non-orthonormal) DFT of $x$ (as an intermediate object) by

$$
X_k := \sum_{n=0}^{N-1} x_n\,\omega^{-kn},\qquad k=0,\dots,N-1.
$$

> We will show:
$$
R_k = |X_k|^2.
$$

**Proof:**

Using $(R^m x)_n=x_{n-m\ (\mathrm{mod}\ N)}$, the autocorrelation is
$$
r_m = \langle x,R^m x\rangle
= \sum_{n=0}^{N-1} \overline{x_n}\,x_{n-m}.
$$

Now take the DFT in $m$:
$$
R_k = \sum_{m=0}^{N-1} r_m\,\omega^{km}
= \sum_{m=0}^{N-1}\left(\sum_{n=0}^{N-1}\overline{x_n}\,x_{n-m}\right)\omega^{km}.
$$
Swap sums:
$$
R_k = \sum_{n=0}^{N-1}\overline{x_n}\sum_{m=0}^{N-1} x_{n-m}\,\omega^{km}.
$$
Change variable $p=n-m$ (mod $N$); as $m$ runs over $0,\dots,N-1$, so does $p$:
$$
\sum_{m=0}^{N-1} x_{n-m}\,\omega^{km}
= \sum_{p=0}^{N-1} x_p\,\omega^{k(n-p)}
= \omega^{kn}\sum_{p=0}^{N-1} x_p\,\omega^{-kp}.
$$
Thus
$$
R_k
= \left(\sum_{n=0}^{N-1}\overline{x_n}\,\omega^{kn}\right)\left(\sum_{p=0}^{N-1} x_p\,\omega^{-kp}\right).
$$
The first factor is the complex conjugate of the second:
$$
\sum_{n=0}^{N-1}\overline{x_n}\,\omega^{kn}=\overline{\sum_{n=0}^{N-1} x_n\,\omega^{-kn}}=\overline{X_k}.
$$
Therefore
$$
R_k = \overline{X_k}\,X_k = |X_k|^2.
$$
$\square$

**Immediate consequences:**
- Each $R_k$ is **real** and $R_k\ge 0$.
- The data $\{r_m\}$ determine the **power spectrum** $\{|X_k|^2\}$.

---

### 3.3 Mapping $|X_k|^2$ to real Fourier-subspace energies $E_k$

We now connect the real-basis coefficients $\alpha$ (Section 1.3) to the complex DFT coefficients $X_k$.

#### DC and Nyquist (1D subspaces)

- DC:
  $$
  \alpha_0=\langle v_0,x\rangle=\frac{1}{\sqrt{N}}\sum_{n=0}^{N-1}x_n=\frac{X_0}{\sqrt{N}},
  \quad\Rightarrow\quad
  E_0=|\alpha_0|^2=\frac{|X_0|^2}{N}=\frac{R_0}{N}.
  $$

- If $N$ is even (Nyquist):
  note $\omega^{- (N/2)n}=e^{+i\pi n}=(-1)^n$, so
  $$
  X_{N/2}=\sum_{n=0}^{N-1}x_n(-1)^n,\qquad
  \alpha_{N/2}=\langle v_{N/2},x\rangle=\frac{1}{\sqrt{N}}X_{N/2}.
  $$
  Hence
  $$
  E_{N/2}=|\alpha_{N/2}|^2=\frac{|X_{N/2}|^2}{N}=\frac{R_{N/2}}{N}.
  $$

#### The $(c_k,s_k)$ plane (2D subspaces)

For $k=1,\dots,\left\lfloor\frac{N-1}{2}\right\rfloor$, use
$$
\cos(k\theta_n)=\frac{1}{2}\left(\omega^{-kn}+\omega^{kn}\right),\qquad
\sin(k\theta_n)=\frac{1}{2i}\left(\omega^{-kn}-\omega^{kn}\right).
$$

Since $c_k$ and $s_k$ are real, $\alpha_{c,k}=\langle c_k,x\rangle=\sum_n c_{k,n}x_n$ and similarly for $s_k$.
A direct substitution yields:

$$
\alpha_{c,k}
= \sqrt{\frac{2}{N}}\sum_{n=0}^{N-1} x_n\cos(k\theta_n)
= \sqrt{\frac{2}{N}}\cdot\frac{1}{2}\left(X_k + X_{N-k}\right),
$$
$$
\alpha_{s,k}
= \sqrt{\frac{2}{N}}\sum_{n=0}^{N-1} x_n\sin(k\theta_n)
= \sqrt{\frac{2}{N}}\cdot\frac{1}{2i}\left(X_k - X_{N-k}\right),
$$
where we used $\omega^{kn}=\omega^{-(N-k)n}$, hence $\sum_n x_n\omega^{kn}=X_{N-k}$.

Now compute the plane energy:
$$
|\alpha_{c,k}|^2+|\alpha_{s,k}|^2
=
\frac{2}{N}\left(\frac{1}{4}|X_k+X_{N-k}|^2+\frac{1}{4}|X_k-X_{N-k}|^2\right).
$$
Use the identity (valid for any complex $a,b$):
$$
|a+b|^2+|a-b|^2 = 2(|a|^2+|b|^2).
$$
With $a=X_k$, $b=X_{N-k}$, we get:
$$
|\alpha_{c,k}|^2+|\alpha_{s,k}|^2
=
\frac{2}{N}\cdot\frac{1}{4}\cdot 2\left(|X_k|^2+|X_{N-k}|^2\right)
=
\frac{|X_k|^2+|X_{N-k}|^2}{N}.
$$
Therefore, using $R_k=|X_k|^2$,
$$
E_k = |\alpha_{c,k}|^2+|\alpha_{s,k}|^2
= \frac{R_k+R_{N-k}}{N}.
$$

---

## 4) Full algorithm (uses only allowed measurements + arithmetic)

### Inputs
- $x\in\mathbb{C}^N$
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
(Complex arithmetic allowed internally; by Section 3.2, output satisfies $R_k\in\mathbb{R}_{\ge 0}$.)

### Step 3 — output real Fourier-subspace energies
- DC:
  $$
  E_0 = \frac{R_0}{N}.
  $$
- For $k=1,\dots,\left\lfloor\frac{N-1}{2}\right\rfloor$:
  $$
  E_k = \frac{R_k + R_{N-k}}{N}.
  $$
- If $N$ even:
  $$
  E_{N/2} = \frac{R_{N/2}}{N}.
  $$

These $E_k$ are exactly:
- $E_0=|\alpha_0|^2$
- $E_k=|\alpha_{c,k}|^2+|\alpha_{s,k}|^2$
- $E_{N/2}=|\alpha_{N/2}|^2$ (if applicable)

---

## 5) Example: $N=6$

Let
$$
x=(1,\ i,\ 0,\ 0,\ 0,\ 0)^\top,\qquad N=6,
$$
where $i=\sqrt{-1}$.
Use the right shift $R(x_0,\dots,x_5)=(x_5,x_0,\dots,x_4)$.

### 5.1 Autocorrelation measurements
Compute:
$$
r_m=\langle x,R^m x\rangle,\quad m=0,\dots,5.
$$

Direct evaluation gives:
- $r_0=\langle x,x\rangle = |1|^2+|i|^2 = 2$
- $r_1=\langle x,Rx\rangle = \overline{1}\cdot 0 + \overline{i}\cdot 1 = (-i)\cdot 1 = -i$
- $r_2=\langle x,R^2x\rangle = 0$
- $r_3=\langle x,R^3x\rangle = 0$
- $r_4=\langle x,R^4x\rangle = 0$
- $r_5=\langle x,R^5x\rangle = \overline{1}\cdot i + \overline{i}\cdot 0 = i$

So:
$$
(r_0,r_1,r_2,r_3,r_4,r_5)=(2,\ -i,\ 0,\ 0,\ 0,\ i),
$$
and indeed $r_{6-1}=r_5=\overline{r_1}$.

### 5.2 DFT of $r$
Let $\omega=e^{-i2\pi/6}$. Compute
$$
R_k=\sum_{m=0}^{5} r_m\,\omega^{km},\quad k=0,\dots,5.
$$

Because only $m=0,1,5$ are nonzero:
$$
R_k = 2 + (-i)\omega^{k} + (i)\omega^{5k}.
$$
But $\omega^{5k}=\overline{\omega^{k}}$, so $R_k$ is real and simplifies to:
$$
R_k = 2 + 2\,\mathrm{Im}(\omega^k).
$$

Evaluating for $k=0,\dots,5$ yields:
$$
(R_0,R_1,R_2,R_3,R_4,R_5)=(2,\ 2-\sqrt{3},\ 2-\sqrt{3},\ 2,\ 2+\sqrt{3},\ 2+\sqrt{3}).
$$
All values are real and nonnegative, consistent with $R_k=|X_k|^2$.

### 5.3 Convert to real Fourier-subspace energies
For $N=6$, the real subspaces are:
- $E_0$ (DC)
- $E_1$ for $\mathrm{span}\{c_1,s_1\}$
- $E_2$ for $\mathrm{span}\{c_2,s_2\}$
- $E_3$ (Nyquist)

Use:
$$
E_0=\frac{R_0}{6},\qquad
E_1=\frac{R_1+R_5}{6},\qquad
E_2=\frac{R_2+R_4}{6},\qquad
E_3=\frac{R_3}{6}.
$$

Compute:
- DC:
  $$
  E_0=\frac{2}{6}=\frac{1}{3}.
  $$
- $k=1$ plane:
  $$
  E_1=\frac{(2-\sqrt{3})+(2+\sqrt{3})}{6}=\frac{4}{6}=\frac{2}{3}.
  $$
- $k=2$ plane:
  $$
  E_2=\frac{(2-\sqrt{3})+(2+\sqrt{3})}{6}=\frac{4}{6}=\frac{2}{3}.
  $$
- Nyquist:
  $$
  E_3=\frac{2}{6}=\frac{1}{3}.
  $$

These are exactly the energies in the real Fourier subspaces, computed using only $\langle x,R^m x\rangle$ measurements.

---

## Summary

- The allowed measurements $\{r_m=\langle x,R^m x\rangle\}$ do **not** uniquely determine the complex Fourier coefficients (phases are lost; global phase and shifts preserve $\{r_m\}$).
- They **do** uniquely determine the **power spectrum** $|X_k|^2$ via the DFT of $\{r_m\}$.
- From that power spectrum, you can compute the **energy in each real Fourier subspace**:
  DC, each $(c_k,s_k)$ plane, and (if $N$ even) Nyquist.