# Symmetry-Adapted Orthonormal Basis in **N** Dimensions (Real, 0-Indexed)

This document defines a systematic way to construct an orthonormal basis of $\mathbb{R}^N$ that is adapted to the **cyclic shift symmetry** of coordinates indexed by $0,1,\dots,N-1$. The basis is built from sampled cosine/sine patterns on the $N$-point circle and has three key properties:

- It is an **orthonormal** and **complete** basis of $\mathbb{R}^N$.
- It decomposes $\mathbb{R}^N$ into **independent frequency subspaces**.
- The **cyclic shift operator** becomes block-diagonal: a 1D invariant part, several 2D rotation blocks, and (if $N$ is even) a 1D sign-flip part.

---

## 1) Construction (for any $N$, real vectors, 0-start index)

### 1.1 Notation

Vectors are real column vectors
$$
x = (x_0, x_1, \dots, x_{N-1})^\top \in \mathbb{R}^N,
$$
with standard inner product
$$
\langle x, y\rangle = \sum_{n=0}^{N-1} x_n y_n.
$$

Define the angle grid
$$
\theta_n = \frac{2\pi n}{N}, \qquad n=0,1,\dots,N-1.
$$

### 1.2 Basis vectors

#### Constant (DC) vector
$$
v_0 = \frac{1}{\sqrt{N}}(1,1,\dots,1)^\top.
$$

#### Cosine/Sine paired vectors
For each integer
$$
k = 1,2,\dots,\left\lfloor \frac{N-1}{2}\right\rfloor,
$$
define two vectors in $\mathbb{R}^N$:
$$
c_k = \sqrt{\frac{2}{N}}
\begin{pmatrix}
\cos(k\theta_0)\\
\cos(k\theta_1)\\
\vdots\\
\cos(k\theta_{N-1})
\end{pmatrix},
\qquad
s_k = \sqrt{\frac{2}{N}}
\begin{pmatrix}
\sin(k\theta_0)\\
\sin(k\theta_1)\\
\vdots\\
\sin(k\theta_{N-1})
\end{pmatrix}.
$$

#### Alternating (Nyquist) vector (only when $N$ is even)
If $N$ is even, include
$$
v_{N/2} = \frac{1}{\sqrt{N}}
\begin{pmatrix}
1\\
-1\\
1\\
-1\\
\vdots
\end{pmatrix}
\quad\text{(entry } n \text{ is }(-1)^n\text{)}.
$$

### 1.3 The basis set and a recommended ordering

- If $N$ is odd, the basis set is
$$
\mathcal{B}=\{v_0\}\ \cup\ \{c_k,s_k:\ k=1,\dots,(N-1)/2\}.
$$

- If $N$ is even, the basis set is
$$
\mathcal{B}=\{v_0\}\ \cup\ \{c_k,s_k:\ k=1,\dots,N/2-1\}\ \cup\ \{v_{N/2}\}.
$$

A common column ordering for the change-of-basis matrix $U$ is:

- If $N$ is odd:
$$
U = [\,v_0,\ c_1,\ s_1,\ c_2,\ s_2,\ \dots,\ c_{(N-1)/2},\ s_{(N-1)/2}\,].
$$

- If $N$ is even:
$$
U = [\,v_0,\ c_1,\ s_1,\ \dots,\ c_{N/2-1},\ s_{N/2-1},\ v_{N/2}\,].
$$

(Orthogonality will imply $U^\top U = I$, so $U^{-1}=U^\top$.)

### 1.4 Example: $N=6$ (Even case)

For $N=6$, the basis contains 6 orthonormal vectors. Here they are listed with indices $0\dots5$ corresponding to the columns of $U$. We extract scalar factors so the vector entries are integers.

1.  **Index 0 (DC, $v_0$):**
    $$
    v_0 = \frac{1}{\sqrt{6}}
    \begin{pmatrix} 1 \\ 1 \\ 1 \\ 1 \\ 1 \\ 1 \end{pmatrix}
    $$

2.  **Index 1 (Cosine $k=1$, $c_1$):**
    $$
    c_1 = \frac{1}{2\sqrt{3}}
    \begin{pmatrix} 2 \\ 1 \\ -1 \\ -2 \\ -1 \\ 1 \end{pmatrix}
    $$

3.  **Index 2 (Sine $k=1$, $s_1$):**
    $$
    s_1 = \frac{1}{2}
    \begin{pmatrix} 0 \\ 1 \\ 1 \\ 0 \\ -1 \\ -1 \end{pmatrix}
    $$

4.  **Index 3 (Cosine $k=2$, $c_2$):**
    $$
    c_2 = \frac{1}{2\sqrt{3}}
    \begin{pmatrix} 2 \\ -1 \\ -1 \\ 2 \\ -1 \\ -1 \end{pmatrix}
    $$

5.  **Index 4 (Sine $k=2$, $s_2$):**
    $$
    s_2 = \frac{1}{2}
    \begin{pmatrix} 0 \\ 1 \\ -1 \\ 0 \\ 1 \\ -1 \end{pmatrix}
    $$

6.  **Index 5 (Nyquist, $v_3$):**
    $$
    v_3 = \frac{1}{\sqrt{6}}
    \begin{pmatrix} 1 \\ -1 \\ 1 \\ -1 \\ 1 \\ -1 \end{pmatrix}
    $$

---

## 2) Cheat sheet: converting representations (vanilla $\leftrightarrow$ new basis)

### 2.1 What “coordinates in the new basis” mean

Let the new-basis coefficient vector be $\alpha \in \mathbb{R}^N$. Using the column ordering of $U$ above:

- $x$ is the **vanilla-coordinate** vector in $\mathbb{R}^N$.
- $\alpha$ is the **new-basis-coordinate** vector such that
$$
x = U\alpha.
$$

Because the basis is orthonormal (proved later), we also have
$$
\alpha = U^\top x.
$$

### 2.2 Vanilla → new basis (forward)

Coefficients are dot products with the basis vectors:

- DC:
$$
\alpha_0 = \langle v_0, x\rangle
= \frac{1}{\sqrt{N}}\sum_{n=0}^{N-1} x_n.
$$

- For each $k=1,\dots,\left\lfloor\frac{N-1}{2}\right\rfloor$:
$$
\alpha_{c,k} = \langle c_k, x\rangle
= \sqrt{\frac{2}{N}} \sum_{n=0}^{N-1} x_n \cos\!\left(\frac{2\pi k n}{N}\right),
$$
$$
\alpha_{s,k} = \langle s_k, x\rangle
= \sqrt{\frac{2}{N}} \sum_{n=0}^{N-1} x_n \sin\!\left(\frac{2\pi k n}{N}\right).
$$

- If $N$ is even (alternating mode):
$$
\alpha_{N/2} = \langle v_{N/2}, x\rangle
= \frac{1}{\sqrt{N}} \sum_{n=0}^{N-1} x_n (-1)^n.
$$

### 2.3 New basis → vanilla (inverse)

Reconstruct each coordinate $x_n$ from coefficients:

- If $N$ is odd:
$$
x_n
=
\frac{\alpha_0}{\sqrt{N}}
+
\sum_{k=1}^{(N-1)/2}
\sqrt{\frac{2}{N}}
\left[
\alpha_{c,k}\cos\!\left(\frac{2\pi k n}{N}\right)
+
\alpha_{s,k}\sin\!\left(\frac{2\pi k n}{N}\right)
\right].
$$

- If $N$ is even:
$$
x_n
=
\frac{\alpha_0}{\sqrt{N}}
+
\sum_{k=1}^{N/2-1}
\sqrt{\frac{2}{N}}
\left[
\alpha_{c,k}\cos\!\left(\frac{2\pi k n}{N}\right)
+
\alpha_{s,k}\sin\!\left(\frac{2\pi k n}{N}\right)
\right]
+
\frac{\alpha_{N/2}}{\sqrt{N}}(-1)^n.
$$

### 2.4 Matrix one-liners

With $U$ defined in Section 1.3:

- Forward:
$$
\alpha = U^\top x.
$$

- Inverse:
$$
x = U\alpha.
$$

### 2.5 Example: $N=6$

Given data vector $x=(x_0, \dots, x_5)^\top \in \mathbb{R}^6$ and coefficients $\alpha = (\alpha_0, \alpha_{c,1}, \alpha_{s,1}, \alpha_{c,2}, \alpha_{s,2}, \alpha_3)^\top$.

**Forward Map ($x \to \alpha$):**
$$ \alpha_0 = \tfrac{1}{\sqrt{6}}(x_0+x_1+x_2+x_3+x_4+x_5) $$
$$ \alpha_{c,1} = \tfrac{1}{2\sqrt{3}}(2x_0 + x_1 - x_2 - 2x_3 - x_4 + x_5) $$
$$ \alpha_{s,1} = \tfrac{1}{2}(x_1 + x_2 - x_4 - x_5) $$
$$ \alpha_{c,2} = \tfrac{1}{2\sqrt{3}}(2x_0 - x_1 - x_2 + 2x_3 - x_4 - x_5) $$
$$ \alpha_{s,2} = \tfrac{1}{2}(x_1 - x_2 + x_4 - x_5) $$
$$ \alpha_{3} = \tfrac{1}{\sqrt{6}}(x_0 - x_1 + x_2 - x_3 + x_4 - x_5) $$

**Inverse Map ($\alpha \to x$):**
Reconstructing all components $x_0 \dots x_5$:

$$
x_0 = \frac{\alpha_0}{\sqrt{6}} + \frac{2\alpha_{c,1}}{2\sqrt{3}} + 0 + \frac{2\alpha_{c,2}}{2\sqrt{3}} + 0 + \frac{\alpha_3}{\sqrt{6}}
$$
$$
x_1 = \frac{\alpha_0}{\sqrt{6}} + \frac{\alpha_{c,1}}{2\sqrt{3}} + \frac{\alpha_{s,1}}{2} - \frac{\alpha_{c,2}}{2\sqrt{3}} + \frac{\alpha_{s,2}}{2} - \frac{\alpha_3}{\sqrt{6}}
$$
$$
x_2 = \frac{\alpha_0}{\sqrt{6}} - \frac{\alpha_{c,1}}{2\sqrt{3}} + \frac{\alpha_{s,1}}{2} - \frac{\alpha_{c,2}}{2\sqrt{3}} - \frac{\alpha_{s,2}}{2} + \frac{\alpha_3}{\sqrt{6}}
$$
$$
x_3 = \frac{\alpha_0}{\sqrt{6}} - \frac{2\alpha_{c,1}}{2\sqrt{3}} + 0 + \frac{2\alpha_{c,2}}{2\sqrt{3}} + 0 - \frac{\alpha_3}{\sqrt{6}}
$$
$$
x_4 = \frac{\alpha_0}{\sqrt{6}} - \frac{\alpha_{c,1}}{2\sqrt{3}} - \frac{\alpha_{s,1}}{2} - \frac{\alpha_{c,2}}{2\sqrt{3}} + \frac{\alpha_{s,2}}{2} + \frac{\alpha_3}{\sqrt{6}}
$$
$$
x_5 = \frac{\alpha_0}{\sqrt{6}} + \frac{\alpha_{c,1}}{2\sqrt{3}} - \frac{\alpha_{s,1}}{2} - \frac{\alpha_{c,2}}{2\sqrt{3}} - \frac{\alpha_{s,2}}{2} - \frac{\alpha_3}{\sqrt{6}}
$$

---

## 3) Proof of orthogonality, normalization, completeness

This section proves that the vectors in Section 1 form an **orthonormal basis** of $\mathbb{R}^N$.

### 3.1 A core summation identity

For any integer $m$, define the complex number
$$
\omega = e^{i2\pi/N}.
$$
Then
$$
\sum_{n=0}^{N-1}\omega^{mn}
=
\begin{cases}
N, & m \equiv 0 \ (\mathrm{mod}\ N),\\
0, & m \not\equiv 0 \ (\mathrm{mod}\ N).
\end{cases}
$$
This is a geometric series: if $\omega^m \neq 1$, then
$$
\sum_{n=0}^{N-1}(\omega^m)^n = \frac{1-(\omega^m)^N}{1-\omega^m}=\frac{1-1}{1-\omega^m}=0.
$$

From this, one also obtains the real consequences (for $m\not\equiv 0 \pmod N$):
$$
\sum_{n=0}^{N-1}\cos(m\theta_n)=0,\qquad
\sum_{n=0}^{N-1}\sin(m\theta_n)=0.
$$

### 3.2 Pairwise orthogonality of cosine/sine families

Use product-to-sum identities:
$$
\cos a\cos b=\tfrac12(\cos(a-b)+\cos(a+b)),
$$
$$
\sin a\sin b=\tfrac12(\cos(a-b)-\cos(a+b)),
$$
$$
\cos a\sin b=\tfrac12(\sin(a+b)+\sin(b-a)).
$$

Let $k,\ell \in \{1,\dots,\lfloor (N-1)/2\rfloor\}$. Summing over $n$ and applying Section 3.1 to the resulting $\cos(m\theta_n)$ and $\sin(m\theta_n)$ sums yields:

- Cos–cos:
$$
\sum_{n=0}^{N-1}\cos(k\theta_n)\cos(\ell\theta_n)
=
\begin{cases}
\frac{N}{2}, & k=\ell,\\
0, & k\neq \ell,
\end{cases}
$$

- Sin–sin:
$$
\sum_{n=0}^{N-1}\sin(k\theta_n)\sin(\ell\theta_n)
=
\begin{cases}
\frac{N}{2}, & k=\ell,\\
0, & k\neq \ell,
\end{cases}
$$

- Cos–sin:
$$
\sum_{n=0}^{N-1}\cos(k\theta_n)\sin(\ell\theta_n)=0.
$$

After multiplying by the normalization factor $\sqrt{2/N}$, these become:

$$
\langle c_k, c_\ell\rangle = \delta_{k\ell},\qquad
\langle s_k, s_\ell\rangle = \delta_{k\ell},\qquad
\langle c_k, s_\ell\rangle = 0.
$$

### 3.3 Orthogonality with the constant vector $v_0$

For any $k\ge1$,
$$
\langle v_0, c_k\rangle
=
\frac{1}{\sqrt{N}}\sqrt{\frac{2}{N}}
\sum_{n=0}^{N-1}\cos(k\theta_n)
=0,
$$
$$
\langle v_0, s_k\rangle
=
\frac{1}{\sqrt{N}}\sqrt{\frac{2}{N}}
\sum_{n=0}^{N-1}\sin(k\theta_n)
=0.
$$
Thus $v_0$ is orthogonal to every $c_k$ and $s_k$.

### 3.4 The even-$N$ alternating vector $v_{N/2}$

Assume $N$ is even. Then $v_{N/2}$ has entries $(-1)^n$.

- Normalization:
$$
\|v_{N/2}\|^2 = \frac{1}{N}\sum_{n=0}^{N-1}1=1.
$$

- Orthogonality with $v_0$:
$$
\langle v_0, v_{N/2}\rangle
=
\frac{1}{N}\sum_{n=0}^{N-1}(-1)^n = 0.
$$

- Orthogonality with $c_k, s_k$ for $k=1,\dots,N/2-1$:  
write $(-1)^n=\cos(\pi n)$ and note $\pi n = (N/2)\theta_n$. Then products like $(-1)^n\cos(k\theta_n)$ become combinations of $\cos((N/2\pm k)\theta_n)$, whose sums vanish by Section 3.1 because $N/2\pm k\not\equiv 0 \pmod N$ for those $k$. Similarly for sine. Hence
$$
\langle v_{N/2}, c_k\rangle = 0,\qquad \langle v_{N/2}, s_k\rangle = 0.
$$

### 3.5 Completeness

Count the vectors:

- If $N$ is odd: $1 + 2\cdot\frac{N-1}{2} = N$.
- If $N$ is even: $1 + 2\cdot(\frac{N}{2}-1) + 1 = N$.

You have $N$ pairwise orthonormal vectors in $\mathbb{R}^N$. Therefore they form a complete orthonormal basis.

Equivalently, with $U$ as the column matrix, orthonormality implies
$$
U^\top U = I,
$$
so $U$ is invertible and $\mathrm{span}(\mathcal{B})=\mathbb{R}^N$.

---

## 4) Proof of the cyclic shift structure (block-diagonal action)

### 4.1 Define the cyclic shift operator

Define $R:\mathbb{R}^N\to\mathbb{R}^N$ by shifting indices forward by 1 (mod $N$):

$$
(Rx)_n = x_{n-1\ (\mathrm{mod}\ N)}.
$$

This is an orthogonal linear operator (it permutes coordinates), so it preserves inner products and norms.

### 4.2 Action on the constant and alternating modes

- Constant mode: for all $n$,
$$
(Rv_0)_n = (v_0)_{n-1} = \frac{1}{\sqrt{N}},
$$
so
$$
Rv_0 = v_0.
$$
Thus $\mathrm{span}\{v_0\}$ is a 1D invariant subspace and $R$ acts as the identity there.

- Alternating mode (only if $N$ even): since $(v_{N/2})_n = \frac{1}{\sqrt{N}}(-1)^n$,
$$
(Rv_{N/2})_n = (v_{N/2})_{n-1} = \frac{1}{\sqrt{N}}(-1)^{n-1} = -\frac{1}{\sqrt{N}}(-1)^n,
$$
so
$$
Rv_{N/2} = -v_{N/2}.
$$
Thus $\mathrm{span}\{v_{N/2}\}$ is invariant and $R$ acts as multiplication by $-1$.

### 4.3 Action on each cosine/sine pair: a 2D rotation

Fix $k \in \{1,\dots,\lfloor (N-1)/2\rfloor\}$. Consider the pair $(c_k, s_k)$.

Compute $(Rc_k)_n$:
$$
(Rc_k)_n = (c_k)_{n-1}
= \sqrt{\frac{2}{N}}\cos\!\big(k\theta_{n-1}\big).
$$

Since $\theta_{n-1} = \theta_n - \frac{2\pi}{N}$, we have
$$
k\theta_{n-1} = k\theta_n - \frac{2\pi k}{N}.
$$

Using $\cos(A-B)=\cos A\cos B+\sin A\sin B$,
$$
\cos\!\left(k\theta_n - \frac{2\pi k}{N}\right)
=
\cos(k\theta_n)\cos\!\left(\frac{2\pi k}{N}\right)
+
\sin(k\theta_n)\sin\!\left(\frac{2\pi k}{N}\right).
$$

Therefore
$$
Rc_k
=
\cos\!\left(\frac{2\pi k}{N}\right)c_k
+
\sin\!\left(\frac{2\pi k}{N}\right)s_k.
$$

Similarly, for $s_k$:
$$
(Rs_k)_n = (s_k)_{n-1}
= \sqrt{\frac{2}{N}}\sin\!\left(k\theta_n - \frac{2\pi k}{N}\right),
$$
and using $\sin(A-B)=\sin A\cos B-\cos A\sin B$,
$$
Rs_k
=
-\sin\!\left(\frac{2\pi k}{N}\right)c_k
+
\cos\!\left(\frac{2\pi k}{N}\right)s_k.
$$

### 4.4 Block-diagonal form in the new basis

Let $\phi_k = \frac{2\pi k}{N}$. In the ordered basis $\{c_k, s_k\}$, the operator $R$ acts as:

$$
R|_{\mathrm{span}\{c_k,s_k\}} \;\equiv\;
\begin{pmatrix}
\cos\phi_k & -\sin\phi_k\\
\sin\phi_k & \cos\phi_k
\end{pmatrix},
$$
i.e. a **2D rotation** by angle $\phi_k$.

Putting everything together, in the full new basis $\mathcal{B}$, the matrix of $R$ is block-diagonal:

- A $1\times 1$ block $[1]$ on $v_0$,
- For each $k$: a $2\times 2$ rotation block on $(c_k,s_k)$,
- If $N$ is even: a final $1\times 1$ block $[-1]$ on $v_{N/2}$.

This is the precise statement that the basis is **symmetry-adapted** to cyclic shifts.