## Numerical Data Files Documentation

### Simulation Data Files

The `Simulation~.mat` files contain all numerical data for classically simulable decompositions of the corresponding quantum state families. These files store critical information about the classical simulability decomposition parameters.

#### Data Structure

- **`nonzeroElemIndex`**: Indices of nonzero elements in the coefficient vector $c_\mu$ (as discussed in the main text), identified within numerical precision.
- **`c_value`**: Numerical values for $c_\mu$.
- **`Tau_value`**: Numerical values for the Hermitian $\tau_{i,\mu}$.
- **`CV`**: Numerical classical visibility.
- **`Dtable`**: Deterministic strategy.
- **`rho`**: State family.
- **`d`**: Dimension of Hilbert space.

#### Numerical Example: Decomposition Analysis

Taking `SimulationForMUBs2.mat` as a concrete example:

For each nonzero coefficient index, the numerical values `c_value` are approximately equal to $1/36$. 

The key constants appearing in the decomposition are:
$$
\begin{align}
\frac{1}{9} &\approx 0.111111111111111 \\
\frac{\sqrt{3}}{6} &\approx 0.288675134594813 \\
\frac{\sqrt{3}}{9} &\approx 0.192450089729875 \\
\frac{\sqrt{3}}{18} &\approx 0.096225044864938
\end{align}
$$

Let `A = Tau_value{1,nonzeroElemIndex(1)} * 36`:

```
A = 
[ 0.111111087338697 + 0.000000000000000i, -0.055555543667109 + 0.096225024273455i,
 -0.055555543667109 - 0.096225024273455i,  0.111111087338925 + 0.000000000000000i,
 -0.222222174653105 + 0.192450048538202i, -0.055555543667019 - 0.288675072803231i,
 
 -0.222222174653105 - 0.192450048538202i;
 -0.055555543667019 + 0.288675072803231i;
  0.777777611187483 + 0.000000000000000i ]
```

Similarly, let`B = Tau_value{2,nonzeroElemIndex(1)} * 36`:

```
B = 
[ 0.111111087338867 + 0.000000000000000i, -0.222222174653166 - 0.192450048538436i,
 -0.222222174653166 + 0.192450048538436i,  0.777777611187347 + 0.000000000000000i,
 -0.055555543667131 - 0.096225024273520i, -0.055555543667177 + 0.288675072803129i,
 
 -0.055555543667131 + 0.096225024273520i;
 -0.055555543667177 - 0.288675072803129i;
  0.111111087338891 + 0.000000000000000i ]
```

And`C = Tau_value{3,nonzeroElemIndex(1)} * 36`:

```
C = 
[ 0.777777611187541 + 0.000000000000000i,  0.277777718320317 + 0.096225024264909i,
  0.277777718320317 - 0.096225024264909i,  0.111111087338833 + 0.000000000000000i,
  0.277777718320278 - 0.096225024264609i,  0.111111087334112 + 0.000000000000102i,
  
  0.277777718320278 + 0.096225024264609i;
  0.111111087334112 - 0.000000000000102i;
  0.111111087338731 + 0.000000000000000i ]
```

In symbolic form, these matrices are approximately:


$$
A \approx \begin{pmatrix}
\frac{1}{9} & -\frac{1}{18}+i\frac{\sqrt{3}}{18} & -\frac{2}{9}-i\frac{\sqrt{3}}{9} \\[4pt]
-\frac{1}{18}-i\frac{\sqrt{3}}{18} & \frac{1}{9} & -\frac{1}{18}+i\frac{\sqrt{3}}{6} \\[4pt]
-\frac{2}{9}+i\frac{\sqrt{3}}{9} & -\frac{1}{18}-i\frac{\sqrt{3}}{6} & \frac{7}{9}
\end{pmatrix}
$$

$$
B \approx \begin{pmatrix}
\frac{1}{9} & -\frac{2}{9}-i\frac{\sqrt{3}}{9} & -\frac{1}{18}+i\frac{\sqrt{3}}{18} \\[4pt]
-\frac{2}{9}+i\frac{\sqrt{3}}{9} & \frac{7}{9} & -\frac{1}{18}-i\frac{\sqrt{3}}{6} \\[4pt]
-\frac{1}{18}-i\frac{\sqrt{3}}{18} & -\frac{1}{18}+i\frac{\sqrt{3}}{6} & \frac{1}{9}
\end{pmatrix}
$$

$$
C \approx \begin{pmatrix}
\frac{7}{9} & \frac{5}{18}+i\frac{\sqrt{3}}{18} & \frac{5}{18}+i\frac{\sqrt{3}}{18} \\[4pt]
\frac{5}{18}-i\frac{\sqrt{3}}{18} & \frac{1}{9} & \frac{1}{9} \\[4pt]
\frac{5}{18}-i\frac{\sqrt{3}}{18} & \frac{1}{9} & \frac{1}{9}
\end{pmatrix}
$$

#### Key Properties: Commutativity and Simultaneous Diagonalization

All three matrices $A$, $B$, and $C$ are Hermitian and pairwise commute:

$$
AB=BA,\quad AC=CA,\quad BC=CB.
$$

This commutativity property ensures that they can be simultaneously diagonalized by a unitary transformation.

#### Simultaneous Diagonalization

A common orthonormal eigenbasis can be chosen as:

$$
u_1=\frac{1}{6}
\begin{pmatrix}
2\\
-1-i\sqrt{3}\\
-4+2i\sqrt{3}
\end{pmatrix},\quad
u_2=\frac{1}{6}
\begin{pmatrix}
2\\
-4+2i\sqrt{3}\\
-1-i\sqrt{3}
\end{pmatrix},\quad
u_3=\frac{1}{6}
\begin{pmatrix}
5+i\sqrt{3}\\
2\\
2
\end{pmatrix}.
$$

Define the unitary matrix:

$$
U= \begin{pmatrix} |&|&|\\ u_1&u_2&u_3\\ |&|&| \end{pmatrix}=
\frac{1}{6}
\begin{pmatrix}
2&2&5+i\sqrt{3}\\
-1-i\sqrt{3}&-4+2i\sqrt{3}&2\\
-4+2i\sqrt{3}&-1-i\sqrt{3}&2
\end{pmatrix}.
$$

This unitary transformation satisfies $U^\dagger U=I$ and diagonalizes all three matrices:

$$
U^\dagger A U=
\begin{pmatrix}
1&0&0\\
0&0&0\\
0&0&0
\end{pmatrix},\quad
U^\dagger B U=
\begin{pmatrix}
0&0&0\\
0&1&0\\
0&0&0
\end{pmatrix},\quad
U^\dagger C U=
\begin{pmatrix}
0&0&0\\
0&0&0\\
0&0&1
\end{pmatrix}.
$$

#### Eigenvalue Decompositions

The matrices admit the following spectral decompositions:

$$
\boxed{
A=
U
\begin{pmatrix}
1&0&0\\
0&0&0\\
0&0&0
\end{pmatrix}
U^\dagger
}
$$

$$
\boxed{
B=
U
\begin{pmatrix}
0&0&0\\
0&1&0\\
0&0&0
\end{pmatrix}
U^\dagger
}
$$

$$
\boxed{
C=
U
\begin{pmatrix}
0&0&0\\
0&0&0\\
0&0&1
\end{pmatrix}
U^\dagger
}
$$

#### General Validity

The matrices $A$, $B$, $C$ can be diagonalized simultaneously by the same unitary transformation. By the same reasoning, the spectral decomposition holds for the `Tau_value` matrices at any nonzero coefficient index position. 

Ultimately, the decomposition always satisfies the classical simulability condition:

$$
\begin{equation}
	\mathrm{CV}\rho_x + (1-\mathrm{CV})\frac{\mathbb{I}}{d} =\sum_{\mu,i}D(i|x,\mu)\,\tau_{i,\mu},\quad \forall x.
\end{equation}
$$

where $\mathrm{CV}$ denotes the classical visibility of the quantum state.



### Witness Data Files

The `Witness~.mat` files contain all numerical witness matrix data for the corresponding quantum state family.

#### Data Structure

- **`W_opt`**: Numerical witness values.
- **`t`**: Numerical threshold value.
- **`rho`**: State family.
- **`d`**: Dimension of Hilbert space.
