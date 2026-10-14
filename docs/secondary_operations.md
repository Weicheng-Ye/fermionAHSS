# The secondary differentials `Tau` and `Psi`

`Tau` is $d_3:E_3^{n,0}\to E_3^{n+3,-2}$ (integral row $0$ to binary row
$-2$) and `Psi` is $d_4:E_4^{n,-1}\to E_4^{n+4,-4}$ (binary row $-1$ to
integral row $-4$); the other $d_3$, from row $-2$ to row $-4$, is the
primary $\operatorname{Dtilde}=\beta_s(\operatorname{Sq}^2+\omega)$. For a
fixed defining cochain $b$ the cochain representatives are $\tau'$ (binary,
degree $n+3$) and $\psi'$ (integral, degree $n+4$). The helper family is
`chi7_tail` with the secondary selectors $\epsilon=(1,0,0)$ and
$\eta=(1,0,1)$. Notation: [conventions.md](conventions.md). This page
transcribes `natural_secondary.gi` and `phase_eval.source_splitting`,
`low_phases.secondary` of the previous implementation; the formulas are
on normalized simplicial cochains (their status on a resolution is in
[natural_evaluation.md](natural_evaluation.md), Section 7).

## 1. Input and the first defining cochain

Start with a binary cocycle $a\in Z^n(X;\mathbb F_2)$; for integral input
$a=\rho A$ with $A\in Z^n(X;\mathbb Z_s)$. Put

$$
B=\tfrac12\,\delta\widetilde a,\qquad e=\rho B,\qquad C_B=\tfrac12\bigl(B+\widetilde e\bigr),\qquad
u=a\smile_{n-2}a,\qquad v=\omega\smile a,\qquad w=s\smile e
\tag{S1}
$$

($|B|=|e|=n+1$, $|u|=|v|=|w|=n+2$; the **plus** sign in $C_B$ is part of
the convention; both divisions are exact). In the interval-cut convention
$e=a\smile_{n-1}a$ literally, so $u+v+w=D(a)$. Choose a binary $b$ of degree $n+1$ with

$$
\delta b=u+v+w=D(a).
\tag{S2}
$$

If no such $b$ exists (the primary obstruction is nonzero) the secondary
operation is not defined on $a$.

## 2. The common binary cochain and integer cochain

$$
\boxed{\;
\begin{aligned}
F={}&Q^2(b)+\omega\smile b+\zeta_{2,n}(\omega,a)+\chi_n(a)
+u\smile_{n+1}v+u\smile_{n+1}w+v\smile_{n+1}w\\
&+\zeta_{1,n+1}(s,e)+(\omega\smile_1s)\smile e+s\smile u+s\smile s\smile\rho C_B,
\end{aligned}\qquad
q=\widetilde\omega\smile B+B\smile_{n-1}B .\;}
\tag{S3, S4}
$$

$F$ is reduced modulo two as a whole ($Q^2(b)=b\smile_{n-1}b+b\smile_n\delta b$);
$q$ is an integer cochain of degree $n+3$ whose products use the integral
cup signs (C5). The uncorrected integral expression is
$\psi_{\rm ref}=\tfrac14\delta_s(2\widetilde F+q)$.

## 3. `Psi`

With the **integer sum** $Z=\widetilde{s\smile u}+\widetilde{\omega\smile e}$,

$$
\boxed{\;\psi'=\tfrac14\,\delta_s\bigl(2\widetilde F+q+2Z\bigr)=\psi_{\rm ref}+\tfrac12\delta_sZ,\;}
\tag{S6, S7}
$$

whose cohomological normalization relative to the reference is
$\operatorname{Psi}_n(a)=\operatorname{Psi}_{{\rm ref},n}(a)+\beta_s(s\operatorname{Sq}^2[a])+\beta_s(\omega\operatorname{Sq}^1[a])$,
i.e. the coefficients of $\beta_s(s\operatorname{Sq}^2a),\beta_s(s^3a),\beta_s(\omega\operatorname{Sq}^1a)$
are $\eta=(1,0,1)$. The division by four is exact for a valid defining
system, and $\delta_s\psi'=0$.

## 4. `Tau`

Now $\delta_sA=0$. Put

$$
K=\tfrac12\bigl(A-\widetilde a\bigr),\qquad t=\rho K,\qquad x=\delta t,\qquad y=s\smile a
\tag{S9}
$$

(so $\delta t=e+s\smile a$ modulo two), and

$$
G=Q^2(t)+\omega\smile t+x\smile_ny+\zeta_{1,n}(s,a)+(\omega\smile_1s)\smile a+s\smile b,\qquad
L=\tfrac12\bigl(q-\delta_s\widetilde G\bigr),\qquad M=\widetilde F+L,\qquad \kappa=s\smile s\smile s\smile a .
\tag{S10}
$$

$G$ is binary of degree $n+2$; $L,M$ are integer cochains of degree $n+3$. Then

$$
\boxed{\;\tau'=\rho M+\kappa=F+\rho\Bigl(\tfrac12\bigl(q-\delta_s\widetilde G\bigr)\Bigr)+s\smile s\smile s\smile a,\;}
\tag{S11}
$$

so $\operatorname{Tau}_n(A)=\operatorname{Tau}_{{\rm ref},n}(A)+s^3\rho A$: the
coefficients on the basis $s^3\rho A,\ s\omega\rho A,\ \operatorname{Sq}^1(\omega)\rho A$
are $\epsilon=(1,0,0)$. This tuple belongs to `chi7_tail` and must not be
transferred to another $\chi$ family. $\delta\tau'=0$.

## 5. The matched integral lift

With the binary degree-$(n+2)$ cochain
$h=s\smile b+s\smile s\smile t+\omega\smile t+(\omega\smile_1s)\smile a$ and
$M'=M+Z-\delta_s\widetilde h$:

$$
\delta h=\rho Z+\kappa,\qquad \rho M'=\tau',\qquad \delta_sM'=2\psi',\qquad \delta\tau'=0,\qquad \delta_s\psi'=0 .
\tag{S14}
$$

For integral input the evaluator keeps $M'$ and obtains the matched
$\psi'$ as $\delta_sM'/2$; for binary input it uses (S7). The square-zero
statement is cohomological: $\operatorname{Dtilde}[\tau']=0$ in
$H^{n+6}(X;\mathbb Z_s)$, not a literally zero cochain.

## 6. The A-only splitting and the form used by the stacking model

The part of $\tau'$ that does not involve $b$ is isolated as follows
(`phase_eval.source_splitting`); every letter is binary except
$A,B,C_B,q_{\rm int},L_A$:

$$
\begin{aligned}
F_A&=\zeta_{2,n}(\omega,a)+\chi_n(a)+u\smile_{n+1}v+u\smile_{n+1}w+v\smile_{n+1}w+\zeta_{1,n+1}(s,e)+(\omega\smile_1s)\smile e+s\smile u+s\smile s\smile\rho C_B,\\
g&=Q^2(t)+\omega\smile t+(\delta t)\smile_n(s\smile a)+\zeta_{1,n}(s,a)+(\omega\smile_1s)\smile a,\\
q_{\rm int}&=\widetilde\omega\smile B+B\smile_{n-1}B,\qquad
L_A=\tfrac12\bigl(q_{\rm int}-\delta_s\widetilde g+\widetilde s\smile\widetilde p\bigr),\qquad p=u+v+w,\\
k_0&=F_A+\rho L_A+s\smile s\smile s\smile a,
\end{aligned}
\tag{S}
$$

so that $F=Q^2(b)+\omega\smile b+F_A=E(b)+F_A$, $G=g+s\smile b$, and
$\delta k_0=D(p)$. With $y=Q_D(b)$ and the **splitting carry**

$$
\lambda=\rho\Bigl(\tfrac12\bigl(\widetilde g+\widetilde s\smile\widetilde b-\widetilde{g+s\smile b}\bigr)\Bigr)=g\cdot(s\smile b)\quad\text{(pointwise product)},
$$

$$
\boxed{\;\tau'_n(A;b)=y+k_0+\delta\lambda=Q_D(b)+k_0(A)+\delta\lambda(A,b).\;}
\tag{S23}
$$

$p$ and $k_0$ depend only on $(A,s,\omega)$ (on $A$ through $a$ and $t$, so
modulo four); $\lambda$ also on $b$. For $n=1$ the tertiary phase uses the
**adjusted** source instead: $e_0=s\smile t$, $q=(\omega+s\smile s)\smile a$,
$k=k_0+Q_D(e_0)+\delta e_0\smile_2p+s\smile(\delta e_0\smile_3p)$ (Section 4
of [tertiary_operations.md](tertiary_operations.md)); for $n=2,3$, $q=p$ and
$k=k_0$. The degree-zero instance is (C12)–(C13) of the conventions.

## 7. The helpers $\zeta$ and $\chi$

The Cartan helpers $\zeta_{1,n},\zeta_{2,n}$ are the words (C4), zero for
$n=0$, with the cochain identities (for binary cocycles, $e_y=\rho(\delta\widetilde y/2)$)

$$
\delta\zeta_{1,n}(x,a)=\operatorname{Sq}^2(xa)+x\operatorname{Sq}^2a+x^2e_a,\qquad
\delta\zeta_{2,n}(x,a)=\operatorname{Sq}^2(xa)+x\operatorname{Sq}^2a+x^2a+e_xe_a
$$

($|x|=1$, resp. $2$). The Adem helper $\chi_n$, $0\le n\le7$, is the
quadratic Boolean polynomial

$$
\chi_n(a)=\sum_{c\in\mathcal C_n}\ \prod_{j\in\operatorname{faces}(c)}a(F_j)\pmod2,
$$

where $F_1,\dots,F_M$ are the $(n+1)$-vertex faces of the $(n+3)$-simplex
in lexicographic order and $\mathcal C_n$ is the `monomials` list of the
degree-$n$ record of the data file `data/chi-calibrated-degree7-anf.g`
(7.5 MB; the new repository must carry it). Decoding: expand the integer
$c$ in radix $2^{\text{bitsPerFace}}$ (7 for $n\le6$, 8 for $n=7$); its
nonzero digits are the one-based face indices of the monomial. The
boundary identity is

$$
\delta\chi_n(a)=\operatorname{Sq}^2\operatorname{Sq}^2a+e_a\smile_{n-2}e_a,\qquad \chi_n(0)=0,
$$

and $\chi_0=\chi_1=0$. Sizes:

| $n$ | output degree | faces | word summands | ANF monomials |
| ---: | ---: | ---: | ---: | ---: |
| 2 | 5 | 20 | 24 | 33 |
| 3 | 6 | 35 | 267 | 496 |
| 4 | 7 | 56 | 2 896 | 5 680 |
| 5 | 8 | 84 | 22 147 | 32 720 |
| 6 | 9 | 120 | 160 350 | 149 704 |
| 7 | 10 | 165 | 1 248 069 | 516 329 |

$\chi_n$ is a sum of surjection words (`tail_words`), i.e. a natural
operation on simplicial cochains, calibrated by a suspension comparison in
degree seven; other families are not interchangeable with it.

## 8. Domains, indeterminacy and the page groups

$$
\operatorname{Tau}_n:\ker\bigl(\operatorname{Dbar}:H^n(\mathbb Z_s)\to H^{n+2}(\mathbb F_2)\bigr)\longrightarrow
\frac{\ker\bigl(\operatorname{Dtilde}:H^{n+3}(\mathbb F_2)\to H^{n+6}(\mathbb Z_s)\bigr)}{DH^{n+1}(\mathbb F_2)},
\qquad
\operatorname{Psi}_n:\frac{\ker\bigl(D:H^n(\mathbb F_2)\to H^{n+2}(\mathbb F_2)\bigr)}{\operatorname{Dbar}H^{n-2}(\mathbb Z_s)}\longrightarrow
\frac{H^{n+4}(\mathbb Z_s)}{\operatorname{Dtilde}H^{n+1}(\mathbb F_2)} .
$$

A change of the allowed $b$ is absorbed by the image quotients. The
underlying relation among primary operations is

$$
D\circ D=\rho\operatorname{Dtilde}\circ(\operatorname{Sq}^1+s),
$$

both sides equal to $\operatorname{Sq}^3\operatorname{Sq}^1+s\operatorname{Sq}^3+s\operatorname{Sq}^2\operatorname{Sq}^1+s^3\operatorname{Sq}^1+(\operatorname{Sq}^1\omega)\operatorname{Sq}^1+s(\operatorname{Sq}^1\omega)$
(untwisted: $\operatorname{Sq}^2\operatorname{Sq}^2=\operatorname{Sq}^3\operatorname{Sq}^1$).
The runtime family covers $0\le n\le7$; the window needs
$\operatorname{Tau}_0,\dots,\operatorname{Tau}_3$ and $\operatorname{Psi}_0,\dots,\operatorname{Psi}_4$.

## 9. Page-class reductions

These hold for the classes on the page, not for the literal cochains
needed by a tertiary defining system.

- **Even integral input.** If $A=2U$ is a cocycle, take $b=0$: then $a=0$, $F=q=0$, $G=E(\rho U)$ and $[\operatorname{Tau}(2U)]=[\rho\operatorname{Dtilde}(\rho U)]$ modulo the $D$ indeterminacy. Only a primary operation is needed.
- **$\operatorname{Psi}_0$ on $BG$.** A closed degree-zero input is constant; for input $1$ the primary condition is $\delta b=\omega$, and $F=b\smile\delta b+\delta b\smile b=\delta(b\smile b)$, so $\psi'=\beta_s$ of a coboundary: $\operatorname{Psi}_0=0$ whenever defined.
- **Odd constant $A$ with $s=0$.** The primary condition is $[\omega]=0$; $B=e=q=0$, $F=\delta(b\smile b)$ and $G=t\,\omega$ with $t=(A-1)/2$ constant, so $[\operatorname{Tau}_0(A)]=0$.
- **Reuse of the matched lift.** For a binary input with a signed integral lift $A$, $\rho M'=\tau'(A;b)$ and $\delta_sM'=2\psi'(a;b)$ with the same $b$: the integral secondary evaluation already yields `Psi`.
