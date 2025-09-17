# Z-Transform Exercises 
**Author:** Juan Rodriguez Esteban

---


## A1 — Finite Sequences → Polynomials

**Answers**

* i) For $x[n]=\{1,2,5\}$ at $n=\{0,1,2\}$:
  $X_1(z)=1+2z^{-1}+5z^{-2}$.
* ii) For $x[n]=\{0,3,0,4\}$ at $n=\{0,1,2,3\}$:
  $X_2(z)=3z^{-1}+4z^{-3}$.

**Quick MATLAB (symbolic or numeric)**

```matlab
% Symbolic (if available)
syms n z; X1 = 1 + 2*z^(-1) + 5*z^(-2); X2 = 3*z^(-1) + 4*z^(-3);
pretty(X1); pretty(X2);

% Numeric check (no toolbox)
x1=[1 2 5]; n1=0:2; X1d=@(zz) sum(x1.*zz.^(-n1)); X1p=@(zz) 1+2./zz+5./zz.^2;
x2=[0 3 0 4]; n2=0:3; X2d=@(zz) sum(x2.*zz.^(-n2)); X2p=@(zz) 3./zz+4./zz.^3;
ztest=[2+3j,-1.2+0.5j,0.7-0.9j];
max(abs(arrayfun(X1d,ztest)-arrayfun(X1p,ztest)))
max(abs(arrayfun(X2d,ztest)-arrayfun(X2p,ztest)))
```


* Finite sequences map to polynomials in $z^{-1}$; there’s no infinite series to converge. The ROC is the entire $z$-plane except possible singularities at $z=0$ (and $z=\infty$ if negative-time samples existed). Verifying numerically avoids dependency on toolboxes.

---

## A2 — Infinite Sequences & ROC

**Answers**

| Case | $x[n]$             | $X(z)$                                  | ROC | Poles |       |          |
| ---- | ------------------ | --------------------------------------- | --- | ----- | ----- | -------- |
| (a)  | $0.6^{n}u[n]$      | $\frac{1}{1-0.6z^{-1}}=\frac{z}{z-0.6}$ | (   | z     | >0.6) | $z=0.6$  |
| (b)  | $({-0.8})^{n}u[n]$ | $\frac{1}{1+0.8z^{-1}}=\frac{z}{z+0.8}$ | (   | z     | >0.8) | $z=-0.8$ |
| (c)  | $-0.9^{n}u[-n-1]$  | $\frac{1}{1-0.9z^{-1}}=\frac{z}{z-0.9}$ | (   | z     | <0.9) | $z=0.9$  |

**MATLAB snippet**

```matlab
% Symbolic route (if available)
syms n z; a = sym(0.6); b = sym(-0.8); r = sym(0.9);
Xa = symsum((a^n)*z^(-n), n, 0, inf);
Xb = symsum((b^n)*z^(-n), n, 0, inf);
% Left-sided (change of variables)
syms m positive; Xc = -symsum((r^(-m))*z^(m), m, 1, inf);
simplify([Xa, Xb, Xc])

% Numeric fallback: function handles
Fa=@(zz) 1./(1-0.6./zz); Fb=@(zz) 1./(1+0.8./zz); Fc=@(zz) 1./(1-0.9./zz);
```


* Right-sided exponentials converge outside the outermost pole radius; left-sided converge inside the innermost. The ROC never includes poles. Visualizing ROC as circles in the $z$-plane makes pole/zero placements intuitive.

---

## A3 — Properties: Linearity & Shifting

**Given** $x_1[n]=0.5^n u[n]$, $x_2[n]=(-0.5)^n u[n]$.

**Answers**

* (a) $Z\{2x_1[n]-3x_2[n]\}=\displaystyle 2\frac{1}{1-0.5z^{-1}}-3\frac{1}{1+0.5z^{-1}}=\frac{-1+2.5z^{-1}}{1-0.25z^{-2}}$; ROC $|z|>0.5$; poles at $\pm 0.5$.
* (b) $Z\{x_1[n-3]\}=z^{-3}X_1(z)=\displaystyle \frac{z^{-3}}{1-0.5z^{-1}}=\frac{1}{z^{2}(z-0.5)}$; ROC $|z|>0.5$.

**MATLAB**

```matlab
% Symbolic (if available)
syms n z
x1=(sym(1)/2)^n*heaviside(n); x2=(-sym(1)/2)^n*heaviside(n);
X_lin = simplify(ztrans(2*x1-3*x2, n, z));
X_shift = simplify(ztrans(subs(x1,n,n-3), n, z));

% Numeric fallback (no symbolic)
X1=@(zz) 1./(1-0.5./zz); X2=@(zz) 1./(1+0.5./zz);
Xlin=@(zz) 2*X1(zz)-3*X2(zz); Xshift=@(zz) (zz).^(-3).*X1(zz);
```


* Linearity and time shifting are fast paths to results that would otherwise require sums. The shift inserts a factor $z^{-k}$, moving poles and adding zeros at the origin as expected.

---

## A4 — Inverse Z-Transform

**Answers**

* (a) $X_a(z)=\frac{1}{1-0.7z^{-1}}\Rightarrow x_a[n]=0.7^n u[n]$ (pole at 0.7; causal ROC $|z|>0.7$).
* (b) $X_b(z)=\frac{1-0.5z^{-1}}{1-0.8z^{-1}}=\tfrac{5}{8}+\tfrac{3}{8}\frac{1}{1-0.8z^{-1}}\Rightarrow x_b[n]=\tfrac{5}{8}\,\delta[n]+\tfrac{3}{8}(0.8)^n u[n]$.
  Equivalent: $x_b[n]=\delta[n]+0.3\,(0.8)^{n-1}u[n-1]$.

**MATLAB**

```matlab
syms n z
Xa = 1/(1 - 0.7*z^(-1));  xa = iztrans(Xa, z, n);
Xb = (1 - 0.5*z^(-1))/(1 - 0.8*z^(-1)); xb = iztrans(Xb, z, n);
```


* Matching to known pairs and quick partial fractions make $\text{iztrans}$ almost redundant. Both forms of $x_b[n]$ are useful depending on whether you prefer an explicit $\delta[n]$ or a single shifted exponential.

---

## A5 — H(z), Poles/Zeros & Frequency Response

**System**
$H(z)=\dfrac{1-2.4z^{-1}+2.88z^{-2}}{1-0.8z^{-1}+0.64z^{-2}}$.

**Poles/Zeros (numeric)**

* Zeros: $z=1.2\pm j1.2$ (|z|≈1.697, ∠≈±45°).
* Poles: $z=0.4\pm j0.6928$ (|z|=0.8, ∠≈±60°).
* Stable (poles inside unit circle), non-minimum-phase (zeros outside).

**Plots**
After running the code below, attach these files:

![Circle.bmp](https://github.com/user-attachments/files/22384350/Circle.bmp)
![Bode.bmp](https://github.com/user-attachments/files/22384352/Bode.bmp)

```matlab
b=[1 -2.4 2.88]; a=[1 -0.8 0.64]; ensure_fig_dir();
% Pole–zero
figure('Color','w'); zplane(b,a); grid on; title('Pole–Zero Plot');
exportgraphics(gcf,'figures/A5_pz.png','Resolution',200);

% Frequency response (512 pts)
[H,w] = freqz(b,a,512);
figure('Color','w');
subplot(2,1,1); plot(w/pi,abs(H)); grid on; ylabel('|H|'); title('Magnitude');
subplot(2,1,2); plot(w/pi,unwrap(angle(H))); grid on; xlabel('\omega/\pi'); ylabel('Phase (rad)');
exportgraphics(gcf,'figures/A5_freqz.png','Resolution',200);

% Optional: test signal
n=0:511; x=sin(0.2*pi*n)+0.5*sin(0.8*pi*n); y=filter(b,a,x);
```

**Key readings from the plots**

* Peak near $\omega/\pi\approx 0.348$ with $|H|\approx 5.51$; DC≈1.76; Nyquist≈2.57.
* Interpretation: resonant **peaking/band-pass** behavior centered near $\omega\approx \pi/3$.
* Zeros at ≈±45° shape the response; poles at ±60° set the resonance.



* The pole radius (0.8) sets a strong resonance near the pole angles on the unit circle. Zeros outside the unit circle confirmed a non-minimum-phase design while still being stable. The response emphasizes mid-band energy relative to DC and Nyquist, consistent with the pole–zero geometry.

