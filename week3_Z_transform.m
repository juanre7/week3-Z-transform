%% Discrete-time basic signals (stem + labels)
% Author: Juan Rodriguez Esteban
clear; close all; clc;

%% Common index axis (finite support for discrete n)
% Task 1: MATLAB Programming
% 
% A1 – Warm‑up: Finite Sequences → Polynomials
% (a) For each sequence, write X(z) explicitly as a polynomial in z^{-1} and verify with MATLAB symbolic tools.
%    i) x[n] = {1, 2, 5} at n = {0,1,2}
%   ii) x[n] = {0, 3, 0, 4} at n = {0,1,2,3}
syms n z
% i) Construct from samples
x1 = [1 2 5]; n1 = 0:2;
X1 = 1 + 2*z^(-1) + 5*z^(-2);
X1_alt = sum(x1 .* z.^(-n1));
simplify(X1 - X1_alt)   % -> 0

pretty(X1)

% ii) Construct from samples (include the zeros explicitly)
x2 = [0 3 0 4]; n2 = 0:3;
X2 = 0 + 3*z^(-1) + 0*z^(-2) + 4*z^(-3);
X2_alt = sum(x2 .* z.^(-n2));
simplify(X2 - X2_alt)   % -> 0

pretty(X2)

% (b) Briefly explain (2–3 sentences): Why do finite sequences behave “nicely” regarding ROC?
% 
% Finite-length sequences have z-transforms that are finite sums of z^{-1} powers (polynomials), 
% so convergence issues vanish except at singularities (z=0 and possibly z=∞).
% Hence, for right-sided finite sequences the ROC is the entire z-plane except z=0; if both 
% positive and negative indices exist, the ROC is the entire plane except possibly z=0 and z=∞.

%% A2 – Infinite Sequences & ROC
% For each sequence, derive X(z) and specify the ROC. Then, verify using symbolic MATLAB where possible.
% (a) x[n] = a^n u[n], with a = 0.6
% (b) x[n] = (−0.8)^n u[n]
% (c) x[n] = −(0.9)^n u[−n−1]   (a left‑sided sequence)

    syms n z real

    % (a) a = 0.6, right-sided
    a = sym(0.6);
    Xa = symsum((a^n)*z^(-n), n, 0, inf);         % -> z/(z - a)
    Xa = simplify(Xa);                             % == 1/(1 - a*z^-1)

    % (b) a = -0.8, right-sided
    b = sym(-0.8);
    Xb = symsum((b^n)*z^(-n), n, 0, inf);
    Xb = simplify(Xb);                             % == 1/(1 - b*z^-1) = 1/(1+0.8 z^-1)

    % (c) left-sided: x[n] = -(0.9)^n u[-n-1]
    r = sym(0.9);
    syms m positive
    % sum_{n=-inf}^{-1} [-(r^n) z^{-n}]  -> change n=-m, m=1..inf
    Xc = -symsum((r^(-m))*z^(m), m, 1, inf);      % -> z/(z - r)
    Xc = simplify(Xc);                             % == 1/(1 - r*z^-1)

    pretty(Xa); pretty(Xb); pretty(Xc);

% Explain in 3–5 sentences how the ROC relates to |z| and to poles.
% ROC: the set of z where X(z)’s power series converges. 
% For one-sided exponentials, % |a/z|<1 (right) or |z/a|<1 (left), i.e., outside or inside a circle. 
% Poles are never in the ROC; for rational X(z): 
% Right-sided → outside the outermost pole 
% Left-sided → inside the innermost
% Two-sided → between poles

%% A3 – Properties: Linearity & Shifting 
% Let x1[n] = (0.5)^n u[n] and x2[n] = (−0.5)^n u[n].
% (a) Use linearity to compute Z{2x1[n] − 3x2[n]}.
% (b) Use the time-shift property to compute Z{x1[n−3]}.
% Verify both with MATLAB symbolic functions and include the outputs in your report.

syms n z
x1 = (sym(1)/2)^n * heaviside(n);
x2 = (-sym(1)/2)^n * heaviside(n);

% Individual transforms
X1 = ztrans(x1, n, z);      % -> 1/(1 - 0.5*z^-1)
X2 = ztrans(x2, n, z);      % -> 1/(1 + 0.5*z^-1)

% (a) Linearity
X_lin  = ztrans(2*x1 - 3*x2, n, z);
X_linS = simplify(X_lin);    % -> (-1 + 2.5*z^-1)/(1 - 0.25*z^-2)

% (b) Time shift
X_shift  = ztrans( subs(x1, n, n-3), n, z );
X_shiftS = simplify(X_shift); % -> z^-3/(1 - 0.5*z^-1)

pretty(X1); pretty(X2); pretty(X_linS); pretty(X_shiftS);



%% A4 – Inverse Z-Transform 
% Find x[n] for each X(z) using MATLAB’s iztrans and (briefly) show how you could get it by inspection of poles/ROC.
% (a) X(z) = 1 / (1 − 0.7 z^{−1})
% (b) X(z) = (1 − 0.5 z^{−1}) / (1 − 0.8 z^{−1})

syms n z
% (a)
Xa = 1/(1 - 0.7*z^(-1));
xa = iztrans(Xa, z, n)        % -> 0.7^n*heaviside(n)

% (b)
Xb = (1 - 0.5*z^(-1)) / (1 - 0.8*z^(-1));
xb = iztrans(Xb, z, n)        % -> dirac(n) + 0.3*(0.8)^(n-1)*heaviside(n-1)
                              % (equivalently  (3/8)*(0.8)^n*heaviside(n) + (5/8)*dirac(n))



%% A5 – H(z), Poles/Zeros & Frequency Response 
% Consider the filter frequency response:
%     H(z) = (1 − 2.4 z^{−1} + 2.88 z^{−2}) / (1 − 0.8 z^{−1} + 0.64 z^{−2})
% (a) Plot poles and zeros with zplane and list their numeric values.

b = [1 -2.4 2.88];
a = [1 -0.8 0.64];

% Numeric zeros/poles
z = roots(b);
p = roots(a);
fprintf('Zeros:\n'); disp(z.')
fprintf('Poles:\n'); disp(p.')
fprintf('Zero magnitudes/angles (deg):\n'); disp([abs(z).', rad2deg(angle(z)).'])
fprintf('Pole magnitudes/angles (deg):\n'); disp([abs(p).', rad2deg(angle(p)).'])

% Pole–zero plot: try zplane if available, else manual
if exist('zplane','file')
    figure; zplane(b,a); grid on; title('Pole–Zero Plot');
else
    th = linspace(0,2*pi,512);
    figure; plot(cos(th),sin(th),'k--'); axis equal; hold on; grid on;
    plot(real(z),imag(z),'o','MarkerSize',8,'LineWidth',1.5); % zeros
    plot(real(p),imag(p),'x','MarkerSize',10,'LineWidth',1.5); % poles
    xlabel('Re\{z\}'); ylabel('Im\{z\}'); title('Pole–Zero Plot (manual)');
    legend('Unit circle','Zeros','Poles');
end


% (b) Plot magnitude and phase responses with freqz (use 256 or 512 points).

N = 512;
if exist('freqz','file')
    [H,w] = freqz(b,a,N);             % uses z^{-1} convention
else
    w = linspace(0,pi,N).';
    kB = (0:numel(b)-1); kA = (0:numel(a)-1);
    Hb = exp(-1j*w*kB) * b.';         % sum_k b(k)*e^{-jwk}
    Ha = exp(-1j*w*kA) * a.';         % sum_k a(k)*e^{-jwk}
    H  = Hb ./ Ha;
end

figure;
subplot(2,1,1);
plot(w/pi, abs(H)); grid on; ylabel('|H(e^{j\omega})|'); title('Magnitude');
subplot(2,1,2);
plot(w/pi, unwrap(angle(H))); grid on; xlabel('\omega/\pi'); ylabel('Phase (rad)');

% (c) Identify the filter type (low-pass, high-pass, band-pass, etc.) and justify in 4–6 sentences by referring to the plots and pole locations.

% Poles at radius 0.8 and angles ±60° give resonance near ω≈π/3. 
% The magnitude peaks around ω/π≈0.35 at |H|≈5.5, while DC≈1.76 and Nyquist≈2.57, so mid 
% frequencies are emphasized. Zeros at 1.2±j1.2 (outside the unit circle) shape the response 
% without notches. Poles inside ⇒ stable; zeros outside ⇒ non-minimum-phase. 
% Overall, it’s a resonant peaking/band-pass filter near ω≈π/3.

%% Task 2: GitHub Submission
% 
% Upload your MATLAB code files to the repository.
% Add a README.md file that includes:
%    • Clear answers, plots, and short explanations (bullet points are fine).
%    • Screenshots or exported figures for pole-zero plots and frequency responses.
%    • Brief reflections (2–4 sentences per task) on what you learned/noticed.
% Deliverables
% GitHub Repository link (MATLAB code + README + plots).
% 