# Supplementary-doc-for-delayed-UIO
This repository contains the simulation parameters of "Finite-delayed State Reconstruction with Unknown Inputs via Geometric Quotient Outputs" and the MATLAB function of the numerical calculation algorithms 
# delayed_uio_design

**Design a fixed-delay functional unknown-input observer using geometric quotient spaces.**

`delayed_uio_design` is an **offline design function**. It computes geometric subspaces, quotient-coordinate maps, a stabilizing output injection (a *friend*), and target recovery maps. It returns the matrices needed to implement an observer; it does **not** acquire measurements, run a simulation, or generate plots itself.

Source: [`delayed_uio_design.m`](delayed_uio_design.m). Place this README next to that file in the repository. All helper functions required by the designer are local functions in the same `.m` file.

## Contents

[Syntax](#syntax) · [Description](#description) · [Requirements](#requirements) · [Examples](#examples) · [Input Arguments](#input-arguments) · [Output Arguments](#output-arguments) · [Algorithms](#algorithms) · [Numerical Considerations](#numerical-considerations) · [Troubleshooting](#troubleshooting)

## Syntax

```matlab
g = delayed_uio_design(A,B,Bbar,C,Q,r)
g = delayed_uio_design(A,B,Bbar,C,Q,r,opts)
```

The first syntax uses the default options. The second supplies options in a **scalar structure**; the function does not accept a list of name–value arguments.

## Description

The function considers the discrete-time linear time-invariant plant

$$
\begin{aligned}
x(k+1)&=Ax(k)+Bu(k)+\bar B d(k),\\
y(k)&=Cx(k),
\end{aligned}
$$

and designs an estimate of the prescribed linear quantity

$$
    \lim_{k\to\infty}
    \|Qx(k-r)-\widehat q_r(k)\|
    =0,\qquad k\geq r.
$$

Here, `u(k)` is known, `d(k)` is unknown, and `r` is a nonnegative integer delay **in samples**. An estimate produced at sample `k` concerns the state at sample `k-r`. The function has no sampling-time argument; with sampling period `Ts`, the physical delay is `r*Ts`.

The design tests the geometric condition

$$
\mathscr{I}_r\subseteq\ker Q,\qquad
\mathscr{I}_r=\mathscr W_{g,r}^{\ast}\cap\ker\bar C_r.
$$

`g.I` contains a basis of the ambiguity subspace $\mathscr I_r$. For a feasible target, the returned recovery maps satisfy

$$
Q=E_{Q,r}P_{g,r}+F_{Q,r}\bar C_r.
$$

With the returned stabilizing friend, the quotient-state error has Schur dynamics. In the exact, noise-free model, these identities yield asymptotic reconstruction independently of the unknown-input sequence.

> **Two flags serve different purposes.** `g.feasible` reports the numerical target-feasibility test. `g.canRun` reports whether a feasible target **and** a verified stable observer realization are available. Check `g.canRun` before using `g.observer`.

The output equation must have the supplied form `y = C*x`: **unknown-input direct feedthrough is not supported**. The function does not design a noise-optimal filter, reconstruct `d(k)`, or handle time-varying matrices directly.

## Requirements

Save `delayed_uio_design.m` in the current folder or add its folder to the MATLAB path. In the folder containing the file, use:

```matlab
addpath(pwd)
which delayed_uio_design -all
help delayed_uio_design
```

| Gain method | Dependency and behavior |
| --- | --- |
| `'riccati'` | Uses the deterministic dual Riccati iteration implemented in this file. No Control System Toolbox is needed for this branch. |
| `'place'` | Uses `place` from Control System Toolbox for the assignable quotient. |
| `'auto'` | Uses `place` when the function is found on the MATLAB path; otherwise uses `'riccati'` when no poles were specified. |

`'auto'` selects a method based on the availability of `place`; it does **not** catch a subsequent `place` failure and retry with Riccati. Specify `'riccati'` explicitly for an example that does not depend on `place` availability.

The supplied source does not declare a minimum tested MATLAB release or GNU Octave compatibility. It uses numeric MATLAB linear algebra, including `svd`, `schur`, `ordeig`, and `ordschur`; it does not use symbolic computations.

## Examples

### 1. Design a delayed observer

This small example illustrates the interface; it is not a traffic-network benchmark. The first state is driven by an arbitrary unknown input. The measurement reports its one-step-delayed value. The third state has stable, unmeasured dynamics driven by the known input.

```matlab
A = [0 0 0;
     1 0 0;
     0 0 0.6];
B = [1; 0; 1];
Bbar = [1; 0; 0];
C = [0 1 0];
Q = eye(3);                     % Request all three state components.
r = 1;                          % One sample of reconstruction delay.

opts = struct('gainMethod','riccati', 'verbose',true);
g = delayed_uio_design(A,B,Bbar,C,Q,r,opts);

assert(g.canRun, 'A runnable observer was not returned.');

disp(g.dim)
disp(g.design)
disp(g.checks)

% Access the principal design matrices.
Pg   = g.Pg;                    % Quotient-coordinate map.
Cbar = g.Cbar;                  % Horizon output map.
L    = g.L;                     % Stabilizing friend.
Ae   = g.Ae;                    % Quotient error dynamics.
E    = g.E;                     % Dynamic-state recovery map.
F    = g.F;                     % Horizon-output recovery map.
```

For this example, the expected structural results are `g.feasible = true`, `g.canRun = true`, `g.rbar = 1`, `g.dim.maxReconstructible = 3`, `g.dim.observerState = 2`, and `g.dim.assignable = 1`. The unmeasured stable mode at `0.6` remains a fixed mode; it is not a pole that can be reassigned freely.

The individual entries of the returned bases and maps need not match another run or implementation. Their coordinates depend on the numerical basis choices. Compare subspaces, dimensions, spectra, and defining identities instead.

### 2. Run the observer and plot signed errors

Continue from Example 1. First generate noise-free data for illustration. In this example, `X(:,k+1)` stores `x(k)`, `Y(:,k+1)` stores `y(k)`, and `U(:,k+1)` stores `u(k)`.

```matlab
K = 100;                        % State samples are k = 0,...,K.
kInput = 0:K-1;
U = 0.3*sin(0.11*kInput);        % Known inputs u(0),...,u(K-1).
D = 0.4*sin(0.37*kInput) ...     % Unknown inputs: plant generation only.
  + 0.25*cos(0.19*kInput);

X = zeros(g.dim.n,K+1);
X(:,1) = [2; -1; 3];
for k = 0:K-1
    X(:,k+2) = g.A*X(:,k+1) ...
             + g.B*U(:,k+1) ...
             + g.Bbar*D(:,k+1);
end
Y = g.C*X;
```

Now run the observer. This loop reads **only** the design, `Y`, `U`, and the observer's own state. It does not read `X` or `D`.

```matlab
assert(g.canRun, 'Check g.feasible and enable observer design first.');
r = g.r;
K = size(Y,2)-1;
assert(K >= r, 'At least r+1 output samples are required.');
assert(size(Y,1) == g.dim.measuredOutputs, 'Unexpected output dimension.');
assert(size(U,1) == g.dim.knownInputs && size(U,2) == K, ...
       'U must contain the known input samples u(0),...,u(K-1).');

obs = g.observer;
z = zeros(g.dim.observerState,1);  % z_r(r), not an estimate of full x(0).
Qhat = nan(g.dim.targetRows,K+1);  % Leave unavailable initial estimates as NaN.

for k = r:K
    Yblock = Y(:,k-r+1:k+1);      % y(k-r),...,y(k): r+1 samples.
    Ublock = U(:,k-r+1:k);        % u(k-r),...,u(k-1): r samples.
    Ywin = Yblock(:);            % Time-ordered stacking, oldest sample first.
    Uwin = Ublock(:);

    yXi = obs.Hy*Ywin + obs.Hu*Uwin;

    % Produce the estimate BEFORE advancing the observer state.
    Qhat(:,k+1) = obs.E*z + obs.F*yXi;

    if k < K
        uTarget = U(:,k-r+1);    % u(k-r), used only in the next-state update.
        z = obs.Az*z + obs.Bu*uTarget + obs.By*yXi;
    end
end
```

For `r = 0`, `Ublock` is empty. The update uses `u(k)` only **after** producing the estimate at `k`; in an online implementation, perform that update when `u(k)` is available, before producing the next estimate. For `B = []`, store the input history as `U = zeros(0,K)`.

Compare against the correctly delayed target. The true state is used here only to evaluate the simulation error.

```matlab
kEstimate = r:K;
kTarget = kEstimate-r;
qTrue = g.Q*X(:,kTarget+1);
errorSigned = qTrue-Qhat(:,kEstimate+1);

figure;
plot(kEstimate,errorSigned.','LineWidth',1.2);
grid on;
set(gca,'FontSize',9, 'TickLabelInterpreter','latex');
xlabel('Estimation time $k$ (samples)','Interpreter','latex');
ylabel('$Qx(k-r)-\widehat q(k)$','Interpreter','latex');
legend({'$x_1$','$x_2$','$x_3$'}, ...
       'Interpreter','latex','Location','best','FontSize',9);
drawnow;
set(gca,'LooseInset',get(gca,'TightInset'));
```

This is an ordinary **linear-axis plot of signed component errors**. The legend names are specific to Example 1; adjust them when changing the target. The first valid point is at estimation time `k = r`, and it estimates `Q*x(0)`. It is not an estimate available at time zero when `r > 0`.

### 3. Request the maximal reconstructible quotient

Append an unmeasured, unknown-initialized unstable state to Example 1. At delay one, the original three-state target remains reconstructible, but the appended state does not.

```matlab
A4 = blkdiag(A,1.1);
B4 = [B; 0];
Bbar4 = [Bbar; 0];
C4 = [C 0];

% Test full-state feasibility without attempting gain design.
gFull = delayed_uio_design(A4,B4,Bbar4,C4,eye(4),1, ...
                          struct('designObserver',false));
disp(gFull.feasible)             % Expected: false.

% Inspect an algebraic direction that violates the target condition.
v = gFull.witness;
disp(norm(v))                   % Approximately one when infeasible.
disp(norm(gFull.Q*v))            % Nonzero when infeasible.

% Q=[] asks the function to choose coordinates of X/I_r.
gMax = delayed_uio_design(A4,B4,Bbar4,C4,[],1, ...
                         struct('gainMethod','riccati'));
assert(gMax.canRun);
Qmax = gMax.Q;                   % This is gMax.PI, not eye(4).
disp(gMax.dim.maxReconstructible) % Expected: 3.
```

`Qmax*x` represents the maximal reconstructible **linear information** at the selected delay. Its coordinates need not be named original state components. Use `gMax.Q` to define the true target when scoring an estimate.

`gFull.witness` is a state direction in the computed ambiguity subspace. It is **not** a complete zero-output trajectory or an unknown-input sequence; constructing either requires additional steps.

> `Q = []` does **not** request a zero target. Use `Q = zeros(1,n)` for the identically zero scalar target. Any empty `Q`, including a zero-row matrix, takes the maximal-target branch.

### 4. Find the minimum feasible delay for a fixed target

Use a **fixed, explicitly supplied target** while scanning delays. The function reports the horizon saturation index `g.rbar`; it does not search for the minimum target delay automatically.

```matlab
% Use the three-state A,B,Bbar,C,Q from Example 1.
scanOpts = struct('designObserver',false, 'verbose',false);
g0 = delayed_uio_design(A,B,Bbar,C,Q,0,scanOpts);

delays = (0:g0.rbar).';
feasibleAtDelay = false(size(delays));
maxInformation = zeros(size(delays));

for j = 1:numel(delays)
    gj = delayed_uio_design(A,B,Bbar,C,Q,delays(j),scanOpts);
    feasibleAtDelay(j) = gj.feasible;
    maxInformation(j) = gj.dim.maxReconstructible;
end

firstFeasible = find(feasibleAtDelay,1,'first');
rStar = Inf;
if ~isempty(firstFeasible)
    rStar = delays(firstFeasible);
end

disp(table(delays,feasibleAtDelay,maxInformation))
fprintf('Minimum feasible delay for this target: %g samples\n',rStar);
```

For Example 1, the expected minimum delay is `1`. The maximal information dimensions at delays `0` and `1` are `2` and `3`, respectively.

`g.rbar` is the saturation index of the horizon-kernel recursion, **not** generally the minimum delay for `Q`. In the underlying geometric construction, a fixed target that fails at every delay from zero through saturation cannot become feasible at a later delay. This implementation evaluates that construction numerically, so inspect the numerical diagnostics when interpreting a scan.

Do not use `Q = []` to find the minimum delay of a fixed target: it selects a potentially different target at each delay. Also, an exception during a scan is a numerical/design failure, not a returned infeasibility result.

### 5. Select assignable poles explicitly

Only the quotient $\mathscr X/\mathscr S_r^{\ast}$ has freely assignable poles. Its dimension is `g.dim.assignable`, which can differ from `g.dim.observerState`.

```matlab
geom = delayed_uio_design(A,B,Bbar,C,Q,1, ...
                          struct('designObserver',false));
na = geom.dim.assignable;

poles = [];
if na == 1
    poles = 0.3;
elseif na > 1
    poles = linspace(0.2,0.6,na);
end

optsPlace = struct('gainMethod','place', 'assignablePoles',poles);
gPlace = delayed_uio_design(A,B,Bbar,C,Q,1,optsPlace);

disp(gPlace.spectra.assignable)
disp(gPlace.spectra.good)         % Retained fixed Schur modes.
disp(gPlace.spectra.observer)     % Complete quotient-observer spectrum.
```

A nonempty assignable quotient requires `place` for this example. Specified poles must lie strictly inside the unit circle, and complex poles must appear in conjugate pairs. When there are no assignable modes, leave `assignablePoles` empty; no call to `place` is needed.

### 6. Print and save the design

```matlab
format long g
names = {'Cbar','PW','PS','Pg','PI','Lbase','L','Af','Ae','E','F'};
for j = 1:numel(names)
    fprintf('\ng.%s =\n',names{j});
    disp(g.(names{j}));
end

disp(g.dim)
disp(g.design)
disp(g.spectra)
disp(g.diagnostics)
disp(g.checks)
save('delayed_uio_design_result.mat','g');
```

Save the whole `g` structure rather than combining a projection from one design call with gains from another. These matrices share one coordinate choice within a design result.

## Input Arguments

Let `n` be the state dimension, `ell` the number of known inputs, `m` the number of supplied unknown-input columns, `p` the number of measured outputs, and `q` the number of rows of the actual target.

| Argument | Size / accepted value | Description |
| --- | --- | --- |
| `A` | Nonempty `n`-by-`n` matrix | Discrete-time state-transition matrix. |
| `B` | `n`-by-`ell` matrix or `[]` | Known-input matrix. `[]` becomes `zeros(n,0)`. |
| `Bbar` | `n`-by-`m` matrix or `[]` | Unknown-input matrix. `[]` becomes `zeros(n,0)`. Redundant columns are allowed. |
| `C` | `p`-by-`n` matrix or `[]` | Measured-output matrix. `[]` becomes `zeros(0,n)`, meaning no measurements. |
| `Q` | Nonempty `q`-by-`n` matrix or `[]` | Prescribed target map. `eye(n)` requests the full state; `[]` selects `g.PI` as the target. |
| `r` | Nonnegative integer scalar | Reconstruction delay in samples, including zero. |
| `opts` | Scalar structure, `[]`, or omitted | Design and numerical options. Empty or omitted uses all defaults. |

`A`, `B`, `Bbar`, `C`, and a nonempty `Q` must be finite, real, numeric, two-dimensional matrices. Internally they are converted to **full double-precision** matrices. Sparse storage is therefore not preserved. The implementation does not require `Bbar` to have full column rank or require `Q` to have independent rows. `g.dim.targetRows` counts target rows, not their rank.

The geometry uses the column space of `Bbar`, but `g.Tunknown` retains the **original supplied columns**. Thus `g.dim.unknownInputColumns` and `g.dim.unknownInputRank` can differ.

### opts — Design options

Supply only the fields to change. Names are case-sensitive structure field names; unrecognized fields raise an error.

```matlab
opts = struct;
opts.gainMethod = 'riccati';
opts.designObserver = true;
opts.verbose = true;
```

| Field | Default | Description |
| --- | --- | --- |
| `rankTol` | `1e-10` | Relative/absolute-floor singular-value cutoff for bases, ranks, null spaces, and pseudoinverses: `rankTol*max(1,sigmaMax)`. |
| `checkTol` | `1e-7` | Tolerance for subspace agreement and normalized defining-identity residuals. |
| `targetTol` | `1e-8` | Threshold for `g.targetResidual`, the target inclusion test. |
| `unitCircleTol` | `1e-8` | Fixed eigenvalues with modulus at least `1-unitCircleTol` are conservatively classified as bad. |
| `designObserver` | `true` | Design and verify stabilized quotient dynamics. `false` skips this stage but still computes geometry and feasible recovery maps. |
| `gainMethod` | `'auto'` | `'auto'`, `'place'`, or `'riccati'`. |
| `assignablePoles` | `[]` | Exactly `g.dim.assignable` desired poles for pole placement. Empty uses `0.4` for one assignable mode, or `linspace(0.2,0.7,na)` for `na > 1`. |
| `riccatiTol` | `1e-12` | Relative successive-iterate stopping tolerance for the built-in Riccati iteration. |
| `maxRiccatiIter` | `30000` | Maximum Riccati iteration count. |
| `verbose` | `false` | Print dimensions, target feasibility, run readiness, selected gain method, and maximum verification residual. |
| `warn` | `true` | Enable the fixed-pole unit-circle and ill-conditioned horizon-map warnings. Does not disable validation errors. |

`rankTol`, `checkTol`, `targetTol`, and `riccatiTol` must be finite real scalars strictly between zero and one. `unitCircleTol` must be in `[0,0.1)`. `maxRiccatiIter` must be a positive integer. Boolean options accept logical scalars or numeric `0`/`1`.

`assignablePoles` must be a finite numeric vector strictly inside the unit circle. Pole count and conjugate-pair requirements are checked when pole placement is performed. Nonempty specified poles cannot be used with the Riccati branch. With `designObserver = false`, pole assignment is not performed.

The Riccati branch is an optional deterministic stabilization method with identity weights. It does **not** interpret the unknown input as white noise, take input covariance data, or optimize the final target-error variance. There are no user-supplied Riccati weight options in this implementation.

## Output Arguments

### g — Design result

A scalar structure containing the converted plant data, the actual target, the requested delay, filled-in options, geometric objects, gains, and diagnostics.

#### Feasibility and run readiness

| Field | Meaning |
| --- | --- |
| `g.feasible` | Whether the numerical test `g.targetResidual <= g.options.targetTol` passes. |
| `g.canRun` | `g.feasible && g.options.designObserver && g.design.stable`. |
| `g.Q` | The supplied nonempty target, or `g.PI` when maximal-target mode was requested. |
| `g.maximalTargetRequested` | Whether the input `Q` was empty. |
| `g.targetResidual` | `norm(g.Q*g.I,'fro')/norm(g.Q,'fro')`; defined as zero when `norm(g.Q,'fro')` is zero. |
| `g.witness` | An `n`-by-`1` unit ambiguity direction violating the target condition when infeasible; otherwise `n`-by-`0`. |
| `g.rbar` | Numerically detected saturation index of the horizon-kernel recursion. |

For a successful return, the main cases are:

| Target feasible? | `designObserver` | Result |
| --- | --- | --- |
| Yes | `true` | Stabilized dynamics, recovery maps, and `g.observer` are available; `g.canRun = true`. |
| Yes | `false` | Geometry and recovery maps are available. `g.L` and `g.Ae` are not designed; `g.observer = []` and `g.canRun = false`. |
| No | Either | `g.E = []`, `g.F = []`, and `g.observer = []`; a witness is returned and `g.canRun = false`. Stable quotient dynamics may still have been designed when requested. |

Gain design occurs **before** the target-feasibility test. Use `designObserver = false` for a feasibility-only query, so a gain-design failure does not prevent that query from returning.

A zero-dimensional quotient is legitimate. Some feasible maps can consequently be empty matrices. Use the flags, not `isempty(g.E)` or an assumed positive observer dimension, to decide whether the design is usable.

#### Principal matrices

Write `pr = g.dim.horizonOutput`, `nz = g.dim.observerState`, and `hY = (r+1)*p`.

| Field | Size | Meaning |
| --- | --- | --- |
| `g.Cbar` | `pr`-by-`n` | Coordinates of the canonical horizon projection $P_{\mathscr V_r}$. |
| `g.Pg` | `nz`-by-`n` | Coordinates of $P_{\mathscr W_{g,r}^{\ast}}$. |
| `g.PI` | `(n-g.dim.I)`-by-`n` | Coordinates of $P_{\mathscr I_r}$; maximal reconstructible target map. |
| `g.PW` | `(n-g.dim.W)`-by-`n` | Coordinates of $P_{\mathscr W_r^{\ast}}$. |
| `g.PS` | `(n-g.dim.S)`-by-`n` | Coordinates of $P_{\mathscr S_r^{\ast}}$. |
| `g.Lbase` | `n`-by-`pr` | Initial friend of $\mathscr W_r^{\ast}$; computed even in geometry-only mode. |
| `g.L` | `n`-by-`pr` | Final stabilizing friend when observer design is enabled. |
| `g.Ae` | `nz`-by-`nz` | Induced quotient dynamics and quotient-error transition matrix. |
| `g.E` | `q`-by-`nz` | Recovery map multiplying the observer state. |
| `g.F` | `q`-by-`pr` | Recovery map multiplying the horizon output. |
| `g.Hwindow` | `pr`-by-`hY` | Converts the corrected output window to the effective horizon output. |
| `g.Psi` | `(nz+pr)`-by-`n` | Stacked map `[g.Pg; g.Cbar]` used to construct recovery maps. |

The following identities are checked, subject to the selected tolerances and the relevant feasibility/design flags:

```matlab
g.Pg*g.Bbar                           % Approximately zero.
g.Pg*(g.A + g.L*g.Cbar) - g.Ae*g.Pg  % Approximately zero when designed.
g.E*g.Pg + g.F*g.Cbar - g.Q          % Approximately zero when feasible.
```

> **A quotient map is not a square projector.** `g.Pg*x` is an `nz`-dimensional quotient coordinate. Its kernel is the column space of `g.Wg`. In this orthonormal realization, `g.Pg'*g.Pg` is the square orthogonal projector onto the orthogonal complement of that space; `g.Wg*g.Wg'` projects onto the space itself. Do not invert `g.Pg` to obtain the full state. Full-state recovery, when feasible, uses both `g.E` and `g.F`.

#### g.observer — Ready-to-use realization

This is a structure only when `g.canRun` is true; otherwise it is `[]`.

```matlab
yXi  = g.observer.Hy*Ywin + g.observer.Hu*Uwin;
qhat = g.observer.E*z + g.observer.F*yXi;
zNext = g.observer.Az*z ...
      + g.observer.Bu*uTarget ...
      + g.observer.By*yXi;
```

Here `uTarget = u(k-r)`, `z = z_r(k)`, `qhat` estimates `Q*x(k-r)`, and `zNext = z_r(k+1)`.

| Field | Definition |
| --- | --- |
| `Az` | `g.Ae` |
| `Bu` | `g.Pg*g.B` |
| `By` | `-g.Pg*g.L` |
| `E`, `F` | `g.E`, `g.F` |
| `Hy` | `g.Hwindow` |
| `Hu` | `-g.Hwindow*g.Tknown` |
| `BYwindow` | `By*g.Hwindow` |
| `BUwindow` | `-By*g.Hwindow*g.Tknown` |
| `DYwindow` | `F*g.Hwindow` |
| `DUwindow` | `-F*g.Hwindow*g.Tknown` |

The last four fields provide an equivalent raw-window realization:

```matlab
obs = g.observer;
qhat = obs.E*z + obs.DYwindow*Ywin + obs.DUwindow*Uwin;
zNext = obs.Az*z + obs.Bu*uTarget ...
      + obs.BYwindow*Ywin + obs.BUwindow*Uwin;
```

Do not omit `obs.Bu*uTarget`: it is separate from the known-input correction inside the measurement window. The sign of `By` is already included; do not negate it a second time.

#### g.dim — Dimension summary

| Fields | Meaning |
| --- | --- |
| `n`, `knownInputs`, `measuredOutputs`, `targetRows` | State, known-input, measurement, and actual target row counts. |
| `unknownInputColumns`, `unknownInputRank` | Supplied unknown-input column count and numerical rank. |
| `horizonOutput` | Number of rows of `g.Cbar`, not the raw-window length `(r+1)*p`. |
| `Vr`, `Vstar`, `W`, `S`, `Wg`, `I` | Dimensions of the correspondingly named geometric subspaces. |
| `fixed`, `badFixed`, `goodFixed` | Dimension of `S/W` and its bad/good spectral parts. |
| `maxReconstructible` | `n - g.dim.I`, the maximal independent linear information dimension. |
| `observerState` | `n - g.dim.Wg`, the dynamic quotient-observer dimension. |
| `assignable` | `n - g.dim.S`, the number of freely assignable poles. |

`observerState` excludes the storage needed for the input/output window and is not a minimum-order guarantee. It is also not the same quantity as `maxReconstructible`.

<details>
<summary><strong>Advanced return fields: subspace bases, horizon matrices, and spectral coordinates</strong></summary>

### Subspace bases and iteration histories

Each listed basis has orthonormal **columns** in state coordinates and size `n`-by-its-subspace-dimension. The zero subspace is represented by an `n`-by-`0` matrix.

| Field | Represented subspace |
| --- | --- |
| `g.Bbasis` | $\operatorname{Im}\bar B$. |
| `g.Vr` | Horizon output-nulling kernel $\mathscr V_r$. |
| `g.Vstar` | Stationary output-nulling subspace $\mathscr V^{\ast}$. |
| `g.W` | Infimal conditioned-invariant subspace $\mathscr W_r^{\ast}$. |
| `g.S` | Infimal unobservability subspace $\mathscr S_r^{\ast}$. |
| `g.Vfix` | $\mathscr S_r^{\ast}\cap(\mathscr W_r^{\ast})^{\perp}$, a choice of representatives for `S/W`; not generally `Vstar`. |
| `g.Xbad`, `g.Xgood` | State-space lifts of the bad and good spectral subspaces of the fixed quotient map. |
| `g.Wg` | $\mathscr W_{g,r}^{\ast}$, obtained by adjoining the bad fixed directions to `W`. |
| `g.I` | $\mathscr I_r=\mathscr W_{g,r}^{\ast}\cap\ker\bar C_r$. |

`g.Vseq`, `g.Wseq`, and `g.Sseq` are cell arrays containing the successive bases, starting with `ker(C)`, the zero subspace, and the full state space, respectively. An unchanged terminal iterate is not appended again. In particular, `g.Vseq{1}` represents $\mathscr V_0$, and the sequence is not extended by duplicate entries when the requested delay exceeds `g.rbar`.

### Horizon matrices

Let `hY = (r+1)*p`, `hU = r*ell`, `hD = r*m`, and `hOmega = size(g.Omega,1)`.

| Field | Size | Meaning |
| --- | --- | --- |
| `g.O` | `hY`-by-`n` | `col(C, C*A, ..., C*A^r)`. |
| `g.Tknown` | `hY`-by-`hU` | Known-input response over the output window. |
| `g.Tunknown` | `hY`-by-`hD` | Unknown-input response using the supplied `Bbar` columns. |
| `g.Rbasis` | `hY`-by-`rank(Tunknown)` | Basis of the unknown-input reachable output subspace. |
| `g.Omega` | `hOmega`-by-`hY` | Quotient-coordinate map annihilating `g.Tunknown`. |
| `g.J` | `hOmega`-by-`pr` | Injective factor `g.Omega*g.O*g.Cbar'`. |
| `g.Jleft` | `pr`-by-`hOmega` | Numerical left inverse of `g.J`. |
| `g.Hwindow` | `pr`-by-`hY` | `g.Jleft*g.Omega`. |

The window contains **`r+1` output samples but only `r` known-input samples**:

$$
\begin{aligned}
Y_{\mathrm{win}}(k)&=\operatorname{col}(y(k-r),\ldots,y(k)),\\
U_{\mathrm{win}}(k)&=\operatorname{col}(u(k-r),\ldots,u(k-1)).
\end{aligned}
$$

Consequently,

$$
y_{\xi_r}(k)=H_{\mathrm{window}}
\bigl(Y_{\mathrm{win}}(k)-T_{\mathrm{known}}U_{\mathrm{win}}(k)\bigr)
=\bar C_r x(k-r).
$$

For `r = 0`, both input-response matrices have zero columns. Increasing `r` beyond `g.rbar` still constructs the full requested window; saturation of geometric information does not truncate the returned horizon matrices.

### Fixed and assignable quotient dynamics

| Field | Meaning |
| --- | --- |
| `g.Abase` | `g.A + g.Lbase*g.Cbar`. |
| `g.Af` | `g.Vfix'*g.Abase*g.Vfix`, the fixed dynamics on `S/W`. |
| `g.Zbad`, `g.Zgood` | Bases in the fixed quotient's coordinates, selected by separate real Schur reorderings. |
| `g.Tbad`, `g.Tgood` | Reduced matrices paired with `Zbad` and `Zgood`, respectively. |
| `g.PY` | Output-coordinate map whose kernel is the image of `g.Cbar*g.S`. |
| `g.Aassign` | `g.PS*g.Abase*g.PS'`. |
| `g.Cassign` | `g.PY*g.Cbar*g.PS'`. |
| `g.Lassign` | Injection designed on the assignable quotient; zero-initialized and not designed in geometry-only mode. |

The reduced spectral matrices obey `Af*Zbad = Zbad*Tbad` and `Af*Zgood = Zgood*Tgood` up to tolerance. `Xbad` and `Xgood` are separately re-orthonormalized state-space lifts; their columns need not retain exactly the same coordinates as these reduced spectral matrices. For a nonnormal fixed map, the good and bad invariant subspaces need not be mutually orthogonal.

</details>

#### Design and numerical diagnostics

| Structure | Useful fields |
| --- | --- |
| `g.design` | `requestedMethod`, `methodUsed`, `nAssignable`, `nEffectiveOutputs`, `desiredPoles`, `iterations`, `riccatiResidual`, `stable`, `rho`. |
| `g.spectra` | `fixed`, `bad`, `good`, `nearUnitCircle`, `assignable`, `observer`. |
| `g.diagnostics` | `maxResidual`, `windowNorm2`, `Jcondition`, `recoveryNorm2`, `boundaryClassification`, `numericalOnly`. |
| `g.checks` | Normalized residuals of the defining identities. |
| `g.absoluteChecks` | Corresponding unnormalized Frobenius residuals. |

`g.design.rho` is the spectral radius of the complete designed quotient matrix `g.Ae`, not of the plant matrix `A`. An empty quotient has radius zero by convention. In geometry-only mode, `g.design.stable` is false, `rho` is `NaN`, and `methodUsed` is `'not-designed'`; this does not mean that the target is infeasible.

`g.design.riccatiResidual` is the relative difference between the last two Riccati iterates, not a separately evaluated algebraic Riccati equation residual. `g.spectra.assignable` and `g.spectra.observer` are empty when observer design is disabled. `g.diagnostics.numericalOnly` is always true.

Important checks include `horizonFactorization`, `windowState`, `windowUnknown`, `leftInverse`, `unknownDecoupling`, `generalDecomposition`, `badInvariant`, `goodInvariant`, and `projectionRows`. Observer design additionally checks `friend`, `inducedMap`, `fixedUnchanged`, and `assignableLift`; target feasibility additionally enables `recovery`.

Each algebraic residual is the Frobenius norm divided by `max(1,scale)`, with the scale chosen for that identity in the source. A residual above `checkTol`, or a nonfinite residual, raises an error instead of returning a successful design.

## Algorithms

The implementation follows six stages.

1. **Horizon output-nulling spaces.** Compute
   
   $$
   \mathscr V_0=\ker C,\qquad
   \mathscr V_{j+1}=\ker C\cap A^{-1}
   (\mathscr V_j+\operatorname{Im}\bar B)
   $$
   
   until stationarity. Here $A^{-1}$ means a subspace preimage, not a matrix inverse. Select $\mathscr V_r$ and construct its quotient-coordinate map `Cbar`.

2. **Window realization.** Build the initial-state, known-input, and unknown-input response matrices. Quotient out the unknown-input output subspace and construct `J`, `Jleft`, and `Hwindow`.

3. **Conditioned-invariant and unobservability spaces.** Compute `W` from the increasing recursion starting at zero, and `S` from the decreasing recursion starting at the full state space. Construct `PW` and `PS` from the final bases.

4. **Fixed spectral decomposition.** Construct `Lbase` and `Af`. Use ordered real Schur decompositions to obtain the bad and good invariant subspaces separately. Lift the bad directions to form `Wg`, and intersect with `ker(Cbar)` to obtain `I`.

5. **Stabilization on the assignable quotient.** Design `Lassign` only on `X/S`, then lift it through `PS'`:
   
   $$
   L=L_{\mathrm{base}}+P_S^{\mathsf T}L_{\mathrm{assign}}P_Y,
   \qquad
   A_e=P_g(A+L\bar C_r)P_g^{\mathsf T}.
   $$
   
   The implementation uses the **plus-sign convention** `A + L*Cbar`. With pole placement, `Lassign = -place(Aassign',Cassign',poles)'`. Stable fixed modes remain fixed; the bad fixed component is removed in the quotient construction.

6. **Target recovery and verification.** Test the inclusion using `Q*I`. For a feasible target, split `Q*pinv([Pg; Cbar])` into `E` and `F`, using the implementation's singular-value threshold. Verify the defining identities before returning the result.

For a runnable design, with $\varepsilon_r(k)=P_gx(k-r)-z_r(k)$, the intended exact-arithmetic identities are

$$
\varepsilon_r(k+1)=A_e\varepsilon_r(k),\qquad
Qx(k-r)-\widehat q(k)=E\varepsilon_r(k).
$$

## Numerical Considerations

**Finite precision and scaling.** This is a numerical implementation, not an exact symbolic decision procedure. Subspace dimensions, target feasibility, and pole classification depend on tolerances. The rank cutoff has an absolute floor through `max(1,sigmaMax)`, so changes of physical units can affect rank decisions. For poorly scaled data, use consistent state/input/output scaling and investigate the singular values rather than simply loosening every tolerance.

**Near-unit-circle fixed poles.** The conservative classification can include a mathematically Schur pole lying within `unitCircleTol` of the unit circle. This can enlarge the computed ambiguity space and reduce the reported reconstructible information. Inspect `g.spectra.fixed`, `g.spectra.nearUnitCircle`, and `g.diagnostics.boundaryClassification` before interpreting a boundary case as a precise impossibility result.

**Conditioning versus guarantees.** `windowNorm2`, `Jcondition`, and `recoveryNorm2` help identify potentially sensitive inversions. They are not closed-loop noise-performance certificates. The warning threshold on `Jcondition` is `1e10`; the code does not separately threshold `windowNorm2` to emit a large-gain warning.

**Large horizons.** The raw input-response matrices grow in both row and column count with `r`, and the source converts matrices to dense double precision. Long windows can therefore require substantial memory and can amplify numerical issues from powers of `A`. The code checks for nonfinite horizon entries and disagreement between horizon and recursive ranks.

**Scope of the result.** The implementation assumes an exact known LTI model and no unknown-input direct feedthrough. It does not add robustness to measurement noise or modeling error, assume unknown-input statistics, impose physical constraints on the unknown inputs, guarantee minimum observer order, or generate unknown-input estimates. When an application is embedded in an unrestricted-input model, feasibility applies to trajectories of that model; a numerical infeasibility result should not automatically be interpreted as impossibility for a smaller, additionally constrained trajectory class.

## Troubleshooting

| Symptom / identifier | Meaning and action |
| --- | --- |
| Function not found | Put the `.m` file in the current folder or add its folder to the path. Use `which delayed_uio_design -all` to check for duplicate versions. |
| `g.feasible == false` | The computed ambiguity space affects the requested target. Inspect `g.targetResidual` and `g.witness`; scan feasible delays or change the target/sensing model. This is not itself a thrown error. |
| `g.feasible == true`, `g.canRun == false` | On a successful return, observer design was disabled. Rerun with `designObserver = true` to obtain `g.observer`. |
| `delayed_uio:Arguments`, `:Matrix`, `:Dimensions`, or `:Delay` | Check the six required inputs, finite real numeric matrices, compatible sizes, and nonnegative integer delay. |
| `delayed_uio:Options` | Check exact structure field names, option types, and admissible ranges. |
| `delayed_uio:Toolbox` | Pole placement is unavailable. Use `'riccati'` with empty `assignablePoles`, or make `place` available. |
| `delayed_uio:PoleCount` | Supply one pole per assignable mode, not one per observer state. Query `g.dim.assignable` in geometry-only mode. |
| `delayed_uio:Poles` or `:ComplexGain` | Check pole locations, conjugate pairs, and compatibility with the selected method. |
| `delayed_uio:HorizonOverflow` or `:HorizonRank` | Inspect state scaling, horizon length, and numerical rank thresholds. |
| `delayed_uio:Rank`, `:Recursion`, `:QuotientDimension`, or `:LiftRank` | A subspace iteration or dimension consistency check failed. Investigate scaling and near-dependent directions. |
| `delayed_uio:Riccati` | The built-in iteration overflowed or did not converge within the limit. Inspect scaling and the assignable pair before changing iteration settings. |
| `delayed_uio:Assignable` or `:Unstable` | The numerical assignable/observer quotient failed a stability or effective-output check. Inspect the design data; do not treat the error as an infeasibility certificate. |
| `delayed_uio:Residual` | A defining identity exceeded `checkTol`. The message identifies the failed check. No verified result is returned. |
| `delayed_uio:UnitCircle` warning | Fixed poles near the unit circle were conservatively retained in the bad component. |
| `delayed_uio:WindowCondition` warning | The computed horizon factor has a large condition number. |
| Errors do not converge in a custom simulation | Check data stacking, known-input subtraction, delayed target indexing, the sign of `By`, and that `qhat` is evaluated before advancing `z`. Also check model/noise assumptions. |

### Common indexing mistakes

Do not compare `Qhat(:,k+1)` with `Q*X(:,k+1)` when `r > 0`; compare it with `Q*X(:,k-r+1)`. Do not stack each output channel's complete time history before the next channel: MATLAB's `Yblock(:)` correctly stacks the full output vector at each successive time. Do not replace the unavailable estimates at `k < r` with artificial zeros.

### Source and documentation scope

This page documents the accompanying `delayed_uio_design.m`, including its implemented defaults, numerical decisions, and return fields. The examples are usage illustrations, not captured MATLAB session output or a claim of cross-release testing. The MATLAB source is unchanged by this documentation.
