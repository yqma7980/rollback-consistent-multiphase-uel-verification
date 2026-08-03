# WP10-I5 Pre-registration: Six-block, Nine-block, and Full-matrix P1 Directional Verification

## 1. Scope and claim boundary

WP10-I5 verifies the integrated M2a Jacobian assembled by the unchanged I4 source at the registered P1 element states. It does not modify the UEL and does not test an Abaqus Standard model, nonlinear convergence order, full-domain behavior, production readiness, or thermodynamic consistency.

The source under test is:

`paper1_methods/verification/WP10/I5/src/UEL_TWOPHASE_UPS_WP10I5_integratedDirectionalVerification.for`

Its required SHA-256 is:

`DCD660846A5BB6498181B94F6F24D5B2D5294FD33D6B0D6C0914983E2B2A85CA`

This hash must equal the I4 source hash. Any source difference blocks I5.

## 2. Discrete map and sign

The primary vector and residual are

\[
x=[u,p_w,S_n], \qquad R=[R_U,R_P,R_S].
\]

The 16 element degrees of freedom are ordered node-wise as

`[u_x,u_y,p_w,S_n]`.

The registered indices are:

| Field | Indices |
|---|---|
| `u` | 1, 2, 5, 6, 9, 10, 13, 14 |
| `p_w` | 3, 7, 11, 15 |
| `S_n` | 4, 8, 12, 16 |

The nine physical blocks are:

| Block | Residual rows | Trial columns |
|---|---|---|
| `KUU` | `R_U` | `u` |
| `KUP` | `R_U` | `p_w` |
| `KUS` | `R_U` | `S_n` |
| `KPU` | `R_P` | `u` |
| `KPP` | `R_P` | `p_w` |
| `KPS` | `R_P` | `S_n` |
| `KSU` | `R_S` | `u` |
| `KSP` | `R_S` | `p_w` |
| `KSS` | `R_S` | `S_n` |

One global sign is frozen for every block and full-matrix direction:

\[
D_h(\mathrm{RHS})[d]\approx -J_{\mathrm{phys}}d.
\]

No block-specific or direction-specific sign selection is permitted.

## 3. Physical and returned matrices

The source returns

\[
J_{\mathrm{returned}}=J_{\mathrm{phys}}+J_{\mathrm{num}},
\qquad J_{\mathrm{num}}=10^{-12}I_{16}.
\]

The I5 physical matrix is formed only in the harness:

\[
J_{\mathrm{phys}}=J_{\mathrm{returned}}-10^{-12}I_{16}.
\]

The UEL source is not changed. The I4 KPP decomposition remains a regression control; the I5 block and full-matrix consistency classifications use `J_phys`.

## 4. Registered states and branches

The 13 I3/I4 states are retained without alteration:

`I0`, `I1`, `I2`, `LIMITER_OFF`, `LIMITER_INACT`, `RATIO_POS`, `RATIO_NEG`, `BETA_MIN_POS`, `BETA_MIN_NEG`, `FINAL_CLIP_LOW`, `FINAL_CLIP_HIGH`, `SE_CLIP_LOW`, and `SE_CLIP_HIGH`.

`I0`, `I1`, and `I2` are the core smooth states. The remaining states are branch-preserving regression states. Exact switching states are excluded from smooth central-difference claims. Every plus/minus evaluation records its branch signature; branch-crossed rows are retained and are not counted as smooth passes.

## 5. Directions and characteristic scales

The field directions are copied from I3/I4:

- `D1`: normalized uniform direction;
- `D2`: normalized alternating-sign direction;
- `D3`: normalized dense direction generated with seed `20260723`.

The characteristic scales are frozen as:

\[
X_u=10^{-5}\ \mathrm{m},\qquad X_p=10^7\ \mathrm{Pa},\qquad X_S=10^{-1}.
\]

Three mixed full-matrix directions are pre-registered in scaled coordinates:

| Mixed direction | `u` field | `p_w` field | `S_n` field |
|---|---|---|---|
| `M1` | `D1` | `D2` | `D3` |
| `M2` | `D2` | `D3` | `D1` |
| `M3` | `D3` | `D1` | `D2` |

Each field contribution is divided by `sqrt(3)`, so every mixed vector has unit Euclidean norm in the concatenated scaled coordinate system. Actual scaled and physical components must be saved; the seed alone is insufficient.

## 6. Perturbation sequence

Central differences retain all eleven non-dimensional step sizes:

`1E-2, 3E-3, 1E-3, 3E-4, 1E-4, 3E-5, 1E-5, 3E-6, 1E-6, 3E-7, 1E-7`.

For every sign, state, direction, block, and step size:

1. Start from the same committed `SVARS` copy.
2. Hold `TIME`, `DTIME`, `PROPS`, `JPROPS`, `COORDS`, and old state fixed.
3. Perturb only the declared trial primary variables.
4. Preserve Abaqus `U`/`DU` semantics.
5. Discard candidate `SVARS` after the call.
6. Never pass a trial result to another sign, step, direction, or state.
7. Scan the used RHS, AMATRX, and physical candidate `SVARS` for finite values.

## 7. Six-block, nine-block, and full-matrix sets

The registered six-block set is the original WP09 flow-coupling set:

`KPU, KPP, KPS, KSU, KSP, KSS`.

The nine-block set additionally contains:

`KUU, KUP, KUS`.

The full-matrix test perturbs all 16 primary degrees of freedom simultaneously with `M1`-`M3` and compares all 16 residual components with `-J_phys d`.

## 8. Norms and scaling

Blockwise non-null comparisons use unscaled block L2 norms and

\[
e_{rel}=\frac{\|D_hR+J_{phys}d\|_2}
{\max(\|D_hR\|_2,\|J_{phys}d\|_2,10^{-30})}.
\]

For mixed full-matrix comparisons, residual rows have different units. A baseline, state-specific row scale is therefore frozen before perturbations:

\[
y_i=\max\left[
\left(\sum_j(J_{phys,ij}X_j)^2\right)^{1/2},10^{-30}
\right],
\]

where `X_j` is `X_u`, `X_p`, or `X_S` according to the column field. Both finite-difference and Jacobian-vector residuals are divided componentwise by the same `y_i`. This scale depends only on the baseline physical matrix and registered variable scales; it is not tuned to finite-difference results.

## 9. Acceptance gates

For non-null block and full-matrix series:

- relative L2 error `<=1E-5`;
- cosine similarity `>=0.99999`;
- at least two adjacent `h` values pass;
- a central-difference descent region is visible before the roundoff floor;
- no branch crossing or non-finite value;
- independent run normalized difference `<=1E-13`, with exact equality reported when available.

For registered nullspace series:

- scaled FD norm `<=1E-10`;
- scaled Jacobian-vector norm `<=1E-10`;
- scaled absolute mismatch L2 and infinity norms `<=1E-10`;
- at least two adjacent `h` values pass.

`KSP-D1` uses the I4 nullspace-aware rule. A near-zero series is never rejected solely because its relative error or cosine is ill-conditioned.

The full-matrix core passes only if all `I0`-`I2`/`M1`-`M3` series pass. Branch-state mixed directions must either pass the same gates or receive a pre-registered, evidence-backed nullspace/non-smooth classification; unexplained branch-state failures block I5.

## 10. Replay, parity, and repeatability

Required exact gates are:

- equal-state replay after mixed rejected-trial ordering;
- RHS, returned AMATRX, physical candidate `SVARS`, closure values, and branch IDs;
- I4/I5 same-state raw parity;
- I4/I5 source identity;
- run1/run2 CSV values and keys when the floating environment is identical.

All six-block, nine-block, and full-matrix records must retain every registered `h`. No best-step filtering or failure-row deletion is allowed.

## 11. Classification and completion

Physical blocks may be classified as `CONSISTENT_NUMERICALLY`, `CONSISTENT_NUMERICALLY_NULLSPACE_AWARE`, `KNOWN_OMISSION_REMAINS`, `NONSMOOTH_BRANCH_ONLY`, `STATE_LIFECYCLE_FAILURE`, `NONFINITE_FAILURE`, or `INCONCLUSIVE`.

The returned KPP remains separately classified as `APPROXIMATE_BOUNDED_BY_JNUM`; this does not downgrade the physical KPP classification.

I5 is `DONE` only if:

1. the source is bitwise identical to I4;
2. replay and I4/I5 parity are exact;
3. all six physical blocks have complete, repeatable classifications;
4. all nine physical blocks have complete, repeatable classifications;
5. all core full-matrix series pass;
6. all branch-state mixed series are passed or explained without post-hoc filtering;
7. all 11 step sizes are retained;
8. all numerical records are finite;
9. no sign, DOF, branch, or lifecycle ambiguity remains.

Passing I5 authorizes only a decision on I6. It does not establish quadratic Newton convergence, full-domain readiness, production readiness, or thermodynamic consistency.

## 12. Hard stops

Stop without entering I6 if any source hash differs, replay or parity fails, a DOF/sign map is ambiguous, an unregistered matrix-only term prevents `J_phys` construction, non-target source content changes, a required state cannot be reproduced, values become non-finite, branch crossings dominate a registered smooth series, full-matrix results are not independently repeatable, or passing would require changing thresholds, directions, states, residuals, or the UEL.
