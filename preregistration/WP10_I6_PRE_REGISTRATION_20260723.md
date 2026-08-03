# WP10-I6 Pre-registration: Patch, Rollback, and Newton Order

Date: 2026-07-23  
Task: `WP10_I6_patchRollback_NewtonConvergenceOrder_P1`  
Status: pre-registered before any I6 numerical result

## 1. Scope and claim boundary

WP10-I6 tests the unchanged M2a residual and branchwise physical Jacobian on
P1 deterministic benchmarks. It does not modify or validate thermodynamic
closure, exact switching points, full-domain behavior, Abaqus production
readiness, long-time integration, or M3 fixed-stress coupling.

The primary Newton operator is

```text
D(RHS)[d] = -J_phys d
J_returned = AMATRX
J_num      = 1E-12 I
J_phys     = J_returned - J_num
g(x)       = RHS(x) - RHS(x_target)
J_phys dx  = g(x)
x_next     = x + dx
```

`M2_PHYS` is the only profile eligible for a near-quadratic claim.
`M2_RETURNED` is a sensitivity result. `M1_PHYS` is a quasi-Newton comparator.

## 2. Frozen inputs

| Input | Required SHA-256 |
|---|---|
| M2/I5 UEL and bitwise-identical I6 copy | `DCD660846A5BB6498181B94F6F24D5B2D5294FD33D6B0D6C0914983E2B2A85CA` |
| I5 integrated harness | `4C5C11B486A4AEECDB1E8E22FEAB097CD19EDC1228AEA990994F2A20AC8BFA06` |
| I5 final summary | `C94343334440EF0EFAD5BE6E88AB4FBADEFFCEE542DC58714E04A1A3534E017B` |
| I5b summary | `39DA427F4AEF480028D9D20065800EFC2787F9917AD4B8FDA3ED879C3DF38759` |
| M1 source | `DF26D71744945ACD2693FCA95A9CB91713773D6989474C4439DEE196685A7675` |
| Frozen I5 mixed directions | `35038E483D7A05DF9FC49F696F732800B8D399C97CE808FCD7478201A7029579` |
| Frozen I5 states | `3610CE0DCE34B4B894844B1B17F3A97D7F7D893B77112CCE21059DFD32F42E37` |

The I6 UEL must remain byte-for-byte identical to I5. Only standalone
harnesses, deterministic assembly/linear algebra, audits, and plots may be new.

## 3. Environment and deterministic profile

Use Intel Fortran Classic with the exact I5 setup and flags
`/Od /traceback /check:bounds /extend-source:132 /fpp`. Set
`OMP_NUM_THREADS=1`, `MKL_NUM_THREADS=1`, and `ABAQUS_CPUS=1`. Reuse the I5
M2a environment: all optional capacity/anchor/diagonal, PC slope, front,
state-dependent filter, and output controls are off; active mobility is gated
past the test increment; CAP_REF force is one; local Picard is off.

## 4. DOF, scales, and finite gates

The local order is `(u_x,u_y,p_w,S_n)` at each of four Q4 nodes. Local index
sets are `u=(1,2,5,6,9,10,13,14)`, `p=(3,7,11,15)`, and
`S=(4,8,12,16)`. Characteristic variable scales are
`X_u=1E-5 m`, `X_p=1E7 Pa`, and `X_s=1E-1`.

All used RHS, AMATRX, physical SVARS, derived norms, and iterates must be
finite with magnitude at most `1E100`. Two-process normalized differences are
at most `1E-13`; exact equality is reported wherever achieved.

## 5. Four-element assembled patch

The patch contains nine nodes at integer coordinates on `[0,2] x [0,2]` and
four counter-clockwise elements `(1,2,5,4)`, `(2,3,6,5)`, `(4,5,8,7)`, and
`(5,6,9,8)`. The central perturbation is `h=1E-4`.

Frozen modes are `TX`, `TY`, `ROT`, `EXX`, `EYY`, `GXY`, `P_UNIFORM`,
`P_AFFINE_X`, `P_AFFINE_Y`, `S_UNIFORM`, `S_AFFINE_X`, `S_AFFINE_Y`, and
`MIXED_AFFINE`. Null/rigid modes use scaled L2 and infinity gates `1E-10`.
Non-null modes require relative L2 `<=1E-8` and cosine `>=0.999999`.
Term/assembly closure is `<=1E-12`. Every plus/minus state must preserve the
baseline branch. Single-element I2/I3 subpatches are regression evidence only.

## 6. Rollback histories

`H0` is clean equal-state replay. `H1` permutes two finite rejected trials.
`H2` commits two half increments. `H3` first evaluates two oversized full-step
trials without commit and then performs the identical two half increments.
`H4` inserts a discarded `KINC=0` terminal call. `H5` replays the declared
state after mixed patch/Newton activity. Candidate SVARS are copied to the
committed source only after an accepted synthetic increment.

For the same declared state, RHS, raw AMATRX, J_phys, closure slots, CAP_REF,
physical SVARS, and branch IDs should be bitwise exact. The fallback normalized
gate is `1E-13`. Any failure is `STATE_LIFECYCLE_FAILURE`.

## 7. M1/M2 parity

Compare all 13 I5 states, all patch baselines, and both Newton targets. RHS,
closure values, K_abs, p_eq, p_c, mobility values, branch IDs, and physical
candidate SVARS require exact parity. Only pre-registered Jacobian entries may
differ. Any value or state drift is a hard stop.

## 8. Newton benchmarks

`N1` uses I5 state `I1`, mixed direction `M2`, primary scaled amplitude
`0.75`, and near-root diagnostic amplitude `0.15`. `N2` uses I5 state `I2`,
mixed direction `M3`, primary amplitude `1.00`, and near-root diagnostic
amplitude `0.20`. These values are frozen before execution.

Each target load is manufactured as `RHS(x_target)`. The old primary state,
TIME, DTIME, PROPS, committed SVARS, and target branch remain fixed during a
solve. Fixed local DOFs are `1,2,3,4,6`; free local DOFs are
`5,7,8,9,10,11,12,13,14,15,16`. All initial offsets are projected onto this
free set and renormalized in scaled coordinates. The initial state is frozen as
`x_0=x_target-a d`, where `a` is the registered primary or near-root amplitude;
the sign is not selected after observing a convergence history.

The reduced scaled matrix uses fixed target row scaling and the frozen variable
scales. It must be full rank with condition number `<=1E12`; each deterministic
pivoted solve must have relative residual `<=1E-12`. No line search, damping,
trust region, PNEWDT, cutback, or matrix switching is allowed. Maximum Newton
updates are 20.

## 9. Newton-order gates

Use the frozen scaled total and RU/RP/RS residual norms. Define

```text
q_k = log(r_(k+1)/r_k) / log(r_k/r_(k-1))
C_k = r_(k+1) / r_k^2
```

The primary N1 and N2 M2_PHYS sequences each require at least five valid
residual states, at least three consecutive asymptotic `q_k>=1.8`, finite
`C_k` with `max(C_k)/min(C_k)<=10` over that sequence, no branch crossing or
globalization, final scaled residual `<=1E-12`, and scaled root error
`<=1E-10`. Values below the `1E-14` scaled residual floor are retained but are
not used to manufacture an order estimate. Insufficient iterations are
`INSUFFICIENT_ASYMPTOTIC_ITERATIONS`, not permission to alter the initial state.

## 10. Completion and hard stops

I6 passes only if all assembled patch classes pass, H0-H5 pass, M1/M2 value
parity is exact, N1 and N2 M2_PHYS pass the registered Newton-order gate, two
independent runs repeat, and all evidence is finite. Stop on source drift,
residual/state drift, replay failure, rank/conditioning failure, branch
crossing, linear-solve failure, missing asymptotic window, nonfinite data, or
repeatability failure. Do not tune states, amplitudes, thresholds, residuals,
or Jacobians after execution.
