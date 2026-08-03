# WP10-I5b Pre-registration: FINAL_CLIP_LOW Active-set Verification

Date: 2026-07-23

## 1. Purpose and prior knowledge

This study is a supplemental, non-blind verification of the only unresolved registered branch in WP10-I5. The original I5 evidence is immutable. In particular, the strict mixed full-matrix result remains 36/39, and the following three series remain `INCONCLUSIVE` before I5b:

| State/direction | Best relative L2 error | Minimum cosine |
|---|---:|---:|
| `FINAL_CLIP_LOW/M1` | 5.02854987766851189E-05 | 5.05571213561891275E-01 |
| `FINAL_CLIP_LOW/M2` | 1.52192591442967598E-05 | 9.03298678855495640E-01 |
| `FINAL_CLIP_LOW/M3` | 1.76650333514179487E-05 | 9.99999999562147579E-01 |

I5b does not reselect a favorable finite-difference step. Its primary evidence is an independently implemented fixed-active-set directional derivative of the residual. High-precision finite differences are auxiliary only.

## 2. Frozen evidence fingerprints

| Evidence | SHA-256 |
|---|---|
| I5 UEL | `DCD660846A5BB6498181B94F6F24D5B2D5294FD33D6B0D6C0914983E2B2A85CA` |
| I5 integrated harness | `4C5C11B486A4AEECDB1E8E22FEAB097CD19EDC1228AEA990994F2A20AC8BFA06` |
| I5 final summary | `C94343334440EF0EFAD5BE6E88AB4FBADEFFCEE542DC58714E04A1A3534E017B` |
| I5 mixed directions | `35038E483D7A05DF9FC49F696F732800B8D399C97CE808FCD7478201A7029579` |
| I5 state manifest | `3610CE0DCE34B4B894844B1B17F3A97D7F7D893B77112CCE21059DFD32F42E37` |
| I5 branch manifest | `F23B787088306C6BA765822A38CAA7F77C353210E04F8C038815B7F74DFC6881` |
| I5 baseline snapshot | `7807D6699C11D7ACA0DC9D68496EAE1BDC2394718D33E9026A0A4A0F472DD411` |
| I5 full-matrix run 1 | `D37372B2B49E1DEC2F1C697D6F4D35EC49BBFA62D35ABD45CB8BF5EA80DFBB01` |
| I5 full-matrix run 2 | `4FD532D1D986C54F75B8D09D304E8B655D9A19183CD55C625CA8204F0AB4F213` |

Any mismatch is a hard stop. The UEL, integrated harness, and all I5 evidence files are read-only.

## 3. Execution environment

- Processor policy: P1 standalone only (`OMP_NUM_THREADS=1`, `ABAQUS_CPUS=1`).
- Compiler: Intel Fortran Compiler Classic 2021.13.0, build 20240602_000000.
- Toolchain: Visual Studio 2017 Developer Command Prompt 15.9.78 and Intel oneAPI x64 environment.
- Registered compile profile: `/Od /traceback /check:bounds /extend-source:132 /fpp`.
- No Abaqus Standard engineering solve is authorized.
- No physical or trial-state `SAVE`, `DATA`, module cache, or file cache may be introduced.

## 4. Frozen state, directions, and scales

The state is copied exactly from `wp10_i5_states_run1.csv` for `FINAL_CLIP_LOW`. Its node-major trial saturation values are all `-0.1`, its committed Gauss-point saturation values are all `0.2`, `TIME=(600,1260) s`, `DTIME=60 s`, `KSTEP=3`, and `KINC=10`.

The actual 16-component `M1`, `M2`, and `M3` vectors are copied from `wp10_i5_mixed_directions_run1.csv`; no component is regenerated. The variable scales remain:

- `Xu=1E-5 m`
- `Xp=1E7 Pa`
- `Xs=1E-1`

The complete central-difference list remains:

`1E-2, 3E-3, 1E-3, 3E-4, 1E-4, 3E-5, 1E-5, 3E-6, 1E-6, 3E-7, 1E-7`.

The global relationship is frozen as

```text
D(RHS)[d] ~= -J_phys d
J_phys = J_returned - 1E-12 I16
```

No direction-specific sign is permitted.

## 5. Active-set object

The registered baseline is provisionally identified as a lower hard-clip active branch, not a switching equality. The distinction between the two saturation paths is frozen:

1. The interpolated raw nodal saturation is below zero.
2. `UPS_KINEM_SN` maps the residual closure input `SN_GP` to zero and sets its derivative with respect to the raw nodal saturation to zero.
3. `GSN` is still evaluated from nodal saturation gradients; it is not reset by `UPS_KINEM_SN`.
4. The storage/strain use path receives the clipped `SN_GP` and must follow the registered `SN_GP_RAW_USE` branch.
5. The closure, relative permeability, capillary pressure, and flux follow the separate clipped `SN_GP` path.
6. A lower storage-use clip does not authorize deleting explicit `u` or `p_w` dependencies, nor any branch-preserved flux dependency.
7. `KSTAB=1E-12 I16` belongs to `J_num` and is excluded from the physical oracle.

The branch ID, all relevant switching distances, and plus/minus branch preservation are measured before the oracle is classified. The active-set derivative is

```text
D_A R(x)[d] = lim(epsilon -> 0+) [R_A(x + epsilon d) - R_A(x)] / epsilon,
```

where membership in the baseline active set `A` is fixed. This is not a Clarke generalized derivative and does not characterize an exact switching point.

## 6. Independent residual mirror and oracle

The value mirror must independently replay the complete 16-row `RHS` value path for `FINAL_CLIP_LOW`. The analytic oracle propagates directional derivatives through the same fixed branches. It must not:

- call or inspect the UEL tangent assembly;
- read `AMATRX`, any named Jacobian block, `J_phys`, or a stored matrix-vector product;
- return a hard-coded `-J_phys*d`;
- omit residual terms merely because their response is small.

The UEL may be linked only in the driver to supply immutable reference `RHS`, `AMATRX`, candidate physical `SVARS`, and replay evidence. Matrix data are passed to the comparison layer, never to the oracle module.

The oracle records all 16 rows and the following term groups: `RU_effective_stress`, `RU_peq_force`, `RU_other`, `RP_storage`, `RP_strain`, `RP_flux`, `RP_paired_numerical`, `RS_storage`, `RS_strain`, `RS_mobility`, `RS_pc_first`, `RS_pc_second`, and `RS_other_paired_terms`.

## 7. Frozen acceptance gates

### 7.1 Baseline reproduction

The three original I5 series must reproduce exact keys, 11 step sizes, finite-difference values, `J_phys*d`, branch IDs, `INCONCLUSIVE` classifications, and run1/run2 equality.

### 7.2 Value parity

For the baseline and all 66 `M1`-`M3`, plus/minus, 11-step perturbation states:

- branch ID: exact;
- normalized maximum `RHS` difference: preferably bitwise exact and otherwise `<=1E-13`;
- normalized maximum term difference: `<=1E-13`;
- no oracle mutation of physical candidate `SVARS`.

Any violation is `ORACLE_VALUE_PARITY_FAILURE` and a hard stop.

### 7.3 Branch-distance gate

Distances are measured to the final lower/upper clip, ratio activation, beta minimum, beta-one clamp, effective-saturation clips, capillary-slope cap, permeability bounds, and every other active limiter/cap. The baseline must not lie on an equality. At least two adjacent registered step sizes must preserve every relevant branch, and each reported branch distance must exceed ten times the corresponding maximum perturbation used for the active-set conclusion. Otherwise the classification is `BRANCH_DISTANCE_INSUFFICIENT` or `NONSMOOTH_SWITCH_ONLY`.

### 7.4 Analytic oracle versus physical Jacobian

With the I5 state-specific row scaling, every `M1`-`M3` full-vector comparison must satisfy:

- normalized L2 error `<=1E-10`;
- normalized infinity error `<=1E-10`;
- cosine similarity `>=0.9999999999` for nonzero vectors;
- scaled absolute row error `<=1E-10` for every nonzero row.

For a near-zero row, the scaled oracle norm, scaled Jacobian norm, and scaled mismatch must each be `<=1E-10`; relative error and cosine are not pass/fail metrics.

### 7.5 Term closure, finite values, replay, and repeatability

- Per-equation term sums must equal the oracle total within a normalized `1E-13` closure gate.
- No `NaN`, `Inf`, or magnitude above `1E100` is permitted.
- Equal-state and rejected-trial ordering replay must be exact before and after each oracle phase.
- Two independent P1 processes must have exact keys, branches, values, derivatives, `J_phys*d`, and classifications; where not bitwise exact, normalized maximum difference must be `<=1E-13`.

## 8. High-precision auxiliary protocol

A compiler capability probe records `selected_real_kind`, `precision`, `range`, and `storage_size`. A type is eligible only when decimal precision is at least 30, exponent range at least 300, and storage size at least 128 bits.

If unavailable, the status is `UNAVAILABLE_ENVIRONMENT`; this does not block a passing analytic oracle. If available, all 11 high-precision central-difference steps are retained and compared with both the analytic oracle and `-J_phys*d`. The auxiliary gate is relative L2 error `<=1E-8`, cosine `>=0.99999999`, and at least two adjacent branch-preserving finite steps. Any conflict with a passing analytic oracle is a hard stop.

The binary64 cancellation diagnostic is descriptive only and cannot by itself close I5.

## 9. Permitted classifications

Each direction can only be classified as one of:

- `ACTIVE_SET_CONSISTENT_ANALYTIC`
- `FD64_RESOLUTION_LIMITED_ACTIVE_SET_ORACLE_CONSISTENT`
- `TRUE_BRANCH_JACOBIAN_MISMATCH`
- `NONSMOOTH_SWITCH_ONLY`
- `BRANCH_DISTANCE_INSUFFICIENT`
- `ORACLE_VALUE_PARITY_FAILURE`
- `HIGH_PRECISION_ORACLE_CONFLICT`
- `STATE_LIFECYCLE_FAILURE`
- `NONFINITE_FAILURE`
- `INCONCLUSIVE`

No engineering-acceptance or approximate-pass category is allowed.

## 10. Hard stops and completion rule

Stop without I6 if any input SHA differs, I5 cannot be exactly reproduced, value parity exceeds `1E-13`, the oracle reads matrix data, a UEL/residual change is required, branch stability is insufficient, any sign/DOF/term-closure ambiguity remains, replay or independent repeatability fails, a non-finite value occurs, a threshold/h/direction must be changed, or high-precision evidence conflicts with the oracle.

I5b passes only if all fingerprints and original evidence reproduce, `FINAL_CLIP_LOW` is a stable active branch, the value mirror passes, all three analytic directions and all 16 rows pass, term closure/replay/repeatability/finite gates pass, and no frozen source or I5 evidence changes. High-precision unavailability is allowed. I6 remains unauthorized until this complete gate is evaluated and explicitly recorded.

