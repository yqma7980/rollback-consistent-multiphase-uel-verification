# WP10-I4 Pre-registration: KPP Stabilization Separation and KSP Nullspace Policy

Date: 2026-07-23  
Task: `WP10_I4_KPPStabilizationSeparation_KSPNullspace_P1`  
Canonical workspace: `F:\Abaqus learning\code learning\papers_rollback_consistent_multiphase_poromechanics_20260722`

## 1. Scope and claim boundary

WP10-I4 is a P1 standalone verification stage. It does not alter the UEL physical residual, physical closure, physical state variables, rollback semantics, or any of the nine physical Jacobian blocks. Its sole purposes are:

1. to separate the returned matrix into a physical residual derivative and a matrix-only numerical guard;
2. to test `KPP=dRP/dp_w` after excluding the unpaired global diagonal guard;
3. to test whether the raw `KPP` discrepancy is quantitatively explained by that guard; and
4. to classify `KSP=dRS/dp_w` using a pre-registered absolute gate for the uniform-pressure gradient nullspace and the existing relative/cosine gate for non-null pressure directions.

Passing WP10-I4 does not establish a complete nine-block Jacobian, quadratic Newton convergence, thermodynamic consistency, full-domain validity, or production readiness. WP10-I5, WP10-I6, WP11, Abaqus Standard, full-domain output, 0.1 d, 30 d, and formal matrix cases are outside this task.

## 2. Frozen sources and fingerprints

| Level | Source | Required SHA-256 |
|---|---|---|
| Frozen M1 | `F:\Abaqus learning\code learning\tmp\paper_method_source\cp189B1\core\UEL_TWOPHASE_UPS_cp189B1_localPicardEngineeringParity_4185_4260_recovery4320.for` | `DF26D71744945ACD2693FCA95A9CB91713773D6989474C4439DEE196685A7675` |
| I1 | `paper1_methods\verification\WP10\I1\src\UEL_TWOPHASE_UPS_WP10I1_KabsDerivative_KPUKSU.for` | `C1FCA99B08507ED8EB43C20EE2C38C82C606F742874853DB6FFA4655848D660E` |
| I2 | `paper1_methods\verification\WP10\I2\src\UEL_TWOPHASE_UPS_WP10I2_fullScaleMechanics_KUUKUPKUS.for` | `EB89144BEA1C361015810D2276F4371F55BE33FE47641FFBEF0D7494F8EC60F9` |
| I3 baseline | `paper1_methods\verification\WP10\I3\src\UEL_TWOPHASE_UPS_WP10I3_completeKSSBranchwise.for` | `DCD660846A5BB6498181B94F6F24D5B2D5294FD33D6B0D6C0914983E2B2A85CA` |

The I4 UEL source is a single byte-for-byte copy of the I3 baseline named `UEL_TWOPHASE_UPS_WP10I4_KPPKSPPolicy.for`. Its initial and final SHA-256 must equal the I3 SHA. Any I3/I4 source difference is a hard stop.

## 3. Execution environment

- Architecture: Windows x64.
- Compiler: Intel Fortran Compiler Classic (`ifort`) initialized through Intel oneAPI.
- C/C++ environment: Visual Studio 2017 Professional Developer Command Prompt, version 15.9.78.
- Build flags: `/nologo /Od /traceback /check:bounds /extend-source:132 /fpp`.
- Parallelism: P1 only; `OMP_NUM_THREADS=1`, `ABAQUS_CPUS=1`.
- Abaqus Standard is not invoked.
- Two complete I4 executions must run as separate OS processes under the same environment.

## 4. Discrete map and global sign

The element has 16 degrees of freedom ordered by node as `[u_x,u_y,p_w,S_n]`.

| Field | Full indices |
|---|---|
| `u` | `1,2,5,6,9,10,13,14` |
| `p_w` | `3,7,11,15` |
| `S_n` | `4,8,12,16` |

The two target blocks are:

| Block | Residual rows | Variable columns | Meaning |
|---|---|---|---|
| `KPP` | `3,7,11,15` | `3,7,11,15` | `dRP/dp_w` |
| `KSP` | `4,8,12,16` | `3,7,11,15` | `dRS/dp_w` |

The single global sign is frozen as:

```text
D_h(RHS)[d] approximately -K d
```

No block-, state-, direction-, or step-specific sign choice is permitted.

## 5. Physical and numerical matrices

For the registered main profile:

```text
J_returned = AMATRX_returned
J_num      = KSTAB I_16
J_phys     = J_returned - J_num
KSTAB      = 1.0E-12
```

The I3 source adds `KSTAB` to every returned diagonal entry but adds no paired residual contribution. Consequently, `J_num` is excluded from the physical residual derivative. The primary consistency relation is:

```text
D_h(RHS)[d] approximately -J_phys d
```

For each KPP finite-difference record:

```text
e_raw  = D_h(RHS_P)[d_p] - (-KPP_raw d_p)
e_phys = D_h(RHS_P)[d_p] - (-KPP_phys d_p)
e_num  = KPP_num d_p
```

and the registered decomposition identity is:

```text
e_raw - e_num approximately e_phys
```

The sign of `e_num` is fixed by this definition and must not be selected after inspecting data.

The source scan must separately ledger: global `KSTAB`; paired pressure/saturation capacity; paired pressure/saturation anchor; paired pressure/saturation diagonal regularization; optional GP-valued pressure diagonal regularization; and any other pressure/saturation numerical controls. The main profile requires all paired or state-dependent optional terms to be off. If an additional unpaired matrix-only term remains and cannot be uniquely separated, I4 is `JPHYS_JNUM_AMBIGUOUS` and stops.

## 6. Main verification profile

The runtime environment is inherited from I3, with these controls fixed:

```text
DBG_FIELD_CAPACITY_REG=0
DBG_P_CAPACITY_REG=0
DBG_S_CAPACITY_REG=0
DBG_FIELD_ANCHOR_REG=0
DBG_P_ANCHOR_REG=0
DBG_S_ANCHOR_REG=0
DBG_REG_P_DIAG=0
DBG_REG_S_DIAG=0
DBG_PC_SLOPE_CAP_REG=0
DBG_CAPRAMP_SN_JUMP_FILTER=0
DBG_CAPRAMP_PC_FRONT_REG=0
DBG_CAPRAMP_FRONT_JUMP_REG=0
DBG_CAPRAMP_FRONT_TFILTER=0
DBG_CAPRAMP_FRONT_RU_REG=0
DBG_CAPRAMP_WELLB_FORCE_REG=0
DBG_CP179B2_TAN_RAMP=0
DBG_MOBILITY_SPLIT_MODE=0
DBG_MOBILITY_MIN_KINC=999999
DBG_FRONT_PNEWDT=0
```

`KSTAB` remains present in `J_returned` and is subtracted only in the standalone analysis layer. This subtraction is not a UEL source correction.

## 7. Registered states

The I3 state definitions are copied exactly:

`I0`, `I1`, `I2`, `LIMITER_OFF`, `LIMITER_INACT`, `RATIO_POS`, `RATIO_NEG`, `BETA_MIN_POS`, `BETA_MIN_NEG`, `FINAL_CLIP_LOW`, `FINAL_CLIP_HIGH`, `SE_CLIP_LOW`, and `SE_CLIP_HIGH`.

The block-level classification requires all three internal states `I0/I1/I2`. The remaining ten states are branch-preserving regressions. Every state must be finite, replay exact, away from switching equality, branch-stable under plus/minus perturbations, and exactly equal to I3 in physical candidate state.

## 8. Registered pressure directions and scaling

The pressure characteristic scale is:

```text
X_p = 1.0E7 Pa
```

The actual I3 direction components are copied and saved, not regenerated:

- `D1`: normalized uniform pressure direction `[0.5,0.5,0.5,0.5]`;
- `D2`: normalized alternating pressure direction `[0.5,-0.5,0.5,-0.5]`;
- `D3`: the saved dense direction produced with seed `20260723`.

All directions are normalized in scaled coordinates; the physical direction is `X_p d`.

## 9. Perturbation levels and state isolation

The complete central-difference sequence is frozen as:

```text
1E-2, 3E-3, 1E-3, 3E-4, 1E-4, 3E-5,
1E-5, 3E-6, 1E-6, 3E-7, 1E-7
```

Each plus/minus trial starts from the same committed `SVARS` copy. `TIME`, `DTIME`, `PROPS`, `JPROPS`, `COORDS`, and old state remain fixed. Only trial `p_w` is perturbed. Candidate `SVARS` are discarded. Trial state may not propagate between signs, `h`, directions, states, or blocks. All 11 levels are retained, including roundoff-dominated and branch-crossed rows.

## 10. Error definitions

For a non-null comparison between `f=D_h(RHS)[d]` and `g=-K d`:

```text
absolute_error_l2 = ||f-g||_2
relative_error_l2 = ||f-g||_2 / max(||f||_2,||g||_2,1E-30)
absolute_error_inf = ||f-g||_inf
relative_error_inf = ||f-g||_inf / max(||f||_inf,||g||_inf,1E-30)
cosine = f dot g / (||f||_2 ||g||_2)
```

For KSP nullspace scaling, the pressure direction is already physical (`X_p d`). The RS residual characteristic scale for each registered state is frozen before classification as:

```text
R_S,scale(state) = max(1, ||-KSP_phys (X_p D2)||_2, ||-KSP_phys (X_p D3)||_2)
```

The dimensionless nullspace values are:

```text
scaled_norm_fd    = ||D_h RHS_S[D1]||_2 / R_S,scale
scaled_norm_kd    = ||-KSP_phys (X_p D1)||_2 / R_S,scale
scaled_error_l2   = ||D_h RHS_S[D1] + KSP_phys (X_p D1)||_2 / R_S,scale
scaled_error_inf  = ||D_h RHS_S[D1] + KSP_phys (X_p D1)||_inf / R_S,scale
```

The direction-gradient norm is computed from the Q4 gradient operator and saved for D1/D2/D3. It is evidence, not a post-hoc scaling parameter.

## 11. KPP gates and classifications

`KPP_phys` is `CONSISTENT_NUMERICALLY` only if, for I0/I1/I2 and D1/D2/D3:

- relative L2 error is at most `1E-5`;
- cosine similarity is at least `0.99999`;
- at least two adjacent `h` values form a stable low-error region;
- the pre-roundoff region shows a central-difference decrease;
- there is no branch crossing or non-finite value; and
- run1/run2 satisfy the repeatability gate.

`KPP_raw` is expected, but not forced, to be `APPROXIMATE_BOUNDED_BY_JNUM`. The KSTAB explanation passes only if:

- `cosine(e_raw,e_num) >= 0.99999` for nonzero vectors;
- `e_raw-e_num` matches `e_phys` within the registered decomposition closure gate;
- I0/I1/I2, all three directions, and both runs support the same explanation; and
- no additional unexplained fixed error platform remains.

The dimensionless decomposition closure gate is `1E-10`, using the denominator `max(||e_raw||,||e_num||,||e_phys||,1E-30)`.

The predicted scale `KSTAB*X_p=1E-5` is contextual only. The result must use the actual matrix-vector product and saved direction.

## 12. KSP gates and classifications

For uniform-pressure `D1`, relative error and cosine are not primary gates. `D1` is `NULLSPACE_CONSISTENT` only when at least two adjacent `h` values simultaneously satisfy:

```text
scaled_norm_fd <= 1E-10
scaled_norm_kd <= 1E-10
scaled_error_l2 <= 1E-10
scaled_error_inf <= 1E-10
```

For non-null D2/D3, `CONSISTENT_NUMERICALLY` requires relative L2 error at most `1E-5`, cosine at least `0.99999`, at least two adjacent stable low-error `h` values, no branch crossing, and finite values.

The KSP block is `CONSISTENT_NUMERICALLY_NULLSPACE_AWARE` only if I0/I1/I2 D1 are nullspace-consistent, I0/I1/I2 D2/D3 are numerically consistent, the remaining branch states show no unexplained failure, and parity/replay/repeatability pass.

## 13. Parity, replay, finite, and repeatability gates

I3/I4 exact parity is required for RHS, closure values, `p_eq`, `K_abs`, physical candidate `SVARS`, branch IDs, raw returned AMATRX, and every raw block entry. `J_phys` and `J_num` are harness-derived and never written back into the UEL.

Replay is exercised before FD, after mixed rejected-trial call order, after KPP FD, and after KSP FD. RHS, raw AMATRX, `J_phys`, `J_num`, physical `SVARS`, closure values, and branch IDs must be exact.

All RHS, AMATRX, physical `SVARS`, derived matrices, vectors, and metrics must be finite and have absolute value no greater than `1E100`.

Two independent processes must be exact where bitwise comparison is meaningful. Otherwise, the normalized maximum difference must not exceed `1E-13` and the non-bitwise source must be explicitly reported.

The seven non-target blocks `KUU/KUP/KUS/KPU/KPS/KSU/KSS` must match I3 exact snapshots and FD records. Raw KPP/KSP AMATRX entries must also remain exact.

## 14. Completion and hard-stop rules

I4 is `DONE` only if source identity, unique `J_phys/J_num` separation, KPP stabilization explanation, physical KPP consistency, KSP null/non-null classifications, exact parity/replay/non-target regression, two-run repeatability, full 11-level curves, finite checks, and unambiguous sign/DOF/branch/lifecycle evidence all pass.

The task stops as `BLOCKED` on any source SHA mismatch; I3/I4 source difference; ambiguous numerical-matrix separation; residual/value/SVARS/raw-AMATRX drift; replay or regression failure; unexplained KPP physical platform; unexplained non-null KSP behavior; non-finite value; branch ambiguity; irreproducibility; need to modify residual/closure/source; need for Standard/full-domain execution; or need to change any pre-registered threshold, state, direction, or `h` range.

No automatic source correction, I5, I6, WP11, engineering run, or Git operation is authorized.

