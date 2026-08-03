# WP08 dynamic lifecycle verification report

Version: v1.0, 2026-07-22  
Decision: `PASS_LIMITED`  
Next package: `WP09 READY_FOR_SEPARATE_P1_DIRECTIONAL_FD_ONLY`

## 1. Scope and frozen configuration

WP08 tests whether the rollback-safe M1 state transaction remains history-independent under controlled trial-call, rejected-increment, output-request, restart and parallel-environment perturbations. It does not modify the physical residual, Jacobian, constitutive closure, mesh or production loading.

The tested UEL is the frozen `SRC-B1` source:

`tmp/paper_method_source/cp189B1/core/UEL_TWOPHASE_UPS_cp189B1_localPicardEngineeringParity_4185_4260_recovery4320.for`

SHA-256:

`DF26D71744945ACD2693FCA95A9CB91713773D6989474C4439DEE196685A7675`

The host benchmark uses a minimal 2 x 2 patch. The result therefore concerns lifecycle semantics on the benchmark, not full-domain production readiness.

## 2. Verification hierarchy

| Gate | Perturbation | Required comparison | Result |
|---|---|---|---|
| H1 | direct frozen-UEL call-history replay | `RHS`, `AMATRX`, physical `SVARS(1:14)` | `PASS` |
| H2 | clean increments versus deterministic rejected attempts | common accepted times, current state, CAP_REF, nodal `U` | `PASS` |
| H3a | output-request versus no-output deck | common accepted current state and CAP_REF | `PASS` |
| H3b | inserted terminal `KINC=0` call | equal-state replay after terminal call | `PASS_STANDALONE`; host call not observed |
| H4 | uninterrupted path versus restart continuation | accepted GP state, CAP_REF and nodal ODB state | `PASS` |
| H5 | P1 versus available parallel environment | accepted physical state inside inherited envelope | `PASS_P1_VS_P14` |
| H5 diagnostic I/O | P14 concurrent diagnostic writer | one header and parse-safe records | `FAIL_DIAGNOSTIC_ONLY` |

## 3. H1 standalone replay

The harness compiles and calls the exact frozen UEL. It executes five lifecycle tests and two control cases:

- `RC-T1_ORDER`: replay after an intervening trial call;
- `RC-T2_REJECT`: replay after a discarded candidate state;
- `RC-T3_TERMINAL`: replay after an inserted `KINC=0` call;
- `RC-T4_SVARS_SERIALIZE`: binary round trip of the declared state payload;
- `RC-T5_REPEAT_PROCESS`: repeated evaluation and an independent process;
- `RC-T2_LEGACY_COUNTEREXAMPLE`: intentionally unsafe SAVE-cache model;
- `RC-T2_SAFE_MICROKERNEL`: declared-state counterpart.

All M1 relative differences in `RHS`, `AMATRX` and physical `SVARS(1:14)` are exactly zero in the harness output. The unsafe M0 counterexample produces a nonzero `RHS` difference of `10`, while the safe microkernel returns zero. Two independent harness processes produce identical CSV evidence.

This gate establishes that the declared transaction can distinguish the known unsafe and corrected mechanisms. It is not, by itself, proof of Abaqus host cutback or restart behavior.

## 4. Abaqus host results

### 4.1 Clean and rejected-attempt replay

The P1 clean path accepts 19 increments. The cutback path also accepts the same 19 physical increments but contains six rejected `.sta` rows, including the controlled `60 -> 30 s` and `45 -> 30 s` attempts. Both jobs complete without a hard-failure signature or a non-finite audit value.

At common accepted states, the largest scaled current-state difference is `5.662407726E-10` in `p_w`, below the `1E-8` gate. CAP_REF differences are exactly zero. The maximum nodal ODB displacement difference is `4.649171993E-40 m` and the scaled displacement difference is `4.645445761E-38`. The ODB reaction-force difference is retained as a diagnostic (`1.376006964E-6` scaled) and is not used to upgrade the lifecycle claim.

### 4.2 Output non-interference and terminal-call boundary

The output-enabled and no-output P1 decks produce exact equality in the compared current-state and CAP_REF fields. Thus, the tested host output request does not perturb the accepted GP state on this patch.

The host logs do not expose a directly auditable `KINC=0` call. Terminal-call non-interference is therefore supported only by the explicit H1 inserted-call replay. This limitation prevents a claim of complete host terminal-path coverage.

### 4.3 Restart continuation

The uninterrupted P1 path and restart continuation have exact equality in all compared current-state fields, CAP_REF fields and nodal ODB `U/RF` values at common accepted times. The restart seed and continuation complete without hard-failure or finite-audit violations.

This is restart parity for the minimal patch and the frozen state payload. It does not establish full-domain restart production readiness.

### 4.4 Parallel environment and diagnostic separation

P16 did not enter analysis because the runtime exposed only 14 CPUs. H5 therefore compares P1 with P14. The largest scaled accepted-state difference is `5.373956729E-8` for `p_w`, below the inherited engineering envelope `2.30E-6`. The maximum nodal displacement difference is `9.313225746E-10 m`, corresponding to `9.305761355E-8` after scaling. CAP_REF differences are zero.

The P14 diagnostic files expose a distinct concurrency defect: `uentry.csv` contains two header lines and begins with a data row, while `finite.csv` contains four header lines. This is a thread-safety failure in diagnostic file initialization. It is not evidence of physical-state divergence, because the physical accepted-state and ODB gates pass independently. Conversely, the physical pass does not authorize P14 diagnostic CSVs as production evidence.

## 5. Decision

WP08 is classified as `PASS_LIMITED` because:

1. the frozen M1 transaction passes nonzero standalone replay;
2. deterministic rejected attempts preserve the accepted patch state within the declared gate;
3. tested output requests do not perturb accepted GP state;
4. restart continuation reproduces the uninterrupted accepted patch state exactly;
5. P1/P14 accepted physical states remain within the frozen engineering envelope;
6. complete host terminal-call coverage and parallel diagnostic I/O safety are not established.

The decision does not establish:

- full-domain rollback, restart or parallel safety;
- P16 behavior;
- thread-safe parallel diagnostic output;
- thermodynamic consistency;
- a consistent Newton Jacobian;
- conservation, time-step convergence or production readiness.

## 6. WP09 handoff boundary

WP09 is unlocked only as a separate single-core directional finite-difference package. It shall quantify the six flow-related blocks traced by WP06:

`KPU`, `KPP`, `KPS`, `KSU`, `KSP`, `KSS`.

WP09 must use P1, the same frozen source and the declared-state benchmark. It must not be combined with parallel-writer repair, physical-source changes, WP10 implementation, full-domain output work or a long-window solve. A WP09 result may quantify Jacobian error; it may not retroactively upgrade WP08 from patch-level lifecycle evidence to production readiness.

## 7. Evidence and reproduction

Primary machine-readable evidence:

- `outputs/wp08_harness_run1.csv`
- `outputs/wp08_harness_run2.csv`
- `outputs/wp08_harness_summary.json`
- `outputs/host/wp08_host_case_audit.csv`
- `outputs/host/wp08_host_comparisons.csv`
- `outputs/host/wp08_host_summary.json`

Canonical audit command:

```powershell
python run/audit_wp08_host_results_fixed.py
```

Expected terminal decision:

```text
WP08_OVERALL=PASS_LIMITED
WP09_STATUS=READY_FOR_SEPARATE_P1_DIRECTIONAL_FD_ONLY
PARALLEL_DIAGNOSTIC_IO=FAIL_HEADER_RACE_CLASSIFIED_DIAGNOSTIC_ONLY
```

No Git staging or commit is part of WP08 completion.
