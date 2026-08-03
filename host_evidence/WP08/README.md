# WP08 dynamic lifecycle benchmark

WP08 is separated from WP09. It tests declared-state transaction semantics without changing any production residual or Jacobian term. The exact frozen source SHA-256 is `DF26D71744945ACD2693FCA95A9CB91713773D6989474C4439DEE196685A7675`.

## Stage H1: standalone frozen-UEL replay

`run/build_and_run_wp08_harness.cmd` compiles the exact frozen UEL with a direct call-history harness and executes:

- `RC-T1_ORDER`: intervening trial call followed by equal-declared-state replay;
- `RC-T2_REJECT`: rejected-trial/cutback-style replay with candidate state discarded;
- `RC-T3_TERMINAL`: inserted `KINC=0` call with committed state restored by the host harness;
- `RC-T4_SVARS_SERIALIZE`: binary round-trip of the declared SVARS payload;
- `RC-T5_REPEAT_PROCESS`: repeated equal-state evaluation and a second independent process;
- an explicit M0 SAVE-cache counterexample and its M1 state-input counterpart.

H1 compares `RHS`, `AMATRX` and physical `SVARS(1:14)` for all four Gauss points. Diagnostic slots are excluded from the physical-state gate. H1 passed with zero M1 differences; the M0 counterexample reproduced a nonzero history-dependent difference.

## Abaqus host gates

- H2: clean fixed-30 s path versus deterministic PNEWDT rejected attempts. Accepted times and CAP_REF align; maximum scaled current-state difference is `5.6624E-10`.
- H3: output-request versus no-output path is exactly equal at accepted GP states. The host did not expose a directly logged `KINC=0`; terminal-call evidence therefore remains the H1 inserted-call result.
- H4: uninterrupted versus restart continuation is exactly equal for accepted current-state GP fields, CAP_REF and nodal ODB `U/RF` on this patch.
- H5: P1 versus P14 physical accepted states pass the inherited engineering envelope (`2.30E-6`). P16 did not start because the runtime exposed only 14 CPUs.
- Diagnostic boundary: P14 `uentry.csv` and `finite.csv` contain concurrent header-write races. This is classified as diagnostic I/O failure, not physical-state failure. Parallel diagnostic CSVs are not production evidence.

## Result and boundary

`outputs/host/wp08_host_summary.json` records `overall_wp08=PASS_LIMITED` and `wp09_status=READY_FOR_SEPARATE_P1_DIRECTIONAL_FD_ONLY`.

This result validates the lifecycle transaction on the 2 x 2 patch only. It does not establish full-domain production readiness, P16 diagnostic safety, complete host terminal-call coverage, full thermodynamic consistency or a consistent Jacobian.

WP09 must run separately with P1 and the frozen source. It must not be combined with diagnostic-writer repair, physical-model changes, WP10 or a long-window solve.

## Reproduction

1. Run `run/build_and_run_wp08_harness.cmd`.
2. Run the host cases through `run/run_wp08_host_case.cmd` using the inputs in `inputs/`.
3. Extract nodal accepted-frame data with `run/extract_wp08_odb_nodal.py` under `abaqus python`.
4. Run `python run/audit_wp08_host_results_fixed.py`.
