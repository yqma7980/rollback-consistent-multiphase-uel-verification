# Rollback-consistent multiphase UEL verification

This repository is the public reproducibility artifact for the manuscript:

"State-conditioned residual-Jacobian verification for a nonlinear multiphase
poromechanics user element"

The archived Zenodo record retains its original artifact title,
"Rollback-consistent multiphase UEL verification artifact."

The manuscript title was revised during development; the archived scientific
evidence and DOI remain unchanged.

## Persistent identifiers

- GitHub: https://github.com/yqma7980/rollback-consistent-multiphase-uel-verification
- Archived release: https://doi.org/10.5281/zenodo.21773371
- DOI: 10.5281/zenodo.21773371

## Rebuild the figures

Python 3.12 with NumPy 2.4.3 and Matplotlib 3.10.8 is the frozen figure
environment. Run:

    python -m pip install -r requirements-figures.txt
    python scripts/smoke_test_v082.py

The smoke test verifies MANIFEST_SHA256.csv, rebuilds six main and twelve
supplementary figures from the distributed CSV files, and compares all
PDF/SVG/PNG outputs with EXPECTED_FIGURE_HASHES.csv.

## Contents

- src/: byte-preserved I5 UEL source and integrated standalone harness.
- inputs/: frozen states, directions, branches, and benchmark manifests.
- evidence/: all CSV inputs consumed by the figure builders.
- scripts/: package-relative v0.8.2 figure builders and smoke test.
- generated/: expected figures and machine-readable layout audits.
- preregistration/: principal frozen verification protocols.
- host_evidence/WP08/: available repaired-build inputs, logs, parsed outputs,
  source material, and evidence-boundary notes. Local launchers are excluded.
- manuscript/: v0.8.2 manuscript and supplementary-information sources/PDFs.
- MANIFEST_SHA256.csv: byte count and SHA-256 for every distributed object
  except the manifest itself.

## Evidence boundary

The package is self-contained for figure reconstruction. Numerical regeneration
is not a one-command cross-platform workflow because it requires licensed
Abaqus 2026 and Intel Fortran. The available host archive does not contain the
complete unsafe-versus-repaired Abaqus replay. The strict mixed finite-
difference record remains 36/39; the three lower-clip series are supplemented,
not overwritten, by active-set and higher-precision evidence. Newton order, a
common basin, full-domain accepted-state conservation, mesh/time convergence,
engineering qualification, and production readiness remain unestablished.

The UEL source contains dormant generic C:\abq_runs diagnostic fallback
strings. They are retained solely to preserve the verified source SHA-256.
Registered diagnostic routes must remain unset for the reported physical runs.

## Authors

Weiji Sun, Yangqi Ma, Bing Liang, Shi He, and Bo Liang.

Affiliations:

1. School of Mechanics and Engineering, Liaoning Technical University,
   Fuxin, Liaoning 123000, China.
2. China Coal Research Institute, Beijing 100013, China.

## License

Code in src/ and scripts/ is released under the MIT License. Data,
documentation, figures, and manuscript-supporting material are released under
Creative Commons Attribution 4.0 International. Third-party trademarks and
licensed Abaqus software are not included or licensed by this repository.
