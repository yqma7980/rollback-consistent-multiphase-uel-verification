# Available Abaqus host evidence and its boundary

This directory preserves the repaired-build WP08 patch inputs, launcher logs, parsed diagnostic/state outputs, and host audits used for the bounded clean/cutback/restart statements in Paper 1 v0.8.

It is not a complete unsafe-versus-repaired host replay. The unsafe M0 source and the original solver-native `.sta`, `.msg`, and `.dat` files were not retained in this canonical evidence tree. The distributed `*_analysis.log` files are launcher summaries, not substitutes for those raw solver files. Accordingly, the causal CAP_REF counterexample remains established by source tracing and the standalone call-order replay; the files here provide complementary repaired-build host behavior only.

The P14 CSVs are retained as negative evidence. Their known header race means they must not be treated as production diagnostic output.
