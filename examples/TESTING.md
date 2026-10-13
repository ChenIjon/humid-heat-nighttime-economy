# Synthetic demonstration verification

This record concerns the public synthetic demonstration only. The study analyses were performed on a server, using the research environment stated in the main README. This test environment is not the study environment.

On 2026-10-13, the demonstration was checked on Windows using R 4.5.1, mgcv 1.9-3 and dplyr 1.1.4. No real transaction input or real aggregated panel was read. Compatibility with the server's R 4.4.1 environment has not been verified by this local test.

Commands, run from the repository root:

```bash
Rscript --vanilla examples/generate_demo_data.R
Rscript --vanilla examples/run_demo.R
```

The generator produced 360 rows across three city labels and 120 dates. TC, TS and AS each completed four successful bootstrap refits, with finite confidence intervals and converged final models. The demonstration's analysis step took approximately four seconds in this test; this is not a production benchmark. Generated inputs, fitted models, summaries and bootstrap caches were written only inside the demonstration directories, with generated outputs excluded by `.gitignore`.

The runtime emitted nonfatal locale startup warnings; they did not prevent the demonstration from completing. Each demonstration run records its actual package/runtime environment in `output_demo/demo_session_info.txt`.

All repository R scripts passed syntax parsing. Git-style ignore matching was checked with the installed Python pathspec library: eight representative permitted paths were retained and sixteen representative input/output paths were excluded, including nested Timing data, exported analysis input, fitted model objects and bootstrap caches. The nested Timing ignore rules add exclusions and contain no inclusion exceptions. Git metadata was not available in this local folder, so this check does not establish that an existing remote repository or its history contains no restricted files.

Running `Main_model.R` without the restricted panel was also checked as an expected failure: it stops before fitting and points to the README input schema and the public synthetic demonstration.
