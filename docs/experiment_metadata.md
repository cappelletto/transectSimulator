# Batch provenance and replication

Run from `process/` so `../data/` resolves correctly. For example, after verifying that the design fits the site:

```octave
rng(1234, "twister");
batchProcess([1 20 1 3 0], 50, ["FARO"; "COTE"], [0 30], "../runs/example_01");
```

The optional output directory must be new. Omitting it creates and prints a unique timestamped run directory under the working directory. Four-argument calls remain valid, but CSVs now go into that directory rather than directly into the working directory. CSV names use `<SITE>_<P|S><LENGTH>_N<COMPONENT_COUNT>_M<MORTALITY>.csv`, while `survey_` column names are retained. For example, `ISLA_P10_N6_M22.csv` means ISLA, parallel 10 m transects, six component transects per survey, and 22% mortality. `N` does not indicate the number of surveys. Width and spacing remain in the metadata; designs differing only in these fields still require separate batches. Different designs that map to the same filename are rejected before generation; split them into separate invocations.

Each site has one `<SITE>_metadata.yaml`, containing the whole set of design × mortality experiments for that invocation. The YAML writer uses quoted strings and JSON flow arrays, valid in YAML 1.2, and requires no external Octave package.

Shared metadata records UTC start/end times, batch designs and mortality rates, survey count, input paths and SHA-256 hashes, input colony row count, Octave version, and Git revision. Source-file inventories and input simulation-ID lists are omitted. CRS is not inferred from the input text files.

Each `outputs` entry records its CSV name, named design parameters, mortality percentage, timestamps, effective retry limit, validation results, and file checksum. The current effective retry limit is 15. Validation retains `status`, `rows`, `columns`, and `surviving_colony_total`. Simulation-ID lists, including the former `omitted_realization_ids`, are not exported.

The batch preserves design → site → mortality execution order and never resets the RNG between cases. Each manifest records the full ordered `site_names` list. RNG seeds and states are not exported. Metadata operations do not draw from the random stream.

## Checks and completion status

Containment retry exhaustion emits a console warning and continues with the last complete placement. Such cases may have outside corners and still be marked `complete` with validation `passed`: the output checks verify table integrity and population totals, not spatial containment. Warnings are not currently persisted in YAML. Invalid geometry and overlong designs still fail.

Site-specific errors name the site code and its `SITE_NAMES` row. Errors during an experiment also identify the design row, mortality percentage, and output filename; the same message is stored in the failed case's metadata. Original error identifiers are preserved.

Before publication, each matrix must have `survey_count + 2` columns with nonnegative finite integers and every input simulation ID exactly once in sorted order. Surviving totals must sum to `input_rows - ceil(input_rows * mortality / 100)` and cannot exceed the input count for an individual ID. Zero-survivor simulations remain as rows with zero survey counts and `total_colonies = 0`. A 100% mortality case contains all original input IDs with zero counts. Populations with survivors but no detected colonies also retain their rows. IDs absent from the original input are not synthesized.

The temporary CSV is read back and compared with the generated matrix and expected header before it is renamed to its final name. SHA-256 is then recorded. CSVs and YAML updates use temporary files in the same output directory followed by rename. A temporary file is not a completed dataset.

Output entries advance from `pending` to `running` to `complete` or `failed`. A site becomes `complete` only when all its entries complete. A generator/check/write failure records the active case's error, marks its site `failed` and unfinished other sites `interrupted`, preserves completed entries, and stops the batch. Check all site statuses before using a run. Preflight input/name failures create no run directory. An abrupt process kill may leave `running` status, which must also be treated as incomplete; automatic restart/resume is not implemented.

## Reproduction limits

The compact sidecar documents inputs, parameters, checksums, and results, but cannot by itself reproduce the random draws. For exact batch replay, retain the initial RNG state separately and restore it with matching inputs, source code, Octave version, and batch order. The example's explicit `rng(1234, "twister")` call is caller-controlled and is not recorded by the exporter.

Metadata contains absolute input paths for traceability and relative CSV filenames for relocation. After moving a run or inputs, resolve the new paths and verify hashes rather than assuming the original absolute paths still exist. Exact replay across different runtimes is not promised.

## Verification

From the repository root:

```bash
octave --quiet --no-gui --eval 'addpath(fullfile(pwd(),"process")); addpath(fullfile(pwd(),"tests")); test_transect_placement(); test_batch_metadata();'
```

The batch test prints `FIXTURE=/tmp/...`. To independently parse the actual YAML, check checksums, completion/failure records, retained validation totals, and absence of the removed fields, use:

```bash
python3 tests/check_experiment_metadata.py /tmp/fixture_path
```

The Python verification script requires PyYAML; this is a test dependency only. The Octave batch exporter has no YAML package dependency.
