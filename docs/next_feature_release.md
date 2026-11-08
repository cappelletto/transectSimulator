# Next feature release guide

The pre-refactor release is complete. This guide records the observed pain points and proposed work for the next cycle; it does not authorize changes to the completed release's sampling behavior. Assessment date: 2026-11-08; inspected baseline: `dc85549`. This is a source snapshot reference, not a new release tag.

## Completed baseline and constraints

- Octave generates datasets for parallel and series surveys, retains every realization ID present in the input, and writes integer CSVs with `realization_id,survey_1,...,survey_N,total_colonies`.
- CSV filenames include component count: `<SITE>_<P|S><LENGTH>_N<COMPONENT_COUNT>_M<MORTALITY>.csv`. Each batch creates a new run directory and one compact YAML manifest per site with input/output hashes, parameters, timestamps, and validation summaries.
- Placement generation attempts bounded seeds and four-corner checks, then warns and uses the final complete candidate after 15 retries. Boundary contact is allowed. Overlong designs and invalid geometry still fail. Four-corner checking does not establish full containment in concave windows.
- Python is a point-sampling proof of concept with explicit rectangular transects, metre-suffixed fields, and a CLI using application-specific JSON. It is not an Octave replacement yet.
- Preserve the existing datasets and label future behavior changes. Use file based storage and GeoJSON for spatial interchange; keep Euclidean computation in a verified projected metre CRS and geographic conversion at adapters.
- Keep default metadata compact. Source-file inventories, bulk RNG state, and realization-ID lists were deliberately removed. Any optional replay artifact needs an explicit interface decision rather than silently restoring those fields.

See [Octave parameters](octave_parameters.md), [experiment metadata](experiment_metadata.md), and the [Python prototype](python_prototype.md) for current usage and the proposed QGIS MVP.

## Prioritized pain points and tasks

Priority means recommended order for the next cycle, not a commitment to ship every task together. Items marked **decision** need agreement before their implementation.

| Priority / item | Evidence and consequence | Next action and completion criterion |
| --- | --- | --- |
| P0 — Persist placement quality | `sampleTransectPlacement` warns on containment exhaustion, but the batch can still report `complete` and validation `passed`. Console output is the only record of outside corners. | Return placement diagnostics separately from count data and persist compact warning/retry/outside-corner summaries per CSV. A run reopened later must distinguish table integrity from spatial quality. Keep the current warning-and-continue policy as the default unless explicitly changed. |
| P0 — Reproduction contract (**decision**) | Manifests omit RNG information by choice. The shared stream advances through design → site → mortality loops, including retries, so case order affects results. | Decide whether the caller retains RNG state externally or a separate opt-in replay artifact is offered; decide whether independent per-case streams are required for future parallel jobs. Document and test exact replay in a matched runtime without expanding the default YAML with bulk state. |
| P0 — Establish comparison fixtures | The Python sampler returns encounters; Octave exports aggregates and does not persist its sampled rectangle coordinates. Random output alone cannot isolate counting differences from placement differences. | Save a small synthetic fixture of explicit points, component rectangles, realization IDs, and expected counts. Both implementations must match on the agreed boundary/overlap rules, including zero detections and zero survivors. Keep generated experiment datasets out of test rewrites. |
| P1 — Define the public sampler API (**decision**) | Python's `sample(distribution, request)` returns observations, but window geometry, grouped surveys, placement diagnostics, and batch identity are not represented. The adapter protocol delegates observation logic to the distribution implementation. | Specify the snapshot, sampling-request, and result contracts before extending them. Define who owns geometry queries, how components belong to a survey, and how zero-observation requests remain identifiable. Version the file schemas and provide matching Python/CLI examples. |
| P1 — Deliver the QGIS MVP | Current Python JSON is not GeoJSON. Neither Python nor Octave exports the actual survey footprints as GIS layers. | Implement the MVP described in `python_prototype.md`: colony points, surveyed footprint polygons with counts, observation points, a descriptive summary, and a compact manifest. Confirm the three spatial layers load together in QGIS and reconcile with the raw counts, including overlaps and zero hits. |
| P1 — Validate CRS and import legacy files | Python compares CRS strings and `units='m'`; this does not prove projection or units. Octave text inputs have no CRS metadata or stable colony IDs. | Verify source CRS per dataset; `window_260.qpj` identifies EPSG:2202 for that shapefile only. Import CSVs with explicit schemas and stable entity/realization IDs; transform standard GeoJSON longitude/latitude into the chosen metre CRS before sampling. Reject ambiguous CRS instead of assuming UTM for every file. |
| P1 — Separate computation from presentation | Octave interleaves mortality, placement, counting, plotting, and printing. Python has a cleaner core, but its CLI writes directly to the requested output path. | Isolate pure geometry/counting and explicit random-state ownership. Add headless operation and controlled output publication. Demonstrate unchanged seeded Octave counts for the same policy and inputs, and useful Python CLI errors without partial published results. |
| P1 — Finish input/failure checks | Batch validation is clearer, but direct helper calls have fewer shape/type checks; polygon topology and general site/output argument types are not comprehensively checked. Metadata writing during error handling can itself fail and obscure the original failure. | Add focused checks for malformed/nonfinite window coordinates, degenerate regions, direct-call parameter shapes, and invalid site/output argument types. Preserve the original experiment error if recording it fails. Tests must identify the offending field/site and leave the run explicitly incomplete. |
| P2 — Remove naming ambiguity safely | `NUMBER_TRANSECTS`/`N_TRANSECTS` count whole surveys; `N_SEGMENTS` counts component transects; other segment variables refer to window regions. `suma`, `temp`, `uu`, and `vv` obscure meaning. | Use the vocabulary below in the new API. Rename internal variables in small changes; retain positional compatibility or introduce a documented wrapper for existing calls. Verify a seeded numerical comparison after each executable rename. |
| P2 — Improve experiment orchestration | Filenames omit width and spacing and can collide; batches stop at the first hard error and have no resume path. Metadata stores absolute input paths. | Decide a full configuration identity and a separate resume/continue-on-error policy before adding parallel execution. Keep readable filenames, input hashes, and relocatable file references. A resumed job must never silently mix incompatible configurations or overwrite completed artifacts. |
| P2 — Define statistical summaries (**decision**) | Raw survey encounters are summed across components; an overlapping colony may be counted more than once. Detection probability and “accuracy” are not defined statistical estimands. | Define encounter counts versus unique-colony counts, denominators, effective sampled area when rectangles escape the window, and effort normalization. Start with descriptive counts; add probabilities and uncertainty only after those definitions are agreed and tested. |
| Later — Distribution and temporal modules | Analytical realization generation, colony footprints/cover, mortality transitions, and fragmentation lineage are absent from the Python pipeline. | Specify generators and transitions as separate modules that output stable snapshots. Preserve entity and parent IDs, time, units, and provenance. Do this after the sampler/import/export interfaces stabilize. |

## Naming and interface decisions

Use these meanings consistently in new APIs; existing executable names remain compatibility aliases until a deliberate migration:

| Concept | Preferred term | Current naming to reconcile |
| --- | --- | --- |
| One spatial point-pattern snapshot | `realization_id` | Octave CSV now uses this; several source comments and validation messages still say simulation |
| One complete sampling placement | `survey_id`, `survey_count` | Python uses `sampling_id`; Octave count arguments are `NUMBER_TRANSECTS` and `N_TRANSECTS` |
| One rectangle in a survey | `transect_id`, `component_count` | Python `Transect.id`; Octave `N_SEGMENTS` |
| One quadrilateral of the survey window | `window_region` | Octave also calls these segments; avoid conflating them with component transects |
| Distance and direction | `_m`, `_deg` suffixes | Python already uses these; Octave's dimension variables have implicit metre units |
| Realization population after mortality | `total_colonies` | Keep the CSV field; includes zero and differs from the number detected |

The interface specification must also settle boundary inclusion, numerical tolerance, grouping and deduplication of encounters, empty realizations supplied explicitly, and diagnostic return values. Python currently includes edges using a floating-point tolerance; Octave uses `inpolygon`. Treat differences as scientific behavior changes requiring explicit fixtures, not incidental cleanup.

Do not lose the all-input-realizations rule: the point input alone cannot represent a realization with no initial entities. Supporting those requires an explicit realization inventory in the future input contract, not fabricated IDs in the current exports.

## Proposed next-cycle sequence and release gate

1. Freeze the compatibility fixtures and name the spatial-quality and reproduction policies.
2. Define versioned snapshot/request/result APIs and implement the CRS/import boundary.
3. Deliver the CLI-to-QGIS MVP with compact diagnostics and descriptive summaries.
4. Address remaining orchestration, performance, and scientific-summary work as separate features.

Before shipping the next feature release:

- Run the existing Octave argument, placement, and batch checks, plus independent YAML/CSV validation. Use temporary fixtures so production experiment files are untouched.
- Run Python tests for rotated rectangles, boundaries, overlaps, zero observations, invalid payloads, and CRS handling. Add GeoJSON round-trip and CLI publication tests as adapters arrive.
- Verify warning fallback still returns complete surveys and that spatial-quality flags survive export; do not interpret `validation.status='passed'` as proof of containment.
- Open a georeferenced fixture in QGIS and reconcile displayed footprints/hits with its summary.
- Document compatibility changes, supported runtime/dependencies, example commands, and historical headerless/`simulation_id` versus current `realization_id` outputs. Add CI for synthetic fixtures; decide dataset licensing/provenance before broader distribution.

The current repository has no CI or declared license. Performance improvements such as spatial indexes should follow a representative benchmark, after correctness and API stability; the current Python scan examines every point for every transect.
