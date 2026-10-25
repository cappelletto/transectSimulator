# Repository guidance and implementation assessment

## Current project phase

The historical Octave implementation has been assessed. A separate Python sampler proof of concept is now in scope at the user's request; a full migration and statistical analysis remain future work. See `docs/python_prototype.md` for its contract and limits.

The pre-refactor decision is to use GeoJSON for file based spatial interchange and QGIS visualization. The prototype CLI currently uses application-specific JSON. Keep metric Euclidean sampling in a verified projected CRS; implement GeoJSON reprojection at an adapter boundary. `docs/python_prototype.md` records the evidence, field naming convention, and remaining work.

When working in the current phase:

- Preserve the behavior, input data, and committed results unless a task explicitly requests a change.
- Prefer evidence-backed documentation, inventory, reproducibility checks, and characterization tests.
- Keep the Octave implementation and committed datasets intact while developing the isolated Python prototype. Do not silently change scientific behavior.
- Record ambiguities and suspected defects separately from confirmed behavior.
- No GitHub issue exists for this assessment. Do not create one unless explicitly requested.

## Scientific purpose

The repository contains a GNU Octave/MATLAB-style spatial simulation used to estimate how likely different transect sampling designs are to detect rare benthic colonies. It is intended to support a later spatial point-pattern analysis.

The original Spanish brief is `Explicacion_Fran_Simulaciones.pdf`. It describes:

- Six irregular reef survey windows following a 5-9 m depth band: Faro, Cote, Salina, Isla Larga, Yonqui, and Noronqui.
- Observed colony coordinates in UTM, datum SIRGAS-REGVEN.
- One hundred simulated spatial point patterns per reef, produced by a hybrid additive interpoint-area interaction model.
- Random placement of transect designs and counting colonies intersected by them.
- Three designs to compare:
  - five parallel 10 m x 1 m AGRRA transects, separated by 5 m;
  - four collinear 20 m x 1 m transects, separated by 10 m;
  - four collinear 50 m x 1 m transects, separated by 10 m.
- The scientific questions: detection probability for each method and which method is most accurate.

The PDF is a high-level requirements record, not proof that every described output was implemented.

## Repository map

- `process/transectGenerator.m`: core random placement, geometry, mortality filtering, point-in-polygon counting, plotting, and matrix assembly.
- `process/stripPolygon.m`: converts a specially ordered strip polygon into a sequence of quadrilateral index sets.
- `process/batchProcess.m`: loops over designs, sites, and mortality rates and saves ASCII matrices.
- `data/Shape/*_SHP.txt`: simplified two-column window boundaries used by the implementation.
- `data/Sim_*.csv`: three comma-separated columns: X, Y, simulation ID. Each file currently contains IDs 1 through 100.
- `data/points_*.csv`: two comma-separated columns containing observed colony X and Y coordinates.
- `data/original/`: original shapefile components and an associated CSV for one source window.
- `process/results/AGGRA1/`: early observed/AGRRA outputs and Octave state.
- `process/results/SinMortalidad/`: outputs without simulated mortality.
- `process/results/Tests01/`: batch outputs for six sites, three designs, and mortality rates 0, 10, 20, 30, and 40 percent.
- `process/results/Tests01/_automate.sh`: historical in-place output formatting script; it is not the simulation driver.

The historical Octave implementation has no dependency manifest or automated test suite. The Python prototype has `pyproject.toml` and a README. There is currently no continuous-integration configuration or declared license.

## Implemented data flow

`transectGenerator(shape_filename, colonies_filename, N_TRANSECTS, TEMPLATE_PARAMS, TRANSECT_TYPE, MORTALITY_RATE)` is the main entry point.

`TEMPLATE_PARAMS` is `[WIDTH, LENGTH, SKIP, N_SEGMENTS]`. `TRANSECT_TYPE` is `0` for the parallel AGRRA-style design and `1` for the collinear series design.

At a high level, the function:

1. Loads a closed, specially ordered strip boundary from a whitespace-delimited text file.
2. Splits it into adjacent quadrilateral regions with `stripPolygon`.
3. Approximates each region with local U/V vectors and lengths.
4. Loads all simulated colonies for a site.
5. Randomly marks a requested percentage of all rows as removed by setting their simulation ID to `-1`.
6. Generates `N_TRANSECTS` random normalized seed positions with `unifrnd`.
7. Constructs and transforms each component rectangle for the selected design.
8. Uses `inpolygon` to count colonies inside each component, grouped by simulation ID.
9. Returns a matrix with one row per surviving simulation ID: simulation ID in column 1 and one aggregate count per random design placement in the remaining columns.

The committed batch results normally have 100 rows and 101 columns: one simulation-ID column plus 100 random placements. This wide matrix is information-equivalent to part of the requested long-form result, but it is not the 10,000-row `simulation x transect placement x count` format illustrated in the PDF.

`batchProcess(TEST_PARAMS, NUMBER_TRANSECTS, SITE_NAMES, MORTALITY_RATES)` calls the generator for each parameter combination and writes a file named `<SITE>_<P|S><LENGTH>_M<MORTALITY>.csv` in the current working directory. Historical usage appears to have run it from `process/results/Tests01`, because its input paths are relative to that location.

## Geometry assumptions

The window files are not treated as arbitrary GIS polygons. Correct operation depends on all of these assumptions:

- The first coordinate is repeated at the end, forming a closed loop.
- Excluding the repeated endpoint, the boundary has an even number of vertices.
- Vertices describe the two sides of a narrow strip in the ordering expected by `stripPolygon`: for successive quadrilaterals, indices are formed as `[i, i+1, N-i, N-i+1]`.
- Vertex orientation and the selected corner convention align the local V direction with the intended east-west survey path.
- Each quadrilateral is close enough to a parallelogram for a two-vector affine approximation.

The code does not read the committed shapefile directly, inspect CRS metadata, validate topology, or use a general centerline/polygon clipping algorithm. The source brief identifies SIRGAS-REGVEN UTM coordinates, but the UTM zone and unit validation are not encoded in the program.

## Assessment of the current implementation

### What is present

- Parameterized transect width, length, separation, segment count, design type, placement count, and mortality rate.
- Segment-wise orientation across a bent survey window.
- Counts separated by spatial simulation ID and aggregated over components of one design placement.
- Batch generation over sites, designs, and mortality scenarios.
- A full committed `Tests01` result set for the configured historical experiment.
- Version tag `v1.0.0`; history shows the principal implementation was developed from March 18 through April 5, 2018.

### Confirmed limitations and behavior mismatches

- The code produces raw colony counts. It does not calculate detection probabilities, compare method accuracy, estimate uncertainty, or perform the downstream spatial point-pattern analysis.
- `data/points_*.csv` observed colonies are not accepted separately or used by the current main function. The current generator consumes the simulated-colony file only; observed-colony outputs in `results/AGGRA1` are historical artifacts.
- Randomness is not seeded or exposed through an API, so a run cannot be reproduced from its output alone.
- Mortality is sampled over all rows in a reef file, not independently at an exact rate within each simulation. A simulation ID vanishes if all of its colonies are removed. This is confirmed by committed mortality outputs with only 98 or 99 rows instead of 100.
- Series seed-area length is hard-coded as `4 * LENGTH + 3 * SKIP`, even though iteration uses the configurable `N_SEGMENTS`. Values other than four are therefore inconsistent for series designs.
- Seed width is sampled from `[0, 0.95]` as a documented quick fix. Longitudinal seed placement is restricted so the modeled design fits within the total strip length. This does not directly implement the source brief's statement that partially out-of-window transects may be acceptable and only intersecting colonies matter.
- Component rectangles are tested directly with `inpolygon`, but they are not clipped against the reef window. A colony inside a rectangle is counted regardless of whether a portion of that rectangle lies outside the window; because colonies are expected to lie within the window this may approximate the desired count, but sampled area is not calculated.
- Boundary points returned by `inpolygon` as `on` are not separately excluded despite the nearby comment; `in` includes points on the boundary in MATLAB/Octave semantics.
- Plotting and console output are embedded in the computational function, making headless batch execution and performance control difficult.
- File paths depend on the caller's working directory, output metadata are encoded only in filenames, and saved files have no headers.
- Input validation is absent for file shape, polygon ordering, parameter ranges, negative seed area, expected simulation IDs, and degenerate regions.
- Several comments are stale or inconsistent with the files, including a statement that simulated data contain an explicit row ID. Current simulation files contain only X, Y, and simulation ID.
- `unifrnd` introduces an Octave Statistics Package/MATLAB Statistics and Machine Learning Toolbox dependency that is not declared or loaded by the repository.

### Items requiring scientific clarification before migration

- Define “detection probability”: for example, probability of at least one colony per complete design placement, per component transect, or per unit sampled area.
- Define “most accurate” and its reference truth, loss function, uncertainty measure, and whether equal effort or equal sampled area is required across methods.
- Confirm whether transects may extend outside the survey window and how effective sampled area should be handled.
- Confirm whether colonies exactly on a transect boundary count as detected.
- Confirm whether mortality represents independent Bernoulli survival, removal of an exact per-simulation proportion, or a single reef-wide thinning process.
- Confirm whether mortality scenarios should be coupled across rates and methods, and whether random streams must be reproducible.
- Reconcile the six-reef statements with the PDF section that names only four reefs for one simulated-colony output.
- Confirm the authoritative site abbreviations. The repository uses `SALI`, `ISLA`, `YONQ`, and `NORO`, while the original brief also uses `SAL`, `IL`, `YNQ`, and `NOR`.
- Identify the authoritative CRS, including UTM zone, and provenance/licensing constraints for input and generated data.
- Decide whether committed result matrices are reference fixtures, exploratory outputs, or disposable generated artifacts.

## Historical context

Useful milestones from Git history:

- `e0ced98`: initial data, shapefiles, observed results, an early Octave algorithm, and an abandoned C/C++ experiment.
- `ca3f8e2`: added the `0.95` width sampling quick fix.
- `3630d9e`: implemented independently oriented series components.
- `32e85db`: first complete simulation-by-transect output matrix.
- `06db96f`: completed the simulated mortality input.
- `4f4a538`: added configurable segment count and batch processing; renamed several site files.
- `654d312`: committed the six-site x five-mortality-rate x three-design `Tests01` output set.
- `58ed4a5` / `v1.0.0`: reformatted numeric outputs to avoid scientific notation.

Use history to distinguish intentional legacy behavior from accidental regressions. Do not treat the `v1.0.0` tag as evidence of scientific validation; no automated validation is committed.

## Working and verification guidance

- Run Octave commands from `process/` only when their relative input/output paths have been checked. `batchProcess` may need to run from its historical result directory or be given a controlled working directory.
- Load the Statistics package before calling `transectGenerator` if the runtime requires it: `pkg load statistics`.
- Avoid overwriting committed result files during exploratory runs. Write to a temporary or explicitly named characterization directory.
- Fix the random seed when creating characterization fixtures and record the Octave/MATLAB version and package versions.
- Check both dimensions and semantics of outputs. In particular, assert the expected simulation-ID set rather than checking only row count.
- Treat plots as diagnostics, not verification. Verify coordinates, design geometry, membership counts, and edge cases numerically.
- MATLAB compatibility is unverified. The code uses Octave-specific syntax such as default function arguments and chained indexing (`size(x)(1)`), so do not describe it as MATLAB-compatible without testing or adaptation.
- The current execution environment may not have Octave installed. If unavailable, report that clearly; do not claim runtime verification based only on static inspection.

## Priorities for the later migration-preparation phase

When migration is explicitly authorized, preserve a trusted behavioral baseline before porting:

1. Resolve the scientific-definition questions above and define canonical input/output schemas.
2. Establish data provenance, CRS metadata, expected site identifiers, and license status.
3. Create deterministic characterization fixtures for geometry, counts, mortality, and boundary behavior using a known-compatible Octave version.
4. Separate computation, random-number generation, plotting, file I/O, and batch orchestration.
5. Decide which legacy behaviors must remain compatible and which confirmed defects should be corrected behind explicit versioned semantics.
6. Compare old and new implementations on synthetic edge cases and selected real datasets, using scientific invariants rather than byte-identical random outputs alone.
7. Add environment/dependency documentation, automated tests, and a reproducible analysis entry point before changing the implementation language or storage format.

Do not select a migration language, library stack, or new architecture solely from this assessment; that decision belongs to the next phase and should follow clarified scientific and operational requirements.
