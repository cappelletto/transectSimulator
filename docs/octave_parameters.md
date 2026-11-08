# Octave batch parameter reference

This describes the current, functional Octave implementation in `process/batchProcess.m` and `process/transectGenerator.m`. It is a naming and input reference for data generation; the scripts retain their existing signatures and behavior.

## The two counts that are easy to confuse

| Existing name | Clear meaning | Where used |
| --- | --- | --- |
| `NUMBER_TRANSECTS` | Number of independent, randomly placed **design placements** per batch case | Argument 2 of `batchProcess`; passed unchanged to the generator |
| `N_TRANSECTS` | The same number of **design placements** | Argument 3 of `transectGenerator`; controls seed draws, the outer loop, and the number of sampled-count columns |
| `N_SEGMENTS` | Number of rectangular **component transects within one placement** | Element 4 of `DESIGN_PARAMS`; controls the inner loop |

One placement consists of `N_SEGMENTS` rectangles. The count for that placement is the sum of colony hits across its rectangles, separately for each spatial simulation ID. Thus `NUMBER_TRANSECTS = 100` and `N_SEGMENTS = 5` produce **100 sampled-count columns**, each representing one five-rectangle design, not 500 columns. The output header uses `survey_1` through `survey_100` to reflect this distinction.

## Batch call and parameter mapping

```octave
batchProcess(DESIGN_CONFIGS, NUMBER_TRANSECTS, SITE_NAMES, MORTALITY_RATES, output_dir)
```

`DESIGN_CONFIGS` is a numeric matrix with **one design per row and exactly five columns**:

Invalid arguments raise `batchProcess:InvalidParameters` with the argument name and failed rule. Design-value errors also show the row, column, field name, and value, for example: `DESIGN_CONFIGS(2,3) [SPACING] = -1: must be finite and >= 0 metres`. Shape/type errors report the received class and dimensions; mortality-value errors identify the vector index.

| Column | Passed to generator as | Meaning and units |
| --- | --- | --- |
| 1 | `DESIGN_PARAMS(1)` → `WIDTH` | Width of each rectangle, metres |
| 2 | `DESIGN_PARAMS(2)` → `LENGTH` | Length of each rectangle, metres |
| 3 | `DESIGN_PARAMS(3)` → `SPACING` | Parallel design: spacing between component centerlines, metres. Series design: gap between consecutive rectangle ends, metres |
| 4 | `DESIGN_PARAMS(4)` → `N_SEGMENTS` | Number of component rectangles per placement |
| 5 | `TRANSECT_TYPE` | `0` = parallel AGGRA-style (`P` in filename); `1` = collinear series (`S`) |

For each row, site, and mortality rate, `batchProcess` calls:

```octave
transectGenerator(shape_file, simulation_file, NUMBER_TRANSECTS, ...
                  DESIGN_CONFIGS(row, 1:4), DESIGN_CONFIGS(row, 5), mortality_rate)
```

`NUMBER_TRANSECTS` is a single count shared by every batch case. `SITE_NAMES` is a character matrix with **one site code per row**; the current six data pairs use the equal-length codes `FARO`, `COTE`, `SALI`, `ISLA`, `YONQ`, and `NORO`. Avoid padded site names: the script uses each entire row in filenames. `MORTALITY_RATES` must be a **row vector** because the script reads its second dimension; each entry is a percentage, such as `[0 10 20 30 40]`.

The generator removes `ceil(number_of_colony_rows * rate / 100)` randomly selected rows from the entire site simulation file by changing their simulation IDs to `-1`. This is a site-wide thinning operation before counting, not an exact percentage removed separately from each simulation. Every simulation ID present in the input is retained in sorted order, including IDs with no surviving colonies. IDs absent from the input are not synthesized.

## File inputs and working directory

For site code `FARO`, the batch script loads `../data/Shape/FARO_SHP.txt` and `../data/Sim_FARO.csv` relative to the **current working directory**. Run from `process/`, or from a controlled directory where `../data/` resolves to the intended inputs. Inputs are resolved before creating the output directory. The optional fifth argument selects a new directory that must not already exist. With four arguments, the batch creates and prints a fresh `run_<UTC timestamp>_<unique suffix>/` directory. It returns `1` on success and raises an error on failure.

- The shape file is whitespace-delimited, with two numeric columns (`X Y`), no header, a repeated first coordinate at the end, and the special ordered strip boundary required by `stripPolygon.m`. It is not an arbitrary polygon format.
- The simulation file is comma-delimited, with three numeric columns (`X,Y,realization_id`) and no header. Its current IDs are in `1..100`; there is no separate colony row-ID column.
- `data/points_*.csv` contains observed point coordinates but is not read by `batchProcess`.
- Coordinates and transect dimensions are treated as projected metres. CRS metadata are not read from these text inputs; verify the CRS before combining datasets.
- The placement generator uses Octave's built-in `rand`; the Statistics package is not required for placement sampling.

The batch validates parameter ranges, site codes, required files, finite three-column colony data with integer IDs in `1..100`, and output filename uniqueness before starting. Polygon topology and CRS still require verification. Filenames omit width and spacing: conflicting configurations must be split into separate batches.

## Example: two sites, parallel 20 m design, 30% mortality

From `process/` in Octave:

```octave
DESIGN_CONFIGS = [1 20 5 5 0];       % width, length, spacing, components, parallel type
NUMBER_TRANSECTS = 100;           % 100 placements of the full five-component design
SITE_NAMES = ["FARO"; "COTE"];
MORTALITY_RATES = [30];           % row vector; one batch case per site
batchProcess(DESIGN_CONFIGS, NUMBER_TRANSECTS, SITE_NAMES, MORTALITY_RATES);
```

If the design fits, this writes `FARO_P20_N5_M30.csv` and `COTE_P20_N5_M30.csv` in the new run directory, with `FARO_metadata.yaml` and `COTE_metadata.yaml`. New exports have integer values and headers:

```text
realization_id,survey_1,survey_2,...,survey_100,total_colonies
```

Each output row corresponds to an original input simulation ID, including those with zero survivors or zero detections. The last column is that ID's total surviving colonies after the 30% site-wide removal. A zero-survivor row has zeros in all survey columns and `total_colonies`; a surviving population with zero detections retains its positive `total_colonies`. The older files committed under `process/results/` predate this CSV header and final-column format.

## Bounded placement generation

U is the cross-window fraction and V is normalized distance along the strip in both modes. Parallel proposals reserve lateral space for `WIDTH + (N_SEGMENTS-1)*SPACING` and longitudinal space for `LENGTH`. Series proposals reserve lateral space for `WIDTH` and longitudinal space for `N_SEGMENTS*LENGTH + (N_SEGMENTS-1)*SPACING`; component starts are mapped into their respective window regions. Lateral bounds use the perpendicular width of those regions to account for skewed U/V axes.

The complete component set is constructed from one seed before any colonies are counted. Every component must pass `transectCornersInside`: all four corners must lie inside or on the window boundary. This deliberately does not test rectangle edges through concavities. A rejected placement is redrawn up to the global constant `R = 15` retries after its initial attempt. Exhaustion prints a warning with the placement index, window filename, and design details, then uses the last complete candidate for counting even when corners are outside. Infeasible lateral bounds use U = 0.5. All requested component transects and surveys are retained. An overlong design fails immediately.

The five-component, 5 m-spacing example above spans 21 m laterally and may not fit the supplied reef windows. It may emit containment warnings and export placements with outside corners. Retry exhaustion does not prove no valid seed exists; it only means the fixed attempt budget was exhausted. The new bounds and retries change random placements compared with earlier runs.

Run focused geometry and generator checks from the repository root with:

```bash
octave --quiet --no-gui --eval 'addpath("process"); addpath("tests"); test_transect_placement();'
```
