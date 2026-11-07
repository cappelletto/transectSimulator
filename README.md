# Transect simulator

This repository contains the historical Octave implementation in `process/` and a separate Python proof of concept in `src/transect_sampler/`. The Python package queries an already specified point distribution with one or more explicitly placed rectangular transects and returns individual observations. See [the prototype design](docs/python_prototype.md) for the pipeline, schemas, semantics, and follow-up work.

For the functional Octave batch scripts, see the [parameter and input reference](docs/octave_parameters.md), including the distinction between placement count and component transect count.

The [pipeline diagram](docs/pipeline_prototype.svg) highlights the implemented sampler and the planned upstream and QGIS modules.

GeoJSON is the selected spatial interchange format for the planned QGIS companion workflow. The current CLI still uses its prototype JSON schema; the GeoJSON adapter and CRS transformation are follow-up work documented in the prototype design.

The planned visualization MVP will export colony points, transect footprints with counts, observed points, a summary, and a run manifest so one sampling result can be inspected in QGIS.

With Python 3.10 or newer, run the included example without installing dependencies:

```bash
PYTHONPATH=src python -m transect_sampler sample \
  --distribution examples/point_distribution.json \
  --transects examples/transects.json \
  --output /tmp/transect_observations.json
```

Or install the local package with `python -m pip install -e .` and use the `transect-sampler` command. The command performs one independent job, suitable for a future batch runner or web service. It does not place transects randomly or calculate summary statistics.
