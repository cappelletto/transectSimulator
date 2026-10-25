# Transect simulator

This repository contains the historical Octave implementation in `process/` and a separate Python proof of concept in `src/transect_sampler/`. The Python package queries an already specified point distribution with one or more explicitly placed rectangular transects and returns individual observations. See [the prototype design](docs/python_prototype.md) for the pipeline, schemas, semantics, and follow-up work.

GeoJSON is the selected spatial interchange format for the planned QGIS companion workflow. The current CLI still uses its prototype JSON schema; the GeoJSON adapter and CRS transformation are follow-up work documented in the prototype design.

With Python 3.10 or newer, run the included example without installing dependencies:

```bash
PYTHONPATH=src python -m transect_sampler sample \
  --distribution examples/point_distribution.json \
  --transects examples/transects.json \
  --output /tmp/transect_observations.json
```

Or install the local package with `python -m pip install -e .` and use the `transect-sampler` command. The command performs one independent job, suitable for a future batch runner or web service. It does not place transects randomly or calculate summary statistics.
