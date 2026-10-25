"""Versioned JSON boundary for CLI input and observation output."""

import json
from dataclasses import asdict
from pathlib import Path
from typing import Any

from .models import Observation, Point, PointDistribution, SamplingRequest, Transect


def _document(path: Path) -> dict[str, Any]:
    data = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(data, dict) or data.get("schema_version") != "1":
        raise ValueError(f"{path}: expected a JSON object with schema_version='1'")
    return data


def load_distribution(path: Path) -> PointDistribution:
    data = _document(path)
    if data.get("kind") != "points":
        raise ValueError("prototype supports only kind='points'; other distribution adapters are future work")
    return PointDistribution(
        id=data["distribution_id"],
        realization_id=data.get("realization_id"),
        crs=data["crs"],
        units=data["units"],
        points=tuple(Point(**item) for item in data["points"]),
    )


def load_request(path: Path) -> SamplingRequest:
    data = _document(path)
    return SamplingRequest(
        id=data["sampling_id"],
        crs=data["crs"],
        units=data["units"],
        transects=tuple(Transect(**item) for item in data["transects"]),
    )


def result_document(
    distribution: PointDistribution,
    request: SamplingRequest,
    observations: tuple[Observation, ...],
) -> dict[str, Any]:
    return {
        "schema_version": "1",
        "distribution_id": distribution.id,
        "realization_id": distribution.realization_id,
        "sampling_id": request.id,
        "crs": distribution.crs,
        "units": distribution.units,
        "transect_ids": [transect.id for transect in request.transects],
        "observations": [asdict(item) for item in observations],
    }


def write_result(path: Path, data: dict[str, Any]) -> None:
    path.write_text(json.dumps(data, indent=2, allow_nan=False) + "\n", encoding="utf-8")
