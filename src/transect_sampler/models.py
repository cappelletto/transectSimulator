"""Small, immutable domain objects for the first point-distribution adapter."""

from dataclasses import dataclass, field
from math import isfinite
from typing import Any, Mapping, Protocol, Sequence


def _nonempty(value: str, label: str) -> None:
    if not isinstance(value, str) or not value.strip():
        raise ValueError(f"{label} must be a non-empty string")


def _finite(value: float, label: str) -> None:
    if isinstance(value, bool) or not isinstance(value, (int, float)) or not isfinite(value):
        raise ValueError(f"{label} must be a finite number")


@dataclass(frozen=True)
class Point:
    id: str
    x_m: float
    y_m: float
    attributes: Mapping[str, Any] = field(default_factory=dict)

    def __post_init__(self) -> None:
        _nonempty(self.id, "point id")
        _finite(self.x_m, "point x_m")
        _finite(self.y_m, "point y_m")


@dataclass(frozen=True)
class Transect:
    id: str
    start_x_m: float
    start_y_m: float
    length_m: float
    width_m: float
    angle_deg: float

    def __post_init__(self) -> None:
        _nonempty(self.id, "transect id")
        for name in ("start_x_m", "start_y_m", "length_m", "width_m", "angle_deg"):
            _finite(getattr(self, name), f"transect {name}")
        if self.length_m <= 0 or self.width_m <= 0:
            raise ValueError("transect length and width must be positive")


@dataclass(frozen=True)
class Observation:
    sampling_id: str
    transect_id: str
    point_id: str
    x_m: float
    y_m: float
    along_m: float
    offset_m: float
    attributes: Mapping[str, Any]


class SpatialDistribution(Protocol):
    """Adapter seam for empirical or generated spatial realizations."""

    id: str
    crs: str
    units: str

    def observe(self, transect: Transect, sampling_id: str) -> Sequence[Observation]: ...


@dataclass(frozen=True)
class PointDistribution:
    id: str
    crs: str
    units: str
    points: tuple[Point, ...]
    realization_id: str | None = None

    def __post_init__(self) -> None:
        _nonempty(self.id, "distribution id")
        _nonempty(self.crs, "CRS")
        if self.units != "m":
            raise ValueError("prototype coordinates and transect dimensions must use metres (units='m')")
        if self.realization_id is not None:
            _nonempty(self.realization_id, "realization id")
        ids = [point.id for point in self.points]
        if len(ids) != len(set(ids)):
            raise ValueError("point ids must be unique within a distribution")

    def observe(self, transect: Transect, sampling_id: str) -> tuple[Observation, ...]:
        from .sampler import observe_points

        return observe_points(self.points, transect, sampling_id)


@dataclass(frozen=True)
class SamplingRequest:
    id: str
    crs: str
    units: str
    transects: tuple[Transect, ...]

    def __post_init__(self) -> None:
        _nonempty(self.id, "sampling id")
        _nonempty(self.crs, "CRS")
        if self.units != "m":
            raise ValueError("prototype coordinates and transect dimensions must use metres (units='m')")
        if not self.transects:
            raise ValueError("sampling request requires at least one transect")
        ids = [transect.id for transect in self.transects]
        if len(ids) != len(set(ids)):
            raise ValueError("transect ids must be unique within a sampling request")
