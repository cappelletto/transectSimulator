"""Pure sampling operations; no random placement, file I/O, or summary statistics."""

from math import cos, radians, sin, ulp

from .models import Observation, Point, SamplingRequest, SpatialDistribution, Transect


def observe_points(
    points: tuple[Point, ...], transect: Transect, sampling_id: str
) -> tuple[Observation, ...]:
    """Return point centers inside a closed, oriented rectangular transect.

    The centerline starts at (start_x_m, start_y_m). Angle is counterclockwise
    from +X. Both end caps and both long edges are included.
    """
    direction_x = cos(radians(transect.angle_deg))
    direction_y = sin(radians(transect.angle_deg))
    hits: list[Observation] = []
    for point in points:
        dx_m = point.x_m - transect.start_x_m
        dy_m = point.y_m - transect.start_y_m
        along_m = dx_m * direction_x + dy_m * direction_y
        offset_m = -dx_m * direction_y + dy_m * direction_x
        tolerance_m = 16 * ulp(max(1.0, abs(dx_m), abs(dy_m), transect.length_m, transect.width_m))
        if -tolerance_m <= along_m <= transect.length_m + tolerance_m and abs(offset_m) <= transect.width_m / 2 + tolerance_m:
            hits.append(
                Observation(
                    sampling_id=sampling_id,
                    transect_id=transect.id,
                    point_id=point.id,
                    x_m=point.x_m,
                    y_m=point.y_m,
                    along_m=along_m,
                    offset_m=offset_m,
                    attributes=point.attributes,
                )
            )
    return tuple(hits)


def sample(
    distribution: SpatialDistribution, request: SamplingRequest
) -> tuple[Observation, ...]:
    """Query each transect independently; overlapping transects retain both hits."""
    if distribution.crs != request.crs or distribution.units != request.units:
        raise ValueError("distribution and sampling request must use the same CRS and units")
    return tuple(
        observation
        for transect in request.transects
        for observation in distribution.observe(transect, request.id)
    )
