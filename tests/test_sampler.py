import unittest

from transect_sampler import Point, PointDistribution, SamplingRequest, Transect, sample


class SamplerTests(unittest.TestCase):
    def test_rotated_rectangle_includes_edges_and_excludes_outside_points(self):
        distribution = PointDistribution(
            id="d",
            crs="LOCAL:METRES",
            units="m",
            points=(
                Point("start", 0, 0),
                Point("edge", -0.5, 2),
                Point("end", 0, 4),
                Point("outside", 0.6, 2),
            ),
        )
        request = SamplingRequest(
            id="s",
            crs="LOCAL:METRES",
            units="m",
            transects=(Transect("north", 0, 0, 4, 1, 90),),
        )
        observations = sample(distribution, request)
        self.assertEqual([item.point_id for item in observations], ["start", "edge", "end"])
        self.assertAlmostEqual(observations[1].along_m, 2)
        self.assertAlmostEqual(observations[1].offset_m, 0.5)

    def test_overlapping_transects_keep_two_observations(self):
        distribution = PointDistribution(
            id="d", crs="LOCAL:METRES", units="m", points=(Point("p", 1, 0),)
        )
        request = SamplingRequest(
            id="s",
            crs="LOCAL:METRES",
            units="m",
            transects=(Transect("a", 0, 0, 2, 1, 0), Transect("b", 0, 0, 2, 1, 0)),
        )
        self.assertEqual([hit.transect_id for hit in sample(distribution, request)], ["a", "b"])

    def test_coordinate_system_mismatch_is_rejected(self):
        distribution = PointDistribution(id="d", crs="A", units="m", points=())
        request = SamplingRequest(
            id="s", crs="B", units="m", transects=(Transect("t", 0, 0, 1, 1, 0),)
        )
        with self.assertRaisesRegex(ValueError, "same CRS"):
            sample(distribution, request)

    def test_projected_utm_coordinates_use_metre_offsets(self):
        distribution = PointDistribution(
            id="d",
            crs="EPSG:2202",
            units="m",
            points=(Point("inside", 762505, 1302000.4), Point("outside", 762505, 1302000.6)),
        )
        request = SamplingRequest(
            id="s",
            crs="EPSG:2202",
            units="m",
            transects=(Transect("east", 762500, 1302000, 10, 1, 0),),
        )
        observations = sample(distribution, request)
        self.assertEqual([item.point_id for item in observations], ["inside"])
        self.assertAlmostEqual(observations[0].along_m, 5)
        self.assertAlmostEqual(observations[0].offset_m, 0.4)


if __name__ == "__main__":
    unittest.main()
