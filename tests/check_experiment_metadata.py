"""Check an Octave test_batch_metadata fixture; requires PyYAML for this test only."""

import csv
import hashlib
from pathlib import Path
import sys

import yaml


fixture = Path(sys.argv[1]).resolve()
for site in ("A", "B"):
    manifest = yaml.safe_load((fixture / "success" / f"{site}_metadata.yaml").read_text())
    assert manifest["status"] == "complete"
    assert len(manifest["outputs"]) == 6
    assert "source_files" not in manifest
    assert "batch_rng_before" not in manifest
    assert "input_realization_ids" not in manifest
    for item in [manifest["window_input"], manifest["colonies_input"]]:
        assert hashlib.sha256(Path(item["path"]).read_bytes()).hexdigest() == item["sha256"]
    for case in manifest["outputs"]:
        assert case["status"] == "complete" and case["retry_limit"] == 15
        content = (fixture / "success" / case["filename"]).read_bytes()
        assert hashlib.sha256(content).hexdigest() == case["file"]["sha256"]
        rows = list(csv.reader(content.decode().splitlines()))
        assert rows[0] == ["realization_id", "survey_1", "survey_2", "total_colonies"]
        assert len(rows) - 1 == case["validation"]["rows"]
        assert "rng_before" not in case and "rng_after" not in case
        assert "realization_ids" not in case["validation"]
        assert "omitted_realization_ids" not in case["validation"]
        assert [int(row[0]) for row in rows[1:]] == [1, 3]
        if case["mortality_percent"] == 100:
            assert all(int(value) == 0 for row in rows[1:] for value in row[1:])
        assert case["validation"]["status"] == "passed"
        assert case["validation"]["columns"] == 4
        assert case["validation"]["surviving_colony_total"] == sum(int(row[-1]) for row in rows[1:])
    failed = yaml.safe_load((fixture / "failed" / f"{site}_metadata.yaml").read_text())
    assert failed["outputs"][0]["status"] == "complete"
    assert failed["status"] == ("failed" if site == "A" else "interrupted")
    assert failed["outputs"][1]["status"] == ("failed" if site == "A" else "pending")
    if site == "A":
        assert "Site A (SITE_NAMES row 1)" in failed["outputs"][1]["error"]["message"]

print("PASS: compact YAML, hashes, statuses, CSV headers, and validation totals")

segments = yaml.safe_load((fixture / "segments" / "ISLA_metadata.yaml").read_text())
assert [case["filename"] for case in segments["outputs"]] == [
    "ISLA_P10_N3_M22.csv", "ISLA_P10_N6_M22.csv"
]
assert all(case["status"] == "complete" for case in segments["outputs"])
print("PASS: component-count filenames and matching metadata")

warned = yaml.safe_load((fixture / "warning" / "A_metadata.yaml").read_text())
assert warned["status"] == "complete"
assert warned["outputs"][0]["status"] == "complete"
assert warned["outputs"][0]["validation"]["status"] == "passed"
print("PASS: containment warnings still produce completed datasets")
