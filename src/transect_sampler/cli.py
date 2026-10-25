"""Command line adapter for one independent sampling job."""

import argparse
from pathlib import Path

from .io import load_distribution, load_request, result_document, write_result
from .sampler import sample


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(prog="transect-sampler")
    commands = parser.add_subparsers(dest="command", required=True)
    sampling = commands.add_parser("sample", help="query one distribution with a set of transects")
    sampling.add_argument("--distribution", type=Path, required=True)
    sampling.add_argument("--transects", type=Path, required=True)
    sampling.add_argument("--output", type=Path, required=True)
    args = parser.parse_args(argv)

    try:
        distribution = load_distribution(args.distribution)
        request = load_request(args.transects)
        observations = sample(distribution, request)
        write_result(args.output, result_document(distribution, request, observations))
    except (OSError, ValueError, KeyError, TypeError) as error:
        parser.exit(2, f"transect-sampler: {error}\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
