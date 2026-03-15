#!/usr/bin/env python3

import argparse
import sys

from scripts.symsolver_c3c4_fun import (
    build_symsolver_c3c4_cache,
    get_symsolver_c3c4_cache_path,
)

SUPPORTED_PATHWAYS = ("Type-I-C3-C4", "NADP-ME-C4")


def _parse_args():
    parser = argparse.ArgumentParser(
        description=(
            "Build and persist the slow symbolic solver cache so later runs "
            "can load the parsed solution directly."
        )
    )
    parser.add_argument(
        "--pathway",
        choices=(*SUPPORTED_PATHWAYS, "all"),
        default="all",
        help="Which pathway cache to build.",
    )
    parser.add_argument(
        "--overwrite",
        action="store_true",
        help="Overwrite an existing cache file.",
    )
    return parser.parse_args()


def _iter_pathways(pathway):
    if pathway == "all":
        return SUPPORTED_PATHWAYS
    return (pathway,)


def main():
    args = _parse_args()

    for pathway_option in _iter_pathways(args.pathway):
        target_path = get_symsolver_c3c4_cache_path(pathway_option)
        print(f"[build] {pathway_option}")
        print(f"[cache] {target_path}")

        try:
            cache_path = build_symsolver_c3c4_cache(
                pathway_option,
                overwrite=args.overwrite,
            )
        except FileExistsError:
            print("[skip] cache already exists; rerun with --overwrite to rebuild")
            continue
        except KeyboardInterrupt:
            print("\n[stop] interrupted while building cache", file=sys.stderr)
            return 130

        print(f"[done] {cache_path}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
