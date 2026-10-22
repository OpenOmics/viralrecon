#!/usr/bin/env python3
"""Reconstruct the collector's layout from Nextflow-staged reporting files only."""
import argparse
import json
import shutil
import sys
from pathlib import Path, PurePosixPath


def stage(manifest, outdir):
    entries = json.loads(Path(manifest).read_text())
    if not isinstance(entries, list) or not entries:
        raise ValueError('Expected a nonempty list of report inputs')
    outdir = Path(outdir)
    if outdir.exists():
        raise ValueError(f'Report input directory already exists: {outdir}')
    validated, seen = [], set()
    for entry in entries:
        if not isinstance(entry, dict):
            raise ValueError('Each report input must be an object')
        destination = entry.get('destination')
        source = entry.get('source')
        if not isinstance(destination, str) or not destination or not isinstance(source, str) or not source:
            raise ValueError('Every report input needs a source and destination')
        relative = PurePosixPath(destination)
        if relative.is_absolute() or any(p in ('', '.', '..') for p in destination.split('/')) or '\\' in destination:
            raise ValueError(f'Unsafe report destination: {destination!r}')
        if destination in seen:
            raise ValueError(f'Duplicate report destination: {destination}')
        seen.add(destination)
        source = Path(source)
        if not source.is_file():
            raise ValueError(f'Missing report input: {source}')
        validated.append((source, outdir / relative))
    # Check all entries before writing, so missing inputs cannot look like valid zero counts.
    outdir.mkdir(parents=True)
    for source, destination in validated:
        destination.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(source, destination)
    return len(validated)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--manifest', required=True)
    parser.add_argument('--outdir', required=True)
    args = parser.parse_args()
    try:
        count = stage(args.manifest, args.outdir)
        print(f'Staged {count} report inputs in {args.outdir}')
    except (ValueError, OSError) as exc:
        print(f'ERROR: {exc}', file=sys.stderr)
        return 1
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
