#!/usr/bin/env python3
"""
proffast_bin_to_hdf5.py
-----------------------
Convert PROFFAST/OPUS solar spectrum .BIN files to HDF5.

Each .BIN file contains a Windows CRLF text header (dollar-sign-delimited
sections) followed by a block of little-endian float32 spectral values.
The binary block has 5 padding floats before the ngrid spectral values.

Usage
-----
# Convert one or more files, one HDF5 per BIN:
    python proffast_bin_to_hdf5.py 260406_123334SM.BIN 260406_123334SN.BIN

# Merge all BIN files from the same scan into a single HDF5:
    python proffast_bin_to_hdf5.py 260406_123334SM.BIN 260406_123334SN.BIN --merge

# Specify output directory:
    python proffast_bin_to_hdf5.py *.BIN --outdir /path/to/hdf5/

# Glob a directory:
    python proffast_bin_to_hdf5.py /data/tccon/20260406/ --outdir ./hdf5/

Output HDF5 layout (single-file mode)
--------------------------------------
/
├── attrs:
│   ├── location         str      e.g. "GSFC"
│   ├── date             str      e.g. "2026-04-06"
│   ├── date_raw         str      original 6-digit string "260406"
│   ├── time_ut_h        float64  effective UT time in decimal hours
│   ├── solar_elevation_deg   float64
│   ├── azimuth_deg      float64
│   ├── duration_s       float64
│   ├── latitude_deg     float64  (+ North)
│   ├── longitude_deg    float64  (+ East)
│   ├── altitude_km      float64
│   ├── opus_name        str
│   ├── filter           str
│   ├── opd_max_cm       float64
│   ├── semi_fov_rad     float64
│   ├── ils_type         int      1=simple, 2=extended
│   ├── modulation_efficiency   float64
│   ├── modulation_efficiency_sigma  float64
│   ├── comments         str
│   ├── source_file      str      original filename
│   └── created_by       str      this script + version
│
└── spectrum/
    ├── wavenumber       float64[ngrid]  cm⁻¹
    └── intensity        float32[ngrid]  normalized transmittance units

Output HDF5 layout (merged mode)
----------------------------------
/
├── attrs: (shared scan metadata, same fields as above minus filter/instrument)
└── filter_10/
│   ├── attrs: filter-specific instrument params
│   ├── spectrum/wavenumber   float64[ngrid]
│   └── spectrum/intensity    float32[ngrid]
└── filter_12/
    ├── attrs: filter-specific instrument params
    ├── spectrum/wavenumber   float64[ngrid]
    └── spectrum/intensity    float32[ngrid]
"""

from __future__ import annotations

import argparse
import re
import sys
from datetime import datetime, timezone
from pathlib import Path

import h5py
import numpy as np

__version__ = "1.0.0"
SCRIPT_NAME = f"proffast_bin_to_hdf5.py v{__version__}"

# ---------------------------------------------------------------------------
# Parsing
# ---------------------------------------------------------------------------

def _parse_date(date_raw: str) -> str:
    """Convert PROFFAST 6-digit date (DDMMYY or YYMMDD) to ISO 8601.

    PROFFAST uses DDMMYY: '260406' → 26 Apr 2006? or 06 Apr 2026?
    The OPUS spectrum name '260406SN245.0000' and the filename '260406_123334'
    suggest the date token is YYMMDD (year=26, month=04, day=06 → 2026-04-06).
    We assume YYMMDD with 2000+ century.
    """
    if len(date_raw) != 6 or not date_raw.isdigit():
        return date_raw  # return as-is if unexpected
    yy, mm, dd = date_raw[:2], date_raw[2:4], date_raw[4:]
    year = 2000 + int(yy)
    return f"{year}-{mm}-{dd}"


def parse_bin(filepath: Path) -> dict:
    """Parse a single PROFFAST .BIN file.

    Returns a dict with keys:
        header      : raw header text (str)
        wn          : wavenumber axis (float64 ndarray)
        spectrum    : spectral intensities (float32 ndarray)
        metadata    : dict of all parsed header fields
    """
    raw = filepath.read_bytes()

    # Locate all '$\r\n' section terminators
    dollar_positions = [m.start() for m in re.finditer(rb"\$\r\n", raw)]
    if not dollar_positions:
        raise ValueError(f"{filepath.name}: no '$\\r\\n' section markers found — not a PROFFAST BIN?")

    binary_start = dollar_positions[-1] + 3  # byte offset where float32 data begins
    header_text = raw[:binary_start].decode("latin-1")
    parts = header_text.split("$\r\n")
    # parts[0] = column labels (geo block)
    # parts[1] = geo values
    # parts[2] = instrument labels + values
    # parts[3] = ILS labels + values
    # parts[4] = spectral grid labels + values
    # parts[5] = comments + trailing

    def section_values(idx: int) -> list[str]:
        if idx >= len(parts):
            return []
        return parts[idx].strip().split()

    # --- Geolocation block ---
    geo = section_values(1)
    if len(geo) < 10:
        raise ValueError(f"{filepath.name}: geolocation section has fewer fields than expected.")

    date_raw = geo[1]
    date_iso = _parse_date(date_raw)

    metadata: dict = {
        "location":              geo[0].strip(),
        "date":                  date_iso,
        "date_raw":              date_raw,
        "time_ut_h":             float(geo[2]),
        "solar_elevation_deg":   float(geo[3]),
        "azimuth_deg":           float(geo[4]),
        "duration_s":            float(geo[5]),
        "latitude_deg":          float(geo[6]),
        "longitude_deg":         float(geo[7]),
        "altitude_km":           float(geo[8]),
        "opus_name":             geo[9],
    }

    # --- Instrument block ---
    inst = section_values(2)
    metadata["filter"]       = inst[0] if inst else ""
    metadata["opd_max_cm"]   = float(inst[1]) if len(inst) > 1 else float("nan")
    metadata["semi_fov_rad"] = float(inst[2]) if len(inst) > 2 else float("nan")

    # --- ILS block ---
    ils_text = parts[3].strip() if len(parts) > 3 else ""
    ils_lines = [l.strip() for l in ils_text.splitlines() if l.strip()]
    metadata["ils_type"] = int(ils_lines[0]) if ils_lines else 0
    if len(ils_lines) > 1:
        # "9.841E-01, 2.999E-03" → two floats
        eff_parts = [x.strip() for x in ils_lines[1].split(",")]
        metadata["modulation_efficiency"]       = float(eff_parts[0]) if eff_parts else float("nan")
        metadata["modulation_efficiency_sigma"] = float(eff_parts[1]) if len(eff_parts) > 1 else float("nan")
    else:
        metadata["modulation_efficiency"]       = float("nan")
        metadata["modulation_efficiency_sigma"] = float("nan")

    # --- Spectral grid block ---
    spec_vals = section_values(4)
    if len(spec_vals) < 4:
        raise ValueError(f"{filepath.name}: spectral grid section incomplete.")
    firstnue = float(spec_vals[0])
    lastnue  = float(spec_vals[1])
    deltanue = float(spec_vals[2])
    ngrid    = int(spec_vals[3])

    metadata["firstnue_cm-1"] = firstnue
    metadata["lastnue_cm-1"]  = lastnue
    metadata["deltanue_cm-1"] = deltanue
    metadata["ngrid"]         = ngrid

    # --- Comments block ---
    comments_text = ""
    if len(parts) > 5:
        comments_text = parts[5].strip()
    metadata["comments"] = comments_text

    # --- Binary spectral data ---
    binary_data = raw[binary_start:]
    n_floats = len(binary_data) // 4
    all_vals = np.frombuffer(binary_data[: n_floats * 4], dtype=np.float32)

    # There are (n_floats - ngrid) padding floats before the actual spectrum.
    padding = n_floats - ngrid
    if padding < 0:
        raise ValueError(
            f"{filepath.name}: binary block ({n_floats} floats) smaller than ngrid={ngrid}."
        )
    spectrum = all_vals[padding:].copy()

    # Wavenumber axis (use header-specified endpoints for accuracy)
    wn = np.linspace(firstnue, lastnue, ngrid, dtype=np.float64)

    return {
        "header":   header_text,
        "wn":       wn,
        "spectrum": spectrum,
        "metadata": metadata,
    }


# ---------------------------------------------------------------------------
# HDF5 writing helpers
# ---------------------------------------------------------------------------

def _write_metadata_attrs(h5obj, metadata: dict, source_file: str) -> None:
    """Write all metadata fields as HDF5 attributes."""
    for key, val in metadata.items():
        if isinstance(val, float) and np.isnan(val):
            h5obj.attrs[key] = "NaN"
        else:
            h5obj.attrs[key] = val
    h5obj.attrs["source_file"] = source_file
    h5obj.attrs["created_by"]  = SCRIPT_NAME
    h5obj.attrs["created_utc"] = datetime.now(timezone.utc).isoformat()


def _write_spectrum_group(parent, wn: np.ndarray, spectrum: np.ndarray) -> None:
    """Write wavenumber + intensity datasets inside a 'spectrum' group."""
    grp = parent.require_group("spectrum")

    ds_wn = grp.create_dataset(
        "wavenumber", data=wn, dtype=np.float64,
        compression="gzip", compression_opts=5,
    )
    ds_wn.attrs["units"]       = "cm-1"
    ds_wn.attrs["long_name"]   = "Wavenumber"
    ds_wn.attrs["description"] = "Spectral wavenumber axis (linearly spaced)"

    ds_sp = grp.create_dataset(
        "intensity", data=spectrum, dtype=np.float32,
        compression="gzip", compression_opts=5,
    )
    ds_sp.attrs["units"]       = "1"
    ds_sp.attrs["long_name"]   = "Spectral intensity"
    ds_sp.attrs["description"] = (
        "Normalised solar transmittance (PROFFAST preprocess5 output). "
        "Values near 1.0 indicate continuum; deep absorption features approach 0."
    )


# ---------------------------------------------------------------------------
# Single-file conversion
# ---------------------------------------------------------------------------

def convert_single(bin_path: Path, out_dir: Path) -> Path:
    """Convert one .BIN → one .h5 file."""
    parsed = parse_bin(bin_path)
    meta   = parsed["metadata"]

    out_path = out_dir / (bin_path.stem + ".h5")
    with h5py.File(out_path, "w") as hf:
        _write_metadata_attrs(hf, meta, bin_path.name)
        _write_spectrum_group(hf, parsed["wn"], parsed["spectrum"])

    return out_path


# ---------------------------------------------------------------------------
# Merged conversion (multiple BIN → one HDF5)
# ---------------------------------------------------------------------------

def convert_merged(bin_paths: list[Path], out_path: Path) -> Path:
    """Merge several .BIN files from the same scan into one .h5."""
    parsed_list = [parse_bin(p) for p in bin_paths]

    # Use geolocation metadata from the first file; instrument params per-filter
    shared_keys = {
        "location", "date", "date_raw", "time_ut_h",
        "solar_elevation_deg", "azimuth_deg", "duration_s",
        "latitude_deg", "longitude_deg", "altitude_km", "opus_name",
    }
    filter_keys = {
        "filter", "opd_max_cm", "semi_fov_rad",
        "ils_type", "modulation_efficiency", "modulation_efficiency_sigma",
        "firstnue_cm-1", "lastnue_cm-1", "deltanue_cm-1", "ngrid",
        "comments",
    }

    base_meta = {k: v for k, v in parsed_list[0]["metadata"].items() if k in shared_keys}
    source_files = ", ".join(p.name for p in bin_paths)

    with h5py.File(out_path, "w") as hf:
        _write_metadata_attrs(hf, base_meta, source_files)

        for parsed, bin_path in zip(parsed_list, bin_paths):
            filt = parsed["metadata"].get("filter", bin_path.stem)
            grp_name = f"filter_{filt}"
            grp = hf.require_group(grp_name)

            # Per-filter instrument metadata as group attributes
            for key in filter_keys:
                val = parsed["metadata"].get(key)
                if val is None:
                    continue
                if isinstance(val, float) and np.isnan(val):
                    grp.attrs[key] = "NaN"
                else:
                    grp.attrs[key] = val
            grp.attrs["source_file"] = bin_path.name

            _write_spectrum_group(grp, parsed["wn"], parsed["spectrum"])

    return out_path


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def _collect_bin_files(inputs: list[str]) -> list[Path]:
    """Expand directories and glob patterns to a list of .BIN paths."""
    paths: list[Path] = []
    for inp in inputs:
        p = Path(inp)
        if p.is_dir():
            paths.extend(sorted(p.glob("*.BIN")))
            paths.extend(sorted(p.glob("*.bin")))
        elif p.exists():
            paths.append(p)
        else:
            # Try glob expansion (shell may not have done it)
            expanded = sorted(Path(".").glob(inp))
            if expanded:
                paths.extend(expanded)
            else:
                print(f"Warning: '{inp}' not found — skipping.", file=sys.stderr)
    return paths


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description="Convert PROFFAST .BIN solar spectra to HDF5.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    parser.add_argument(
        "inputs", nargs="+",
        help="One or more .BIN files or a directory containing them.",
    )
    parser.add_argument(
        "--outdir", "-o", default=".",
        help="Output directory for .h5 files (default: current directory).",
    )
    parser.add_argument(
        "--merge", "-m", action="store_true",
        help=(
            "Merge all input files into a single HDF5. "
            "The output filename is derived from the first input file's stem."
        ),
    )
    parser.add_argument(
        "--merged-name",
        help="Explicit output filename when using --merge (overrides auto-naming).",
    )
    parser.add_argument(
        "--version", action="version", version=SCRIPT_NAME,
    )
    args = parser.parse_args(argv)

    out_dir = Path(args.outdir)
    out_dir.mkdir(parents=True, exist_ok=True)

    bin_files = _collect_bin_files(args.inputs)
    if not bin_files:
        print("Error: no .BIN files found.", file=sys.stderr)
        return 1

    print(f"Found {len(bin_files)} file(s):")
    for f in bin_files:
        print(f"  {f}")

    if args.merge:
        stem = Path(args.merged_name).stem if args.merged_name else bin_files[0].stem
        out_path = out_dir / (stem + "_merged.h5")
        print(f"\nMerging → {out_path}")
        try:
            convert_merged(bin_files, out_path)
            print(f"  OK  ({out_path.stat().st_size / 1024:.1f} KB)")
        except Exception as exc:
            print(f"  ERROR: {exc}", file=sys.stderr)
            return 1
    else:
        for bin_path in bin_files:
            out_path = out_dir / (bin_path.stem + ".h5")
            print(f"\n  {bin_path.name} → {out_path.name}")
            try:
                convert_single(bin_path, out_dir)
                print(f"  OK  ({out_path.stat().st_size / 1024:.1f} KB)")
            except Exception as exc:
                print(f"  ERROR: {exc}", file=sys.stderr)

    return 0


if __name__ == "__main__":
    sys.exit(main())
