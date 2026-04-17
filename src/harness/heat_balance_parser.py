"""Parse OpenFOAM postProcessing output for heat balance analysis."""

from __future__ import annotations

import re
from dataclasses import dataclass, field
from pathlib import Path


@dataclass
class HeatBalance:
    """Energy budget summary from OpenFOAM wallHeatFlux + volAverageT."""

    heater_input_W: float
    wall_loss_W: float
    vent_loss_W: float = 0.0
    vol_avg_T: float = 0.0  # Volume-averaged temperature [K]
    patch_fluxes: dict[str, float] = field(default_factory=dict)

    @property
    def imbalance_W(self) -> float:
        return self.heater_input_W + self.wall_loss_W + self.vent_loss_W

    @property
    def imbalance_pct(self) -> float:
        if abs(self.heater_input_W) < 1e-12:
            return 0.0
        return self.imbalance_W / abs(self.heater_input_W) * 100.0


_HEATER_PATCHES = {"heater_wall"}
_VENT_PATCHES = {"supply_vent", "exhaust_vent"}


def parse_wall_heat_flux(case_dir: Path) -> dict[str, list[tuple[float, float]]]:
    """Parse wallHeatFlux postProcessing output.

    Returns dict mapping patch name to list of (time, integrated_flux_W) tuples.
    Positive flux = heat into domain, negative = heat out.
    """
    pp_dir = case_dir / "postProcessing" / "wallHeatFlux"
    if not pp_dir.exists():
        return {}

    result: dict[str, list[tuple[float, float]]] = {}

    for time_dir in sorted(pp_dir.iterdir()):
        if not time_dir.is_dir():
            continue

        # OpenFOAM v2312 writes wallHeatFlux.dat (combined format)
        combined = time_dir / "wallHeatFlux.dat"
        if combined.exists():
            _parse_single_dat(combined, result)
            continue

        # Older format: surfaceFieldValue.dat
        sv_file = time_dir / "surfaceFieldValue.dat"
        if sv_file.exists():
            _parse_single_dat(sv_file, result)
            continue

        # Per-patch files (e.g. floor_wallHeatFlux.dat)
        _parse_per_patch_files(time_dir, result)

    return result


def _parse_per_patch_files(
    time_dir: Path, result: dict[str, list[tuple[float, float]]]
) -> None:
    """Parse per-patch wallHeatFlux files (OpenFOAM v2312 format).

    OpenFOAM v2312 writes files like ``floor_wallHeatFlux.dat``.
    We strip the ``_wallHeatFlux`` suffix to recover the patch name.
    """
    for f in sorted(time_dir.iterdir()):
        if not f.is_file() or not f.name.endswith(".dat"):
            continue
        patch_name = f.stem.removesuffix("_wallHeatFlux")
        entries = _read_dat_file(f)
        if entries:
            result.setdefault(patch_name, []).extend(entries)


def _parse_single_dat(
    dat_file: Path, result: dict[str, list[tuple[float, float]]]
) -> None:
    """Parse a combined wallHeatFlux.dat file.

    OpenFOAM v2312 writes rows like:
        time \\t patch \\t min \\t max \\t integral
    The integral column (index 4) is the total wall heat flux [W].
    """
    lines = dat_file.read_text(encoding="utf-8").splitlines()

    for line in lines:
        if line.startswith("#") or not line.strip():
            continue
        parts = line.split()
        if len(parts) < 2:
            continue
        try:
            time_val = float(parts[0])
        except ValueError:
            continue

        if len(parts) >= 5:
            # Format: time patch min max integral
            patch_name = parts[1]
            try:
                integral = float(parts[4])
            except (ValueError, IndexError):
                continue
            result.setdefault(patch_name, []).append((time_val, integral))
        elif len(parts) == 2:
            # Single-value format
            result.setdefault("total", []).append((time_val, float(parts[1])))


def _read_dat_file(path: Path) -> list[tuple[float, float]]:
    """Read a simple time-value .dat file (# comments, tab/space separated)."""
    entries: list[tuple[float, float]] = []
    for line in path.read_text(encoding="utf-8").splitlines():
        if line.startswith("#") or not line.strip():
            continue
        parts = line.split()
        if len(parts) >= 2:
            try:
                entries.append((float(parts[0]), float(parts[1])))
            except ValueError:
                continue
    return entries


def parse_vol_average_t(case_dir: Path) -> list[tuple[float, float]]:
    """Parse volAverageT postProcessing output.

    Returns list of (time, volume_averaged_T_in_K) tuples.
    """
    pp_dir = case_dir / "postProcessing" / "volAverageT"
    if not pp_dir.exists():
        return []

    entries: list[tuple[float, float]] = []
    for time_dir in sorted(pp_dir.iterdir()):
        if not time_dir.is_dir():
            continue
        dat_file = time_dir / "volFieldValue.dat"
        if dat_file.exists():
            entries.extend(_read_dat_file(dat_file))

    return entries


def compute_heat_balance(
    wall_fluxes: dict[str, list[tuple[float, float]]],
    vol_avg_t: list[tuple[float, float]] | None = None,
) -> HeatBalance:
    """Compute heat balance from parsed postProcessing data.

    Uses the last time step values from each patch.
    """
    patch_last: dict[str, float] = {}
    for patch, series in wall_fluxes.items():
        if series:
            patch_last[patch] = series[-1][1]

    _skip = _HEATER_PATCHES | _VENT_PATCHES | {"total"}
    heater_input = sum(v for k, v in patch_last.items() if k in _HEATER_PATCHES)
    vent_loss = sum(v for k, v in patch_last.items() if k in _VENT_PATCHES)
    wall_loss = sum(
        v for k, v in patch_last.items()
        if k not in _skip
    )

    avg_t = vol_avg_t[-1][1] if vol_avg_t else 0.0

    return HeatBalance(
        heater_input_W=heater_input,
        wall_loss_W=wall_loss,
        vent_loss_W=vent_loss,
        vol_avg_T=avg_t,
        patch_fluxes=patch_last,
    )
