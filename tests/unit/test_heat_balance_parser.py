"""Tests for heat balance parser module."""

from __future__ import annotations

from pathlib import Path

from harness.heat_balance_parser import (
    HeatBalance,
    compute_heat_balance,
    parse_vol_average_t,
    parse_wall_heat_flux,
)


class TestHeatBalance:
    def test_imbalance_calculation(self) -> None:
        hb = HeatBalance(heater_input_W=18000.0, wall_loss_W=-17500.0, vent_loss_W=-200.0)
        assert hb.imbalance_W == 300.0
        assert abs(hb.imbalance_pct - 1.667) < 0.01

    def test_zero_heater_no_division_error(self) -> None:
        hb = HeatBalance(heater_input_W=0.0, wall_loss_W=0.0)
        assert hb.imbalance_pct == 0.0

    def test_patch_fluxes_stored(self) -> None:
        hb = HeatBalance(
            heater_input_W=18000.0,
            wall_loss_W=-17000.0,
            patch_fluxes={"heater_wall": 18000.0, "floor": -5000.0},
        )
        assert hb.patch_fluxes["heater_wall"] == 18000.0


class TestParseWallHeatFlux:
    def test_empty_when_no_dir(self, tmp_path: Path) -> None:
        result = parse_wall_heat_flux(tmp_path)
        assert result == {}

    def test_per_patch_files(self, tmp_path: Path) -> None:
        pp = tmp_path / "postProcessing" / "wallHeatFlux" / "0"
        pp.mkdir(parents=True)

        (pp / "heater_wall.dat").write_text(
            "# Time wallHeatFlux\n"
            "100\t18000.5\n"
            "200\t18001.0\n",
            encoding="utf-8",
        )
        (pp / "floor.dat").write_text(
            "# Time wallHeatFlux\n"
            "100\t-5000.0\n"
            "200\t-5100.0\n",
            encoding="utf-8",
        )

        result = parse_wall_heat_flux(tmp_path)
        assert "heater_wall" in result
        assert "floor" in result
        assert len(result["heater_wall"]) == 2
        assert result["heater_wall"][-1] == (200.0, 18001.0)
        assert result["floor"][-1] == (200.0, -5100.0)

    def test_skips_comments_and_blanks(self, tmp_path: Path) -> None:
        pp = tmp_path / "postProcessing" / "wallHeatFlux" / "0"
        pp.mkdir(parents=True)

        (pp / "wall.dat").write_text(
            "# Header\n"
            "# Another comment\n"
            "\n"
            "100\t-3000.0\n",
            encoding="utf-8",
        )

        result = parse_wall_heat_flux(tmp_path)
        assert len(result["wall"]) == 1

    def test_strips_wallheatflux_suffix(self, tmp_path: Path) -> None:
        """OpenFOAM v2312 writes files like floor_wallHeatFlux.dat."""
        pp = tmp_path / "postProcessing" / "wallHeatFlux" / "0"
        pp.mkdir(parents=True)

        (pp / "heater_wall_wallHeatFlux.dat").write_text(
            "# Time wallHeatFlux\n200\t18000.0\n", encoding="utf-8"
        )
        (pp / "floor_wallHeatFlux.dat").write_text(
            "# Time wallHeatFlux\n200\t-5000.0\n", encoding="utf-8"
        )

        result = parse_wall_heat_flux(tmp_path)
        assert "heater_wall" in result
        assert "floor" in result
        assert "heater_wall_wallHeatFlux" not in result

    def test_multiple_time_dirs(self, tmp_path: Path) -> None:
        for t in ["0", "100"]:
            pp = tmp_path / "postProcessing" / "wallHeatFlux" / t
            pp.mkdir(parents=True)
            (pp / "wall.dat").write_text(
                f"# Time flux\n{t}\t-{int(t) + 1000}.0\n",
                encoding="utf-8",
            )

        result = parse_wall_heat_flux(tmp_path)
        assert len(result["wall"]) == 2


class TestParseVolAverageT:
    def test_empty_when_no_dir(self, tmp_path: Path) -> None:
        result = parse_vol_average_t(tmp_path)
        assert result == []

    def test_reads_values(self, tmp_path: Path) -> None:
        pp = tmp_path / "postProcessing" / "volAverageT" / "0"
        pp.mkdir(parents=True)

        (pp / "volFieldValue.dat").write_text(
            "# Time volAverage(T)\n"
            "100\t358.2\n"
            "200\t360.1\n",
            encoding="utf-8",
        )

        result = parse_vol_average_t(tmp_path)
        assert len(result) == 2
        assert result[-1] == (200.0, 360.1)


class TestComputeHeatBalance:
    def test_basic_balance(self) -> None:
        fluxes = {
            "heater_wall": [(100.0, 17500.0), (200.0, 18000.0)],
            "floor": [(100.0, -4000.0), (200.0, -5000.0)],
            "ceiling": [(100.0, -6000.0), (200.0, -7000.0)],
            "front": [(100.0, -2000.0), (200.0, -3000.0)],
            "back": [(100.0, -1500.0), (200.0, -2000.0)],
        }
        vol_avg = [(100.0, 355.0), (200.0, 358.2)]

        hb = compute_heat_balance(fluxes, vol_avg)
        assert hb.heater_input_W == 18000.0
        assert hb.wall_loss_W == -17000.0
        assert hb.vent_loss_W == 0.0
        assert hb.imbalance_W == 1000.0
        assert hb.vol_avg_T == 358.2

    def test_with_vent_patches(self) -> None:
        fluxes = {
            "heater_wall": [(200.0, 18000.0)],
            "floor": [(200.0, -5000.0)],
            "supply_vent": [(200.0, -500.0)],
            "exhaust_vent": [(200.0, -1000.0)],
        }

        hb = compute_heat_balance(fluxes)
        assert hb.heater_input_W == 18000.0
        assert hb.wall_loss_W == -5000.0
        assert hb.vent_loss_W == -1500.0

    def test_empty_fluxes(self) -> None:
        hb = compute_heat_balance({})
        assert hb.heater_input_W == 0.0
        assert hb.wall_loss_W == 0.0

    def test_no_vol_avg(self) -> None:
        fluxes = {"heater_wall": [(100.0, 18000.0)]}
        hb = compute_heat_balance(fluxes, vol_avg_t=None)
        assert hb.vol_avg_T == 0.0
