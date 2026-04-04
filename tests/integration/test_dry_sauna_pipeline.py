"""Integration tests for the dry sauna steady-state pipeline.

These tests validate the full pipeline WITHOUT requiring OpenFOAM.
WSL is mocked — only Python-side logic is tested end-to-end.
"""

from __future__ import annotations

from pathlib import Path

from harness.case_builder import build_case
from harness.kpi import evaluate_phase1_kpis
from harness.probe_parser import get_steady_state_values, parse_probe_file
from harness.schema import load_yaml
from harness.solver_runner import check_convergence, parse_residuals

SAMPLE_PROBE_OUTPUT = """\
# Probe 0 (1.5 2.0 1.25)
# Probe 1 (1.5 0.8 1.25)
# Probe 2 (1.5 0.1 1.25)
#        Time            0               1               2
0               293.15          293.15          293.15
100             340.5           315.2           298.7
200             355.8           325.1           302.3
300             358.2           327.5           303.1
"""

SAMPLE_SOLVER_LOG = """\
Time = 998

smoothSolver:  Solving for Ux, Initial residual = 2e-05, Final residual = 1e-07, No Iterations 3
smoothSolver:  Solving for Uy, Initial residual = 3e-05, Final residual = 2e-07, No Iterations 3
smoothSolver:  Solving for Uz, Initial residual = 1e-05, Final residual = 1e-07, No Iterations 3
GAMG:  Solving for p_rgh, Initial residual = 5e-05, Final residual = 3e-07, No Iterations 4
smoothSolver:  Solving for T, Initial residual = 8e-06, Final residual = 5e-08, No Iterations 3
smoothSolver:  Solving for k, Initial residual = 3e-05, Final residual = 2e-07, No Iterations 3
smoothSolver:  Solving for omega, Initial residual = 4e-05, Final residual = 3e-07, No Iterations 3

End
"""


class TestBuildToKpiPipeline:
    """Test the full pipeline: YAML -> build -> (mock solver) -> parse -> KPI."""

    def test_build_creates_valid_case(self, sample_case_path: Path, tmp_path: Path) -> None:
        """YAML -> build produces a valid OpenFOAM case structure."""
        case_dir = build_case(sample_case_path, output_dir=tmp_path / "case")
        assert (case_dir / "system" / "blockMeshDict").exists()
        assert (case_dir / "0" / "T").exists()
        assert (case_dir / "constant" / "thermophysicalProperties").exists()

    def test_probe_to_kpi(self, tmp_path: Path) -> None:
        """Pre-canned probe output -> parse -> KPI produces expected results."""
        probe_file = tmp_path / "T"
        probe_file.write_text(SAMPLE_PROBE_OUTPUT, encoding="utf-8")

        probe_names = ["upper_bench", "lower_bench", "floor_level"]
        data = parse_probe_file(probe_file, probe_names)
        values = get_steady_state_values(data)

        assert values["upper_bench"] == 358.2
        assert values["lower_bench"] == 327.5

        kpis = evaluate_phase1_kpis(values)
        k01 = kpis[0]
        assert k01.kpi_id == "K-01"
        assert k01.value == 30.7
        assert k01.pass_fail == "pass"

        k07 = kpis[1]
        assert k07.kpi_id == "K-07"
        assert k07.value > 0  # positive stratification

    def test_convergence_from_sample_log(self) -> None:
        """Sample solver log -> residual parse -> convergence check."""
        residuals = parse_residuals(SAMPLE_SOLVER_LOG)
        assert len(residuals) == 1
        assert check_convergence(residuals) is True

    def test_full_yaml_roundtrip(self, sample_case_path: Path, tmp_path: Path) -> None:
        """Verify YAML loads, builds, and probe names match."""
        data = load_yaml(sample_case_path)
        probe_names = [p["name"] for p in data["probes"]]
        assert probe_names == ["upper_bench", "lower_bench", "floor_level"]

        case_dir = build_case(sample_case_path, output_dir=tmp_path / "case")
        control = (case_dir / "system" / "controlDict").read_text(encoding="utf-8")
        for name in probe_names:
            assert name in control
