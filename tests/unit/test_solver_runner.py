"""Tests for solver_runner module."""

from __future__ import annotations

from pathlib import Path
from unittest.mock import MagicMock, patch

from harness.solver_runner import check_convergence, parse_residuals, run_solver

SAMPLE_LOG = """\
Create time

Create mesh for time = 0

SIMPLE: No convergence criteria found

Time = 1

smoothSolver:  Solving for Ux, Initial residual = 1, Final residual = 0.1, No Iterations 1
smoothSolver:  Solving for Uy, Initial residual = 1, Final residual = 0.1, No Iterations 1
smoothSolver:  Solving for Uz, Initial residual = 1, Final residual = 0.1, No Iterations 1
GAMG:  Solving for p_rgh, Initial residual = 1, Final residual = 0.01, No Iterations 3
smoothSolver:  Solving for T, Initial residual = 1, Final residual = 0.001, No Iterations 1
smoothSolver:  Solving for k, Initial residual = 1, Final residual = 0.1, No Iterations 1
smoothSolver:  Solving for epsilon, Initial residual = 1, Final residual = 0.1, No Iterations 1

Time = 2

smoothSolver:  Solving for Ux, Initial residual = 0.5, Final residual = 0.05, No Iterations 2
smoothSolver:  Solving for Uy, Initial residual = 0.5, Final residual = 0.05, No Iterations 2
smoothSolver:  Solving for Uz, Initial residual = 0.5, Final residual = 0.05, No Iterations 2
GAMG:  Solving for p_rgh, Initial residual = 0.5, Final residual = 0.005, No Iterations 3
smoothSolver:  Solving for T, Initial residual = 0.5, Final residual = 0.0005, No Iterations 2
smoothSolver:  Solving for k, Initial residual = 0.5, Final residual = 0.05, No Iterations 2
smoothSolver:  Solving for epsilon, Initial residual = 0.5, Final residual = 0.05, No Iterations 2

Time = 3

smoothSolver:  Solving for Ux, Initial residual = 1e-05, Final residual = 1e-07, No Iterations 3
smoothSolver:  Solving for Uy, Initial residual = 1e-05, Final residual = 1e-07, No Iterations 3
smoothSolver:  Solving for Uz, Initial residual = 1e-05, Final residual = 1e-07, No Iterations 3
GAMG:  Solving for p_rgh, Initial residual = 1e-05, Final residual = 1e-07, No Iterations 3
smoothSolver:  Solving for T, Initial residual = 1e-06, Final residual = 1e-08, No Iterations 3
smoothSolver:  Solving for k, Initial residual = 1e-05, Final residual = 1e-07, No Iterations 3
smoothSolver:  Solving for omega, Initial residual = 1e-05, Final residual = 1e-07, No Iterations 3

End
"""


class TestParseResiduals:
    def test_parses_three_iterations(self) -> None:
        residuals = parse_residuals(SAMPLE_LOG)
        assert len(residuals) == 3

    def test_first_iteration_values(self) -> None:
        residuals = parse_residuals(SAMPLE_LOG)
        assert residuals[0]["Ux"] == 0.1
        assert residuals[0]["T"] == 0.001
        assert residuals[0]["p_rgh"] == 0.01

    def test_last_iteration_converged(self) -> None:
        residuals = parse_residuals(SAMPLE_LOG)
        assert residuals[2]["Ux"] == 1e-7
        assert residuals[2]["T"] == 1e-8

    def test_empty_log(self) -> None:
        assert parse_residuals("") == []

    def test_partial_log(self) -> None:
        partial = (
            "smoothSolver:  Solving for Ux, Initial residual = 0.5,"
            " Final residual = 0.05, No Iterations 2\n"
        )
        residuals = parse_residuals(partial)
        assert len(residuals) == 1
        assert residuals[0]["Ux"] == 0.05


class TestCheckConvergence:
    def test_converged(self) -> None:
        residuals = parse_residuals(SAMPLE_LOG)
        assert check_convergence(residuals) is True

    def test_not_converged(self) -> None:
        residuals = [{"Ux": 0.5, "T": 0.1, "p_rgh": 0.5}]
        assert check_convergence(residuals) is False

    def test_empty_residuals(self) -> None:
        assert check_convergence([]) is False

    def test_custom_thresholds(self) -> None:
        residuals = [{"T": 1e-4}]
        assert check_convergence(residuals, thresholds={"T": 1e-3}) is True
        assert check_convergence(residuals, thresholds={"T": 1e-5}) is False

    def test_uses_final_residuals_for_convergence(self) -> None:
        log = (
            "smoothSolver:  Solving for Ux, Initial residual = 1e-2, "
            "Final residual = 1e-7, No Iterations 3\n"
            "smoothSolver:  Solving for Uy, Initial residual = 1e-2, "
            "Final residual = 1e-7, No Iterations 3\n"
            "smoothSolver:  Solving for Uz, Initial residual = 1e-2, "
            "Final residual = 1e-7, No Iterations 3\n"
            "GAMG:  Solving for p_rgh, Initial residual = 1e-2, "
            "Final residual = 1e-7, No Iterations 3\n"
            "smoothSolver:  Solving for T, Initial residual = 1e-3, "
            "Final residual = 1e-8, No Iterations 3\n"
            "smoothSolver:  Solving for k, Initial residual = 1e-2, "
            "Final residual = 1e-7, No Iterations 3\n"
            "smoothSolver:  Solving for epsilon, Initial residual = 1e-2, "
            "Final residual = 1e-7, No Iterations 3\n"
        )
        residuals = parse_residuals(log)
        assert check_convergence(residuals) is True


class TestRunSolver:
    @patch("harness.solver_runner.wsl_exec")
    def test_success(self, mock_wsl, tmp_path: Path) -> None:
        mock_wsl.return_value = MagicMock(stdout=SAMPLE_LOG, returncode=0)
        result = run_solver(tmp_path, solver_name="buoyantSimpleFoam")
        assert result.success is True
        assert result.iterations == 3
        assert result.converged is True
        assert result.log_path == tmp_path / "log.buoyantSimpleFoam"
        assert result.log_path.exists()
