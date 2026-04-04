"""Integration test for batch comparison pipeline."""

from __future__ import annotations

from pathlib import Path

from harness.batch import run_batch


class TestBatchPipeline:
    def test_full_batch_with_report(self, tmp_path: Path) -> None:
        yamls = sorted(Path("configs/cases").glob("*.yaml"))
        report = run_batch(yamls, max_iter=20000, n_profile=40)

        # Generate report
        md = report.summary_table()
        out = tmp_path / "batch_report.md"
        out.write_text(md, encoding="utf-8")

        assert out.exists()
        content = out.read_text(encoding="utf-8")
        assert "dry_sauna_steady" in content
        assert "aufguss_test" in content
        assert "K-01" in content

    def test_aufguss_lower_stratification(self) -> None:
        yamls = [
            Path("configs/cases/dry_sauna_steady.yaml"),
            Path("configs/cases/aufguss_test.yaml"),
        ]
        report = run_batch(yamls, max_iter=20000, n_profile=40)

        dry = next(c for c in report.cases if "dry" in c.case_name)
        aug = next(c for c in report.cases if "aufguss" in c.case_name)

        dry_dt = dry.result.upper_layer_temp - dry.result.lower_layer_temp
        aug_dt = aug.result.upper_layer_temp - aug.result.lower_layer_temp
        assert aug_dt < dry_dt
