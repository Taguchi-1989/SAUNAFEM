"""CLI entry point for SaunaFlow harness."""

import click

from harness import __version__


@click.group()
@click.version_option(version=__version__, prog_name="saunaflow")
def main() -> None:
    """SaunaFlow - Python harness for sauna CFD simulation."""


@main.command()
@click.argument("case_yaml", type=click.Path(exists=True))
def validate(case_yaml: str) -> None:
    """Validate a YAML case definition against the schema."""
    from harness.schema import load_and_validate

    errors = load_and_validate(case_yaml)
    if errors:
        for err in errors:
            click.echo(f"ERROR: {err}", err=True)
        raise SystemExit(1)
    click.echo("Valid.")


@main.command()
@click.argument("case_yaml", type=click.Path(exists=True))
def build(case_yaml: str) -> None:
    """Build OpenFOAM case directory from YAML definition."""
    click.echo(f"[stub] build: {case_yaml}")


@main.command()
@click.argument("case_dir", type=click.Path(exists=True))
def run(case_dir: str) -> None:
    """Run OpenFOAM solver on a built case directory."""
    click.echo(f"[stub] run: {case_dir}")


@main.command()
@click.argument("case_dir", type=click.Path(exists=True))
def report(case_dir: str) -> None:
    """Generate validation report from solver results."""
    click.echo(f"[stub] report: {case_dir}")


if __name__ == "__main__":
    main()
