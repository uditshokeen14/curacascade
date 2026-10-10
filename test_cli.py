"""Test suite for CuraCascade headless CLI interface in run.py."""

from __future__ import annotations

import csv
import json
from pathlib import Path
import subprocess
import sys
import pytest

_ROOT_DIR = Path(__file__).resolve().parent
PYTHON = sys.executable


@pytest.fixture(autouse=True)
def setup_test_env(tmp_path: Path, monkeypatch: pytest.MonkeyPatch):
    """Set up an isolated sandbox environment with required database artifacts."""
    import shutil
    data_dir = tmp_path / "data"
    data_dir.mkdir(parents=True, exist_ok=True)
    
    for db_name in ("curacascade.db", "curacascade_full.db"):
        src_db = _ROOT_DIR / "data" / db_name
        if src_db.exists():
            shutil.copy2(src_db, data_dir / db_name)
            
    monkeypatch.setenv("CURACASCADE_HOME", str(tmp_path))


def run_cli(*args: str) -> subprocess.CompletedProcess:
    """Helper to execute run.py via subprocess."""
    import os
    env = os.environ.copy()
    env["PYTHONIOENCODING"] = "utf-8"
    cmd = [PYTHON, str(_ROOT_DIR / "run.py"), *args]
    return subprocess.run(
        cmd,
        capture_output=True,
        text=True,
        encoding="utf-8",
        errors="replace",
        env=env,
        cwd=str(_ROOT_DIR),
    )


def test_help_messages():
    """Verify top-level and subcommand help messages exit with 0."""
    res_root = run_cli("--help")
    assert res_root.returncode == 0
    assert "CuraCascade" in res_root.stdout
    assert "audit" in res_root.stdout
    assert "verify" in res_root.stdout
    assert "serve" in res_root.stdout

    res_audit = run_cli("audit", "--help")
    assert res_audit.returncode == 0
    assert "--genes" in res_audit.stdout
    assert "--class" in res_audit.stdout
    assert "--format" in res_audit.stdout
    assert "--output" in res_audit.stdout

    res_verify = run_cli("verify", "--help")
    assert res_verify.returncode == 0
    assert "--db" in res_verify.stdout

    res_serve = run_cli("serve", "--help")
    assert res_serve.returncode == 0
    assert "--port" in res_serve.stdout


def test_audit_gene_with_class_a_exits_1():
    """SCN1B has Class A findings in curacascade.db; CLI should exit with code 1."""
    res = run_cli("audit", "--genes", "SCN1B")
    assert res.returncode == 1
    assert "SCN1B" in res.stdout
    assert "Class A" in res.stdout
    assert "CRITICAL AUDIT ALERT" in res.stdout


def test_audit_gene_without_class_a_exits_0():
    """BARD1 has only Class B findings in curacascade.db; CLI should exit with code 0."""
    res = run_cli("audit", "--genes", "BARD1")
    assert res.returncode == 0
    assert "BARD1" in res.stdout
    assert "Class B" in res.stdout
    assert "CLEAN AUDIT" in res.stdout


def test_audit_clean_gene_exits_0():
    """ACTN2 has 0 findings in curacascade.db; CLI should exit with code 0."""
    res = run_cli("audit", "--genes", "ACTN2")
    assert res.returncode == 0
    assert "CLEAN AUDIT" in res.stdout
    assert "No inconsistency findings matched" in res.stdout


def test_audit_class_filter():
    """Filtering by class B on SCN1B returns 0 findings and exits with code 0."""
    res_b = run_cli("audit", "--genes", "SCN1B", "--class", "B")
    assert res_b.returncode == 0
    assert "CLEAN AUDIT" in res_b.stdout

    # Filtering by class A on SCN1B returns findings and exits with code 1
    res_a = run_cli("audit", "--genes", "SCN1B", "--class", "A")
    assert res_a.returncode == 1
    assert "Class A" in res_a.stdout


def test_audit_format_json_stdout():
    """Ensure json format to stdout produces valid parseable JSON."""
    res = run_cli("audit", "--genes", "BARD1", "--format", "json")
    assert res.returncode == 0
    data = json.loads(res.stdout)
    assert isinstance(data, list)
    assert len(data) >= 1
    assert data[0]["gene_symbol"] == "BARD1"


def test_audit_format_csv_stdout():
    """Ensure csv format to stdout produces valid CSV with headers."""
    res = run_cli("audit", "--genes", "BARD1", "--format", "csv")
    assert res.returncode == 0
    reader = csv.DictReader(res.stdout.strip().splitlines())
    rows = list(reader)
    assert len(rows) >= 1
    assert rows[0]["gene_symbol"] == "BARD1"


def test_audit_output_to_file(tmp_path: Path):
    """Ensure --output correctly writes file and reports summary."""
    json_out = tmp_path / "findings.json"
    res = run_cli("audit", "--genes", "SCN1B", "--format", "json", "--output", str(json_out))
    assert res.returncode == 1  # SCN1B has Class A
    assert json_out.is_file()

    content = json.loads(json_out.read_text(encoding="utf-8"))
    assert isinstance(content, list)
    assert any(r["audit_class"] == "Class A" for r in content)

    csv_out = tmp_path / "findings.csv"
    res_csv = run_cli("audit", "--genes", "BARD1", "--format", "csv", "--output", str(csv_out))
    assert res_csv.returncode == 0
    assert csv_out.is_file()
    reader = csv.DictReader(csv_out.read_text(encoding="utf-8").strip().splitlines())
    rows = list(reader)
    assert len(rows) >= 1


def test_verify_success():
    """Verify suite passes on valid benchmark database and exits with code 0."""
    res = run_cli("verify")
    assert res.returncode == 0
    assert "ALL INTEGRITY & SECURITY CHECKS PASSED" in res.stdout


def test_verify_failure_exits_2():
    """Verify suite fails on missing/corrupt database and exits with code 2."""
    res = run_cli("verify", "--db", "scratch/non_existent.db")
    assert res.returncode == 2
    assert "ONE OR MORE INTEGRITY CHECKS FAILED" in res.stdout or "FAIL" in res.stdout
