"""Artifact roles for new outputs and explicit historical-layout adapters."""
import os
import re
from pathlib import Path


ARTIFACT_FOLDERS = {
    "model": "00_model", "data": "01_results", "png": "10_overview",
    "pdf": "11_pdf", "report": "12_reports", "svg": "13_svg",
    "log": "80_logs", "history": "90_history", "config": "99_config",
    "html": "",
}


def artifact_path(root, name, *, role=None, group=None, create_parent=True):
    """Classify by purpose; JSON requires an explicit role."""
    name = str(name)
    if not name or Path(name).name != name or name in (".", "..") or "\\" in name or ":" in name:
        raise ValueError("Expected a single artifact filename")
    suffix = Path(name).suffix.lower()
    if role is None:
        role = {".npz": "data", ".parquet": "data", ".csv": "data",
                ".png": "png", ".pdf": "pdf", ".svg": "svg", ".md": "report",
                ".log": "log", ".zip": "history", ".html": "html"}.get(suffix)
    if role not in ARTIFACT_FOLDERS:
        raise ValueError("Supply a supported artifact role (JSON needs log, data or config)")
    path = Path(root) / ARTIFACT_FOLDERS[role]
    if group is not None:
        group = Path(group)
        if group.is_absolute() or group.drive or ".." in group.parts:
            raise ValueError("Artifact group must stay inside its format directory")
        path /= group
    path /= name
    if create_parent:
        path.parent.mkdir(parents=True, exist_ok=True)
    return path


def table_path(log_directory, name, *, existing=False):
    """Route numeric tables, retaining batch groups and read fallback to v1."""
    directory = Path(log_directory)
    folder = next((p for p in (directory, *directory.parents) if p.name == "80_logs"), None)
    if folder is None:
        raise ValueError("Expected a directory under 80_logs")
    unified = os.environ.get("COMSOL_WORKFLOW_OUTPUT_LAYOUT") == "unified"
    roles = ("data", "log") if unified else ("log", "data")
    group = directory.relative_to(folder)
    paths = [artifact_path(folder.parent, name, role=role, group=group,
                           create_parent=False) for role in roles]
    if existing:
        return next((p for p in paths if p.is_file()), paths[0])
    paths[0].parent.mkdir(parents=True, exist_ok=True)
    return paths[0]


def table_links(text):
    """Keep generated report links consistent with this run's table layout."""
    if os.environ.get("COMSOL_WORKFLOW_OUTPUT_LAYOUT") == "unified":
        return re.sub(r"80_logs/([^\s`\"'<>)]*\.csv)", r"01_results/\1", text)
    return text


def output_path(root, name):
    """Keep v1 paths; public workflow commands opt in to role-based data paths."""
    name = str(name)
    suffix = Path(name).suffix.lower()
    if os.environ.get("COMSOL_WORKFLOW_OUTPUT_LAYOUT") == "unified":
        role = "config" if name.endswith("_preflight.json") else "log" if suffix == ".json" else None
        return artifact_path(root, name, role=role)
    role = "config" if name.endswith("_preflight.json") else {
        ".npz": "data", ".png": "png", ".pdf": "pdf", ".md": "report",
        ".csv": "log", ".json": "log", ".log": "log", ".zip": "history",
        ".html": "html",
    }.get(suffix)
    return artifact_path(root, name, role=role)


def data_directory(case_dir):
    """S4 compact and legacy data layouts; preserve compact-folder precedence."""
    case_dir = Path(case_dir)
    compact = case_dir / "80_logs"
    return compact if compact.is_dir() else case_dir / "01_results"
