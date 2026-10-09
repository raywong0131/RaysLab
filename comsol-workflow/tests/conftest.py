"""Keep validation offline and generated artifacts outside archive/source trees."""
import builtins
from datetime import datetime, timezone
import os
from pathlib import Path
import uuid


def pytest_configure(config):
    if not config.option.basetemp:
        tag = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ") + "_" + uuid.uuid4().hex[:6]
        config.option.basetemp = str(Path(config.rootpath) / "tests/.work" / tag)
    base = Path(config.option.basetemp).resolve()
    previous_import = builtins.__import__
    config._workflow_restore = (previous_import, {key: os.environ.get(key) for key in
                              ("COMSOL_WORKFLOW_TEST_OUTPUT_ROOT", "MPLCONFIGDIR")})
    os.environ["COMSOL_WORKFLOW_TEST_OUTPUT_ROOT"] = str(base / "exports")
    os.environ.setdefault("MPLCONFIGDIR", str(base / "matplotlib"))

    def forbidden_start(*_args, **_kwargs):
        raise AssertionError("COMSOL startup forbidden during offline validation")

    def guarded_import(name, *args, **kwargs):
        module = previous_import(name, *args, **kwargs)
        if name == "mph" and getattr(module, "__file__", None):
            module.start = forbidden_start
        return module

    builtins.__import__ = guarded_import


def pytest_unconfigure(config):
    if not hasattr(config, "_workflow_restore"):
        return
    builtins.__import__, environment = config._workflow_restore
    for key, value in environment.items():
        if value is None:
            os.environ.pop(key, None)
        else:
            os.environ[key] = value
