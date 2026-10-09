"""Keep historical module imports and executable entrypoints on one owner."""
from importlib import import_module
from runpy import run_module
import sys


def redirect(module_name, target):
    if module_name == "__main__":
        run_module(target, run_name="__main__", alter_sys=True)
    else:
        sys.modules[module_name] = import_module(target)
