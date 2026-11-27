#!/usr/bin/env python3
"""Install stub packages so transformers will load the Evo 2 remote code.

transformers' check_imports() statically scans trust_remote_code modelling
files and refuses to load unless every imported package resolves -- even when
every one of those imports is lazy and guarded. The Evo 2 port imports:

  flash_attn          attention.py (inside a function), rotary.py (try/except),
                      layers.py (behind use_flash_rmsnorm, default False)
  transformer_engine  layers.py, inside `try: ... except ImportError: HAS_TE = False`

Neither is reachable when running with attn_implementation="sdpa" and
use_fp8_input_projections=False. Installing the real packages would mean a long
CUDA source build for flash-attn and a conda-only dependency for TE, to satisfy
a static scan for code that never executes.

So: create packages that EXIST (find_spec succeeds, check_imports passes) but
whose contents are unavailable. Importing a submodule raises ModuleNotFoundError,
a subclass of ImportError, which the port's own try/except already handles by
falling back to the pure-PyTorch path. Any attribute access on the top level
raises with an explanatory message, so a genuine dependency can never be
satisfied silently by a stub.

This trades bit-exactness for runnability: sdpa is not bit-identical to FA2,
and bf16 projections are not identical to FP8. For the swap control that is
acceptable, because both arms of the paired contrast run the same machinery and
the drift cancels in the delta. It would NOT be acceptable for reproducing
published absolute Evo 2 numbers.

Usage:
  python install_import_stubs.py            # into the running interpreter's site-packages
  python install_import_stubs.py --check    # verify only, change nothing
"""

from __future__ import annotations

import argparse
import importlib
import importlib.util
import site
import sys
from pathlib import Path

STUBS = {
    "flash_attn": "attn_implementation='sdpa' avoids every real use of this package.",
    "transformer_engine": "use_fp8_input_projections=False avoids every real use of this package.",
}

TEMPLATE = '''"""Stub package: present so transformers check_imports() passes.

{why}

Installed by scripts/install_import_stubs.py. Nothing here is ever called on
the intended code path. Submodule imports raise ModuleNotFoundError, which the
Evo 2 port's own try/except handles by falling back to pure PyTorch.
"""

__version__ = "0.0.0-stub"
__all__: list = []


def __getattr__(name):
    raise ImportError(
        "{pkg}." + name + " was actually requested, but this is a STUB "
        "installed only to satisfy transformers' static import check. "
        "{why} "
        "If you need the real package, install it and remove this stub; do not "
        "ignore this error, because the numerics would otherwise be wrong."
    )
'''


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--check", action="store_true", help="verify only")
    ap.add_argument("--site-packages", type=Path, default=None)
    args = ap.parse_args()

    sp = args.site_packages or Path(site.getsitepackages()[0])
    print(f"site-packages: {sp}")

    rc = 0
    for pkg, why in STUBS.items():
        spec = importlib.util.find_spec(pkg)
        target = sp / pkg
        is_stub = (target / "__init__.py").is_file() and "0.0.0-stub" in (
            target / "__init__.py").read_text() if (target / "__init__.py").is_file() else False

        if spec is not None and not is_stub:
            print(f"  {pkg:<20} REAL package present at {spec.origin} -- leaving alone")
            continue
        if args.check:
            state = "stub present" if is_stub else "MISSING"
            print(f"  {pkg:<20} {state}")
            if not is_stub:
                rc = 1
            continue

        target.mkdir(parents=True, exist_ok=True)
        (target / "__init__.py").write_text(TEMPLATE.format(pkg=pkg, why=why))
        print(f"  {pkg:<20} stub written to {target}")

    if not args.check:
        for pkg in STUBS:
            importlib.invalidate_caches()
            m = importlib.import_module(pkg)
            print(f"  verify: import {pkg} -> {getattr(m, '__version__', '?')}")
    return rc


if __name__ == "__main__":
    raise SystemExit(main())
