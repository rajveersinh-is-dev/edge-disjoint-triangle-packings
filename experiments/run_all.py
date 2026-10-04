"""Run every experiment in order and report a single pass/fail summary.

Usage
-----
    python experiments/run_all.py            # everything (a few minutes)
    python experiments/run_all.py --quick    # skip the exhaustive scans

Each experiment writes its outputs to ``experiments/results/``.  A non-zero exit
status means that at least one experiment reported a failure; the individual
scripts print exactly which check failed.
"""

from __future__ import annotations

import os
import subprocess
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)

EXPERIMENTS = [
    ("exp1_convex", "convex extremal function by DP and by full enumeration", False),
    ("exp3_weakdual", "structural lemmas (3-cycles are faces; all subcubic trees occur)", False),
    ("exp2_falsify", "Conjecture 7.1: search for a configuration below convex position", False),
    ("exp4_general", "sharp bound n-2, Eulerian constructions, all maximal planar graphs", True),
]

FIGURES = "make_figures"


def run(script: str, env_extra: dict | None = None) -> tuple[int, float]:
    env = dict(os.environ)
    if env_extra:
        env.update(env_extra)
    t0 = time.time()
    proc = subprocess.run([sys.executable, os.path.join(HERE, script)],
                          cwd=ROOT, env=env)
    return proc.returncode, time.time() - t0


def main() -> int:
    quick = "--quick" in sys.argv
    print("=" * 78)
    print("planetri: full reproduction run")
    print("=" * 78)
    failures: list[str] = []
    for name, description, expensive in EXPERIMENTS:
        if quick and expensive:
            print(f"\n--- SKIPPED (--quick): {name} ---")
            continue
        print(f"\n--- {name}: {description} ---")
        extra = {"PLANETRI_MAX_ENUM_N": "6"} if (quick and expensive) else None
        code, dt = run(f"{name}.py", extra)
        print(f"--- {name}: exit={code} in {dt:.1f}s ---")
        if code != 0:
            failures.append(name)

    print(f"\n--- {FIGURES}: figures ---")
    code, dt = run(f"{FIGURES}.py")
    print(f"--- {FIGURES}: exit={code} in {dt:.1f}s ---")
    if code != 0:
        failures.append(FIGURES)

    print("\n" + "=" * 78)
    if failures:
        print("FAILED experiments:", ", ".join(failures))
        return 1
    print("ALL EXPERIMENTS PASSED")
    print("outputs in experiments/results/ and figures/")
    print("=" * 78)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())