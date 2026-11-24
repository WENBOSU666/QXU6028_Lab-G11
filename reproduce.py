#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
reproduce.py — regenerate every reported number and figure with one command.

    python reproduce.py
    python reproduce.py --seed 0 --out output

Runs L1 -> L2 -> L3 from the raw files in data/, then merges the stage results
into results.json (schema QXU6028-results-v1) at the repository root.
"""
from __future__ import annotations

import argparse
import json
import random
import runpy
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parent
SRC = ROOT / "src"
STAGES = [("L1", "l1_forward.py"), ("L2", "l2_inverse.py"), ("L3", "l3_map.py"),
          ("FIG", "report_figures.py")]


def pipeline_commit() -> str:
    """Current HEAD, read straight from .git (no subprocess needed)."""
    try:
        head = (ROOT / ".git" / "HEAD").read_text(encoding="utf-8").strip()
        if head.startswith("ref:"):
            ref = ROOT / ".git" / head.split(" ", 1)[1].strip()
            if ref.exists():
                return ref.read_text(encoding="utf-8").strip()
            packed = ROOT / ".git" / "packed-refs"
            if packed.exists():
                for line in packed.read_text(encoding="utf-8").splitlines():
                    if line.endswith(head.split(" ", 1)[1].strip()):
                        return line.split()[0]
        return head
    except Exception:
        return ""


def run_stage(key: str, script: str) -> dict:
    path = SRC / script
    if not path.exists():
        print(f"[skip] {key}: src/{script} missing")
        return {}
    print(f"\n=== {key}: src/{script} ===")
    t0 = time.time()
    sys.path.insert(0, str(SRC))
    try:
        ns = runpy.run_path(str(path), run_name="__main__")
    except SystemExit:
        ns = {}
    except Exception as exc:                       # keep going, report clearly
        print(f"[fail] {key} raised {type(exc).__name__}: {exc}")
        return {}
    print(f"[done] {key} in {time.time() - t0:.1f}s")

    out_file = ROOT / "output" / f"{key}_results.json"
    if key == "FIG":
        return {}
    if out_file.exists():
        return json.loads(out_file.read_text(encoding="utf-8"))
    if isinstance(ns.get("RESULT"), dict):
        return ns["RESULT"]
    return {}


def merge(stage: dict) -> dict:
    L1 = stage.get("L1") or {}
    L2 = stage.get("L2") or {}
    L3 = stage.get("L3") or {}

    l1 = {
        "correction_factor": L1.get("correction_factor"),
        "n_fermi_op_cm2": L1.get("n_fermi_op_cm2"),
        "mobility_corrected_cm2_Vs": L1.get("mobility_corrected_cm2_Vs"),
    }
    grid = 0.25
    half = max(float(L2.get("Eg_uncertainty_kT") or 0.0), grid / 2)
    l2 = {
        "model": L2.get("model"),
        "Eg_kT": L2.get("Eg_kT"),
        "Eg_uncertainty_kT": (round(half, 4) if L2 else None),
        "EF_trajectory_eV": L2.get("EF_trajectory_eV", []),
    }
    l3 = {
        "map": L3.get("map", []),
        "features": L3.get("features", []),
        "resolution_probe_pitch": L3.get("resolution_probe_pitch"),
        "detection_limit_fraction": L3.get("detection_limit_fraction"),
        "pass": L3.get("pass"),
        "confidence": L3.get("confidence"),
    }
    return {
        "schema_version": "QXU6028-results-v1",
        "group": "G11",
        "dataset_set": "regular",
        "pipeline_commit": pipeline_commit(),
        "L1": l1,
        "L2": l2,
        "L3": l3,
        "diagnostics": {
            "C_at_010eV": L1.get("C_at_010eV"),
            "C_at_020eV": L1.get("C_at_020eV"),
            "rms_cone": L2.get("rms_cone"),
            "rms_selected": L2.get("rms_selected"),
            "bootstrap_gap_kT": [L2.get("bootstrap_gap_lo"), L2.get("bootstrap_gap_hi")],
            "hard_gap_fraction": L2.get("hard_gap_fraction"),
            "l3_available": bool(L3),
        },
    }


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--out", default="output")
    args = ap.parse_args()

    random.seed(args.seed)
    try:
        import numpy as np
        np.random.seed(args.seed)
    except ImportError:
        pass

    print(f"repo  : {ROOT}")
    print(f"seed  : {args.seed}")
    (ROOT / args.out).mkdir(parents=True, exist_ok=True)

    stage = {key: run_stage(key, script) for key, script in STAGES}
    payload = merge(stage)
    text = json.dumps(payload, indent=2, ensure_ascii=False, allow_nan=False) + "\n"
    (ROOT / args.out / "results.json").write_text(text, encoding="utf-8")
    (ROOT / "results.json").write_text(text, encoding="utf-8")

    print("\n" + "=" * 60)
    print(f"L1 correction factor : {payload['L1']['correction_factor']}")
    print(f"L2 model / Eg (kBT)  : {payload['L2']['model']} / {payload['L2']['Eg_kT']}")
    print(f"L3 available         : {payload['diagnostics']['l3_available']}")
    print(f"results.json written to {ROOT / 'results.json'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
