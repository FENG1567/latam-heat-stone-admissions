#!/usr/bin/env python3
"""Portable launcher for the preserved LatAm stone analysis scripts.

This wrapper copies source templates into an isolated workspace and delegates
to their original CLIs. It does not reimplement or modify model algorithms.
"""

from __future__ import annotations

import argparse
import json
import os
import shutil
import subprocess
import sys
from pathlib import Path
from typing import Any


THREAD_ENV = {
    "OMP_NUM_THREADS": "1",
    "OPENBLAS_NUM_THREADS": "1",
    "MKL_NUM_THREADS": "1",
    "BLAS_NUM_THREADS": "1",
    "VECLIB_MAXIMUM_THREADS": "1",
    "NUMEXPR_NUM_THREADS": "1",
}


def package_paths() -> tuple[Path, Path, Path]:
    """Return (package root, analysis tree, manifest) in staged or final layout."""
    script = Path(__file__).resolve()
    analysis_tree = script.parent.parent
    if (analysis_tree / "config" / "module_manifest.json").is_file():
        package_root = analysis_tree.parent if analysis_tree.name == "analysis_code" else analysis_tree
        manifest = analysis_tree / "config" / "module_manifest.json"
    else:
        package_root = analysis_tree.parent
        manifest = package_root / "config" / "module_manifest.json"
    return package_root, analysis_tree, manifest


def load_manifest(path: Path) -> dict[str, Any]:
    try:
        config = json.loads(path.read_text(encoding="utf-8"))
    except FileNotFoundError as exc:
        raise RuntimeError(f"module manifest is missing: {path}") from exc
    if config.get("schema_version") != 1 or not isinstance(config.get("modules"), dict):
        raise RuntimeError(f"unsupported module manifest schema: {path}")
    return config


def workspace_path(package_root: Path, configured: str | None) -> Path:
    value = configured or os.environ.get("ANALYSIS_WORKSPACE") or "build/analysis_workspace"
    path = Path(value).expanduser()
    if not path.is_absolute():
        path = package_root / path
    return path.resolve()


def safe_relative(root: Path, value: str) -> Path:
    rel = Path(value)
    if rel.is_absolute() or ".." in rel.parts:
        raise RuntimeError(f"configured path must stay under its root: {value}")
    resolved = (root / rel).resolve()
    if resolved != root.resolve() and root.resolve() not in resolved.parents:
        raise RuntimeError(f"configured path escapes its root: {value}")
    return resolved


def relocation_map(config: dict[str, Any], workspace: Path) -> list[tuple[str, str]]:
    replacements = []
    for item in config.get("historical_root_replacements", []):
        old = str(item["from"])
        new = str(item["to"]).replace("{workspace}", workspace.as_posix())
        replacements.append((old, new))
    return sorted(replacements, key=lambda pair: len(pair[0]), reverse=True)


def materialize(package_root: Path, config: dict[str, Any], workspace: Path) -> None:
    """Copy templates once; refuse to replace any differing workspace file."""
    template_root = package_root / config.get("template_root", "analysis/templates")
    if not template_root.is_dir():
        raise RuntimeError(
            f"analysis templates are missing: {template_root}; expected the preserved source tree under analysis/templates"
        )
    replacements = relocation_map(config, workspace)
    copied = 0
    for group in config.get("template_groups", []):
        source_root = safe_relative(package_root, group["source"])
        target_root = safe_relative(workspace, group["target"])
        if not source_root.is_dir():
            raise RuntimeError(f"template group is missing: {source_root}")
        for source in sorted(source_root.rglob("*")):
            if not source.is_file():
                continue
            rel = source.relative_to(source_root)
            target = target_root / rel
            if source.suffix.lower() in {".r", ".py", ".sh"}:
                text = source.read_text(encoding="utf-8")
                for old, new in replacements:
                    text = text.replace(old, new)
                content = text.encode("utf-8")
            else:
                content = source.read_bytes()
            target.parent.mkdir(parents=True, exist_ok=True)
            if target.exists():
                if target.read_bytes() != content:
                    raise RuntimeError(
                        f"refusing to overwrite a different workspace file: {target}; choose a fresh --workspace"
                    )
                continue
            target.write_bytes(content)
            copied += 1
    print(f"Materialized {copied} source file(s) under {workspace}")


def source_script(package_root: Path, config: dict[str, Any], module: dict[str, Any]) -> Path:
    return safe_relative(package_root, module["script"])


def workspace_script(workspace: Path, config: dict[str, Any], module: dict[str, Any]) -> Path:
    rel = module.get("workspace_script")
    if rel:
        return safe_relative(workspace, rel)
    source = Path(module["script"])
    for group in config.get("template_groups", []):
        source_root = Path(group["source"])
        try:
            rel = source.relative_to(source_root)
        except ValueError:
            continue
        return safe_relative(workspace, str(Path(group["target"]) / rel))
    raise RuntimeError(f"no template group maps script path: {module['script']}")


def check_r_packages(packages: list[str], rscript: str) -> list[str]:
    if not packages:
        return []
    expr = "missing <- commandArgs(TRUE); bad <- missing[!vapply(missing, requireNamespace, logical(1), quietly=TRUE)]; if(length(bad)) {cat(paste(bad, collapse=', ')); quit(status=2)}"
    proc = subprocess.run(
        [rscript, "--vanilla", "-e", expr, *packages],
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        check=False,
        env={**os.environ, **THREAD_ENV},
    )
    if proc.returncode:
        return [f"required R package(s) unavailable: {proc.stdout.strip() or ', '.join(packages)}"]
    return []


def check_parse(script: Path, rscript: str) -> list[str]:
    if script.suffix.lower() != ".r":
        return []
    proc = subprocess.run(
        [rscript, "--vanilla", "-e", "invisible(parse(file=commandArgs(TRUE)[1]))", str(script)],
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        check=False,
        env={**os.environ, **THREAD_ENV},
    )
    return [] if proc.returncode == 0 else [f"R parse failed for {script}: {proc.stdout.strip()}"]


def preflight_one(
    package_root: Path,
    config: dict[str, Any],
    module_id: str,
    module: dict[str, Any],
    workspace: Path,
    rscript: str,
    use_workspace_script: bool = False,
) -> list[str]:
    issues: list[str] = []
    status = module.get("status", "runnable")
    if status != "runnable":
        issues.append(f"manual/external prerequisite: {module.get('manual_reason', status)}")
    script = workspace_script(workspace, config, module) if use_workspace_script else source_script(package_root, config, module)
    if not script.is_file():
        issues.append(f"source script missing: {script}")
    else:
        issues.extend(check_parse(script, rscript))
    for rel in module.get("requires_inputs", []):
        path = safe_relative(workspace, rel)
        if not path.is_file():
            issues.append(f"required input missing: {path}")
    for rel in module.get("requires_outputs", []):
        path = safe_relative(workspace, rel)
        if not path.is_file():
            issues.append(f"required prior output missing: {path}")
    if module.get("r_packages"):
        issues.extend(check_r_packages(module["r_packages"], rscript))
    return issues


def format_command(module: dict[str, Any], workspace: Path, rscript: str, workers: int) -> list[str]:
    executable = module.get("executable", "Rscript")
    resolved_executable = rscript if executable == "Rscript" else executable
    values = {
        "workspace": workspace.as_posix(),
        "workers": str(workers),
        "python": sys.executable,
    }
    return [resolved_executable, *[str(arg).format(**values) for arg in module["args"]]]


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--workspace", help="isolated analysis workspace (default: build/analysis_workspace)")
    parser.add_argument("--rscript", default=os.environ.get("RSCRIPT", "Rscript"), help="Rscript executable")
    sub = parser.add_subparsers(dest="command", required=True)
    sub.add_parser("list", help="list analysis modules and prerequisites")
    prep = sub.add_parser("prepare", help="assemble input shards using scripts/prepare_inputs.py")
    prep.add_argument("--verify-only", action="store_true", help="verify input manifest and assembled targets only")
    pre = sub.add_parser("preflight", help="check templates, input files, R packages, and R syntax")
    pre.add_argument("module", nargs="?", help="module id; omit to check every module")
    run = sub.add_parser("run", help="run one configured module using its original CLI")
    run.add_argument("module", help="module id from `list`")
    run.add_argument("--workers", type=int, default=1, help="worker count passed where supported")
    return parser


def main(argv: list[str] | None = None) -> int:
    package_root, _analysis_tree, manifest_path = package_paths()
    config = load_manifest(manifest_path)
    args = build_parser().parse_args(argv)
    workspace = workspace_path(package_root, args.workspace)
    modules = config["modules"]

    if args.command == "list":
        print("Analysis modules (commands delegate to preserved scripts):")
        for module_id, module in modules.items():
            status = module.get("status", "runnable")
            needs = ", ".join(module.get("requires_outputs", []) + module.get("requires_inputs", [])) or "none declared"
            dependencies = ", ".join(module.get("depends_on", [])) or "none"
            executable = module.get("executable", "Rscript")
            command = " ".join([executable, *module["args"]])
            print(f"- {module_id} [{status}]: {module['label']}")
            print(f"  depends_on: {dependencies}")
            print(f"  requires: {needs}")
            print(f"  command: {command}")
        print(f"Workspace: {workspace}")
        return 0

    if args.command == "prepare":
        helper = package_root / "scripts" / "prepare_inputs.py"
        manifest = package_root / "data" / "analysis_inputs" / "manifest.json"
        if not helper.is_file():
            print(f"ERROR: input assembly helper is missing: {helper}", file=sys.stderr)
            return 2
        if not manifest.is_file():
            print(f"ERROR: input manifest is missing: {manifest}", file=sys.stderr)
            return 2
        workspace.mkdir(parents=True, exist_ok=True)
        command = [sys.executable, str(helper), "--manifest", str(manifest), "--target-root", str(workspace)]
        if args.verify_only:
            command.append("--verify-only")
        print(" ".join(command))
        code = subprocess.run(command, check=False).returncode
        if code == 0 and not args.verify_only:
            materialize(package_root, config, workspace)
        return code

    if args.command == "preflight":
        selected = [args.module] if args.module else list(modules)
        unknown = [module_id for module_id in selected if module_id not in modules]
        if unknown:
            print(f"ERROR: unknown module id(s): {', '.join(unknown)}", file=sys.stderr)
            return 2
        failed = False
        for module_id in selected:
            module = modules[module_id]
            issues = preflight_one(package_root, config, module_id, module, workspace, args.rscript)
            if issues:
                failed = True
                print(f"{module_id}: BLOCKED")
                for issue in issues:
                    print(f"  - {issue}")
            else:
                print(f"{module_id}: ready (syntax/package checks only; no analysis executed)")
        print(f"Preflight scope workspace: {workspace}")
        return 1 if failed else 0

    if args.command == "run":
        if args.module not in modules:
            print(f"ERROR: unknown module id: {args.module}", file=sys.stderr)
            return 2
        if args.workers < 1:
            print("ERROR: --workers must be positive", file=sys.stderr)
            return 2
        module = modules[args.module]
        if module.get("status", "runnable") != "runnable":
            print(f"ERROR: {module.get('manual_reason', 'module has manual/external prerequisites')}", file=sys.stderr)
            return 2
        try:
            materialize(package_root, config, workspace)
            issues = preflight_one(package_root, config, args.module, module, workspace, args.rscript, use_workspace_script=True)
        except RuntimeError as exc:
            print(f"ERROR: {exc}", file=sys.stderr)
            return 2
        if issues:
            print(f"ERROR: preflight blocked {args.module}:", file=sys.stderr)
            for issue in issues:
                print(f"  - {issue}", file=sys.stderr)
            return 2
        command = format_command(module, workspace, args.rscript, args.workers)
        env = {**os.environ, **THREAD_ENV}
        print(f"Running {args.module}; outputs remain under {workspace}")
        return subprocess.run(command, cwd=workspace, env=env, check=False).returncode

    return 2


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except (RuntimeError, OSError, json.JSONDecodeError) as error:
        print(f"ERROR: {error}", file=sys.stderr)
        raise SystemExit(2)
