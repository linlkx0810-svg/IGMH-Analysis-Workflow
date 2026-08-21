#!/usr/bin/env python3
"""Semi-automatic Multiwfn -> VMD/Tachyon IGMH workflow."""

from __future__ import annotations

import argparse
import hashlib
import tempfile
import shutil
import subprocess
import sys
from pathlib import Path
from typing import Any

from igmh_fragments import (
    FragmentConfig,
    fragment_to_vmd_indices,
    read_fchk_atom_count,
    read_fragment_file,
    validate_fragments,
)
from igmh_comparison import build_comparison_plan
from igmh_pair_parser import parse_pair_value_file, read_pair_outputs, write_pair_csv
from igmh_validation import (
    StageStatus,
    final_status,
    status_dict,
    inspect_png_framing,
    validate_cube_outputs,
    validate_pair_outputs,
    validate_png,
    write_manifest,
)

ROOT = Path(__file__).resolve().parent
DEFAULT_SYSTEMS = {"TS3a": "TS3a_SP_PCM.fchk", "TS3b": "TS3b_SP_PCM.fchk"}
REQUIRED_CUBES = ("sl2r.cub", "dg_inter.cub")
TS3_PANEL_CROPS = {"TS3a": (0, 194, 2260, 2294), "TS3b": (44, 25, 2304, 2125)}


def posix_path(path: Path) -> str:
    return path.resolve().as_posix()

def is_ascii_path(path: Path) -> bool:
    try:
        str(path).encode("ascii")
    except UnicodeEncodeError:
        return False
    return True


def vmd_stage_root() -> Path:
    root = Path(tempfile.gettempdir()) / "igmh_vmd_ascii_stage"
    root.mkdir(parents=True, exist_ok=True)
    return root


def stage_directory_for_vmd(path: Path) -> Path:
    resolved = path.resolve()
    if is_ascii_path(resolved):
        return resolved
    digest = hashlib.sha1(str(resolved).encode("utf-8")).hexdigest()[:12]
    link = vmd_stage_root() / f"dir_{digest}"
    if not link.exists():
        proc = subprocess.run(["cmd", "/c", "mklink", "/J", str(link), str(resolved)], text=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, check=False)
        if proc.returncode != 0:
            raise RuntimeError(f"Could not create ASCII VMD junction for {resolved}: {proc.stdout}")
    return link


def stage_file_for_vmd(path: Path) -> Path:
    resolved = path.resolve()
    if is_ascii_path(resolved):
        return resolved
    digest = hashlib.sha1(str(resolved).encode("utf-8")).hexdigest()[:12]
    staged = vmd_stage_root() / f"file_{digest}_{resolved.name}"
    shutil.copyfile(resolved, staged)
    return staged


def stage_output_for_vmd(path: Path) -> Path:
    resolved = path.resolve()
    if is_ascii_path(resolved):
        return resolved
    digest = hashlib.sha1(str(resolved).encode("utf-8")).hexdigest()[:12]
    staged = vmd_stage_root() / f"out_{digest}_{resolved.name}"
    if staged.exists():
        staged.unlink()
    return staged


def collect_vmd_output(staged: Path, final: Path) -> None:
    if staged.resolve() == final.resolve():
        return
    final.parent.mkdir(parents=True, exist_ok=True)
    shutil.copyfile(staged, final)


def find_executable(name_or_path: str) -> str | None:
    candidate = Path(name_or_path)
    if candidate.exists():
        return str(candidate)
    return shutil.which(name_or_path)


def multiwfn_input(fragment1: str, fragment2: str) -> str:
    """Default menu sequence for Multiwfn IGMH interfragment analysis.

    This sequence was audited against the TS3b example output log generated
    by Multiwfn 2026.1.12. Other Multiwfn versions/builds may require a
    different menu sequence.
    """
    lines = [
        "20", "11", "2", fragment1, fragment2, "11", fragment2, "3 A",
        "0.15", "3", "2", "6", "2", "y", "0", "0", "q",
    ]
    return "\n".join(lines) + "\n"


def read_template(path: Path, fragment1: str, fragment2: str) -> str:
    text = path.read_text(encoding="utf-8")
    return text.replace("{fragment1}", fragment1).replace("{fragment2}", fragment2)


def default_fragment_configs(systems: dict[str, str], fragment1: str, fragment2: str) -> dict[str, FragmentConfig]:
    return {system: FragmentConfig(fragment1=fragment1, fragment2=fragment2) for system in systems}


def resolve_input_file(input_dir: Path, filename: str) -> Path:
    direct = Path(filename)
    candidates = [direct] if direct.is_absolute() else [input_dir / filename, ROOT / filename]
    for candidate in candidates:
        if candidate.exists():
            return candidate.resolve()
    searched = ", ".join(str(p) for p in candidates)
    raise FileNotFoundError(f"Could not find input file {filename!r}. Searched: {searched}")


def ensure_cube_outputs(system: str, cubedir: Path) -> None:
    failures = [status for status in validate_cube_outputs(cubedir, REQUIRED_CUBES) if status.status == "FAIL"]
    if failures:
        details = "; ".join(status.message for status in failures)
        raise FileNotFoundError(f"{system}: required Multiwfn cube output validation failed: {details}")


def print_status(status: StageStatus) -> None:
    suffix = f" - {status.message}" if status.message else ""
    print(f"  {status.name:<18} {status.status}{suffix}")


def run_multiwfn(system: str, fchk: Path, cubedir: Path, executable: str, menu_text: str, overwrite: bool, dry_run: bool) -> int | None:
    existing = [cubedir / name for name in REQUIRED_CUBES if (cubedir / name).exists()]
    if existing and not overwrite:
        print(f"[{system}] Existing cube files found; skipping Multiwfn. Use --overwrite to rerun.")
        ensure_cube_outputs(system, cubedir)
        return None

    cmd = [executable, str(fchk)]
    print(f"[{system}] Running Multiwfn in {cubedir}")
    print("  " + " ".join(cmd))
    if dry_run:
        return None

    cubedir.mkdir(parents=True, exist_ok=True)
    inp = cubedir / f"{system}_Multiwfn_input.txt"
    log = cubedir / f"{system}_Multiwfn_output.log"
    inp.write_text(menu_text, encoding="utf-8")
    with log.open("w", encoding="utf-8", errors="replace") as handle:
        proc = subprocess.run(cmd, input=menu_text, text=True, cwd=cubedir, stdout=handle, stderr=subprocess.STDOUT, check=False)
    if proc.returncode != 0:
        raise RuntimeError(f"{system}: Multiwfn exited with code {proc.returncode}. See {log}")
    ensure_cube_outputs(system, cubedir)
    return proc.returncode


def run_vmd(system: str, cubedir: Path, out_tga: Path, vmd: str, center_selection: str | None, dry_run: bool) -> None:
    if not dry_run:
        out_tga.parent.mkdir(parents=True, exist_ok=True)
    script = stage_file_for_vmd(ROOT / "scripts" / "render_IGMH.tcl")
    vmd_cubedir = stage_directory_for_vmd(cubedir)
    vmd_out_tga = stage_output_for_vmd(out_tga)
    cmd = [vmd, "-dispdev", "text", "-e", str(script), "-args", system, posix_path(vmd_cubedir), posix_path(vmd_out_tga)]
    if center_selection:
        cmd.append(center_selection)
    print(f"[{system}] Rendering with VMD/Tachyon")
    print("  " + " ".join(cmd))
    if dry_run:
        return
    proc = subprocess.run(cmd, cwd=ROOT, check=False)
    if not vmd_out_tga.exists():
        raise FileNotFoundError(f"{system}: VMD did not create {vmd_out_tga}")
    if proc.returncode != 0:
        print(f"[{system}] Warning: VMD returned code {proc.returncode}, but the requested render output exists.")
    collect_vmd_output(vmd_out_tga, out_tga)



def run_vmd_comparison(system: str, cubedir: Path, out_tga: Path, vmd: str, transform_tcl: Path, dry_run: bool) -> None:
    if not dry_run:
        out_tga.parent.mkdir(parents=True, exist_ok=True)
    script = stage_file_for_vmd(ROOT / "scripts" / "render_IGMH_comparison.tcl")
    staged_transform = stage_file_for_vmd(transform_tcl)
    vmd_cubedir = stage_directory_for_vmd(cubedir)
    vmd_out_tga = stage_output_for_vmd(out_tga)
    cmd = [vmd, "-dispdev", "text", "-e", str(script), "-args", system, posix_path(vmd_cubedir), posix_path(vmd_out_tga), posix_path(staged_transform)]
    print(f"[{system}] Rendering aligned comparison with VMD/Tachyon")
    print("  " + " ".join(cmd))
    if dry_run:
        return
    proc = subprocess.run(cmd, cwd=ROOT, check=False)
    if not vmd_out_tga.exists():
        raise FileNotFoundError(f"{system}: VMD did not create {vmd_out_tga}")
    if proc.returncode != 0:
        print(f"[{system}] Warning: VMD returned code {proc.returncode}, but the requested render output exists.")
    collect_vmd_output(vmd_out_tga, out_tga)
def convert_image(src: Path, dst: Path, dry_run: bool) -> None:
    print(f"[image] {src.name} -> {dst.name}")
    if dry_run:
        return
    try:
        from PIL import Image
    except ImportError:
        magick = shutil.which("magick")
        if not magick:
            raise RuntimeError("Pillow is not installed and ImageMagick 'magick' was not found. Install Pillow with 'pip install -r requirements.txt'.")
        subprocess.run([magick, str(src), str(dst)], check=True)
        return
    with Image.open(src) as img:
        img.convert("RGB").save(dst)


def combine_panels(left_png: Path, right_png: Path, combined_png: Path, left_label: str, right_label: str, add_labels: bool, dry_run: bool) -> None:
    print(f"[image] Combining panels -> {combined_png.name}")
    if dry_run:
        return
    try:
        from PIL import Image, ImageDraw
    except ImportError as exc:
        raise RuntimeError("Pillow is required to combine panels. Install with 'pip install -r requirements.txt'.") from exc

    left = Image.open(left_png).convert("RGB")
    right = Image.open(right_png).convert("RGB")
    if left_label in TS3_PANEL_CROPS and right_label in TS3_PANEL_CROPS:
        left = left.crop(TS3_PANEL_CROPS[left_label])
        right = right.crop(TS3_PANEL_CROPS[right_label])
        panel_w, panel_h = left.size
        gap = 50
        title_h = 92 if add_labels else 0
        title_font_size = 48
    else:
        size = min(left.width, left.height, right.width, right.height)
        left = left.crop(((left.width - size) // 2, (left.height - size) // 2, (left.width + size) // 2, (left.height + size) // 2))
        right = right.crop(((right.width - size) // 2, (right.height - size) // 2, (right.width + size) // 2, (right.height + size) // 2))
        panel_w, panel_h = left.size
        gap = max(24, size // 45)
        title_h = max(0, size // 15 if add_labels else 0)
        title_font_size = max(28, size // 38)

    canvas = Image.new("RGB", (panel_w * 2 + gap, panel_h + title_h), "white")
    canvas.paste(left, (0, title_h))
    canvas.paste(right, (panel_w + gap, title_h))
    if add_labels:
        draw = ImageDraw.Draw(canvas)
        font = load_font(title_font_size)
        draw_centered(draw, left_label, panel_w // 2, max(10, title_h // 4), font)
        draw_centered(draw, right_label, panel_w + gap + panel_w // 2, max(10, title_h // 4), font)
    canvas.save(combined_png)


def load_font(size: int):
    from PIL import ImageFont
    for candidate in (r"C:\Windows\Fonts\arial.ttf", r"C:\Windows\Fonts\calibri.ttf"):
        path = Path(candidate)
        if path.exists():
            return ImageFont.truetype(str(path), size)
    return ImageFont.load_default()


def draw_centered(draw, text: str, cx: int, y: int, font) -> None:
    bbox = draw.textbbox((0, 0), text, font=font)
    draw.text((cx - (bbox[2] - bbox[0]) / 2, y), text, fill=(0, 0, 0), font=font)


def parse_systems(items: list[str]) -> dict[str, str]:
    systems = {}
    for item in items:
        if "=" not in item:
            raise ValueError(f"System entry must be LABEL=filename, got {item!r}")
        label, filename = item.split("=", 1)
        systems[label.strip()] = filename.strip()
    return systems


def ordered_labels(rendered_pngs: dict[str, Path], panel_order: list[str] | None) -> list[str]:
    labels = panel_order or list(rendered_pngs)
    missing = [label for label in labels if label not in rendered_pngs]
    if missing:
        raise ValueError("Panel order references missing rendered image(s): " + ", ".join(missing))
    if len(labels) < 2:
        raise ValueError("At least two systems are required to combine panels")
    return labels


def existing_rendered_pngs(figures_dir: Path, labels: list[str]) -> dict[str, Path]:
    rendered: dict[str, Path] = {}
    for label in labels:
        path = figures_dir / f"{label}_IGMH_clean.png"
        if path.exists():
            rendered[label] = path
    return rendered


def manifest_system_entry(system: str, fchk: Path | None, fragments: FragmentConfig) -> dict[str, Any]:
    return {
        "system": system,
        "input_file": "" if fchk is None else str(fchk),
        "atom_count": None,
        "fragment1": fragments.fragment1,
        "fragment2": fragments.fragment2,
        "alignment_used": False,
        "alignment_RMSD": None,
        "camera_shared": False,
        "scale_shared": False,
        "statuses": [],
        "warnings": [],
        "final_status": "FAIL",
    }


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description="Run a semi-automatic Multiwfn IGMH analysis and VMD/Tachyon clean rendering workflow."
    )
    parser.add_argument("--input-dir", default="input", type=Path, help="Directory containing local .fchk files.")
    parser.add_argument("--output-dir", default="output", type=Path, help="Directory for Multiwfn cube/log output.")
    parser.add_argument("--figures-dir", default="figures", type=Path, help="Directory for rendered PNG figures.")
    parser.add_argument("--systems", nargs="+", default=[f"{k}={v}" for k, v in DEFAULT_SYSTEMS.items()], help="Systems as LABEL=filename. Default: TS3a=TS3a_SP_PCM.fchk TS3b=TS3b_SP_PCM.fchk")
    parser.add_argument("--fragment1", default="1-72", help="Fragment 1 atom range passed to Multiwfn.")
    parser.add_argument("--fragment2", default="73-94", help="Fragment 2 atom range passed to Multiwfn.")
    parser.add_argument("--fragments-file", type=Path, help="CSV with per-system fragment definitions: system,fragment1,fragment2[,center_selection].")
    parser.add_argument("--auto-center-fragments", action="store_true", help="If center_selection is blank, center VMD rendering on fragment1+fragment2.")
    parser.add_argument("--multiwfn", default="Multiwfn", help="Multiwfn executable name or path.")
    parser.add_argument("--vmd", default="vmd", help="VMD executable name or path.")
    parser.add_argument("--multiwfn-input-template", type=Path, help="Optional Multiwfn menu input template. Use {fragment1} and {fragment2} placeholders.")
    parser.add_argument("--skip-multiwfn", action="store_true", help="Use existing cube files in output/*_IGMH_files.")
    parser.add_argument("--skip-render", action="store_true", help="Skip VMD rendering and only check/generate cubes.")
    parser.add_argument("--overwrite", action="store_true", help="Rerun Multiwfn even if cube files already exist.")
    parser.add_argument("--no-panel-labels", action="store_true", help="Only write the no-text clean comparison; skip the labeled comparison.")
    parser.add_argument("--validate-only", action="store_true", help="Validate inputs/fragments/existing outputs and do not run Multiwfn or VMD.")
    parser.add_argument("--pair-analysis", action="store_true", help="Require completed atom-pair delta-g/IBSIW outputs and export parsed CSV files.")
    parser.add_argument("--panel-order", nargs="+", help="Explicit panel order for combined figures, e.g. --panel-order TS5R TS5S.")
    parser.add_argument("--comparison-config", type=Path, help="JSON config for Kabsch-aligned shared-camera comparison rendering.")
    parser.add_argument("--strict", action="store_true", help="Return non-zero when any requested validation stage fails.")
    parser.add_argument("--manifest", default="output/run_manifest.json", type=Path, help="Path for machine-readable run manifest JSON.")
    parser.add_argument("--dry-run", action="store_true", help="Print planned commands without running external programs.")
    args = parser.parse_args(argv)

    input_dir = (ROOT / args.input_dir).resolve() if not args.input_dir.is_absolute() else args.input_dir
    output_dir = (ROOT / args.output_dir).resolve() if not args.output_dir.is_absolute() else args.output_dir
    figures_dir = (ROOT / args.figures_dir).resolve() if not args.figures_dir.is_absolute() else args.figures_dir
    manifest_path = (ROOT / args.manifest).resolve() if not args.manifest.is_absolute() else args.manifest
    systems = parse_systems(args.systems)

    combine_only = args.skip_render and args.panel_order and not args.validate_only and not args.pair_analysis
    if combine_only:
        if not args.dry_run:
            figures_dir.mkdir(parents=True, exist_ok=True)
        rendered_pngs = existing_rendered_pngs(figures_dir, args.panel_order)
        labels = ordered_labels(rendered_pngs, args.panel_order)
        left_label, right_label = labels[0], labels[1]
        clean_combined = figures_dir / f"{left_label}_vs_{right_label}_IGMH_clean.png"
        combine_panels(rendered_pngs[left_label], rendered_pngs[right_label], clean_combined, left_label, right_label, add_labels=False, dry_run=args.dry_run)
        manifest = {
            "systems": [],
            "comparison": {
                "panel_order": [left_label, right_label],
                "clean_figure": str(clean_combined),
                "camera_shared": False,
                "scale_shared": False,
            },
            "final_status": "PASS",
        }
        if not args.no_panel_labels:
            labeled_combined = figures_dir / f"{left_label}_vs_{right_label}_IGMH_labeled.png"
            combine_panels(rendered_pngs[left_label], rendered_pngs[right_label], labeled_combined, left_label, right_label, add_labels=True, dry_run=args.dry_run)
            manifest["comparison"]["labeled_figure"] = str(labeled_combined)
        if not args.dry_run:
            write_manifest(manifest_path, manifest)
            print(f"[manifest] {manifest_path}")
        print("[done] Panel assembly completed with FINAL_STATUS=PASS.")
        return 0

    fragment_configs = default_fragment_configs(systems, args.fragment1, args.fragment2)
    if args.fragments_file:
        fragments_path = args.fragments_file if args.fragments_file.is_absolute() else ROOT / args.fragments_file
        file_configs = read_fragment_file(fragments_path)
        for system in systems:
            if system not in file_configs:
                raise ValueError(f"{fragments_path} has no fragment row for system {system!r}")
        fragment_configs = {system: file_configs[system] for system in systems}

    if not args.validate_only:
        multiwfn = find_executable(args.multiwfn) or args.multiwfn
        vmd = find_executable(args.vmd) or args.vmd
        if not args.skip_multiwfn and not args.dry_run and not find_executable(args.multiwfn):
            raise FileNotFoundError("Multiwfn executable was not found. Add it to PATH or pass --multiwfn PATH.")
        if not args.skip_render and not args.dry_run and not find_executable(args.vmd):
            raise FileNotFoundError("VMD executable was not found. Add it to PATH or pass --vmd PATH.")
    else:
        multiwfn = args.multiwfn
        vmd = args.vmd

    if not args.dry_run:
        figures_dir.mkdir(parents=True, exist_ok=True)
    rendered_pngs: dict[str, Path] = {}
    cube_dirs: dict[str, Path] = {}
    comparison_statuses: list[StageStatus] = []
    comparison_framing: dict[str, Any] = {}
    manifest: dict[str, Any] = {"systems": [], "comparison": {}, "final_status": "FAIL"}
    all_statuses: list[StageStatus] = []

    for system, filename in systems.items():
        fragments = fragment_configs[system]
        fchk: Path | None = None
        entry = manifest_system_entry(system, None, fragments)
        statuses: list[StageStatus] = []
        print(f"[{system}] Fragment 1: {fragments.fragment1}; Fragment 2: {fragments.fragment2}")
        try:
            fchk = resolve_input_file(input_dir, filename)
            entry["input_file"] = str(fchk)
            atom_count = read_fchk_atom_count(fchk)
            entry["atom_count"] = atom_count
            validation = validate_fragments(fragments.fragment1, fragments.fragment2, atom_count)
            statuses.append(StageStatus("FRAGMENTS", "PASS", f"{len(validation.uncovered_atoms)} uncovered atom(s)"))
            entry["warnings"].extend(validation.warnings)
        except Exception as exc:
            statuses.append(StageStatus("FRAGMENTS", "FAIL", str(exc)))

        cubedir = output_dir / f"{system}_IGMH_files"
        cube_dirs[system] = cubedir
        if args.validate_only:
            statuses.extend(validate_cube_outputs(cubedir, REQUIRED_CUBES))
        elif not any(status.status == "FAIL" for status in statuses):
            if args.multiwfn_input_template:
                template_path = args.multiwfn_input_template if args.multiwfn_input_template.is_absolute() else ROOT / args.multiwfn_input_template
                menu_text = read_template(template_path, fragments.fragment1, fragments.fragment2)
            else:
                menu_text = multiwfn_input(fragments.fragment1, fragments.fragment2)
            center_selection = fragments.center_selection
            if args.auto_center_fragments and not center_selection:
                center_selection = fragment_to_vmd_indices(f"{fragments.fragment1},{fragments.fragment2}")

            if args.skip_multiwfn:
                print(f"[{system}] Skipping Multiwfn; checking existing cube files in {cubedir}")
                if not args.dry_run:
                    statuses.extend(validate_cube_outputs(cubedir, REQUIRED_CUBES))
                    ensure_cube_outputs(system, cubedir)
            else:
                assert fchk is not None
                run_multiwfn(system, fchk, cubedir, multiwfn, menu_text, args.overwrite, args.dry_run)
                if not args.dry_run:
                    statuses.extend(validate_cube_outputs(cubedir, REQUIRED_CUBES))

            if not args.skip_render and not args.comparison_config:
                out_tga = figures_dir / f"{system}_IGMH_clean.tga"
                out_png = figures_dir / f"{system}_IGMH_clean.png"
                run_vmd(system, cubedir, out_tga, vmd, center_selection, args.dry_run)
                convert_image(out_tga, out_png, args.dry_run)
                rendered_pngs[system] = out_png
                if not args.dry_run:
                    statuses.extend(validate_png(out_png, expected_size=(2400, 2400)))

        if args.pair_analysis and not args.dry_run:
            log_path = cubedir / f"{system}_Multiwfn_output.log"
            pair_statuses = validate_pair_outputs(cubedir, log_path, parse_pair_value_file)
            statuses.extend(pair_statuses)
            if fchk is not None and all(status.status == "PASS" for status in pair_statuses):
                records = read_pair_outputs(cubedir, fchk)
                csv_path = cubedir / f"{system}_atom_pair_IGMH.csv"
                write_pair_csv(csv_path, system, records)
                statuses.append(StageStatus("PAIR_CSV", "PASS", str(csv_path)))

        for status in statuses:
            print_status(status)
        entry["statuses"] = [status_dict(status) for status in statuses]
        entry["final_status"] = final_status(statuses)
        manifest["systems"].append(entry)
        all_statuses.extend(statuses)

    if args.comparison_config and not args.skip_render:
        config_path = args.comparison_config if args.comparison_config.is_absolute() else ROOT / args.comparison_config
        plan = build_comparison_plan(config_path, cube_dirs, output_dir)
        manifest["comparison"] = {
            "alignment_used": True,
            "alignment_atom_mapping": plan.alignment_atom_mapping,
            "alignment_RMSD": plan.alignment_rmsd_A,
            "alignment_RMSD_units": "angstrom",
            "volume_registration": plan.volume_registration,
            "camera_rotation_shared": True,
            "camera_shared": True,
            "scale_shared": True,
            "screen_translation_shared": False,
            "projection_shared": plan.projection == "Orthographic",
            "projection": plan.projection,
            "scale": plan.scale,
            "screen_translate": plan.screen_translate,
            "per_panel_screen_translation": plan.per_panel_screen_translate,
            "panel_order": plan.panel_order,
            "framing_status": {},
            "final_status": "FAIL",
        }
        for label in plan.panel_order:
            cubedir = cube_dirs[label]
            out_tga = figures_dir / f"{label}_IGMH_clean.tga"
            out_png = figures_dir / f"{label}_IGMH_clean.png"
            run_vmd_comparison(label, cubedir, out_tga, vmd, plan.transform_tcl, args.dry_run)
            convert_image(out_tga, out_png, args.dry_run)
            rendered_pngs[label] = out_png
            if not args.dry_run:
                render_statuses, framing = inspect_png_framing(out_png, expected_size=(2400, 2400))
                comparison_statuses.extend(render_statuses)
                comparison_framing[label] = {
                    "statuses": [status_dict(status) for status in render_statuses],
                    "framing": framing,
                }
                for status in render_statuses:
                    print_status(status)
        manifest["comparison"]["framing_status"] = comparison_framing
        comparison_statuses.extend([
            StageStatus("ALIGNMENT", "PASS", ", ".join(f"{k}={v:.6f} A" for k, v in plan.alignment_rmsd_A.items())),
            StageStatus("VOLUME_REGISTRATION", "PASS", plan.volume_registration),
            StageStatus("SHARED_CAMERA", "PASS", "identical rotate_matrix used for all panels"),
            StageStatus("SHARED_SCALE", "PASS", f"scale={plan.scale:.10f}"),
        ])
        manifest["comparison"]["final_status"] = final_status(comparison_statuses)
        all_statuses.extend(comparison_statuses)
    if args.skip_render and args.panel_order:
        rendered_pngs.update(existing_rendered_pngs(figures_dir, args.panel_order))

    should_combine = (not args.skip_render and len(rendered_pngs) >= 2) or (args.skip_render and args.panel_order and len(rendered_pngs) >= 2)
    if should_combine:
        labels = ordered_labels(rendered_pngs, args.panel_order)
        left_label, right_label = labels[0], labels[1]
        clean_combined = figures_dir / f"{left_label}_vs_{right_label}_IGMH_clean.png"
        combine_panels(rendered_pngs[left_label], rendered_pngs[right_label], clean_combined, left_label, right_label, add_labels=False, dry_run=args.dry_run)
        comparison = manifest.setdefault("comparison", {})
        comparison.update({
            "panel_order": [left_label, right_label],
            "clean_figure": str(clean_combined),
        })
        comparison.setdefault("camera_shared", False)
        comparison.setdefault("scale_shared", False)
        if not args.no_panel_labels:
            labeled_combined = figures_dir / f"{left_label}_vs_{right_label}_IGMH_labeled.png"
            combine_panels(rendered_pngs[left_label], rendered_pngs[right_label], labeled_combined, left_label, right_label, add_labels=True, dry_run=args.dry_run)
            comparison["labeled_figure"] = str(labeled_combined)

    manifest["final_status"] = final_status(all_statuses)
    if not args.dry_run:
        write_manifest(manifest_path, manifest)
        print(f"[manifest] {manifest_path}")

    print(f"[done] Workflow completed with FINAL_STATUS={manifest['final_status']}.")
    if args.strict and manifest["final_status"] != "PASS":
        return 1
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception as exc:
        print(f"[error] {exc}", file=sys.stderr)
        raise SystemExit(1)
