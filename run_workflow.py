#!/usr/bin/env python3
"""Semi-automatic Multiwfn -> VMD/Tachyon IGMH workflow.

The script keeps the scientific workflow explicit: Multiwfn generates the
numerical IGMH cube files, VMD/Tachyon renders the surfaces, and Python only
checks files, launches programs, converts images, and combines panels.
"""

from __future__ import annotations

import argparse
import csv
import shutil
import subprocess
import sys
from dataclasses import dataclass
from pathlib import Path


ROOT = Path(__file__).resolve().parent
DEFAULT_SYSTEMS = {
    "TS3a": "TS3a_SP_PCM.fchk",
    "TS3b": "TS3b_SP_PCM.fchk",
}
REQUIRED_CUBES = ("sl2r.cub", "dg_inter.cub")
TS3_PANEL_CROPS = {
    "TS3a": (0, 194, 2260, 2294),
    "TS3b": (44, 25, 2304, 2125),
}


@dataclass(frozen=True)
class FragmentConfig:
    fragment1: str
    fragment2: str
    center_selection: str | None = None


def posix_path(path: Path) -> str:
    return path.resolve().as_posix()


def find_executable(name_or_path: str) -> str | None:
    candidate = Path(name_or_path)
    if candidate.exists():
        return str(candidate)
    return shutil.which(name_or_path)


def multiwfn_input(fragment1: str, fragment2: str) -> str:
    """Default menu sequence for Multiwfn IGMH interfragment analysis.

    This sequence was verified with the TS3a/TS3b example workflow. Record the
    actual Multiwfn version used from your local Multiwfn banner or generated
    logs. Other Multiwfn versions/builds may require a different menu sequence.
    """
    lines = [
        "20",
        "11",
        "2",
        fragment1,
        fragment2,
        "11",
        fragment2,
        "3 A",
        "0.15",
        "3",
        "2",
        "6",
        "2",
        "y",
        "0",
        "0",
        "q",
    ]
    return "\n".join(lines) + "\n"


def read_template(path: Path, fragment1: str, fragment2: str) -> str:
    text = path.read_text(encoding="utf-8")
    return text.replace("{fragment1}", fragment1).replace("{fragment2}", fragment2)


def read_fragment_file(path: Path) -> dict[str, FragmentConfig]:
    """Read per-system fragment definitions from CSV.

    Required columns: system, fragment1, fragment2
    Optional column: center_selection
    """
    configs: dict[str, FragmentConfig] = {}
    with path.open(newline="", encoding="utf-8-sig") as handle:
        reader = csv.DictReader(handle)
        required = {"system", "fragment1", "fragment2"}
        if not reader.fieldnames or not required.issubset(reader.fieldnames):
            raise ValueError(f"{path} must contain columns: system, fragment1, fragment2")
        for row in reader:
            system = (row.get("system") or "").strip()
            if not system:
                continue
            center_selection = (row.get("center_selection") or "").strip() or None
            configs[system] = FragmentConfig(
                fragment1=(row["fragment1"] or "").strip(),
                fragment2=(row["fragment2"] or "").strip(),
                center_selection=center_selection,
            )
    return configs


def default_fragment_configs(systems: dict[str, str], fragment1: str, fragment2: str) -> dict[str, FragmentConfig]:
    return {system: FragmentConfig(fragment1=fragment1, fragment2=fragment2) for system in systems}


def fragment_to_vmd_indices(fragment: str) -> str:
    """Convert one-based Multiwfn atom ranges such as 1-72,80,85 to VMD indices."""
    tokens = fragment.replace("–", "-").replace("—", "-").replace(",", " ").split()
    parts: list[str] = []
    for token in tokens:
        if "-" in token:
            start_text, end_text = token.split("-", 1)
            start = int(start_text)
            end = int(end_text)
            if start <= 0 or end <= 0 or end < start:
                raise ValueError(f"Invalid atom range: {token}")
            parts.append(f"{start - 1} to {end - 1}")
        else:
            atom = int(token)
            if atom <= 0:
                raise ValueError(f"Invalid atom index: {token}")
            parts.append(str(atom - 1))
    return "index " + " ".join(parts)


def resolve_input_file(input_dir: Path, filename: str) -> Path:
    direct = Path(filename)
    candidates = []
    if direct.is_absolute():
        candidates.append(direct)
    else:
        candidates.extend([input_dir / filename, ROOT / filename])
    for candidate in candidates:
        if candidate.exists():
            return candidate.resolve()
    searched = ", ".join(str(p) for p in candidates)
    raise FileNotFoundError(f"Could not find input file {filename!r}. Searched: {searched}")


def ensure_cube_outputs(system: str, cubedir: Path) -> None:
    missing = [name for name in REQUIRED_CUBES if not (cubedir / name).exists()]
    if missing:
        raise FileNotFoundError(
            f"{system}: missing required Multiwfn output file(s) in {cubedir}: "
            + ", ".join(missing)
        )


def run_multiwfn(
    system: str,
    fchk: Path,
    cubedir: Path,
    executable: str,
    menu_text: str,
    overwrite: bool,
    dry_run: bool,
) -> None:
    existing = [cubedir / name for name in REQUIRED_CUBES if (cubedir / name).exists()]
    if existing and not overwrite:
        print(f"[{system}] Existing cube files found; skipping Multiwfn. Use --overwrite to rerun.")
        ensure_cube_outputs(system, cubedir)
        return

    cmd = [executable, str(fchk)]
    print(f"[{system}] Running Multiwfn in {cubedir}")
    print("  " + " ".join(cmd))
    if dry_run:
        return

    cubedir.mkdir(parents=True, exist_ok=True)
    inp = cubedir / f"{system}_Multiwfn_input.txt"
    log = cubedir / f"{system}_Multiwfn_output.log"
    inp.write_text(menu_text, encoding="utf-8")

    with log.open("w", encoding="utf-8", errors="replace") as handle:
        proc = subprocess.run(
            cmd,
            input=menu_text,
            text=True,
            cwd=cubedir,
            stdout=handle,
            stderr=subprocess.STDOUT,
            check=False,
        )
    if proc.returncode != 0:
        raise RuntimeError(f"{system}: Multiwfn exited with code {proc.returncode}. See {log}")
    ensure_cube_outputs(system, cubedir)


def run_vmd(
    system: str,
    cubedir: Path,
    out_tga: Path,
    vmd: str,
    center_selection: str | None,
    dry_run: bool,
) -> None:
    if not dry_run:
        out_tga.parent.mkdir(parents=True, exist_ok=True)
    script = ROOT / "scripts" / "render_IGMH.tcl"
    cmd = [
        vmd,
        "-dispdev",
        "text",
        "-e",
        str(script),
        "-args",
        system,
        posix_path(cubedir),
        posix_path(out_tga),
    ]
    if center_selection:
        cmd.append(center_selection)
    print(f"[{system}] Rendering with VMD/Tachyon")
    print("  " + " ".join(cmd))
    if dry_run:
        return
    proc = subprocess.run(cmd, cwd=ROOT, check=False)
    if proc.returncode != 0:
        raise RuntimeError(f"{system}: VMD exited with code {proc.returncode}")
    if not out_tga.exists():
        raise FileNotFoundError(f"{system}: VMD did not create {out_tga}")


def convert_image(src: Path, dst: Path, dry_run: bool) -> None:
    print(f"[image] {src.name} -> {dst.name}")
    if dry_run:
        return
    try:
        from PIL import Image
    except ImportError:
        magick = shutil.which("magick")
        if not magick:
            raise RuntimeError(
                "Pillow is not installed and ImageMagick 'magick' was not found. "
                "Install Pillow with 'pip install -r requirements.txt' or convert the TGA manually."
            )
        subprocess.run([magick, str(src), str(dst)], check=True)
        return
    with Image.open(src) as img:
        img.convert("RGB").save(dst)


def combine_panels(
    left_png: Path,
    right_png: Path,
    combined_png: Path,
    left_label: str,
    right_label: str,
    add_labels: bool,
    dry_run: bool,
) -> None:
    print(f"[image] Combining panels -> {combined_png.name}")
    if dry_run:
        return
    try:
        from PIL import Image, ImageDraw, ImageFont
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


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description="Run a semi-automatic Multiwfn IGMH analysis and VMD/Tachyon clean rendering workflow."
    )
    parser.add_argument("--input-dir", default="input", type=Path, help="Directory containing local .fchk files.")
    parser.add_argument("--output-dir", default="output", type=Path, help="Directory for Multiwfn cube/log output.")
    parser.add_argument("--figures-dir", default="figures", type=Path, help="Directory for rendered PNG figures.")
    parser.add_argument(
        "--systems",
        nargs="+",
        default=[f"{k}={v}" for k, v in DEFAULT_SYSTEMS.items()],
        help="Systems as LABEL=filename. Default: TS3a=TS3a_SP_PCM.fchk TS3b=TS3b_SP_PCM.fchk",
    )
    parser.add_argument("--fragment1", default="1-72", help="Fragment 1 atom range passed to Multiwfn.")
    parser.add_argument("--fragment2", default="73-94", help="Fragment 2 atom range passed to Multiwfn.")
    parser.add_argument(
        "--fragments-file",
        type=Path,
        help="CSV with per-system fragment definitions: system,fragment1,fragment2[,center_selection].",
    )
    parser.add_argument(
        "--auto-center-fragments",
        action="store_true",
        help="If center_selection is blank, center VMD rendering on fragment1+fragment2.",
    )
    parser.add_argument("--multiwfn", default="Multiwfn", help="Multiwfn executable name or path.")
    parser.add_argument("--vmd", default="vmd", help="VMD executable name or path.")
    parser.add_argument(
        "--multiwfn-input-template",
        type=Path,
        help="Optional Multiwfn menu input template. Use {fragment1} and {fragment2} placeholders.",
    )
    parser.add_argument("--skip-multiwfn", action="store_true", help="Use existing cube files in output/*_IGMH_files.")
    parser.add_argument("--skip-render", action="store_true", help="Skip VMD rendering and only check/generate cubes.")
    parser.add_argument("--overwrite", action="store_true", help="Rerun Multiwfn even if cube files already exist.")
    parser.add_argument("--no-panel-labels", action="store_true", help="Only write the no-text clean comparison; skip the TS3a/TS3b labeled comparison.")
    parser.add_argument("--dry-run", action="store_true", help="Print planned commands without running external programs.")
    args = parser.parse_args(argv)

    input_dir = (ROOT / args.input_dir).resolve() if not args.input_dir.is_absolute() else args.input_dir
    output_dir = (ROOT / args.output_dir).resolve() if not args.output_dir.is_absolute() else args.output_dir
    figures_dir = (ROOT / args.figures_dir).resolve() if not args.figures_dir.is_absolute() else args.figures_dir
    systems = parse_systems(args.systems)
    fragment_configs = default_fragment_configs(systems, args.fragment1, args.fragment2)
    if args.fragments_file:
        fragments_path = args.fragments_file if args.fragments_file.is_absolute() else ROOT / args.fragments_file
        file_configs = read_fragment_file(fragments_path)
        for system in systems:
            if system not in file_configs:
                raise ValueError(f"{fragments_path} has no fragment row for system {system!r}")
        fragment_configs = {system: file_configs[system] for system in systems}

    multiwfn = find_executable(args.multiwfn) or args.multiwfn
    vmd = find_executable(args.vmd) or args.vmd

    if not args.skip_multiwfn and not args.dry_run and not find_executable(args.multiwfn):
        raise FileNotFoundError("Multiwfn executable was not found. Add it to PATH or pass --multiwfn PATH.")
    if not args.skip_render and not args.dry_run and not find_executable(args.vmd):
        raise FileNotFoundError("VMD executable was not found. Add it to PATH or pass --vmd PATH.")

    if not args.dry_run:
        figures_dir.mkdir(parents=True, exist_ok=True)
    rendered_pngs: dict[str, Path] = {}

    for system, filename in systems.items():
        fragments = fragment_configs[system]
        if args.multiwfn_input_template:
            template_path = (
                args.multiwfn_input_template
                if args.multiwfn_input_template.is_absolute()
                else ROOT / args.multiwfn_input_template
            )
            menu_text = read_template(template_path, fragments.fragment1, fragments.fragment2)
        else:
            menu_text = multiwfn_input(fragments.fragment1, fragments.fragment2)
        center_selection = fragments.center_selection
        if args.auto_center_fragments and not center_selection:
            center_selection = fragment_to_vmd_indices(f"{fragments.fragment1},{fragments.fragment2}")

        print(f"[{system}] Fragment 1: {fragments.fragment1}; Fragment 2: {fragments.fragment2}")
        cubedir = output_dir / f"{system}_IGMH_files"
        if args.skip_multiwfn:
            print(f"[{system}] Skipping Multiwfn; checking existing cube files in {cubedir}")
            if not args.dry_run:
                ensure_cube_outputs(system, cubedir)
        else:
            fchk = resolve_input_file(input_dir, filename)
            run_multiwfn(system, fchk, cubedir, multiwfn, menu_text, args.overwrite, args.dry_run)

        if not args.skip_render:
            out_tga = figures_dir / f"{system}_IGMH_clean.tga"
            out_png = figures_dir / f"{system}_IGMH_clean.png"
            run_vmd(system, cubedir, out_tga, vmd, center_selection, args.dry_run)
            convert_image(out_tga, out_png, args.dry_run)
            rendered_pngs[system] = out_png

    if not args.skip_render and len(rendered_pngs) >= 2:
        labels = list(rendered_pngs)
        clean_combined = figures_dir / f"{labels[0]}_vs_{labels[1]}_IGMH_clean.png"
        combine_panels(
            rendered_pngs[labels[0]],
            rendered_pngs[labels[1]],
            clean_combined,
            labels[0],
            labels[1],
            add_labels=False,
            dry_run=args.dry_run,
        )
        if not args.no_panel_labels:
            labeled_combined = figures_dir / f"{labels[0]}_vs_{labels[1]}_IGMH_labeled.png"
            combine_panels(
                rendered_pngs[labels[0]],
                rendered_pngs[labels[1]],
                labeled_combined,
                labels[0],
                labels[1],
                add_labels=True,
                dry_run=args.dry_run,
            )

    print("[done] Workflow completed.")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception as exc:
        print(f"[error] {exc}", file=sys.stderr)
        raise SystemExit(1)
