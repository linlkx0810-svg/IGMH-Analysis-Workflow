"""Workflow status checks for Multiwfn/VMD IGMH runs."""

from __future__ import annotations

from dataclasses import dataclass, asdict
from pathlib import Path
import json
import re
from typing import Any

REQUIRED_CUBES = ("sl2r.cub", "dg_inter.cub")


@dataclass(frozen=True)
class StageStatus:
    name: str
    status: str
    message: str = ""

    @property
    def ok(self) -> bool:
        return self.status == "PASS"


def pass_status(name: str, message: str = "") -> StageStatus:
    return StageStatus(name, "PASS", message)


def fail_status(name: str, message: str = "") -> StageStatus:
    return StageStatus(name, "FAIL", message)


def warn_status(name: str, message: str = "") -> StageStatus:
    return StageStatus(name, "WARN", message)


def validate_cube_outputs(cubedir: Path, required: tuple[str, ...] = REQUIRED_CUBES) -> list[StageStatus]:
    statuses = []
    for name in required:
        path = cubedir / name
        if not path.exists():
            statuses.append(fail_status("IGMH_CUBES", f"missing {name}"))
        elif path.stat().st_size == 0:
            statuses.append(fail_status("IGMH_CUBES", f"empty {name}"))
        else:
            statuses.append(pass_status("IGMH_CUBES", f"{name} present"))
    return statuses


def pair_integration_status(log_path: Path) -> StageStatus:
    if not log_path.exists():
        return warn_status("PAIR_INTEGRATION", f"log not found: {log_path}")
    text = log_path.read_text(encoding="utf-8", errors="replace")
    percentages = [float(match.group(1)) for match in re.finditer(r"(\d+(?:\.\d+)?)\s*%", text)]
    if any(value >= 99.999 for value in percentages):
        return pass_status("PAIR_INTEGRATION", "integration reached 100%")
    if percentages:
        return fail_status("PAIR_INTEGRATION", f"highest reported progress was {max(percentages):.1f}%")
    return warn_status("PAIR_INTEGRATION", "no percentage progress found in log")


def _pair_output_status_name(filename: str) -> str:
    if filename.lower() == "atmdg.txt":
        return "ATMDG_OUTPUT"
    if filename.lower() == "ibsiw.txt":
        return "IBSIW_OUTPUT"
    return filename.upper().replace(".TXT", "_OUTPUT")


def validate_text_pair_output(directory: Path, filename: str, parser) -> StageStatus:
    name = _pair_output_status_name(filename)
    path = directory / filename
    if not path.exists():
        return fail_status(name, f"missing {filename}")
    if path.stat().st_size == 0:
        return fail_status(name, f"empty {filename}")
    try:
        parsed = parser(path)
    except Exception as exc:  # pragma: no cover
        return fail_status(name, f"could not parse {filename}: {exc}")
    if not parsed:
        return fail_status(name, f"no pair rows parsed from {filename}")
    return pass_status(name, f"parsed {len(parsed)} pair rows")


def validate_pair_outputs(directory: Path, log_path: Path | None, parser) -> list[StageStatus]:
    statuses = []
    if log_path:
        statuses.append(pair_integration_status(log_path))
    statuses.append(validate_text_pair_output(directory, "atmdg.txt", parser))
    statuses.append(validate_text_pair_output(directory, "IBSIW.txt", parser))
    return statuses


def inspect_png_framing(path: Path, expected_size: tuple[int, int] | None = None, white_threshold: int = 248) -> tuple[list[StageStatus], dict[str, Any]]:
    statuses: list[StageStatus] = []
    summary: dict[str, Any] = {"path": str(path)}
    if not path.exists():
        return [fail_status("RENDERING", f"missing {path.name}")], summary
    if path.stat().st_size == 0:
        return [fail_status("RENDERING", f"empty {path.name}")], summary
    try:
        from PIL import Image
    except ImportError:
        return [warn_status("RENDERING", "Pillow unavailable; image content not inspected")], summary

    with Image.open(path) as img:
        rgb = img.convert("RGB")
        summary["resolution"] = [rgb.width, rgb.height]
        if expected_size and rgb.size != expected_size:
            statuses.append(fail_status("RENDER_RESOLUTION", f"expected {expected_size}, got {rgb.size}"))
        else:
            statuses.append(pass_status("RENDER_RESOLUTION", f"{rgb.width}x{rgb.height}"))

        pixels = rgb.load()
        xs: list[int] = []
        ys: list[int] = []
        for y in range(rgb.height):
            for x in range(rgb.width):
                r, g, b = pixels[x, y]
                if r < white_threshold or g < white_threshold or b < white_threshold:
                    xs.append(x)
                    ys.append(y)
        if not xs:
            statuses.append(fail_status("RENDER_CONTENT", "image appears blank against white background"))
            return statuses, summary
        bbox = [min(xs), min(ys), max(xs), max(ys)]
        margins_px = [bbox[0], bbox[1], rgb.width - 1 - bbox[2], rgb.height - 1 - bbox[3]]
        margins_percent = [
            100.0 * margins_px[0] / rgb.width,
            100.0 * margins_px[1] / rgb.height,
            100.0 * margins_px[2] / rgb.width,
            100.0 * margins_px[3] / rgb.height,
        ]
        min_margin = min(margins_percent)
        summary.update({"bbox": bbox, "margins_px": margins_px, "margins_percent": margins_percent, "min_margin_percent": min_margin})
        message = "outer margins %: " + ", ".join(f"{value:.2f}" for value in margins_percent)
        if min_margin < 1.0:
            statuses.append(fail_status("RENDER_FRAMING", message))
        elif min_margin < 3.0:
            statuses.append(warn_status("RENDER_FRAMING", message))
        else:
            statuses.append(pass_status("RENDER_FRAMING", message))
    return statuses, summary


def validate_png(path: Path, expected_size: tuple[int, int] | None = None, margin_px: int = 1) -> list[StageStatus]:
    statuses, _summary = inspect_png_framing(path, expected_size=expected_size)
    return statuses


def final_status(statuses: list[StageStatus], strict_warnings: bool = False) -> str:
    if any(status.status == "FAIL" for status in statuses):
        return "FAIL"
    if strict_warnings and any(status.status == "WARN" for status in statuses):
        return "FAIL"
    return "PASS"


def write_manifest(path: Path, payload: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2, ensure_ascii=False), encoding="utf-8")


def status_dict(status: StageStatus) -> dict[str, str]:
    return asdict(status)
