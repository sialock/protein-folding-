"""Export simple vector charts without requiring a plotting GUI/backend."""
from __future__ import annotations

from html import escape
from pathlib import Path


def write_trace_svg(path: str | Path, trace: list[dict], metric: str, title: str, unit: str) -> None:
    if not trace:
        raise ValueError("Trace must not be empty")
    values = [float(row[metric]) for row in trace]
    steps = [int(row["step"]) for row in trace]
    low, high = min(values), max(values)
    if high <= low:
        high = low + 1
    left, top, width, height = 74, 65, 680, 300
    max_step = max(steps) or 1
    coords = [
        f"{left + step / max_step * width:.1f},{top + (high - value) / (high - low) * height:.1f}"
        for step, value in zip(steps, values)
    ]
    labels = "".join(
        f'<text x="{left - 9}" y="{top + height - i * height / 4 + 4:.1f}" '
        f'font-size="12" text-anchor="end" fill="#455569">'
        f'{low + (high-low)*i/4:.2f}</text>'
        for i in range(5)
    )
    svg = f'''<svg xmlns="http://www.w3.org/2000/svg" width="820" height="440" viewBox="0 0 820 440">
<rect width="820" height="440" rx="16" fill="#f9fbff"/>
<text x="74" y="36" fill="#14253b" font-family="Arial,sans-serif" font-size="22" font-weight="bold">{escape(title)}</text>
<line x1="{left}" y1="{top}" x2="{left}" y2="{top + height}" stroke="#9da9b8"/>
<line x1="{left}" y1="{top + height}" x2="{left + width}" y2="{top + height}" stroke="#9da9b8"/>
{labels}
<polyline points="{' '.join(coords)}" fill="none" stroke="#1762c4" stroke-width="2.5" stroke-linejoin="round"/>
<text x="{left+width/2}" y="415" font-size="14" text-anchor="middle" fill="#455569">Monte Carlo step</text>
<text x="{left}" y="54" font-size="12" fill="#455569">{escape(unit)}</text>
<text x="{left+width}" y="390" font-size="12" text-anchor="end" fill="#455569">{max_step} steps</text>
</svg>'''
    Path(path).write_text(svg + "\n", encoding="utf-8")