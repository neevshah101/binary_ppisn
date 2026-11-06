#!/usr/bin/env python3
"""
HR diagram of the single-star grid in runs/single_star/.

2x2 panels: rows = non-rotating / rotating (_w_ suffix), columns = metallicity.
Tracks are coloured by initial mass; runs without LOGS/history.data are skipped.

Usage:
    python plot_single_star_hr.py [-o hr_single_star.pdf]
"""

import argparse
import re
from pathlib import Path

import matplotlib.pyplot as plt
import mesa_reader as mr
import numpy as np

REPO = Path(__file__).parent.resolve()
RUNS_DIR = REPO / "runs" / "single_star"

# Sequential blue ramp (light -> dark = low -> high mass)
RAMP = ["#86b6ef", "#6da7ec", "#5598e7", "#3987e5", "#2a78d6",
        "#256abf", "#1c5cab", "#184f95", "#104281", "#0d366b"]

RUN_RE = re.compile(r"^m_(?P<m>[\d.]+)(?:_Pextra[\dp]+)?(?P<w>_w_[\d.]+)?$")


def z_value(zdir):
    """'z_1d-4' -> 1e-4."""
    return float(zdir.name[2:].replace("d", "e"))


def load_track(run_dir):
    h = mr.MesaData(str(run_dir / "LOGS" / "history.data"))
    # drop models superseded by retries/restarts: keep the last occurrence of each model_number
    mn = h.data("model_number")
    keep = np.ones(len(mn), bool)
    keep[:-1] = mn[:-1] < np.minimum.accumulate(mn[::-1])[::-1][1:]
    return h.data("log_Teff")[keep], h.data("log_L")[keep]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("-o", "--output", default=str(RUNS_DIR / "hr_single_star.pdf"))
    args = ap.parse_args()

    zdirs = sorted(RUNS_DIR.glob("z_*"), key=z_value)
    runs = []  # (z, rotating, mass, run_dir)
    for zd in zdirs:
        for rd in zd.iterdir():
            m = RUN_RE.match(rd.name)
            if m and (rd / "LOGS" / "history.data").exists():
                runs.append((z_value(zd), bool(m["w"]), float(m["m"]), rd, m["w"]))
    if not runs:
        raise SystemExit(f"No history.data found under {RUNS_DIR}")

    masses = sorted({r[2] for r in runs})
    colors = {m: RAMP[round(i * (len(RAMP) - 1) / max(1, len(masses) - 1))]
              for i, m in enumerate(masses)}

    fig, axes = plt.subplots(2, len(zdirs), figsize=(5 * len(zdirs), 9),
                             sharex=True, sharey=True, squeeze=False)
    for j, zd in enumerate(zdirs):
        z = z_value(zd)
        for i, rot in enumerate([False, True]):
            ax = axes[i, j]
            sel = sorted((r for r in runs if r[0] == z and r[1] == rot), key=lambda r: r[2])
            for _, _, mass, rd, w in sel:
                teff, lum = load_track(rd)
                ax.plot(teff, lum, color=colors[mass], lw=1.2, label=f"{mass:g}")
                ax.plot(teff[-1], lum[-1], "o", ms=4, color=colors[mass],
                        mec="white", mew=0.8)
            wlabel = ("rotating, $\\omega/\\omega_{crit}$ = " + w.split("_")[-1]) if sel and rot \
                else ("rotating" if rot else "non-rotating")
            ax.set_title(f"Z = {z:g}, {wlabel}", fontsize=11)
            if not sel:
                ax.text(0.5, 0.5, "no models", transform=ax.transAxes,
                        ha="center", va="center", color="0.5")
            ax.grid(color="0.9", lw=0.6)
            ax.tick_params(direction="in", top=True, right=True)
            if i == 1:
                ax.set_xlabel(r"$\log_{10}(T_\mathrm{eff}/\mathrm{K})$")
            if j == 0:
                ax.set_ylabel(r"$\log_{10}(L/L_\odot)$")
    axes[0, 0].invert_xaxis()  # shared x: inverts all panels

    handles = [plt.Line2D([], [], color=colors[m], lw=2) for m in masses]
    fig.legend(handles, [f"{m:g}" for m in masses], title=r"$M_\mathrm{ini}/M_\odot$",
               loc="center right", frameon=False)
    fig.text(0.5, 0.005, "dot = last model in history.data", ha="center",
             fontsize=8, color="0.4")
    fig.tight_layout(rect=(0, 0.02, 0.9, 1))
    fig.savefig(args.output)
    print(f"Saved {args.output} ({len(runs)} tracks)")


if __name__ == "__main__":
    main()
