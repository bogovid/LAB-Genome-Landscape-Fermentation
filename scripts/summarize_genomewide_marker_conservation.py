#!/usr/bin/env python3
"""Summarize pairwise genome-wide distance/ANI versus 28-marker discordance.

Designed for the unique-pair CSV tables archived with the v1.1.0 Zenodo data package.
Expected columns include either `Mash_distance` or `ANI`, plus a marker discordance
column named `Marker_discordance` or `Marker_hamming`.
"""

from __future__ import annotations

import argparse
import pandas as pd


def find_col(df: pd.DataFrame, names: list[str]) -> str:
    for name in names:
        if name in df.columns:
            return name
    raise KeyError(f"None of the expected columns were found: {names}")


def spearman_rho(x: pd.Series, y: pd.Series) -> float:
    """Compute Spearman rho as Pearson correlation of average ranks."""
    return float(x.rank(method="average").corr(y.rank(method="average")))


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("csv", help="Unique-pair CSV or CSV.GZ")
    parser.add_argument("--kind", choices=["mash", "fastani"], required=True)
    args = parser.parse_args()

    df = pd.read_csv(args.csv)
    discord_col = find_col(df, ["Marker_discordance", "Marker_hamming", "Marker_differences", "N_marker_differences"])
    identical = df[discord_col] == 0

    print(f"Unique pairs: {len(df)}")
    print(f"Identical 28-marker profiles: {identical.sum()} ({identical.mean()*100:.6f}%)")

    if args.kind == "mash":
        metric = find_col(df, ["Mash_distance", "mash_distance"])
        rho = spearman_rho(df[metric], df[discord_col])
        subset = df[df[metric] <= 0.10]
        print(f"Mash distance range: {df[metric].min():.6f} to {df[metric].max():.6f}")
        print(f"Spearman rho, Mash distance vs marker discordance: {rho:.6f}")
        print(f"Pairs Mash <=0.10: {len(subset)}; identical profiles: {(subset[discord_col] == 0).sum()}")
    else:
        metric = find_col(df, ["ANI", "ANI_pct", "FastANI"])
        print(f"ANI range: {df[metric].min():.4f} to {df[metric].max():.4f}")
        print(f"Pairs differing at exactly one marker: {(df[discord_col] == 1).sum()}")


if __name__ == "__main__":
    main()
