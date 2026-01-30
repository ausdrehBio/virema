# compare results

from __future__ import annotations

import math
import sys

import matplotlib.pyplot as plt
import pandas as pd


VALID_PATH = "data/valid_data/Alnaji2021_SRR14352106.csv"
OWN_PATHs = ["output/PR8/SRR14352105/Virus_Recombination_Results.csv",
             "output/PR8/SRR14352106/Virus_Recombination_Results.csv",
             "output/PR8/SRR14352107/Virus_Recombination_Results.csv",
             "output/PR8/SRR14352108/Virus_Recombination_Results.csv",
             "output/PR8/SRR14352109/Virus_Recombination_Results.csv",
             "output/PR8/SRR14352110/Virus_Recombination_Results.csv",
             "output/PR8/SRR14352111/Virus_Recombination_Results.csv",
             "output/PR8/SRR14352112/Virus_Recombination_Results.csv",
             "output/PR8/SRR14352113/Virus_Recombination_Results.csv",
             "output/PR8/SRR14352114/Virus_Recombination_Results.csv",
             "output/PR8/SRR14352115/Virus_Recombination_Results.csv",
             "output/PR8/SRR14352116/Virus_Recombination_Results.csv",
             "output/PR8/SRR14352117/Virus_Recombination_Results.csv"]


OWN_PATHs = ["output/PR8/SRR14352106/Virus_Recombination_Results.csv"]

VENN_PATH = "output/venn_comparison.png"

KEY_COLUMNS = ["Reference", "Start", "End"]


def _normalize_valid(df: pd.DataFrame) -> pd.DataFrame:
    df = df.rename(columns={"Segment": "Reference"})
    df = df[KEY_COLUMNS + ["NGS_read_count"]].copy()
    return _normalize_types(df)


def _normalize_own(df: pd.DataFrame) -> pd.DataFrame:
    df = df[KEY_COLUMNS + ["NGS_read_count"]].copy()
    return _normalize_types(df)


def _normalize_types(df: pd.DataFrame) -> pd.DataFrame:
    df["Reference"] = df["Reference"].astype(str).str.strip()
    df["Start"] = pd.to_numeric(df["Start"], errors="raise").astype(int)
    df["End"] = pd.to_numeric(df["End"], errors="raise").astype(int)
    df["NGS_read_count"] = pd.to_numeric(
        df["NGS_read_count"], errors="raise"
    ).astype(int)
    return df


def _aggregate(df: pd.DataFrame) -> pd.DataFrame:
    return (
        df.groupby(KEY_COLUMNS, as_index=False)["NGS_read_count"]
        .sum()
        .sort_values(KEY_COLUMNS, kind="mergesort")
        .reset_index(drop=True)
    )


def compare_results(valid_path: str, own_paths: list[str]) -> int:
    valid_df = _aggregate(
        _normalize_valid(pd.read_csv(valid_path, keep_default_na=False))
    )
    own_frames = [
        _normalize_own(pd.read_csv(path, keep_default_na=False))
        for path in own_paths
    ]
    own_df = _aggregate(pd.concat(own_frames, ignore_index=True))

    key_only = valid_df[KEY_COLUMNS].merge(
        own_df[KEY_COLUMNS], on=KEY_COLUMNS, how="outer", indicator=True
    )
    missing_in_own = key_only[key_only["_merge"] == "left_only"][KEY_COLUMNS]
    extra_in_own = key_only[key_only["_merge"] == "right_only"][KEY_COLUMNS]

    merged_counts = valid_df.merge(
        own_df, on=KEY_COLUMNS, how="inner", suffixes=("_valid", "_own")
    )
    count_mismatch = merged_counts[
        merged_counts["NGS_read_count_valid"] != merged_counts["NGS_read_count_own"]
    ]

    valid_set = set(tuple(row) for row in valid_df[KEY_COLUMNS].to_numpy())
    own_set = set(tuple(row) for row in own_df[KEY_COLUMNS].to_numpy())
    _plot_venn(valid_set, own_set, VENN_PATH)

    if missing_in_own.empty and extra_in_own.empty and count_mismatch.empty:
        print("OK: Beide Dataframes enthalten die selben Ergebnisse.")
        return 0

    if not missing_in_own.empty:
        print("Fehlend in Virus_Recombination_Results.csv (nur in Alnaji2021.csv):")
        missing_rows = valid_df.merge(missing_in_own, on=KEY_COLUMNS, how="inner")
        print(missing_rows.to_string(index=False))
    if not extra_in_own.empty:
        print("Zusätzlich in Virus_Recombination_Results.csv (nicht in Alnaji2021.csv):")
        extra_rows = own_df.merge(extra_in_own, on=KEY_COLUMNS, how="inner")
        print(extra_rows.to_string(index=False))
    if not count_mismatch.empty:
        print("NGS_read_count Unterschiede:")
        print(
            count_mismatch[
                KEY_COLUMNS + ["NGS_read_count_valid", "NGS_read_count_own"]
            ].to_string(index=False)
        )

    return 1


def _circle_intersection_area(r1: float, r2: float, distance: float) -> float:
    if distance >= r1 + r2:
        return 0.0
    if distance <= abs(r1 - r2):
        return math.pi * min(r1, r2) ** 2
    r1_sq = r1 * r1
    r2_sq = r2 * r2
    angle1 = math.acos((distance * distance + r1_sq - r2_sq) / (2 * distance * r1))
    angle2 = math.acos((distance * distance + r2_sq - r1_sq) / (2 * distance * r2))
    overlap = (
        r1_sq * angle1
        + r2_sq * angle2
        - 0.5
        * math.sqrt(
            max(
                0.0,
                (-distance + r1 + r2)
                * (distance + r1 - r2)
                * (distance - r1 + r2)
                * (distance + r1 + r2),
            )
        )
    )
    return overlap


def _solve_circle_distance(r1: float, r2: float, target_area: float) -> float:
    max_area = math.pi * min(r1, r2) ** 2
    target = min(max(target_area, 0.0), max_area)
    if target == 0.0:
        return r1 + r2
    if target == max_area:
        return abs(r1 - r2)
    low = abs(r1 - r2)
    high = r1 + r2
    for _ in range(60):
        mid = (low + high) / 2.0
        area = _circle_intersection_area(r1, r2, mid)
        if area > target:
            low = mid
        else:
            high = mid
    return (low + high) / 2.0


def _plot_venn(valid_set: set[tuple], own_set: set[tuple], out_path: str) -> None:
    only_valid = len(valid_set - own_set)
    only_own = len(own_set - valid_set)
    both = len(valid_set & own_set)
    total_valid = only_valid + both
    total_own = only_own + both

    r1 = math.sqrt(total_valid / math.pi) if total_valid > 0 else 0.0
    r2 = math.sqrt(total_own / math.pi) if total_own > 0 else 0.0
    distance = _solve_circle_distance(r1, r2, both)
    left_center = (0.0, 0.0)
    right_center = (distance, 0.0)

    fig, ax = plt.subplots(figsize=(6, 4))
    ax.set_aspect("equal")
    ax.axis("off")

    if r1 > 0:
        ax.add_patch(
            plt.Circle(left_center, r1, color="#4C78A8", alpha=0.4, linewidth=0)
        )
    if r2 > 0:
        ax.add_patch(
            plt.Circle(right_center, r2, color="#F58518", alpha=0.4, linewidth=0)
        )

    left_x = left_center[0] - (r1 * 0.5 if r1 > 0 else 0.2)
    right_x = right_center[0] + (r2 * 0.5 if r2 > 0 else 0.2)
    overlap_x = (left_center[0] + right_center[0]) / 2.0

    ax.text(left_x, 0.0, str(only_valid), ha="center", va="center", fontsize=12)
    ax.text(right_x, 0.0, str(only_own), ha="center", va="center", fontsize=12)
    ax.text(overlap_x, 0.0, str(both), ha="center", va="center", fontsize=12)

    max_r = max(r1, r2)
    label_y = -max_r * 1.15 if max_r > 0 else -0.5
    ax.text(
        left_center[0],
        label_y,
        "Jens",
        ha="center",
        va="center",
        fontsize=10,
    )
    ax.text(
        right_center[0],
        label_y,
        "Own",
        ha="center",
        va="center",
        fontsize=10,
    )

    padding = max_r * 0.25 if max_r > 0 else 0.5
    x_min = min(left_center[0] - r1, right_center[0] - r2) - padding
    x_max = max(left_center[0] + r1, right_center[0] + r2) + padding
    y_min = -max_r - padding
    y_max = max_r + padding
    ax.set_xlim(x_min, x_max)
    ax.set_ylim(y_min, y_max)

    fig.tight_layout()
    fig.savefig(out_path, dpi=200)
    plt.close(fig)


if __name__ == "__main__":
    sys.exit(compare_results(VALID_PATH, OWN_PATHs))
