"""Converts validate.py results into weighted accuracy and signed counting error per class.

validate.py writes one row per (event, algo, param set). This script can be
changed to change the analysis without re-running the recon or validation(long part)
"""

import os
import pandas as pd

HERE = os.path.dirname(
    os.path.abspath(__file__)
)  # reads the CSV here regardless of cwd


def summarize(csv_path):
    df = pd.read_csv(csv_path)
    df["error"] = df["reco_nbub"] - df["scan_nbub"]  # +over-counts, -under-counts

    grouped = df.groupby(["algo", "params"])
    per_nbub = (
        df.groupby(["algo", "params", "scan_nbub"])["agree"].mean().mul(100).unstack()
    )

    overall = (
        pd.DataFrame(
            {
                "accuracy": grouped["agree"].mean().mul(100),  # every event equal
                "weighted": per_nbub.mean(axis=1),  # nbub=1,2,3... weighted equally
            }
        )
        .round(1)
        .sort_values("weighted", ascending=False)
    )

    bias = (
        df.groupby(["algo", "params", "scan_nbub"])["error"].mean().unstack().round(2)
    )
    bias = bias.reindex(overall.index)

    print(overall.to_string())
    print("\nsigned error per bubble count (>0 over-counts, <0 under-counts):")
    print(bias.to_string())
    return overall, bias


if __name__ == "__main__":
    summarize(os.path.join(HERE, "validation_results.csv"))
