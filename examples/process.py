import re
import argparse
import pandas as pd
from tqdm.auto import tqdm
import json
import numpy as np


parser = argparse.ArgumentParser()

parser.add_argument("--folder", "-f")

args = parser.parse_args()
with open(f"{args.folder}/time.out", "r") as f:
    times = list(
        map(
            lambda x: x.replace("\n", ","),
            filter(lambda x: len(x) > 0, f.read().split("\n\n")),
        )
    )

time_df = pd.DataFrame(
    columns=["name", "parameters", "memory", "cpu_time", "elapsed_time"]
)

for idx, each in enumerate(tqdm(times)):
    b = list(filter(lambda x: len(x) > 0, each.split(",")))
    if b[0].split("/")[-1] in ["entropies.mpl", "no_transcendence.mpl"]:
        continue
    time_df.loc[idx, "name"] = b[0]
    with open(f"{args.folder}/{b[0].split('/')[-1]}", "r") as f:
        file = f.readlines()[4:]
        deg = file[0].split(" ")[-1]
        params = re.findall(r"\[.*\]", file[0])
        length = file[1].split(", ")[0].strip(" #")
        counter = ",".join(
            map(
                lambda x: x.split(" = ")[1],
                re.findall(r"\[.*\]", file[1])[0].strip("[]").split(", "),
            )
        )
        total_count, entr = list(
            map(lambda x: float(x.strip("# ")), file[2].split(", "))
        )
        occurrences = re.findall(r"\[.*\]", file[4])[0]
        sum_occ = sum(
            list(
                map(lambda x: int(x[2:]), re.findall(r"=\s*\d+", occurrences))
            )
        )

        trbasis_degree_counts = sorted(
            list(
                map(
                    lambda x: np.sum(-x * np.log2(x)),
                    map(
                        lambda x: np.array(json.loads(x))
                        / sum(np.array(json.loads(x))),
                        re.findall(r"\[[0-9, ]+\]", file[5]),
                    ),
                )
            ),
            reverse=False,
        )
        entropy_sum = sum(trbasis_degree_counts)

    if "maple: fatal error" in b:
        time_df.loc[idx, "is_error"] = 1
        b = b[6:]
        time_df.loc[idx, "parameters"] = params
        time_df.loc[idx, "memory"] = int(b[-1])
        time_df.loc[idx, "cpu_time"] = sum(
            map(lambda x: float(x), b[1].strip("=").split("+"))
        )
        time_df.loc[idx, "elapsed_time"] = float(b[0])
        time_df.loc[idx, "sum_occurrences"] = sum_occ
        time_df.loc[idx, "occurrences"] = str(occurrences)
        time_df.loc[idx, "lex_entropy"] = str(trbasis_degree_counts)
        time_df.loc[idx, "entropy_sum"] = entropy_sum
        time_df.loc[idx, "degree_of_monomials"] = (
            int(deg.strip()) if deg else "NaN"
        )
        time_df.loc[idx, "length_of_counter"] = (
            int(length.strip("# ")) if length else "NaN"
        )
        time_df.loc[idx, "total_count"] = (
            int(total_count) if total_count else "NaN"
        )
        time_df.loc[idx, "entropy"] = entr if total_count else "NaN"
        time_df.loc[idx, "counter_string"] = counter
    else:
        b = b[4:]
        time_df.loc[idx, "parameters"] = params
        time_df.loc[idx, "memory"] = int(b[-1])
        time_df.loc[idx, "cpu_time"] = sum(
            map(lambda x: float(x), b[1].strip("=").split("+"))
        )
        time_df.loc[idx, "elapsed_time"] = float(b[0])
        time_df.loc[idx, "sum_occurrences"] = sum_occ
        time_df.loc[idx, "occurrences"] = str(occurrences)
        time_df.loc[idx, "lex_entropy"] = str(trbasis_degree_counts)
        time_df.loc[idx, "entropy_sum"] = entropy_sum
        time_df.loc[idx, "degree_of_monomials"] = (
            int(deg.strip()) if deg else "NaN"
        )
        time_df.loc[idx, "length_of_counter"] = (
            int(length.strip("# ")) if length else "NaN"
        )
        time_df.loc[idx, "total_count"] = (
            int(total_count) if total_count else "NaN"
        )
        time_df.loc[idx, "entropy"] = entr if total_count else "NaN"
        time_df.loc[idx, "counter_string"] = counter
time_df["cpu_time"] = time_df["cpu_time"].map(lambda x: f"{x:.2f}")
time_df.sort_values(by="lex_entropy", ascending=True).to_csv(
    f"{args.folder}/times.csv", index=False, float_format="%.6f"
)
