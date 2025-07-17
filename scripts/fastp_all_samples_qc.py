import sys
from pathlib import Path
import json
import pandas as pd
import numpy as np


def fastp_all_samples_qc(files_fastp_json, out_tsv, out_excel):
    title = ["Sample", "RawReads", "RawBases", "CleanReads", "CleanBases", "RawQ20",
             "RawQ30", "CleanQ20", "CleanQ30", "CleanAverageLength", "GC"]
    df = pd.DataFrame(columns=title)
    for js_path in files_fastp_json:
        js_data = json.loads(open(js_path, "r").read())
        sample = Path(js_path).stem
        mean_lengths = np.array(
            [v for k, v in js_data["summary"]["after_filtering"].items() if k.endswith("mean_length")])
        out = [
            sample,
            js_data["summary"]["before_filtering"]["total_reads"],
            js_data["summary"]["before_filtering"]["total_bases"],
            js_data["summary"]["after_filtering"]["total_reads"],
            js_data["summary"]["after_filtering"]["total_bases"],
            js_data["summary"]["before_filtering"]["q20_rate"],
            js_data["summary"]["before_filtering"]["q30_rate"],
            js_data["summary"]["after_filtering"]["q20_rate"],
            js_data["summary"]["after_filtering"]["q30_rate"],
            int(mean_lengths.mean()),
            js_data["summary"]["after_filtering"]["gc_content"],
        ]
        df.loc[len(df)] = out
    df.to_csv(out_tsv, index=False, sep="\t")
    df.to_excel(out_excel, index=False)


if __name__ == "__main__":
    fastp_all_samples_qc(sys.argv[3:], sys.argv[1], sys.argv[2])
