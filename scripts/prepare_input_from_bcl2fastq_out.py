import click
from pathlib import Path


@click.command()
@click.option("--input-dir", type=click.Path(exists=True), required=True,
              help="bcl2fastq 结果 fastq 文件夹, 按照 samplesheet 中的 sample_id 来定义样本名和分配 fastq")
@click.option("--output-file", default="input.tsv", show_default=True, help="慢病毒分析流程输入文件")
@click.help_option("--help", help="显示帮助信息并退出")
def main(input_dir: str, output_file: str):
    fqdict = {}
    for fq in Path(input_dir).glob("*.fastq.gz"):
        if "Undetermined" in fq.name:
            continue
        sampleid = fq.name.split("_S")[0]
        fqdict.setdefault(sampleid, {})
        if "R1" in fq.name:
            fqdict[sampleid].update({"read1": str(fq)})
        elif "R2" in fq.name:
            fqdict[sampleid].update({"read2": str(fq)})
    with open(output_file, "w") as f:
        for k, v in fqdict.items():
            f.write(f"{k}\t{v["read1"]}\t{v["read2"]}\n")


if __name__ == "__main__":
    main()
