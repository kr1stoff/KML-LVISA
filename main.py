import logging
import click
from pathlib import Path
from kml_lvisa import prepare_fastq_by_samptab
from kml_lvisa import get_sample_names_by_samptab
from kml_lvisa import create_snakemake_configfile
from kml_lvisa import run_snakemake


logging.basicConfig(level=logging.DEBUG,
                    format="%(asctime)s - %(levelname)s - %(filename)s - %(message)s",
                    datefmt="%Y-%m-%d %H:%M:%S")


@click.command()
@click.option('--sample-table', type=click.Path(exists=True), required=True, help='样本信息表.')
@click.option('--work-dir', type=str, default='lvisa_result', help='结果生成目录. [default: lvisa_result]')
@click.help_option(help="查看帮助信息.")
def main(work_dir, sample_table):
    """慢病毒插入位点分析流程."""
    logging.info(f'开始分析!')
    sample_table = Path(sample_table).resolve()
    work_dir = Path(work_dir).resolve()

    # fastq
    sample_names = get_sample_names_by_samptab(sample_table)
    prepare_fastq_by_samptab(work_dir, sample_table)

    # snakemake
    create_snakemake_configfile(sample_names, work_dir)
    run_snakemake(work_dir)

    logging.info(f'分析完成!')


if __name__ == '__main__':
    main()
