import logging
import click
from pathlib import Path
from kml_lvisa import prepare_fastq_by_samptab
from kml_lvisa import get_sample_names_by_samptab
from kml_lvisa import get_threads_dict
from kml_lvisa import get_conda_env_dict
from kml_lvisa import get_database_dict
from kml_lvisa import get_my_scripts_path
from kml_lvisa import get_software_dict
from kml_lvisa import create_snakemake_configfile
from kml_lvisa import run_snakemake


logging.basicConfig(level=logging.DEBUG,
                    format="%(asctime)s - %(levelname)s - %(filename)s - %(message)s",
                    datefmt="%Y-%m-%d %H:%M:%S")


@click.command()
@click.option('--sample_table', '-s', type=click.Path(exists=True), required=True, help='样本信息表.')
@click.option('--work_dir', '-w', type=str, default='lvisa_result', help='结果生成目录.')
@click.help_option('-h', '--help')
def main(work_dir, sample_table):
    logging.info(f'开始分析!')
    sample_table = Path(sample_table).resolve()
    work_dir = Path(work_dir).resolve()
    # work_dir = '/data/mengxf/Project/KML240924_lvis_pipeline/result/241031'
    # sample_table = '/data/mengxf/GitHub/KML-LVISA/templates/input.tsv'

    # config
    threads_dict = get_threads_dict()
    conda_env_dict = get_conda_env_dict()
    database_dict = get_database_dict()
    software_dict = get_software_dict()
    my_scripts_path = get_my_scripts_path()
    # fastq
    sample_names = get_sample_names_by_samptab(sample_table)
    prepare_fastq_by_samptab(work_dir, sample_table)

    # snakemake
    create_snakemake_configfile(sample_names, work_dir)
    run_snakemake(work_dir)

    logging.info(f'分析完成!')
