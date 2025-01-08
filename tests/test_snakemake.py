from kml_lvisa import create_snakemake_configfile
from kml_lvisa import get_sample_names_by_samptab
from kml_lvisa import run_snakemake


work_dir = '/data/mengxf/Project/KML241220_lvis_YuShiYan5/results/25010701'
sample_table = '/data/mengxf/Project/KML241220_lvis_YuShiYan5/FASTQ/input_tsv/input250107.tsv'


def test_create():
    sample_names = get_sample_names_by_samptab(sample_table)
    create_snakemake_configfile(sample_names, work_dir)


# def test_run():
#     run_snakemake(work_dir)
