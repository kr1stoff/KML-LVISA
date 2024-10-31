from kml_lvisa import create_snakemake_configfile
from kml_lvisa import get_sample_names_by_samptab


work_dir = '/data/mengxf/Project/KML240924_lvis_pipeline/result/241031'
sample_table = '/data/mengxf/GitHub/KML-LVISA/templates/input.tsv'


def test_create():
    sample_names = get_sample_names_by_samptab(sample_table)
    create_snakemake_configfile(sample_names, work_dir)
