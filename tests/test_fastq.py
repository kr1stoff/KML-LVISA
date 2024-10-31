from kml_lvisa import prepare_fastq_by_samptab
from kml_lvisa import get_sample_names_by_samptab

work_dir = '/data/mengxf/Project/KML240924_lvis_pipeline/result/241031'
sample_table = '/data/mengxf/GitHub/KML-LVISA/templates/input.tsv'


# def test_prepare():
#     prepare_fastq_by_samptab(work_dir, sample_table)


def test_get():
    get_sample_names_by_samptab(sample_table)
