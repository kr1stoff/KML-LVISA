from kml_lvisa import prepare_fastq_by_samptab
from kml_lvisa import get_sample_names_by_samptab

work_dir = '/data/mengxf/Project/KML241220_lvis_YuShiYan5/results/25011001'
sample_table = '/data/mengxf/Project/KML241220_lvis_YuShiYan5/input/input.tsv'


def test_prepare():
    prepare_fastq_by_samptab(work_dir, sample_table)


def test_get():
    names = get_sample_names_by_samptab(sample_table)
    print(names)
