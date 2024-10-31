import logging
from kml_lvisa import prepare_fastq_by_samptab
from kml_lvisa import get_sample_names_by_samptab
from kml_lvisa import get_threads_dict
from kml_lvisa import get_conda_env_dict
from kml_lvisa import get_database_dict
from kml_lvisa import get_my_scripts_path
from kml_lvisa import get_software_dict


logging.basicConfig(level=logging.DEBUG,
                    format="%(asctime)s - %(levelname)s - %(filename)s - %(message)s",
                    datefmt="%Y-%m-%d %H:%M:%S")

work_dir = '/data/mengxf/Project/KML240924_lvis_pipeline/result/241031'
sample_table = '/data/mengxf/GitHub/KML-LVISA/templates/input.tsv'

# main
# config
threads_dict = get_threads_dict()
conda_env_dict = get_conda_env_dict()
database_dict = get_database_dict()
software_dict = get_software_dict()
my_scripts_path = get_my_scripts_path()
# fastq
sample_names = get_sample_names_by_samptab(sample_table)
prepare_fastq_by_samptab(work_dir, sample_table)
