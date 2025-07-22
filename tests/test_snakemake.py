from kml_lvisa import create_snakemake_configfile
from kml_lvisa import get_sample_names_by_samptab
from kml_lvisa import run_snakemake
from pathlib import Path


work_dir = "/data/mengxf/Project/KML250709-lvis-jiance-run1-2/result/250714"
sample_table = "/data/mengxf/Project/KML250709-lvis-jiance-run1-2/input.tsv"
sample_table = Path(sample_table).resolve()


sample_names = get_sample_names_by_samptab(sample_table)
create_snakemake_configfile(sample_names, work_dir)

# run_snakemake(work_dir)
