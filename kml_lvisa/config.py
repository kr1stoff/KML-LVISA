from pathlib import Path
import yaml
import os
import math


def get_conda_env_dict() -> dict:
    """获取 Conda 环境字典"""
    yaml_conda_env = Path(__file__).resolve().parents[1].joinpath('config/conda_env.yaml')
    with open(yaml_conda_env) as f:
        dict_conda_env = yaml.safe_load(f)
    return dict_conda_env


def get_database_dict() -> dict:
    """获取数据库字典"""
    yaml_db = Path(__file__).resolve().parents[1].joinpath('config/database.yaml')
    with open(yaml_db) as f:
        dict_db = yaml.safe_load(f)
    # assets 一些小型数据库文件
    # 3' LTR 识别序列
    dict_db['3ltr'] = str(Path(__file__).resolve().parents[1].joinpath('assets/3LTR/3LTR.fa'))
    # UMI 参考文件
    dict_db['umi'] = str(Path(__file__).resolve().parents[1].joinpath('assets/UMI/ajtk_umi.json'))
    # [250715 mxfa] 阳控参考. 按照阴控逻辑列所有位点，无需预设参考
    # dict_db['pc'] = str(Path(__file__).resolve().parents[1].joinpath('assets/PC/CBPL0002.txt'))
    return dict_db


def get_software_dict() -> dict:
    """获取软件字典"""
    yaml_software = Path(__file__).resolve().parents[1].joinpath('config/software.yaml')
    with open(yaml_software) as f:
        dict_soft = yaml.safe_load(f)
    return dict_soft


def get_threads_dict() -> dict:
    """获取最大线程数, 高线程数为 max cpu count / 2, 低线程数为 high / 4"""
    # high: 当前总线程 / 2; low: high / 4
    max_threads = os.cpu_count()
    if max_threads is None:
        max_threads = 8
    high_threads = math.floor(max_threads / 2)
    low_threads = math.floor(high_threads / 4)
    dict_thr = {'high': high_threads, 'low': low_threads, 'max': max_threads}
    return dict_thr


def get_my_scripts_path() -> str:
    """获取 scripts 目录地址"""
    return str(Path(__file__).resolve().parents[1].joinpath('scripts'))
