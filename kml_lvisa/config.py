from pathlib import Path
import yaml
import logging
import os
import math


def get_conda_env_dict() -> dict:
    """
    获取 Conda 环境字典
    :return:    环境字典
    """
    logging.info('获取环境字典')
    yaml_conda_env = Path(__file__).resolve().parents[1].joinpath('config/conda_env.yaml')

    with open(yaml_conda_env) as f:
        dict_conda_env = yaml.safe_load(f)

    return dict_conda_env


def get_database_dict() -> dict:
    logging.info('获取数据库字典')
    yaml_database = Path(__file__).resolve().parents[1].joinpath('config/database.yaml')

    with open(yaml_database) as f:
        dict_database = yaml.safe_load(f)

    # assets 一些小型数据库文件
    ltr_fa = Path(__file__).resolve().parents[1].joinpath('assets/3LTR/3LTR.fa')
    dict_database['3ltr'] = str(ltr_fa)

    return dict_database


def get_software_dict() -> dict:
    logging.info('获取软件字典')
    yaml_software = Path(__file__).resolve().parents[1].joinpath('config/software.yaml')

    with open(yaml_software) as f:
        dict_soft = yaml.safe_load(f)

    return dict_soft


def get_threads_dict() -> dict:
    """
    获取最大线程数, 高线程分配为 max cpu count, 低线程为 high / 4
    :return dict_thr:   high_threads, low_threads 高/低线程分配数
    """
    logging.info('获取线程数字典')
    max_threads = os.cpu_count()
    high_threads = max_threads
    low_threads = math.floor(high_threads / 4)
    dict_thr = {'high': high_threads, 'low': low_threads}

    return dict_thr


def get_my_scripts_path() -> str:
    logging.info('获取 scripts 目录地址')
    return str(Path(__file__).resolve().parents[1].joinpath('scripts'))
