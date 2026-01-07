from pathlib import Path
import pandas as pd
from subprocess import run
import re
import logging
from multiprocessing import Pool

from kml_lvisa.config import get_threads_dict


def prepare_fastq_by_samptab(workdir, samptab: Path) -> None:
    """
    在项目目录下面准备 fastq 文件
    - 如果未压缩 link, 如果压缩 zcat
    - 支持 .tsv 和 .xlsx 格式
    :param workdir:         工作目录
    :param samptab:     样本信息表 samptab
    :return:
    """
    logging.info('在项目目录下面准备 fastq 文件')
    # 创建 {workdir}/.rawdata
    Path(workdir).joinpath('.rawdata').mkdir(exist_ok=True, parents=True)
    df = samptab2dataframe(samptab)
    # * 多线程软链接或解压
    para_num = get_threads_dict()['max']
    pargs = [(workdir, row[1][0], row[1][1], row[1][2]) for row in df.iterrows()]
    with Pool(para_num) as pool:
        pool.map(copy_fastq, pargs)


def copy_fastq(para_args: tuple) -> None:
    """
    多线程复制或压缩 fastq 到 .rawdata 目, 按照指定格式命名
    :param para_args:   (workdir, name, fq1, fq2)
    :subparam workdir:     分析解雇目录
    :subparam name:        样本名
    :subparam fq1:         fastq1
    :subparam fq2:         fastq2
    """
    workdir, name, fq1, fq2 = para_args
    if fq1.endswith('.gz'):
        cml = f"""
        cp {fq1} {workdir}/.rawdata/{name}_1.fastq.gz
        cp {fq2} {workdir}/.rawdata/{name}_2.fastq.gz
        """
    else:
        cml = f"""
        gzip -c {fq1} > {workdir}/.rawdata/{name}_1.fastq.gz
        gzip -c {fq2} > {workdir}/.rawdata/{name}_2.fastq.gz
        """
    # logging.debug(cml)
    run(cml, shell=True, executable='/bin/bash', capture_output=True)


def get_sample_names_by_samptab(samptab: Path) -> list:
    """
    获取样本名列表
    :param samptab:
    :return sample_names:   样本名列表
    """
    logging.info('获取样本名列表')
    df = samptab2dataframe(samptab)
    return df.iloc[:, 0].to_list()


def samptab2dataframe(samptab: Path) -> pd.DataFrame:
    """
    输入 sample table 转成 DataFrame 格式

    :param samptab:
    :return df: sample table 转的 DataFrame
    """
    samptabstr = str(samptab)
    if samptabstr.endswith('.xlsx'):
        df = pd.read_excel(samptab, header=None)
    elif samptabstr.endswith('.tsv') or samptabstr.endswith('.txt'):
        df = pd.read_table(samptab, sep='\t', header=None)
    else:
        raise ValueError(f'sample table 扩展名必须是 .xlsx or .tsv or .txt : {samptab}')
    # 检查 sample table
    check_samptab(df)
    return df


def check_samptab(df) -> None:
    """
    检查 sample table 文件, 输入 sample table 转的 DataFrame
    :param df:
    :return:
    """
    for row in df.iterrows():
        name, fastq1, fastq2 = row[1]
        # 检查名称
        pattern = r'[\\/:*?"<>| ]'
        assert not re.search(pattern, name), f'样本名称含有非法字符 (\\/:*?"<>| ) : {name}'
        # 检查 fastq 是否存在
        assert Path(fastq1).exists(), f'fastq1 不存在 : {fastq1}'
        assert Path(fastq2).exists(), f'fastq2 不存在 : {fastq2}'
        # 检查 fastq1 和 fastq2 是否相同
        assert fastq1 != fastq2, f'fastq1 和 fastq2 相同 : {fastq1} - {fastq2}'
