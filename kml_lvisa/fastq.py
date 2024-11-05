from pathlib import Path
import pandas as pd
from subprocess import run
import re
import logging


def prepare_fastq_by_samptab(workdir, samptab: str) -> None:
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

    # 软链接或解压
    for row in df.iterrows():
        name, fastq1, fastq2 = row[1]
        copy_fastq(workdir, name, fastq1, fastq2)


def get_sample_names_by_samptab(samptab: str) -> list:
    """
    获取样本名列表
    :param samptab:
    :return sample_names:   样本名列表
    """
    logging.info('获取样本名列表')
    df = samptab2dataframe(samptab)
    return df.iloc[:, 0].to_list()


def samptab2dataframe(samptab: str) -> pd.DataFrame:
    """
    输入 sample table 转成 DataFrame 格式

    :param samptab:
    :return df: sample table 转的 DataFrame
    """
    logging.info('输入 sample table 转成 DataFrame 格式')
    if samptab.endswith('.xlsx'):
        df = pd.read_excel(samptab, header=None)
    elif samptab.endswith('.tsv'):
        df = pd.read_table(samptab, sep='\t', header=None)
    else:
        raise ValueError(f'sample table 扩展名必须是 .xlsx or .tsv : {samptab}')

    # 检查 sample table
    check_samptab(df)

    return df


def copy_fastq(workdir, name, fq1: str, fq2) -> None:
    """
    复制或压缩 fastq 到 .rawdata 目, 按照指定格式命名
    :param workdir:     分析解雇目录
    :param name:        样本名
    :param fq1:         fastq1
    :param fq2:         fastq2
    """
    logging.info('复制或压缩 fastq 到 .rawdata 目, 按照指定格式命名')
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
    logging.debug(cml)
    run(cml, shell=True, executable='/bin/bash', capture_output=True)


def check_samptab(df) -> None:
    """
    检查 sample table 文件, 输入 sample table 转的 DataFrame
    :param df:
    :return:
    """
    logging.info('检查 sample table 文件, 输入 sample table 转的 DataFrame')
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
