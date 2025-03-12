#!/usr/bin/env python
# -*- coding: utf-8 -*-

import gzip
import sys
import math
import multiprocessing
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

def phred_to_error_prob(q):
    """将Phred质量分数转换为错误率"""
    return 10 ** (-q / 10)

def error_prob_to_phred(p):
    """将错误率转换为Phred质量分数"""
    return -10 * math.log10(p) if p > 0 else 0

def read_fastq(fastq_file):
    if fastq_file.endswith('.gz'):
        open_func = gzip.open
        mode = 'rt'
    else:
        open_func = open
        mode = 'r'
    with open_func(fastq_file, mode) as f:
        record_lines = []
        for line in f:
            record_lines.append(line.strip())
            if len(record_lines) == 4:
                yield record_lines
                record_lines = []

def calculate_average_quality(quality_string):
    """计算碱基质量字符串的平均质量"""
    # 将每个碱基的质量分数转换为错误率
    error_probs = [phred_to_error_prob(ord(char) - 33) for char in quality_string]
    # 计算错误率的均值
    avg_error_prob = sum(error_probs) / len(error_probs)
    # 将错误率均值转换回Phred质量分数
    return error_prob_to_phred(avg_error_prob)

def process_fastq(fastq_file, outfile_preffix, threads):
    """处理FASTQ文件，计算每条read的长度和平均碱基质量"""
    results = []
    with multiprocessing.Pool(threads) as pool:
        for record_lines in read_fastq(fastq_file):
            results.append([record_lines[0], len(record_lines[1]) , pool.apply_async(calculate_average_quality, (record_lines[-1],))])
        pool.close()
        pool.join()
    with open(f"{outfile_preffix}.stat", 'w') as stat_file:
        for items in results:
            stat_file.write(f"{items[0]}\t{items[1]}\t{items[2].get():.2f}\n")

def plot_readQC_distribution(stat_file, tag_name, outfile_preffix):
    df = pd.read_table(stat_file, usecols=[1, 2, 3], names=['read_length', 'gc_content', 'read_quality'], dtype={'read_length':'int', 'gc_content':'float', 'read_quality':'float'}, header=None, compression="gzip")
    read_len, read_gc, read_qual = df.read_length.tolist(), df.gc_content.tolist(), df.read_quality.tolist()
    # for read length
    fig, ax = plt.subplots(figsize=(4, 3), constrained_layout=True)
    xmin, xmax = math.floor(np.quantile(read_len, 0.01)), math.ceil(np.quantile(read_len, 0.99))
    ax.hist(
        read_len,
        range=(xmin, xmax),
        bins=50,
        color="#006400",
        alpha=0.7,
    )
    ax.set_xlabel("Read length")
    ax.set_title(f"Read length distribution\n({tag_name})")
    ax.grid(axis="both", color="k", alpha=0.1, which="both")
    ax.spines[["right", "top"]].set_visible(False)
    fig.savefig(f"{outfile_preffix}.length_distribution.png", dpi=300)
    # for read gc content
    fig, ax = plt.subplots(figsize=(4, 3), constrained_layout=True)
    xmin, xmax = 0, 100
    ax.hist(
        read_gc,
        range=(xmin, xmax),
        bins=50,
        color="#006400",
        alpha=0.7,
    )
    ax.set_xlabel("GC content")
    ax.set_title(f"Read GC content distribution\n({tag_name})")
    ax.grid(axis="both", color="k", alpha=0.1, which="both")
    ax.spines[["right", "top"]].set_visible(False)
    fig.savefig(f"{outfile_preffix}.gc_content.png", dpi=300)
    # for read quality
    fig, ax = plt.subplots(figsize=(4, 3), constrained_layout=True)
    xmin, xmax = math.floor(np.quantile(read_qual, 0.01)), math.ceil(np.quantile(read_qual, 0.99))
    ax.hist(
        read_qual,
        range=(xmin, xmax),
        bins=50,
        color="#006400",
        alpha=0.7,
    )
    ax.set_xlabel("Read quality")
    ax.set_title(f"Read quality distribution\n({tag_name})")
    ax.grid(axis="both", color="k", alpha=0.1, which="both")
    ax.spines[["right", "top"]].set_visible(False)
    fig.savefig(f"{outfile_preffix}.quality_distribution.png", dpi=300)

if __name__ == "__main__":

    stat_file = sys.argv[1]
    tag_name = sys.argv[2]
    outfile_preffix = sys.argv[3]
    # threads = int(sys.argv[3])
    # process_fastq(fastq_file_path, outfile_preffix, threads)
    plot_readQC_distribution(stat_file, tag_name, outfile_preffix)
