#!/usr/bin/env python
# -*- coding: utf-8 -*-

import pandas as pd
import altair as alt
import json
import base64
from jinja2 import Environment, FileSystemLoader
import numpy as np
import os
import matplotlib.pyplot as plt
import io


#创建最终目录
final_directory = 'report_files'
if not os.path.exists(final_directory):
    # 如果目录不存在，则创建目录
    os.makedirs(final_directory)
    print(f"Directory '{final_directory}' has been created.")
else:
    print(f"Directory '{final_directory}' already exists.")
    
def check_and_create_subdirectory(parent_directory, sub_directory_name):
    # 拼接父目录和子目录路径
    sub_directory = os.path.join(parent_directory, sub_directory_name)
    
    # 检查子目录是否存在，如果不存在则创建
    if not os.path.exists(sub_directory):
        os.makedirs(sub_directory)
        print(f"Subdirectory '{sub_directory}' has been created.")
    else:
        print(f"Subdirectory '{sub_directory}' already exists.")

check_and_create_subdirectory('report_files', 'fl_reads')
check_and_create_subdirectory('report_files', 'genome_alignments')
check_and_create_subdirectory('report_files', 'sqanti_qc')
check_and_create_subdirectory('report_files', 'transcriptome_alignments')
check_and_create_subdirectory('report_files', 'diff_exp')
check_and_create_subdirectory('report_files', 'diff_splice')
check_and_create_subdirectory('report_files', 'jaffal_fusion')
check_and_create_subdirectory('report_files', 'merge_counts')


def image_to_base64(image_path):
    try:
        with open(image_path, 'rb') as image_file:
            encoded_string = base64.b64encode(image_file.read()).decode('utf-8')
        return f"data:image/png;base64,{encoded_string}"
    except FileNotFoundError:
        print(f"错误：未找到文件 {image_path}")
        return None

#1.2 table

root_dir_1_2 = 'fl_reads'
dfs2 = []
# 用来存储文件路径的列表
file_paths_1_2 = []

# 遍历根目录下的所有子目录
for subdir in os.listdir(root_dir_1_2):
    subdir_path = os.path.join(root_dir_1_2, subdir)
    
    # 确保它是一个目录
    if os.path.isdir(subdir_path):
        # 查找该子目录下以 .identifying_statistic.txt 结尾的文件
        for file_name in os.listdir(subdir_path):
            if file_name.endswith('.identifying_statistic.txt'):
                # 构建文件的相对路径并添加到列表中
                file_paths_1_2.append(f'{root_dir_1_2}/{subdir}/{file_name}')

# 遍历每个文件路径
for file_path in file_paths_1_2:
    sample_name = file_path.split('/')[1]  # 获取Mix1_1, Mix1_2等
    # 读取文件内容
    with open(file_path, 'r') as file:
        lines = file.readlines()
    # 第一个数据框：第2-10行
    part1 = pd.read_csv(io.StringIO(''.join(lines[1:10])), sep='\t', header=None)
    
    part1.columns = ['Total_base_count', 'Valid_base_count', 'Valid_base_proportion(%)']
    result = {
    "Total bases": part1['Total_base_count'][1],
    "Valid bases": part1['Valid_base_count'][1],
    "Valid base %": part1['Valid_base_proportion(%)'][1],
    "Length-filtered reads %": part1["Valid_base_proportion(%)"][4],
    "QC-filtered reads %": part1["Valid_base_proportion(%)"][5],
    "Chimeric reads %": part1["Valid_base_proportion(%)"][6],
    "Full-length reads %": part1["Valid_base_proportion(%)"][8]
    }
    df = pd.DataFrame([result])
    df['Sample'] = sample_name
    dfs2.append(df)

merged_df3 = pd.concat(dfs2, ignore_index=True)
sample_col = merged_df3.pop('Sample')  # 移除Sample列
merged_df3.insert(0, 'Sample', sample_col)  # 在第一列插入Sample列
merged_df3.to_csv("report_files/fl_reads/ReadSummary.csv",index=False)
glycine_df_html = merged_df3.to_html(index=False, classes='data-table no-pagination', border=0)


#2.1 genome alignment table
root_dir_2_1 = 'genome_alignments'
dfs = []
# 用来存储文件路径的列表
file_paths_2_1 = []

# 遍历根目录下的所有子目录
for subdir in os.listdir(root_dir_2_1):
    subdir_path = os.path.join(root_dir_2_1, subdir)
    
    # 确保它是一个目录
    if os.path.isdir(subdir_path):
        # 查找该子目录下以 .identifying_statistic.txt 结尾的文件
        for file_name in os.listdir(subdir_path):
            if file_name == 'read_aln_stats.tsv':
                # 构建文件的相对路径并添加到列表中
                file_paths_2_1.append(f'{root_dir_2_1}/{subdir}/{file_name}')

# 遍历每个文件路径
for file_path in file_paths_2_1:
    # 读取当前文件
    df = pd.read_csv(file_path, sep='\t')
    
    # 获取文件名并将其作为新列添加到DataFrame中
    sample_name = file_path.split('/')[1]  # 获取Mix1_1, Mix1_2等
    df['Sample'] = sample_name
    
    # 将当前DataFrame添加到列表中
    dfs.append(df)

# 合并所有的DataFrame
merged_df = pd.concat(dfs, ignore_index=True)
sample_col = merged_df.pop('Sample')  # 移除Sample列
merged_df.insert(0, 'Sample', sample_col)  # 在第一列插入Sample列
merged_df.to_csv("report_files/genome_alignments/read_aln_stats.csv",index=False)
genome_alignments_table_html  = merged_df.to_html(index=False, classes='data-table no-pagination', border=0)
#print(genome_alignments_table_html)

# 2.3表格
OUT_classification_df = pd.read_csv("sqanti_qc/OUT.transcript_models_classification.txt", sep = "\t")
OUT_classification_df_part = OUT_classification_df[["isoform", "length", "exons", "structural_category", "associated_gene", "associated_transcript", "diff_to_TSS", "diff_to_TTS", "dist_to_CAGE_peak", "within_CAGE_peak"]].head(100)
OUT_classification_df.to_csv("report_files/sqanti_qc/OUT.transcript_models_classification.txt",index=False)
OUT_classification_df_part_html  = OUT_classification_df_part.to_html(index=False, classes='data-table', border=0)

# 3.1 Transcript Quantification and Differential Expression Analysis
root_dir_3_1 = 'transcriptome_alignments'
dfs = []
# 用来存储文件路径的列表
file_paths_3_1 = []

# 遍历根目录下的所有子目录
for subdir in os.listdir(root_dir_3_1):
    subdir_path = os.path.join(root_dir_3_1, subdir)
    
    # 确保它是一个目录
    if os.path.isdir(subdir_path):
        # 查找该子目录下以 .identifying_statistic.txt 结尾的文件
        for file_name in os.listdir(subdir_path):
            if file_name == 'read_aln_stats.tsv':
                # 构建文件的相对路径并添加到列表中
                file_paths_3_1.append(f'{root_dir_3_1}/{subdir}/{file_name}')

# 遍历每个文件路径
for file_path in file_paths_3_1:
    # 读取当前文件
    df = pd.read_csv(file_path, sep='\t')
    
    # 获取文件名并将其作为新列添加到DataFrame中
    sample_name = file_path.split('/')[1]  # 获取Mix1_1, Mix1_2等
    df['Sample'] = sample_name
    
    # 将当前DataFrame添加到列表中
    dfs.append(df)
    
merged_df2 = pd.concat(dfs, ignore_index=True)
sample_col2 = merged_df2.pop('Sample')  # 移除Sample列
merged_df2.insert(0, 'Sample', sample_col2)  # 在第一列插入Sample列
merged_df2.to_csv("report_files/transcriptome_alignments/read_aln_stats.csv",index=False)
genome_alignments_table_html2  = merged_df2.to_html(index=False, classes='data-table no-pagination', border=0)


#3.4 DGE Table
test_results_dge = pd.read_csv("diff_exp/results_dge.tsv",sep = "\t")
test_results_dge.to_csv("report_files/diff_exp/results_dge.tsv",index=False)
test_results_dge_1k = test_results_dge.head(100)
test_results_dge_html = test_results_dge_1k.to_html(index=False, classes='data-table', border=0)

#3.4 DTE Table
test_results_dte = pd.read_csv("diff_exp/results_dte.tsv",sep = "\t")
test_results_dte.to_csv("report_files/diff_exp/results_dte.tsv",index=False)
test_results_dte_1k = test_results_dte.head(100)
test_results_dte_html = test_results_dte_1k.to_html(index=False, classes='data-table', border=0)

#3.5 events.add_genename.psi
events_add_genename_df = pd.read_csv("diff_splice/events.add_genename.psi", sep= "\t")
events_add_genename_df.to_csv("report_files/diff_splice/events.add_genename.psi",index=False)
events_add_genename_df_1k = events_add_genename_df.head(100)
events_add_genename_df_html = events_add_genename_df_1k.to_html(index=False, classes='data-table', border=0)

#3.5 suppa_diffSplice_event.add_genename.dpsi
suppa_diffSplice_event_df = pd.read_csv("diff_splice/suppa_diffSplice_event.add_genename.dpsi", sep= "\t")
suppa_diffSplice_event_df.to_csv("report_files/diff_splice/suppa_diffSplice_event.add_genename.dpsi",index=False)
suppa_diffSplice_event_df_1k = suppa_diffSplice_event_df.head(100)
suppa_diffSplice_event_df_html = suppa_diffSplice_event_df_1k.to_html(index=False, classes='data-table', border=0)

#3.5 iso_isoform.add_genename.psi
iso_isoform_df = pd.read_csv("diff_splice/iso_isoform.add_genename.psi", sep= "\t")
iso_isoform_df.to_csv("report_files/diff_splice/iso_isoform.add_genename.psi",index=False)
iso_isoform_df_1k = iso_isoform_df.head(100)
iso_isoform_df_html = iso_isoform_df_1k.to_html(index=False, classes='data-table', border=0)

#3.5 suppa_diffSplice_iso.add_genename.dpsi
suppa_diffSplice_iso_df = pd.read_csv("diff_splice/suppa_diffSplice_iso.add_genename.dpsi", sep= "\t")
suppa_diffSplice_iso_df.to_csv("report_files/diff_splice/suppa_diffSplice_iso.add_genename.dpsi",index=False)
suppa_diffSplice_iso_df_1k = suppa_diffSplice_iso_df.head(100)
suppa_diffSplice_iso_df_html = suppa_diffSplice_iso_df_1k.to_html(index=False, classes='data-table', border=0)

#jaffa_results
jaffa_results_df = pd.read_csv("jaffal_fusion/jaffa_results.csv", sep= ",")
suppa_diffSplice_iso_df.to_csv("report_files/jaffal_fusion/jaffa_results.csv",index=False)
jaffa_results_df_1k = jaffa_results_df.head(100)
jaffa_results_df_html = jaffa_results_df_1k.to_html(index=False, classes='data-table', border=0)


# -----------------------------
# 4. 准备图片路径
# -----------------------------

# 为第一部分准备图片路径，基于实际文件命名模式
root_dir = 'fl_reads'
mix_categories = []
for subdir in os.listdir(root_dir):
    subdir_path = os.path.join(root_dir, subdir)
    if os.path.isdir(subdir_path) and subdir != '.ipynb_checkpoints':
        mix_categories.append(subdir)
        
fl_reads_images = {}

# 基于文件浏览器中看到的命名模式，每个Mix有6种类型的图片
image_types = [
    "passed.length_distribution.png",
    "passed.quality_distribution.png",
    "passed.gc_content.png",
    "usable.length_distribution.png",
    "usable.quality_distribution.png",
    "usable.gc_content.png"
]

for mix in mix_categories:
    # 将图片转换为 Base64 编码
    fl_reads_images[mix] = [image_to_base64(f"fl_reads/{mix}/{mix}.{img_type}") for img_type in image_types]



# 为第二部分准备图片路径
#transcript_depth_image = "known_transcripts_depth/Profile.png"
transcript_depth_image = image_to_base64("known_transcripts_depth/Profile.png")

#2.3 转录本结构
#transcript_structure_image = "sqanti_qc/OUT.transcript_models.evidence.png"
transcript_structure_image = image_to_base64("sqanti_qc/OUT.transcript_models.evidence.png")

#3.2 基因水平和转录本水平的覆盖率图
gene_transcript_cov_image = image_to_base64("known_transcripts_coverage/coverage.png")


#3.3 TPM table
all_tpm_df = pd.read_csv("merge_counts/all_TPM.add_genename.tsv", sep = "\t")
all_tpm_df.to_csv("report_files/merge_counts/all_TPM.add_genename.tsv",index=False)
all_tpm_df_1k = all_tpm_df.head(100)
all_tpm_df_html = all_tpm_df_1k.to_html(index=False, classes='data-table', border=0)

#3.4  DGE results_dge.maplot.png
results_dge_image = image_to_base64("diff_exp/results_dge.maplot.png")

#3.4 DTE results_dge.maplot.png
results_dte_image = image_to_base64("diff_exp/results_dte.maplot.png")


#3.5 sashimi plot 判断展示
plot_sashimi_dir = 'diff_splice/plot_sashimi'
which_to_show_path = os.path.join(plot_sashimi_dir, 'which_to_show')
which_most_sig_path = os.path.join(plot_sashimi_dir, 'which_most_sig')
# 读取文件内容
sashimi_filename = None
if os.path.exists(which_to_show_path):
    with open(which_to_show_path, 'r') as file:
        sashimi_filename = file.readline().strip()  # 读取文件的第一行
elif os.path.exists(which_most_sig_path):
    with open(which_most_sig_path, 'r') as file:
        sashimi_filename = file.readline().strip()  # 读取文件的第一行

# 如果两个文件都不存在，设置默认值

if sashimi_filename:
    Most_significance_sashimi_plot = image_to_base64(f"diff_splice/plot_sashimi/{sashimi_filename}")
else:
    Most_significance_sashimi_plot = None

#4 gene_fusion_circos_plot.png

fusion_images = {} 
root_dir1 = 'jaffal_fusion'
mix_categories1 = []
for subdir in os.listdir(root_dir1):
    subdir_path = os.path.join(root_dir1, subdir)
    if os.path.isdir(subdir_path) and subdir != '.ipynb_checkpoints':
        mix_categories1.append(subdir)

        
for mix in mix_categories1:
    fusion_images[mix] = image_to_base64(f"jaffal_fusion/{mix}/gene_fusion_circos_plot.png")

#4 which_to_show
fusion_plot_dir = './jaffal_fusion'
fusion_which_to_show_path = os.path.join(fusion_plot_dir, 'which_to_show')

fusion_filename = None
if os.path.exists(fusion_which_to_show_path):
    with open(fusion_which_to_show_path, 'r') as file:
        fusion_filename = file.readline().strip()  # 读取文件的第一行

if fusion_filename:
    fusion_plot = image_to_base64(f"jaffal_fusion/{fusion_filename}")
else:
    fusion_plot = None
    
html_template = """
<!DOCTYPE html>
<html lang="zh">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1">
    <title>Full-Length Transcriptome Analysis Report</title>
    <!-- Vega & Vega-Lite & Vega-Embed的依赖库 -->
    <script src="https://cdn.jsdelivr.net/npm/vega@5"></script>
    <script src="https://cdn.jsdelivr.net/npm/vega-lite@5"></script>
    <script src="https://cdn.jsdelivr.net/npm/vega-embed@6"></script>
    <!-- DataTables CSS -->
    <link rel="stylesheet" type="text/css" href="https://cdn.datatables.net/1.11.5/css/jquery.dataTables.css">
    <!-- jQuery -->
    <script src="https://code.jquery.com/jquery-3.5.1.js"></script>
    <!-- DataTables JS -->
    <script type="text/javascript" charset="utf8" src="https://cdn.datatables.net/1.11.5/js/jquery.dataTables.js"></script>
    <style>
        body {
            font-family: 'Open Sans', "Helvetica Neue", Helvetica, Arial, "PingFang SC", "Hiragino Sans GB", "Heiti SC", "Microsoft YaHei", "WenQuanYi Micro Hei", sans-serif;
            margin: 0;
            display: flex;
        }
        
        p {
        line-height: 1.5;  /* 使用倍数设置行间距 */
        }
        
        br {
            margin-bottom: 1.5em;  /* 设置与行间距相同的值 */
        }
        
        /* 固定导航栏 */
        .navigation-bar {
            height: 80px;  /* 固定高度，防止被Logo撑大 */
            background-color: #077C79;
            width: 100%;
            position: fixed; /* 固定在页面顶部 */
            top: 0;
            left: 0;
            z-index: 1000; /* 确保导航栏在最上层 */
            display: flex;
            justify-content: space-between; /* 将左侧和右侧内容分开 */
            align-items: center; /* 垂直居中对齐 */
            padding: 10px 40px; /* 左右内边距 */
            box-sizing: border-box;
        }
        
        .nav-left {
            color: white;
            font-size: 30px;
            font-weight: bold;
        }
        .nav-right {
            display: flex;
            justify-content: flex-end; /* 右对齐 */
            align-items: center; /* 垂直居中 */
            margin-right: 10px; /* 确保 Logo 不会完全贴到右边 */
        }

        .nav-right img {
            height: auto;
            width: 250px;
        }

        /* 侧边栏 */
        #toc-container {
            position: fixed;
            top: 80px;
            left: 0;
            height: calc(100vh - 100px); /* 侧边栏的高度与导航栏高度相关 */
            width: 250px;
            padding: 10px;
            background-color: white;  /* 设置白色背景 */
            box-shadow: 4px 0 8px rgba(0, 0, 0, 0.1);
            z-index: 1100;
            overflow-y: auto;
            overflow-x: hidden; /* 防止侧边栏出现横向滚动条 */
            display: flex;
            flex-direction: column; /* 让内容纵向排列 */
        }
        
        
        /* 侧边栏的整体布局 */
        #toc {
            min-width: 200px; /* 让菜单比容器宽 */
            padding: 0;
            display: block;
            flex-direction: column;
            white-space: nowrap; /* 避免换行 */
        }
        

        #toc ul {
            list-style-type:  none ;
            padding: 0;
            margin: 0;
            overflow-x: auto;  /* 启用水平滚动条 */
        }

        #toc li {
            margin: 5px 0;
            position: relative;
            list-style-type: none;
        }
        
        .list-unstyled {
            list-style: none !important;
            padding-left: 0 !important;
            margin-left: 0 !important;
        }
        .list-unstyled ul {
            margin-left: 20px !important;
        }

        #toc a {
            font-size: 12px;
            color: black;
            text-decoration: none;
            padding: 8px 8px;
            margin: 5px 0;
            display: block;
            white-space: nowrap;
            overflow-x: auto;
            
            transition: background-color 0.3s;
            width: 100%;
            background-color: transparent; /* 设置为透明 */
            border-radius: 10px;  /* 为标题添加圆角 */
            /*min-width: 100px; 设置最小宽度 */
            max-width: 300px; /* 设置最大宽度 */
        }

        #toc a.active {
            background-color: transparent;
            border-radius: 10px;
        }

        #toc a:hover {
            background-color: transparent;
            border-radius: 8px;
        }
        
        /* 正文部分 */
        .content {
            margin-top: 80px; /* 距离顶部导航栏 80px */
            margin-left: 280px; /* 避免与侧边栏重叠，确保内容从侧边栏右边开始 */
            padding: 20px;
            width: calc(100% - 280px); /* 使正文宽度适应侧边栏的宽度 */
            display: flex; /* 采用flex布局 */
            flex-direction: column; /* 将子项设置为纵向排列 */
            max-width: 100%; /* 防止超出容器宽度 */
            box-sizing: border-box; /* 确保 padding 不会增加总宽度 */
        }

        .content h2, .content h3 {
            margin-top: 20px;
        }
        
        /* 目录样式 */
        .toc-links {
            margin-bottom: 10px;
        }

        .toc-links ul {
            list-style-type: none;
            padding: 0;
            margin: 0;
        }

        .toc-links li::marker {
            content: "";
        }
        

        /* 一级标题样式 */
        .main-link {
            font-size: 14px;
            color: black;
            padding: 8px;
            display: block;
            cursor: pointer;
            text-decoration: none;
            white-space: nowrap;
            overflow: hidden;
            border-radius: 10px; /* 为一级标题添加圆角 */
        }
        
        .accordion-header {
            margin-bottom: 5px; /* 调整一级标题与内容之间的间距 */
        }
        
        .main-link:hover {
            background-color: transparent;
            border-radius: 8px;
        }

        /* 一级标题鼠标悬停效果 */
        .main-link:hover {
            background-color: transparent;
            border-radius: 8px;
        }

        /* 侧边栏容器样式 */
        .accordion-button {
            display: inline-flex; /* 让按钮和符号保持在同一行 */
            align-items: center; /* 垂直居中对齐符号和文本 */
            background-color: #eee;
            color: #444;
            cursor: pointer;
            padding: 8px;
            width: 100%;
            border: none;
            text-align: left;
            outline: none;
            font-size: 15px;
            font-weight: bold;
            transition: 0.4s;
            justify-content: flex-start; /* 确保文本和符号从左侧对齐 */
        }
        .submenu-3 {
            display: none; /* 默认隐藏三级标题 */
            padding-left: 15px; /* 缩进 */
        }
        
        .accordion-button.active, .accordion-button:hover {
          background-color: #ccc;
        }


        .accordion-button:after {
          content: '+'; /* 加号 */
          color: #777;
          font-weight: bold;
          float: left;
          margin-left: 0;
          margin-right: 3px;
          order: -1; /* 确保符号显示在文本的前面 */
        }

        .accordion-button.active:after {
          content: "-"; /* 减号 */
        }
        
        /* 精准调整3.4和3.5按钮，使其文字与普通li对齐 */
        .accordion-button.subsection-aligned { 
            position: relative;
            padding-left: 8px;
            font-weight: normal;            /* 去掉加粗 */
            background-color: transparent;  /* 默认背景透明 */
            box-shadow: none;               /* 防止额外的阴影效果 */
            color: #000000;
        }

        .accordion-button.subsection-aligned:after {
            position: absolute;
            left: -5px;
            top: 50%;
            transform: translateY(-50%);
        }
        
        .accordion-button.subsection-aligned:hover,
        .accordion-button.subsection-aligned.active {
            background-color: transparent;  /* 悬停和激活时也保持无背景 */
            font-weight: normal;            /* 防止悬停时加粗 */
        }
        
        .subsection-li {
            margin: 5px 0;
        }

        /* 手风琴面板样式 */
        .submenu {
          padding: 0 18px;
          background-color: white;
          display: block; /* 默认隐藏 */
          overflow: hidden;
          transition: all 0.3s ease;
          margin-bottom: 10px;
        }

        .submenu ul {
          list-style-type: none;
          padding: 0;
          margin: 0;
        }

        .submenu a {
          font-size: 14px;
          color: black;
          text-decoration: none;
          padding: 8px;
          display: block;
        }
        
        /* 表格部分，铺满整个宽度 */
        .table-section {
            width: 100%; /* 使表格容器铺满内容区域 */
            margin: 0 0 50px 0; /* 添加下边距 */
        }
        .table-section h2, .table-section h3 {
            text-align: center; /* 小标题左对齐 */
            margin-bottom: 20px;
        }
        .data-table {
            width: 100%; /* 确保表格铺满容器 */
            border-collapse: collapse;
            margin-bottom: 50px;
        }
        /* 表头样式 */
        .data-table th {
            background-color: #f2f2f2;
            text-align: center;
            padding: 8px;
            /* 移除边框 */
            border: none;
            /* 添加下边框以区分表头和表体 */
            border-bottom: 2px solid #dddddd;
        }
        /* 表格数据单元格样式 */
        .data-table td {
            text-align: center;
            padding: 8px;
            /* 移除边框 */
            border: none;
            /* 添加底部边框 */
            border-bottom: 1px solid #dddddd;
        }
        /* 交替行颜色 */
        .data-table tr:nth-child(even) {
            background-color: #f9f9f9; /* 偶数行灰色背景 */
        }

        .data-table tr:nth-child(odd) {
            background-color: #ffffff; /* 奇数行白色背景 */
        }

        /* 鼠标悬停行高亮（可选） */
        .data-table tr:hover {
            background-color: #eaeaea;
        }
        
        /* 图表部分 */
        .chart-section {
            max-width: 100%; /* 修改为100%以适应图片网格 */
            margin: 0 auto 50px; /* 居中并添加下边距 */
        }
        .chart-section h2 {
            text-align: left; /* 小标题左对齐 */
            margin-bottom: 20px;
        }
        .chart {
            display: flex;
            justify-content: center; /* 居中图表 */
            margin-bottom: 30px;
        }
        .chart div {
            width: 100%;
        }
        /* 图片容器 */
        .image-container {
            width: 100%;
            text-align: center;
            margin-bottom: 30px;
        }
        .image-container img {
            max-width: 100%;
            height: auto;
        }
        /* 图片网格样式 */
        .image-grid {
            display: flex;
            flex-wrap: wrap;
            justify-content: space-between;
            gap: 10px;
            padding: 0 20px; /* 增加左右内边距，调整图片和两边的距离 */
            max-width: 100%;
            align-items: center;
        }
        .image-row {
            display: flex;
            justify-content: center;
            margin-bottom: 20px;
        }
        
        /* 右移图片的效果（第一行） */
        .image-row:first-child .image-grid {
            margin-left: 40px; /* 向右偏移第一行的图片 */
        }
        
        .image-label {
            font-size: 16px;  /* 标签字体大小 */
            font-weight: bold;
            margin-right: 10px;  /* 标签与图片之间的间距 */
            white-space: nowrap; /* 防止标签换行 */
            display: flex;
            align-items: center; /* 垂直居中标签 */
        }
        
        .image-row h3 {
            width: 100%;
            text-align: left;
            margin-bottom: 10px;
        }
        .image-grid-container {
            display: flex;
            flex-wrap: wrap;
            justify-content: space-between;
        }
        .image-cell {
            flex: 1 1 calc(33.33% - 10px); /* 每个图片单元格占一行的50%，留有间隔 */
            max-width: calc(33.33% - 10px); /* 设置最大宽度 */
        }
        .image-cell img {
            width: 90%; /* 图片大小调小至容器的90% */
            height: auto;
            border: 1px solid #ddd;
            margin: 0 auto; /* 居中显示图片 */
            display: block; /* 确保图片是块级元素 */
        }
        
        .image-text-container {
            display: flex; /* 启用 Flexbox 布局 */
            align-items: center; /* 垂直居中 */
            justify-content: space-between; /* 左右对齐 */
            margin: 20px 0; /* 设置外边距 */
        }

        .image-half {
            flex: 1; /* 让图片占据左半部分 */
            padding-right: 20px; /* 右侧留白 */
            text-align: center; /* 图片居中 */
        }

        .image-half img {
            max-width: 100%; /* 图片最大宽度 */
            height: auto; /* 保持原始比例 */
            border: 1px solid #ddd; /* 添加边框 */
        }

        .text-half {
            flex: 1; /* 让文字占据右半部分 */
            text-align: left; /* 文字左对齐 */
        }
        
        /* 单个图片样式 */
        .single-image {
            width: 100%;
            text-align: center;
            margin: 20px 0;
        }
        .single-image img {
            max-width: 100%;
            height: auto;
            border: 1px solid #ddd;
        }
        
        .download-section {
            margin-top: 20px;
        }

        .download-btn {
            display: inline-block;
            padding: 10px 20px;
            background-color: #007bff;
            color: white;
            text-decoration: none;
            border-radius: 5px;
        }

        .download-btn:hover {
            background-color: #0056b3;
        }
        
        /* 给目标元素添加上边距 */
        #section-1, #section-2, #section-3, #section-4 {
            scroll-margin-top: 80px; /* 设置与导航栏高度相同的值 */
        }
        
        /* 给目标元素添加上边距，避免被导航栏遮挡 */
        #a0, #b0, #c0, #e0, 
        #a1, #a2, #a3, #b1, #b2, #b3, #b4, #b5, #b6, #b7, #b8, #b9, #c1, #e1, #e2 {
            scroll-margin-top: 90px;  /* 设置为导航栏的高度 */
        }
        
        
        /* 字符居中 */
        .center-bold-text {
            text-align: center;
            font-weight: bold;
            font-size: 18px; /* 设置字体大小 */
        }
        
        /* 边框 */
        .bordered {
            border: 1px solid #ddd;
            padding: 5px;
            border-radius: 10px; /* 设置圆角 */
        }

    </style>
</head>

<body>
    <!-- 导航栏 -->
    <div class="navigation-bar">
        <div class="nav-left">Full-Length Transcriptome Analysis Report</div>
        <div class="nav-right">
            <a href="https://www.cyclone-seq.com/" target="_blank" rel="noopener noreferrer">
                <img src="{{ logo_path }}" alt="Logo">
            </a>
        </div>
    </div>

    <!-- 侧边栏 -->
    <div id="toc-container">
        <div class="accordion-item">
            <button class="accordion-button">
                <a href="#e0" style="text-decoration: none; color: inherit;">Read QC</a>
            </button>
            <div class="submenu">
                <ul>
                    <li><a href="#e1">1.1 Read summary: passed reads & full-length reads</a></li>
                    <li><a href="#e2">1.2 Summary of full-length reads identification</a></li>
                </ul>
            </div>
        </div>

        <!-- Transcriptome Reconstruction -->
        <div class="accordion-item">
            <button class="accordion-button">
                <a href="#a0" style="text-decoration: none; color: inherit;">Transcriptome Reconstruction </a>
            </button>
            <div class="submenu">
                <ul>
                    <li><a href="#a1" class="toc-link">2.1 Alignment summary: genome alignment</a></li>
                    <li><a href="#a2" class="toc-link">2.2 Coverage depth on known transcripts</a></li>
                    <li><a href="#a3" class="toc-link">2.3 Transcriptome QC & annotation</a></li>
                </ul>
            </div>
        </div>

        <!-- Transcript Quantification -->
        <div class="accordion-item">
            <button class="accordion-button">
                <a href="#b0" style="text-decoration: none; color: inherit;">Transcript Quantification </a> 
            </button>
            <div class="submenu">
                <ul>
                    <li><a href="#b1" class="toc-link">3.1 Alignment summary: transcriptome alignment</a></li>
                    <li><a href="#b2" class="toc-link">3.2 Coverage on known transcripts and genes</a></li>
                    <li><a href="#b3" class="toc-link">3.3 TPM(Transcripts Per Million)</a></li>

                    <li class="subsection-li">
                        <button class="accordion-button subsection-aligned">
                            <span>3.4 Differential expression analysis</span>
                        </button>
                        <div class="submenu-3">
                            <ul class="list-unstyled">
                                <li><a href="#b5" class="toc-link">3.4.1 Differential gene expression</a></li>
                                <li><a href="#b6" class="toc-link">3.4.2 Differential transcript expression</a></li>
                            </ul>
                        </div>
                    </li>

                    <!-- 三级标题 -->
                    <li class="subsection-li">
                        <button class="accordion-button subsection-aligned">
                            <span>3.5 Differential splicing analysis</span>
                        </button>
                        <div class="submenu-3">
                            <ul class="list-unstyled">
                                <li><a href="#b8" class="toc-link">3.5.1 Differential splicing with local events</a></li>
                                <li><a href="#b9" class="toc-link">3.5.2 Differential transcript usage</a></li>
                            </ul>
                        </div>
                    </li>
                </ul>
            </div>
        </div>

        <!-- Gene Fusion Detection -->
        <div class="accordion-item">
            <button class="accordion-button">
                <a href="#c0" style="text-decoration: none; color: inherit;">Gene Fusion Detection </a>
            </button>
            <div class="submenu">
                <ul>
                     <li><a href="#c1" class="toc-link">4.1 Gene Fusions</a></li>
                </ul>
            </div>
        </div>
    </div>
    
    <!-- 正文 -->
    <div class="content">
        <div id="e0">
            <h2>1. Read QC</h2>
            <p>In third-generation full-length transcriptome sequencing, Read QC (Reads Quality Control) is a critical step to ensure data quality. By rigorously filtering low-quality, chimeric, and non-full-length reads, the accuracy and reliability of subsequent analyses can be significantly improved. </p>
             <div id="e1">
                <h3>1.1 Read summary: passed reads & full-length reads</h3>
             </div>
             
             <p>For both passed reads and full-length reads, length distribution charts, quality score distribution charts, and GC content distribution charts are displayed. The x-axis represents read length, read quality score, and GC content, respectively, while the y-axis represents the number of reads. </p>
            
            
        <!-- 下拉菜单 -->
        <label for="sample-select"> </label>
        <select id="sample-select" onchange="showImages(this)">
            <option value="">--Please choose a sample--</option>
            {% for mix in fl_reads_images.keys() | sort %}
            <option value="{{ mix }}" {% if loop.index == 1 %} selected {% endif %}>{{ mix }}</option>
            {% endfor %}
        </select>

        <!-- 图片显示容器 -->
        <div id="image-grid" style="margin-top: 20px;">
            <!-- 动态加载图片 -->
            {% for mix, images in fl_reads_images.items() %}
            <div id="images-{{ mix }}" class="image-group bordered" style="display: none;">
                <h3>{{ mix }}</h3>

                <!-- 第一行 - 图片和标签 "passed reads" -->
                <div class="image-row">
                    <div class="image-label">passed reads</div> <!-- 标签 -->
                    <div class="image-grid">
                        <div class="image-cell">
                            <img src="{{ images[0] }}" alt="{{ mix }} - {{ images[0].split('/')[-1] }}">
                        </div>
                        <div class="image-cell">
                            <img src="{{ images[1] }}" alt="{{ mix }} - {{ images[1].split('/')[-1] }}">
                        </div>
                        <div class="image-cell">
                            <img src="{{ images[2] }}" alt="{{ mix }} - {{ images[2].split('/')[-1] }}">
                        </div>
                    </div>
                </div>

                <!-- 第二行 - 图片和标签 "full-length reads" -->
                <div class="image-row">
                    <div class="image-label">full-length reads</div> <!-- 标签 -->
                    <div class="image-grid">
                        <div class="image-cell">
                            <img src="{{ images[3] }}" alt="{{ mix }} - {{ images[3].split('/')[-1] }}">
                        </div>
                        <div class="image-cell">
                            <img src="{{ images[4] }}" alt="{{ mix }} - {{ images[4].split('/')[-1] }}">
                        </div>
                        <div class="image-cell">
                            <img src="{{ images[5] }}" alt="{{ mix }} - {{ images[5].split('/')[-1] }}">
                        </div>
                    </div>
                </div>
            </div>
            {% endfor %}
        </div>
            <div id="e2">
                <h3>1.2 Summary of full-length reads identification</h3>
                <p>
                    <strong>Sample</strong>: Sample name. <br>
                    <strong>Total bases</strong>: Total number of sequenced bases. <br>
                    <strong>Valid bases</strong>: The number of bases that passed quality control. <br>
                    <strong>Valid bases (%)</strong>: Percentage of valid bases, calculated as Valid bases / Total bases × 100%. <br>
                    <strong>Length-filtered reads (%)</strong>: Percentage of reads that failed length filtering. <br>
                    <strong>QC-filtered reads (%)</strong>: Percentage of reads that failed quality filtering. <br>
                    <strong>Chimeric reads (%)</strong>: Percentage of chimeric reads. <br>
                    <strong>Full-length reads (%)</strong>: Percentage of full-length reads. <br>
                </p>
                <!-- glycine 表格 -->
                <div class="table-section bordered">
                    <!-- 添加跳转按钮 -->
                    <button onclick="alert('Please open the path in the file manager: fl_reads/ReadSummary.csv');" style="margin-bottom: 5px; padding: 5px; background-color: #077C79; color: white; border: none; border-radius: 5px; cursor: pointer;">
                        File link
                    </button>
                    {{ glycine_df_html | safe }}
                </div>
             </div>

    <!-- JavaScript 部分 -->
    <script type="text/javascript">

        function showImages(selectElement) {
            // 隐藏所有的图片组
            const imageGroups = document.querySelectorAll('.image-group');
            imageGroups.forEach(group => group.style.display = 'none');

            // 获取用户选择的样本
            const selectedSample = selectElement.value;

            // 如果选择了某个样本，显示对应的图片
            if (selectedSample) {
                const selectedImageGroup = document.getElementById(`images-${selectedSample}`);
                selectedImageGroup.style.display = 'block';
            }
        }
    </script>
            
    </div>

        <div id="a0">
            <h2>2. Transcriptome Reconstruction</h2>
            <p>Transcriptome reconstruction is a core step in the analysis of third-generation full-length transcriptome sequencing. It aims to assemble complete and accurate transcript structures from sequencing data, enabling the characterization of alternative splicing, transcription start sites (TSSs), transcription termination sites (TTSs), and more. This process facilitates a deeper understanding of gene expression complexity and provides a foundation for functional studies and molecular mechanism investigations. </p>

            <div id="a1">
                <h3>2.1 Alignment summary: genome alignment</h3>
                
                <p>
                    This module presents the alignment of transcripts to the reference genome to assess the quality of the sequencing data and the efficiency of the alignment.  <br>
                    <strong>Sample</strong>: Sample name.  <br>
                    <strong>PrimAlnPerc</strong>: Percentage of primary aligned reads, used to assess the quality of genome alignment.  <br>
                    <strong>MapPerc</strong>: Percentage of mapped reads, including primary, secondary, and supplementary alignments.  <br>
                    <strong>PrimAln</strong>: The number of primary aligned reads.  <br>
                    <strong>SecAln</strong>: The number of secondary aligned reads.  <br>
                    <strong>SupAln</strong>: The number of supplementary aligned reads.  <br>
                    <strong>Unmapped</strong>: The number of reads that were not aligned to the reference genome, which may be caused by low-quality sequences, sequencing errors, or incomplete genome assembly.  <br>
                    <strong>TotalReads</strong>: The number of full-length reads.  <br>
                    <strong>TotalRecords</strong>: The number of reads retained after alignments, including primary, secondary, and supplementary alignments.  <br>
                </p>
                
                <!-- genome_alignments_table_html 表格 -->
                <div class="table-section bordered">
                    <!-- 添加跳转按钮 -->
                    <button onclick="alert('Please open the path in the file manager: genome_alignments/read_aln_stats.csv');" style="margin-bottom: 5px; padding: 5px; background-color: #077C79; color: white; border: none; border-radius: 5px; cursor: pointer;">
                        File link
                    </button>
                    {{ genome_alignments_table_html | safe }}
                </div>
                
            </div>
            <div id="a2">
                <h3>2.2 Coverage depth on known transcripts</h3>
                
                
                <!-- 添加单个图片：Transcripts depth -->
                <div class="image-text-container bordered">
                    <div class="image-half">
                        <img src="{{ transcript_depth_image }}" alt="transcript_depth_image" style="width: 70%; height: auto;">
                    </div>
                    <div class="text-half">
                        <h3>This module displays the coverage depth distribution.</h3>
                        <p>
                            - <strong>X-axis</strong>: The region from 1kb upstream of the transcription start site (TSS) to 1kb downstream of the transcription termination site (TES). <br>
                            - <strong>Y-axis</strong>: The mean counts per million (CPM) mapped reads. <br>
                        </p>
                    </div> 
                </div>
            </div>
            
            <div id="a3">
                <h3>2.3 Transcriptome QC & annotation</h3>
                
                <div class="image-text-container bordered">
                    <div class="image-half">
                        <img src="{{ transcript_structure_image }}" alt="transcript_structure_image" style="width: 70%; height: auto;">
                    </div>
                    <div class="text-half">
                        <h3>This module annotates and performs quality control of the reconstructed transcripts to assess the structural features of the transcripts and their association with known genes.</h3>
                        <p>
                            The x-axis of the bar chart represents the structural category of the transcripts, including the following four types: <br>
                            - <strong>FSM (Full-splice match)</strong>: Fully matches the structure of known transcripts. <br>
                            - <strong>ISM (Incomplete-splice match)</strong>: Partially matches known transcripts, possibly with differences in splicing sites. <br>
                            - <strong>NIC (Novel in catalog)</strong>: Does not match known transcripts but overlaps with the genomic annotation, which may be a new transcript. <br>
                            - <strong>NNIC (Novel not in catalog)</strong>: Does not match known transcripts or genomic annotations, which may be a new transcript. <br>
                            The y-axis represents the number of transcripts for each structural category. Typically, FSM transcripts have the largest number, indicating that most transcripts are consistent with the known transcript structures. The number of ISM, NIC, and NNIC is relatively small, indicating the presence of new transcripts or differences in transcript structure. <br>
                        </p>
                    </div> 
                </div>
                
                
                <p>
                    The table lists the annotation information of the transcripts in detail.
                </p>
                <p>
                    <strong>isoform</strong>: The ID of the transcript. <br>
                    <strong>length</strong>: The length of the transcript (base pairs). <br>
                    <strong>exons</strong>: The number of exons in the transcript. <br>
                    <strong>structural_category</strong>: The structural category of the transcript. <br>
                    <strong>associated_gene</strong>: The gene associated with the transcript. <br>
                    <strong>associated_transcript</strong>: The known transcript associated with the transcript. <br>
                    <strong>diff_to_TSS</strong>: The distance between the TSS of this transcript and the TSS of the nearest known transcript. <br>
                    <strong>diff_to_TTS</strong>: The distance between the TES of this transcript and the TES of the nearest known transcript. <br>
                    <strong>dist_to_CAGE_peak</strong>: The distance of this transcript from the CAGE peak. CAGE (Cap Analysis of Gene Expression) is used to identify transcription start sites. <br>
                    <strong>within_CAGE_peak</strong>: Whether the TSS of this transcript is located inside the CAGE peak.
                </p>
                
                <!-- OUT.transcript_models_classification 表格 -->
                <div class="table-section bordered" style="overflow-x: auto; max-width: 100%;">
                    <!-- 添加跳转按钮 -->
                    <button onclick="alert('Please open the path in the file manager: sqanti_qc/OUT.transcript_models_classification.txt');" style="margin-bottom: 5px; padding: 5px; background-color: #077C79; color: white; border: none; border-radius: 5px; cursor: pointer;">
                        File link
                    </button>
                    {{ OUT_classification_df_part_html | safe }}
                </div>
            </div>
        </div>

        <div id="b0">
            <h2>3. Transcript Quantification</h2>
            
            <p>
                This section shows differential gene/transcript expression. Salmon was used to assign reads to annotated isoforms defined by the GTF-format annotation. These counts were used to perform a statistical analysis to identify the genes and isoforms that show differences in abundance between the experimental conditions.
            </p>
            
            <div id="b1">
                <h3>3.1 Alignment summary: transcriptome alignment</h3>
                
                <p>
                    Reference transcriptome alignment statistics. For the explanation of each statistical indicator, please refer to section 2.1.
                </p>
                
                <!-- genome_alignments_table 表格 -->
                <div class="table-section bordered">
                    <!-- 添加跳转按钮 -->
                    <button onclick="alert('Please open the path in the file manager: transcriptome_alignments/read_aln_stats.csv');" style="margin-bottom: 5px; padding: 5px; background-color: #077C79; color: white; border: none; border-radius: 5px; cursor: pointer;">
                        File link
                    </button>
                    {{ genome_alignments_table_html2 | safe }}
                </div> 
            </div>
            
            <div id="b2">
                <h3>3.2 Coverage on known transcripts and genes</h3>
                
                <div class="image-text-container bordered">
                    <div class="image-half">
                        <img src="{{ gene_transcript_cov_image }}" alt="gene_transcript_cov_image" style="width: 70%; height: auto;">
                    </div>
                    <div class="text-half">
                        <p>
                            - The statistics of gene coverage and transcript coverage for each sample against known genes (only protein-coding genes and lncRNA genes are counted). <br>
                            - Transcript coverage is the proportion of covered transcripts among all known transcripts; gene coverage is the proportion of covered genes among all known genes. <br>
                            - A transcript is considered covered if it has at least 80% coverage, and a gene is considered covered if at least one of its transcripts is covered. <br>
                        </p>
                    </div> 
                </div>

            </div>
            
            <div id="b3">
                <h3>3.3 TPM(Transcripts Per Million)</h3>
                <p>
                    Table showing the annotated Transcripts Per Million identified by Minimap2 mapping and Salmon transcript detection.
                </p>
                <!-- TPM.tsv 表格 -->
                <div class="table-section bordered">
                    <!-- 添加跳转按钮 -->
                    <button onclick="alert('Please open the path in the file manager: merge_counts/all_TPM.add_genename.tsv');" style="margin-bottom: 5px; padding: 5px; background-color: #077C79; color: white; border: none; border-radius: 5px; cursor: pointer;">
                        File link
                    </button>
                    {{ all_tpm_df_html | safe }}
                </div>
            </div>
            
            <div id="b4">
                <h3>3.4 Differential expression anlysis</h3>
                <div id="b5">
                    <h3>3.4.1 Differential gene expression</h3>
                    
                    <p>
                        Table showing the genes from the edgeR analysis. Information shown includes the log2 fold change between experimental conditions, the log-scaled counts per million measure of abundance and the FDR-corrected p-value (False discovery rate - Benjamini-Hochberg).
                    </p>

                    <p>
                        This table has not been filtered for genes that satisfy statistical or magnitudinal thresholds.
                    </p>
                    
                    <!-- DGE.tsv 表格 -->
                    <div class="table-section bordered">
                        <!-- 添加跳转按钮 -->
                        <button onclick="alert('Please open the path in the file manager: diff_exp/results_dge.tsv');" style="margin-bottom: 5px; padding: 5px; background-color: #077C79; color: white; border: none; border-radius: 5px; cursor: pointer;">
                            File link
                        </button>
                        {{ test_results_dge_html | safe }}
                    </div>
                    
                    <!-- 添加单个图片：results_dge.maplot.png -->
                    <div class="image-text-container bordered">
                        <div class="image-half">
                            <img src="{{ results_dge_image }}" alt="results_dge_image" style="width: 70%; height: auto;">
                        </div>
                        <div class="text-half">
                            <p>
                                The figure below presents the MA figure from this edgeR analysis, which visualises differences in measurements between the two experimental conditions. <br>
                                The y axis (M) is the log2 ratio of gene expression calculated between the conditions. The x axis (A) is a log2 transformed mean expression value.  <br>
                                Genes that satisfy the logFC and FDR-corrected (False discovery rate - Benjamini-Hochberg) p-value thresholds defined are shaded as 'Up-' or 'Down-' regulated. <br>
                            </p>
                        </div> 
                    </div>


                </div>
                
                <div id="b6">
                    <h3>3.4.2 Differential transcript expression</h3>
                    <p>
                        The table presents the results of differential transcript expression analysis obtained using edgeR. As previously mentioned, further elaboration is omitted here.
                    </p>
                    <!-- DTE.tsv 表格 -->
                    <div class="table-section bordered">
                        <!-- 添加跳转按钮 -->
                        <button onclick="alert('Please open the path in the file manager: diff_exp/results_dte.tsv');" style="margin-bottom: 5px; padding: 5px; background-color: #077C79; color: white; border: none; border-radius: 5px; cursor: pointer;">
                            File link
                        </button>
                        {{ test_results_dte_html | safe }}
                    </div>
                    
                    
                    <!-- 添加单个图片：results_dte.maplot.png -->
                    <div class="image-text-container bordered">
                        <div class="image-half">
                            <img src="{{ results_dte_image }}" alt="results_dte_image" style="width: 70%; height: auto;">
                        </div>
                        <div class="text-half">
                            <p>
                                The figure below presents the MA figure from this edgeR analysis, which visualises differences in measurements between the two experimental conditions. <br>
                                As previously mentioned, further elaboration is omitted here.
                            </p>
                        </div> 
                    </div>
                    
                </div>
            </div>
            
            <div id="b7">
                <h3>3.5 Differential splicing analysis</h3>
                
                <p>
                    To study splicing at the transcript isoform or at the local alternative splicing event level, across multiple conditions.  <br>
                    The local alternative splicing events are standard local splicing variations, refer to specific splicing patterns within localized regions of a gene, such as skipped exons (SE), retained introns (RI), etc., which involve only a small segment of the gene (e.g., a single exon or intron); whereas a transcript event is an isoform-centric approach, where each isoform in a gene is described separately, focusing on expression changes of entire transcripts (isoforms).  <br>
                </p>
                <div id="b8">
                    <h3>3.5.1 Differential splicing with local events (differential exon usage)</h3>
                    <p>
                        PSI file for local AS events. The first column is the event ID (transcript event or local AS event ID), the second column is gene name, and the following columns are the PSI values in each sample.   <br>
                        Percent Spliced In (PSI) is a quantitative measure used in SUPPA2 to estimate the inclusion level of alternative splicing events (e.g., skipped exons, retained introns, or alternative splice sites). It represents the proportion of transcripts that include a specific splicing event relative to all transcripts spanning that locus, which are normalized between 0 and 1 ([0,1]). <br>
                        PSI values range from 0 to 1, where: <br>
                        <strong>0</strong> indicates the splicing event is completely excluded (e.g., an exon is never included in mature transcripts).  <br>
                        <strong>1</strong> indicates the splicing event is fully included (e.g., an exon is always retained).  <br>
                    </p>
                    
                    <!-- events.add_genename.tsv 表格 -->
                    <div class="table-section bordered">
                        <!-- 添加跳转按钮 -->
                        <button onclick="alert('Please open the path in the file manager: diff_splice/events.add_genename.psi');" style="margin-bottom: 5px; padding: 5px; background-color: #077C79; color: white; border: none; border-radius: 5px; cursor: pointer;">
                            File link
                        </button>
                        {{ events_add_genename_df_html | safe }}
                    </div>
                    
                    <p>
                        The dpsi file (differential Percent Spliced In file) is a key output generated by SUPPA2 when performing differential alternative splicing analysis. It quantifies changes in splicing patterns between experimental conditions (e.g., control vs. treatment) by calculating the difference in the inclusion levels of alternative splicing events. <br>
                        The dpsi file contains the following columns: <br>
                        <strong>feature:</strong> Event ID, identifier for the alternative splicing event <br>
                        <strong>gene_name:</strong> Gene associated with the splicing event <br>
                        <strong>evnets.events_dPSI:</strong> Event PSI difference (ΔPSI) between Cond1 and Cond2 (ΔPSI = PSI_2 - PSI_1) <br>
                        <strong>evnets.events_pvalue:</strong> Significance of the difference of PSI between Cond1 and Cond2 <br>
                    </p>

                    <!-- suppa_diffSplice_event.add_genename.tsv 表格 -->
                    <div class="table-section bordered">
                        <!-- 添加跳转按钮 -->
                        <button onclick="alert('Please open the path in the file manager: diff_splice/suppa_diffSplice_event.add_genename.dpsi');" style="margin-bottom: 5px; padding: 5px; background-color: #077C79; color: white; border: none; border-radius: 5px; cursor: pointer;">
                            File link
                        </button>
                        {{ suppa_diffSplice_event_df_html | safe }}
                    </div>
                    

                    <!-- 添加单个图片：Most_significance_sashimi_plot -->
                    {% if Most_significance_sashimi_plot != None %}
                        <div class="image-text-container bordered">
                            <div class="image-half">
                                <img src="{{ Most_significance_sashimi_plot }}" alt="Most_significance_sashimi_plot" style="width: 70%; height: auto;">
                            </div>
                            <div class="text-half">
                                <p>
                                    The ggsashimi tool generates visualizations of sequencing data to display splicing patterns, read coverage, and junction reads across genomic regions. <br>
                                    Its output image integrates multiple layers of information to help researchers interpret alternative splicing events, transcript structures, and differences between experimental conditions.  <br>
                                </p>
                            </div> 
                        </div>
                    {% else %}
                        <p class="center-bold-text">There is no significant AS.</p>
                    {% endif %}
                </div>
                
                
                <div id="b9">
                    <h3>3.5.2 Differential transcript usage</h3>
                    
                    <p>
                        PSI file for transcript "events". As previously mentioned, further elaboration is omitted here.
                    </p>
                    
                    <!-- iso_isoform.add_genename.psi 表格 -->
                    <div class="table-section bordered">
                        <!-- 添加跳转按钮 -->
                        <button onclick="alert('Please open the path in the file manager: diff_splice/iso_isoform.add_genename.psi');" style="margin-bottom: 5px; padding: 5px; background-color: #077C79; color: white; border: none; border-radius: 5px; cursor: pointer;">
                            File link
                        </button>
                        {{ iso_isoform_df_html | safe }}
                    </div>

                    <p>
                        The dpsi file for transcript "events". As previously mentioned, further elaboration is omitted here.
                    </p>
                    
                    <!-- suppa_diffSplice_iso.add_genename.dpsi 表格 -->
                    <div class="table-section bordered">
                        <!-- 添加跳转按钮 -->
                        <button onclick="alert('Please open the path in the file manager: diff_splice/suppa_diffSplice_iso.add_genename.dpsi');" style="margin-bottom: 5px; padding: 5px; background-color: #077C79; color: white; border: none; border-radius: 5px; cursor: pointer;">
                            File link
                        </button>
                        {{ suppa_diffSplice_iso_df_html | safe }}
                    </div>
                </div>
                
            </div>
            
        </div>

        <div id="c0">
            <h2>4. Gene Fusion Detection</h2>
            <div id="c1">
                <h3>4.1 Gene Fusions</h3>
             </div>
            
            <p>Please refer to the contents in JAFFAL <a href="https://github.com/Oshlack/JAFFA/wiki/OutputDescription" target="_blank"><strong>OutputDescription</strong></a>.</p>
            
            <!-- jaffa_results 表格 -->
            <div class="table-section bordered" style="overflow-x: auto; max-width: 100%;">
                <!-- 添加跳转按钮 -->
                <button onclick="alert('Please open the path in the file manager: jaffal_fusion/jaffa_results.csv');" style="margin-bottom: 5px; padding: 5px; background-color: #077C79; color: white; border: none; border-radius: 5px; cursor: pointer;">
                    File link
                </button>
                {{ jaffa_results_df_html | safe }}
            </div>
            
            
            <!-- 下拉菜单 -->
            <label for="fusion-sample-select"></label>
            <select id="fusion-sample-select" onchange="showFusionImage(this)">
                <option value="">--Please choose a sample--</option>
                {% for mix in fusion_images.keys()  | sort %}
                <option value="{{ mix }}">{{ mix }}</option>
                {% endfor %}
            </select>
            
            <!-- 图片+文本排版容器 -->
            <div id="fusion-image-container" style="margin-top: 20px;">
                {% for mix, img_data in fusion_images.items()  %}
                <div id="fusion-{{ mix }}" class="fusion-image" style="display: none;">
                    <!-- 使用原来排版方式，图片+文本部分 -->
                    <div class="image-text-container bordered">
                        <div class="image-half">
                            <img src="{{ img_data }}" alt="{{ mix }} gene fusion plot" style="width: 70%; height: auto;">
                        </div>
                        <div class="text-half">
                            <p>
                                The Circos plot displays detected gene fusion events. Arcs connecting genomic regions represent fusions: Red arcs indicate potentially novel fusion events not previously reported (unknown); Blue arcs represent known, previously documented fusion events. The thickness of each arc is proportional to the number of supporting spanning reads identified by the JAFFAL analysis, with more reads resulting in a thicker arc. <br>
                            </p>
                        </div> 
                    </div>
                </div>
                {% endfor %}
            </div>

            <script type="text/javascript">
                function showFusionImage(selectElement) {
                    // 隐藏所有图像
                    const allImages = document.querySelectorAll('.fusion-image');
                    allImages.forEach(image => {
                        image.style.display = 'none';
                    });

                    // 获取用户选择的样本
                    const selectedSample = selectElement.value;

                    // 如果选择了某个样本，显示对应的图片
                    if (selectedSample) {
                        const selectedImage = document.getElementById('fusion-' + selectedSample);
                        if (selectedImage) {
                            selectedImage.style.display = 'block';
                        }
                    }
                }

                // 页面加载自动显示第一个样本（可选）
                    window.onload = function() {
                    // 1.1 默认选择第一个样本
                    const qcSelect = document.getElementById("sample-select");
                    if (qcSelect && qcSelect.options.length > 1) {
                        qcSelect.selectedIndex = 1;
                        showImages(qcSelect);
                    }

                    // 4.1 默认选择第一个样本
                    const fusionSelect = document.getElementById("fusion-sample-select");
                    if (fusionSelect && fusionSelect.options.length > 1) {
                        fusionSelect.selectedIndex = 1;
                        showFusionImage(fusionSelect);
                    }
                };
            </script>
            
            
            
            
            
            <p>
                
            </p>
            <!-- 添加单个图片：single fusion_plot -->
            {% if fusion_plot != None %}
                <div class="image-text-container bordered">
                    <div class="image-half">
                        <img src="{{ fusion_plot }}" alt="single fusion_plot" style="width: 70%; height: auto;">
                    </div>
                    <div class="text-half">
                        <p>
                            The figure below shows the position of the most significant fusion gene pair on the chromosome, along with the number of supporting reads. <br>
                        </p>
                    </div> 
                </div>
            {% else %}
                <p class="center-bold-text">There is no significant gene fusion enent.</p>
            {% endif %}
            
            
        </div>
    </div>
    
    <script type="text/javascript">
        // 页面加载完成后执行
        $(document).ready(function() {
            // 初始化所有 class 为 "data-table" 的表格
            $('.data-table').each(function() {
                // 如果当前表格有 no-pagination 类，则禁用分页
                if ($(this).hasClass('no-pagination')) {
                    $(this).DataTable({
                        "paging": false,  // 禁用分页
                        "info": false,  // 禁用显示 "Showing x to y of z entries" 等信息
                        "lengthChange": false  // 禁用显示 "Show entries" 控件
                    });
                } else {
                    // 否则启用默认分页
                    $(this).DataTable({
                        "paging": true  // 启用分页
                    });
                }
            });
        });
    </script>
    
    <script>
        // 获取所有折叠按钮
        var acc = document.getElementsByClassName("accordion-button");

        // 遍历所有按钮，绑定点击事件
        for (var i = 0; i < acc.length; i++) {
            acc[i].addEventListener("click", function() {
                this.classList.toggle("active"); // 切换按钮状态（加号/减号）

                // 获取当前点击按钮后面的面板
                var submenu = this.nextElementSibling;

                // 切换面板的显示/隐藏
                if (submenu.style.display === "block") {
                    submenu.style.display = "none"; // 隐藏面板
                } else {
                    submenu.style.display = "block"; // 展开面板
                }
            });
        }
    </script>

    
</body>
</html>
"""

# 选择是否使用Data URI
use_data_uri = True  # 设置为True以使用Data URI嵌入Logo

# 定义Logo的文件路径（使用SVG）
logo_path_file = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'logo/logoEN.svg')  # 使用指定的Logo文件路径

if use_data_uri:
    try:
        with open(logo_path_file, 'rb') as image_file:
            encoded_string = base64.b64encode(image_file.read()).decode('utf-8')

        # 生成Data URI，确保MIME类型与Logo格式匹配
        # 对于SVG格式
        logo_data_uri = f"data:image/svg+xml;base64,{encoded_string}"

    except FileNotFoundError:
        print(f"错误：未找到Logo文件 '{logo_path_file}'。请确保文件存在并路径正确。")
        logo_data_uri = logo_path_file  # 使用文件路径作为备选
else:
    # 直接使用文件路径
    logo_data_uri = logo_path_file  # 使用指定的Logo文件路径
    
    
env = Environment(loader=FileSystemLoader('.'))
template = env.from_string(html_template)
final_html = template.render(logo_path=logo_data_uri,
fl_reads_images=fl_reads_images,
glycine_df_html=glycine_df_html,
genome_alignments_table_html=genome_alignments_table_html,
transcript_depth_image=transcript_depth_image,
transcript_structure_image=transcript_structure_image,
OUT_classification_df_part_html=OUT_classification_df_part_html,
genome_alignments_table_html2=genome_alignments_table_html2,
gene_transcript_cov_image=gene_transcript_cov_image,
all_tpm_df_html=all_tpm_df_html,
test_results_dge_html=test_results_dge_html,
results_dge_image=results_dge_image,
test_results_dte_html=test_results_dte_html,
results_dte_image=results_dte_image,
events_add_genename_df_html=events_add_genename_df_html,
suppa_diffSplice_event_df_html=suppa_diffSplice_event_df_html,
iso_isoform_df_html=iso_isoform_df_html,
suppa_diffSplice_iso_df_html=suppa_diffSplice_iso_df_html,
Most_significance_sashimi_plot=Most_significance_sashimi_plot,
jaffa_results_df_html=jaffa_results_df_html,
fusion_plot=fusion_plot,
fusion_images=fusion_images)


output_file = 'report_files/Report_EN.html'  # 修改为您想要的文件名
with open(output_file, 'w', encoding='utf-8') as f:
    f.write(final_html)
