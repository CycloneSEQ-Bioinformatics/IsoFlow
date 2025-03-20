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
genome_alignments_table_html  = merged_df.to_html(index=False, classes='data-table no-pagination', border=0)
#print(genome_alignments_table_html)

# 2.3表格
OUT_classification_df = pd.read_csv("sqanti_qc/OUT.transcript_models_classification.txt", sep = "\t")
OUT_classification_df_part = OUT_classification_df[["isoform", "length", "exons", "structural_category", "associated_gene", "associated_transcript", "diff_to_TSS", "diff_to_TTS", "dist_to_CAGE_peak", "within_CAGE_peak"]].head(100)
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
genome_alignments_table_html2  = merged_df2.to_html(index=False, classes='data-table no-pagination', border=0)


#3.4 DGE Table
test_results_dge = pd.read_csv("diff_exp/results_dge.tsv",sep = "\t")
test_results_dge_1k = test_results_dge.head(100)
test_results_dge_html = test_results_dge_1k.to_html(index=False, classes='data-table', border=0)

#3.4 DTE Table
test_results_dte = pd.read_csv("diff_exp/results_dte.tsv",sep = "\t")
test_results_dte_1k = test_results_dte.head(100)
test_results_dte_html = test_results_dte_1k.to_html(index=False, classes='data-table', border=0)

#3.5 events.add_genename.psi
events_add_genename_df = pd.read_csv("diff_splice/events.add_genename.psi", sep= "\t").head(100)
events_add_genename_df_html = events_add_genename_df.to_html(index=False, classes='data-table', border=0)

#3.5 suppa_diffSplice_event.add_genename.dpsi
suppa_diffSplice_event_df = pd.read_csv("diff_splice/suppa_diffSplice_event.add_genename.dpsi", sep= "\t").head(100)
suppa_diffSplice_event_df_html = suppa_diffSplice_event_df.to_html(index=False, classes='data-table', border=0)

#3.5 iso_isoform.add_genename.psi
iso_isoform_df = pd.read_csv("diff_splice/iso_isoform.add_genename.psi", sep= "\t").head(100)
iso_isoform_df_html = iso_isoform_df.to_html(index=False, classes='data-table', border=0)

#3.5 suppa_diffSplice_iso.add_genename.dpsi
suppa_diffSplice_iso_df = pd.read_csv("diff_splice/suppa_diffSplice_iso.add_genename.dpsi", sep= "\t").head(100)
suppa_diffSplice_iso_df_html = suppa_diffSplice_iso_df.to_html(index=False, classes='data-table', border=0)

#jaffa_results
jaffa_results_df = pd.read_csv("jaffal_fusion/jaffa_results.csv", sep= ",").head(100)
jaffa_results_df_html = jaffa_results_df.to_html(index=False, classes='data-table', border=0)


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
all_tpm_df = pd.read_csv("merge_counts/all_TPM.add_genename.tsv", sep = "\t").head(100)
all_tpm_df_html = all_tpm_df.to_html(index=False, classes='data-table', border=0)

#3.4  DGE results_dge.maplot.png
results_dge_image = image_to_base64("diff_exp/results_dge.maplot.png")


#3.4 DTE results_dge.maplot.png
results_dte_image = image_to_base64("diff_exp/results_dte.maplot.png")


#3.5 sashimi plot 判断展示
plot_sashimi_dir = './diff_splice/plot_sashimi'
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

#4 span_reads_heatmap.png

fusion_heatmap = image_to_base64("jaffal_fusion/span_reads_heatmap.png")

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
    <title>全长转录组生信分析报告</title>
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
            height: 80px;
            width: auto;
            max-height: 100%; /* 确保不会超过导航栏 */
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
            justify-content: flex-start;
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
        <div class="nav-left">全长转录组生信分析报告</div>
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
                <a href="#e0" style="text-decoration: none; color: inherit;">下机序列质控</a>
            </button>
            <div class="submenu">
                <ul>
                    <li><a href="#e1">1.1 通过序列和全长序列</a></li>
                    <li><a href="#e2">1.2 全长序列鉴定</a></li>
                </ul>
            </div>
        </div>

        <!-- Transcriptome Reconstruction -->
        <div class="accordion-item">
            <button class="accordion-button">
                <a href="#a0" style="text-decoration: none; color: inherit;">转录组重建 </a>
            </button>
            <div class="submenu">
                <ul>
                    <li><a href="#a1" class="toc-link">2.1 基因组比对</a></li>
                    <li><a href="#a2" class="toc-link">2.2 已知转录本覆盖深度</a></li>
                    <li><a href="#a3" class="toc-link">2.3 转录组注释及质控</a></li>
                </ul>
            </div>
        </div>

        <!-- Transcript Quantification -->
        <div class="accordion-item">
            <button class="accordion-button">
                <a href="#b0" style="text-decoration: none; color: inherit;">转录本定量 </a> 
            </button>
            <div class="submenu">
                <ul>
                    <li><a href="#b1" class="toc-link">3.1 转录组比对</a></li>
                    <li><a href="#b2" class="toc-link">3.2 已知转录本及基因覆盖度</a></li>
                    <li><a href="#b3" class="toc-link">3.3 TPM(Transcripts Per Million)</a></li>
                    <!-- 三级标题 -->
                    <button class="accordion-button">
                        <a href="#b4" style="text-decoration: none; color: inherit;">3.4 差异表达分析 </a>
                    </button>
                    <div class="submenu-3">
                        <ul class="list-unstyled">
                            <li><a href="#b5" class="toc-link">3.4.1 差异基因表达</a></li>
                            <li><a href="#b6" class="toc-link">3.4.2 差异转录本表达</a></li>
                        </ul>
                    </div>

                    <!-- 三级标题 -->
                    <button class="accordion-button">
                        <a href="#b7" style="text-decoration: none; color: inherit;">3.5 可变剪切分析 </a>
                    </button>
                    <div class="submenu-3">
                        <ul class="list-unstyled">
                            <li><a href="#b8" class="toc-link">3.5.1 局部事件差异剪接</a></li>
                            <li><a href="#b9" class="toc-link">3.5.2 差异转录本使用</a></li>
                        </ul>
                    </div>
                </ul>
            </div>
        </div>

        <!-- Gene Fusion Detection -->
        <div class="accordion-item">
            <button class="accordion-button">
                <a href="#c0" style="text-decoration: none; color: inherit;">融合基因检测 </a>
            </button>
            <div class="submenu">
                <ul>
                     <li><a href="#c1" class="toc-link">4.1 基因融合</a></li>
                </ul>
            </div>
        </div>
    </div>
    
    <!-- 正文 -->
    <div class="content">
        <div id="e0">
            <h2>1.下机序列质控</h2>
            <p>在三代全长转录本测序中，Read QC（Reads质量控制）是确保数据质量的关键步骤。通过对测序reads进行严格的质量控制，可以有效过滤低质量、嵌合及非全长的reads，从而提高后续分析的准确性和可靠性。 </p>
             <div id="e1">
                <h3>1.1 通过序列和全长序列</h3>
             </div>
             
             <p>对passed reads和full-length reads分别展示了长度分布图、质量值分布图和GC含量分布图。其中横坐标分别表示reads长度、reads质量值和GC含量，纵坐标表示reads数量。 </p>
            
            
            <!-- 下拉菜单 -->
        <label for="sample-select"> </label>
        <select id="sample-select" onchange="showImages(this)">
            <option value="">--请选择样本--</option>
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
                <h3>1.2 全长序列鉴定</h3>
                <p>
                    <strong>Sample</strong>: 样本名称。 <br>
                    <strong>Total bases</strong>: 测序总碱基数。 <br>
                    <strong>Valid bases</strong>: 有效碱基数。 <br>
                    <strong>Valid bases (%)</strong>: 有效碱基占比，计算公式为Valid bases / Total bases × 100%。 <br>
                    <strong>Length-filtered reads (%)</strong>: 因长度不符合标准而被过滤的reads占总reads的比例。 <br>
                    <strong>QC-filtered reads (%)</strong>: 因质量值不符合标准而被过滤的reads占总reads的比例。 <br>
                    <strong>Chimeric reads (%)</strong>: 嵌合reads占总reads的比例。 <br>
                    <strong>Full-length reads (%)</strong>: 全长reads占总reads的比例。 <br>
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
        // 页面加载时自动显示第一个样本的图片
        window.onload = function() {
            // 默认显示第一个样本的信息
            const selectElement = document.getElementById("sample-select");
            selectElement.value = selectElement.options[1].value; // 设置为第一个样本的值
            showImages(selectElement);  // 调用函数显示对应的图片
        };

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
            <h2>2.转录组重建</h2>
            <p>转录组重建（Transcriptome Reconstruction）是三代全长转录组测序数据分析的核心环节之一，旨在从测序数据中组装出完整、准确的转录本结构，进而解析基因的可变剪接、转录起始位点（TSS）、转录终止位点（TTS）等信息。通过转录本重建，可以深入了解基因表达的复杂性，为功能研究和分子机制解析提供基础。 </p>


            <div id="a1">
                <h3>2.1 基因组比对</h3>
                
                <p>
                    该模块展示了转录本比对到参考基因组的情况，以评估测序数据的质量和比对的效率。 <br>
                    <strong>Sample</strong>: 样本名称.  <br>
                    <strong>PrimAlnPerc</strong>: 主要比对reads占总reads的比例，反映了测序数据的比对质量.  <br>
                    <strong>MapPerc</strong>: 比对reads占总reads的比例，包含主要比对、次要比对和补充比对.  <br>
                    <strong>PrimAln</strong>: 主要比对reads的数量.  <br>
                    <strong>SecAln</strong>: 次要比对reads的数量.  <br>
                    <strong>SupAln</strong>: 补充比对reads的数量.  <br>
                    <strong>Unmapped</strong>: 未比对到参考基因组的reads数量，可能由低质量序列、测序错误或基因组组装不完整导致.  <br>
                    <strong>TotalReads</strong>: 全长reads的数量.  <br>
                    <strong>TotalRecords</strong>: 比对后保留的reads数，包括主要比对、次要比对和补充比对.  <br>
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
                <h3>2.2 已知转录本覆盖深度</h3>
                
                
                <!-- 添加单个图片：Transcripts depth -->
                <div class="image-text-container bordered">
                    <div class="image-half">
                        <img src="{{ transcript_depth_image }}" alt="transcript_depth_image" style="width: 70%; height: auto;">
                    </div>
                    <div class="text-half">
                        <h3>该模块展示了已知转录本的覆盖深度分布。</h3>
                        <p>
                            - <strong>横坐标</strong>: 表示从转录起始位点（TSS）上游1kb到转录终止位点（TES）下游1kb的区域。 <br>
                            - <strong>纵坐标</strong>: 表示平均每百万比对reads数（Mean CPM，Counts Per Million）。 <br>
                        </p>
                    </div> 
                </div>
            </div>
            
            <div id="a3">
                <h3>2.3 转录组注释及质控</h3>
                
                <div class="image-text-container bordered">
                    <div class="image-half">
                        <img src="{{ transcript_structure_image }}" alt="transcript_structure_image" style="width: 70%; height: auto;">
                    </div>
                    <div class="text-half">
                        <h3>该模块对重建的转录本进行注释和质量控制，以评估转录本的结构特征及其与已知基因的关联。</h3>
                        <p>
                            柱状图的横坐标表示转录本的结构类型（Structural category），包括以下四种： <br>
                            - <strong>FSM (Full-splice match)</strong>: 完全匹配已知转录本的结构。 <br>
                            - <strong>ISM (Incomplete-splice match)</strong>: 与已知转录本部分匹配，可能存在剪接位点差异。 <br>
                            - <strong>NIC (Novel in catalog)</strong>: 与已知转录本不匹配，但与基因组注释有重叠，可能是新的转录本。 <br>
                            - <strong>NNIC (Novel not in catalog)</strong>: 与已知转录本和基因组注释均不匹配，可能是新的转录本。 <br>
                            纵坐标表示每种结构类型的转录本数量。通常情况下，FSM的转录本数量最多，表明大部分转录本与已知转录本相符。ISM、NIC和NNIC的数量相对较少，表示存在新的转录本或者转录本结构差异。 <br>
                        </p>
                    </div> 
                </div>
                
                
                <p>
                    表格详细列出了转录本的注释信息。
                </p>
                <p>
                    <strong>isoform</strong>: 转录本的唯一标识符。 <br>
                    <strong>length</strong>: 转录本的长度（bp）。 <br>
                    <strong>exons</strong>: 转录本的外显子数量。 <br>
                    <strong>structural_category</strong>: 转录本的结构类型。 <br>
                    <strong>associated_gene</strong>: 与该转录本关联的基因。 <br>
                    <strong>associated_transcript</strong>: 与该转录本关联的已知转录本。 <br>
                    <strong>diff_to_TSS</strong>: 该转录本TSS与最近的已知转录本TSS的距离。 <br>
                    <strong>diff_to_TTS</strong>: 该转录本TES与最近的已知转录本TES的距离。 <br>
                    <strong>dist_to_CAGE_peak</strong>: 该转录本与CAGE峰的距离，用于鉴定转录起始位点。 <br>
                    <strong>within_CAGE_peak</strong>: 该转录本的TSS是否位于CAGE峰内部。
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
            <h2>3. 转录本定量</h2>
            
            <p>
                本节展示了差异基因/转录本表达分析结果。使用 Salmon 将测序读段分配到由GTF格式注释定义的已知转录本上。基于这些计数数据，通过统计分析识别出在不同实验条件之间表达量存在显著差异的基因和转录本。
            </p>
            
            <div id="b1">
                <h3>3.1 转录组比对</h3>
                
                <p>
                    比对到参考转录组的比对结果统计。各项统计指标的说明参考2.1部分。
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
                <h3>3.2 已知转录本及基因覆盖度</h3>
                
                <div class="image-text-container bordered">
                    <div class="image-half">
                        <img src="{{ gene_transcript_cov_image }}" alt="gene_transcript_cov_image" style="width: 70%; height: auto;">
                    </div>
                    <div class="text-half">
                        <p>
                            - 各样本对已知基因（仅统计蛋白编码基因和lncRNA基因）的基因覆盖度和转录本覆盖度的统计。 <br>
                            - 转录本覆盖度为，被覆盖的转录本在所有已知转录本的占比；基因覆盖度为，被覆盖的基因在所有已知基因的占比。 <br>
                            - 转录本被覆盖的定义为一个转录本至少有80%的覆盖率，基因被覆盖的定义为该基因至少有一个转录本被覆盖。 <br>
                        </p>
                    </div> 
                </div>

            </div>
            
            <div id="b3">
                <h3>3.3 TPM(Transcripts Per Million)</h3>
                <p>
                    表格中展示了已知转录本的表达定量结果，以TPM作为标准化定量指标。
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
                <h3>3.4 差异表达分析</h3>
                <div id="b5">
                    <h3>3.4.1 差异基因表达</h3>
                    
                    <p>
                        表格展示了使用edgeR得到的差异基因表达分析结果。所展示的信息包括实验条件之间的表达水平的log2倍数变化、按logCPM的丰度度量以及FDR校正后的p值（错误发现率 - Benjamini-Hochberg法）。
                    </p>

                    <p>
                        此表格尚未针对满足统计或量级阈值的基因进行过滤
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
                                此图为本次edgeR分析中的MA图，展示了两组实验条件之间基因表达水平的差异。 <br>
                                纵轴（M）表示两组条件间基因表达的log2比率。  <br>
                                横轴（A）表示经过log2转换的平均表达值。满足定义的logFC和FDR校正p值阈值的基因被标记为“上调”或“下调”。 <br>
                            </p>
                        </div> 
                    </div>


                </div>
                
                <div id="b6">
                    <h3>3.4.2 差异转录本表达</h3>
                    <p>
                        表格展示了使用edgeR得到的差异转录本表达分析结果。与上文相同，不再赘述。
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
                                此图为本次edgeR分析中的MA图，展示了两组实验条件之间基因表达水平的差异。与上文相同，不再赘述。
                            </p>
                        </div> 
                    </div>
                    
                </div>
            </div>
            
            <div id="b7">
                <h3>3.5 可变剪切分析</h3>
                
                <p>
                    差异可变剪切分析在两个维度上进行：局部可变剪接事件，指基因中特定区域的可变剪接模式，如外显子跳跃（Skipped Exons, SE）、内含子保留（Retained Introns, RI）等，仅涉及基因的某一小段区域（如单个外显子或内含子）；转录本水平事件，关注整个转录本（isoform）的表达变化。  <br>
                </p>
                <div id="b8">
                    <h3>3.5.1 局部事件差异剪接</h3>
                    <p>
                        PSI文件（针对局部可变剪接事件）。该文件的第一列为事件ID（转录本事件或局部可变剪接事件ID），第二列为基因名称，后续各列为每个样本中的PSI值。  <br>
                        剪接百分比（PSI）， 是SUPPA2中使用的一种定量指标，用于估计可变剪接事件（例如外显子跳跃、内含子保留或可变剪接位点）的包含水平。它表示包含特定剪接事件的转录本占该基因位点所有转录本的比例，其值经过标准化处理，范围为 0到1（[0,1]）。 <br>
                        具体含义如下：<br>
                        <strong>0</strong> 表示该剪接事件被完全排除（例如，某个外显子从未被包含在成熟转录本中）。  <br>
                        <strong>1</strong> 表示该剪接事件被完全包含（例如，某个外显子始终被保留）。 <br>
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
                        dpsi 文件是 SUPPA2 在进行差异可变剪接分析时生成的关键输出文件。它通过计算可变剪接事件保留水平的差异，量化了实验条件之间（例如对照组 vs. 处理组）剪接模式的变化。<br>
                        dpsi 文件包含以下列： <br>
                        <strong>feature:</strong> 事件ID，可变剪接事件的标识符。 <br>
                        <strong>gene_name:</strong> 与剪接事件相关的基因名称。 <br>
                        <strong>evnets.events_dPSI:</strong> 条件1和条件2之间的事件PSI差异（ΔPSI = PSI_2 - PSI_1）。 <br>
                        <strong>evnets.events_pvalue:</strong> 条件1和条件2之间PSI差异的显著性。 <br>
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
                                    ggsashimi能够通过测序数据生成可视化图像，展示基因组区域内剪接模式、读段覆盖度及连接读段的分布。 <br>
                                    其输出图像整合了多层信息，帮助研究人员解析可变剪接事件、转录本结构以及实验条件之间的差异。  <br>
                                </p>
                            </div> 
                        </div>
                    {% else %}
                        <p class="center-bold-text">There is no significant AS.</p>
                    {% endif %}
                </div>
                
                
                <div id="b9">
                    <h3>3.5.2 差异转录本使用</h3>
                    
                    <p>
                        PSI文件（针对转录本水平事件）。与上文相同，不再赘述。
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
                        dPSI文件（针对转录本水平事件）。与上文相同，不再赘述。
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
            <h2>4. 融合基因检测</h2>
            <div id="c1">
                <h3>4.1 基因融合</h3>
             </div>
            
            <p>请参考JAFFAL官方 <a href="https://github.com/Oshlack/JAFFA/wiki/OutputDescription" target="_blank"><strong>OutputDescription</strong></a> 内容。</p>
            
            <!-- jaffa_results 表格 -->
            <div class="table-section bordered" style="overflow-x: auto; max-width: 100%;">
                <!-- 添加跳转按钮 -->
                <button onclick="alert('Please open the path in the file manager: jaffal_fusion/jaffa_results.csv');" style="margin-bottom: 5px; padding: 5px; background-color: #077C79; color: white; border: none; border-radius: 5px; cursor: pointer;">
                    File link
                </button>
                {{ jaffa_results_df_html | safe }}
            </div>
            
            
            <!-- 添加单个图片：results_dge.maplot.png -->
            <div class="image-text-container bordered">
                <div class="image-half">
                    <img src="{{ fusion_heatmap }}" alt="fusion_heatmap" style="width: 70%; height: auto;">
                </div>
                <div class="text-half">
                    <p>
                        下图展示了高度显著的融合基因的热图,x轴代表样本ID,y轴代表融合基因对,颜色深度反映支持该事件的read数量。 <br>
                    </p>
                </div> 
            </div>
            
            
            
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
                            下图显示了染色体上最显著的融合基因对的位置，以及支持该事件的read数量。 <br>
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
logo_path_file = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'logo/logoCN.svg')  # 使用指定的Logo文件路径

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
fusion_heatmap=fusion_heatmap,
fusion_plot=fusion_plot)


output_file = 'report_files/Report_CN.html'  # 修改为您想要的文件名
with open(output_file, 'w', encoding='utf-8') as f:
    f.write(final_html)
