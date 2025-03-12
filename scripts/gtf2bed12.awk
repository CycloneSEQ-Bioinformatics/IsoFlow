BEGIN {FS="\t";OFS="\t"} 
# 提取外显子信息
$3 == "exon" {
    # 保存转录本的信息（染色体、起始、结束、转录本ID、外显子信息）
    exon_start=length(exon_start)>0 ? exon_start "," $4-1 : $4-1; # 保存外显子起始位置
    exon_end=length(exon_end)>0 ? exon_end "," $5 : $5; # 保存外显子结束位置
}

# 提取转录本信息（仅获取转录本的起始和结束）
$3 == "transcript" {
    if (transcript_id) {
        # 获取外显子信息
        split(exon_start, starts, ",");
        split(exon_end, ends, ",");
        
        block_count=length(starts);  # 外显子数量
        block_sizes="";
        block_starts="";
        # 计算外显子的大小和偏移
        for (i=1; i<=block_count; i++) {
            block_sizes=block_sizes (ends[i]-starts[i])",";  # 外显子的大小
            block_starts=block_starts (starts[i]-start)",";  # 外显子的起始偏移量
        }
        # 输出 BED12 格式
        print chrom, start, end, transcript_id, 0, strand, start, end, "0,0,0", block_count, block_sizes, block_starts;
    }
    if (match($9, /transcript_id "([^"]+)"/)) {
        transcript_id=substr($9, RSTART+15, RLENGTH-16);  # 提取转录本ID
    }
    chrom = $1;
    start = $4 - 1;  # 转录本起始位置（0-based）
    end = $5;  # 转录本结束位置（1-based）
    strand = $7;  # 转录本链方向
    exon_start = "";
    exon_end = "";
}

END{
    split(exon_start, starts, ",");
    split(exon_end, ends, ",");
    block_count=length(starts);
    block_sizes="";
    block_starts="";
    for (i=1; i<=block_count; i++) {
        block_sizes=block_sizes (ends[i]-starts[i])",";
	block_starts=block_starts (starts[i]-start)",";
    }
    print chrom, start, end, transcript_id, 0, strand, start, end, "0,0,0", block_count, block_sizes, block_starts;
}
