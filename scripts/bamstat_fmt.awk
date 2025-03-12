BEGIN {
    FS="\t";
    total_records = 0;
    total_reads = 0;
    primary_mapped = 0;
    primary_perc = 0;
    secondary = 0;
    supplementary = 0;
    mapped = 0;
    map_perc = 0;
}

{
    if ($3 ~ /total/) {total_records = $1;}
    if ($3 ~ /primary$/) {total_reads = $1;}
    if ($3 ~ /primary mapped$/) {primary_mapped = $1;}
    if ($3 ~ /primary mapped %$/) {primary_perc = $1;}
    if ($3 ~ /secondary$/) {secondary = $1;}
    if ($3 ~ /supplementary$/) {supplementary = $1;}
    if ($3 ~ /^mapped$/) {mapped = $1;}
    if ($3 ~ /^mapped %$/) {map_perc = $1;}
}

END {
    unmapped = total_records - mapped;
    printf "PrimAlnPerc\tMapPerc\tPrimAln\tSecAln\tSupAln\tUnmapped\tTotalReads\tTotalRecords\n"
    printf "%.2f\t%.2f\t%d\t%d\t%d\t%d\t%d\t%d\n", primary_perc, map_perc, primary_mapped, secondary, supplementary, unmapped, total_reads, total_records;
}
