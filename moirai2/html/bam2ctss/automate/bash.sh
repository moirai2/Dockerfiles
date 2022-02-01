#$ -q sge
mkdir -p bam2ctss/output
perl moirai2.pl \
-c 'samtools.sif' \
-i '$id->bam->$input,$id->flag->$flag,$id->cutoff->$cutoff' \
-o '$id->ctssbed->$output' \
-f '$output' \
command/bam/bam_to_ctssbed.json \
'$output=bam2ctss/output/$id.ctss.bed'