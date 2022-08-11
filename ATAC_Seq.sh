# !/bin/bash
# Usage: bash ATAC_Seq.sh [SRRID]
# -------------------------------------------------
# 指定所有变量
data=$1
# bwa version: 0.7.17-r1188
# bwa index -p o.sativa ../all.chrs.con
Index=/root/genome/oryza_sativa/bwaIndex/o.sativa
effective_genome_size=373128865
SrrName=$(basename ${data})
SrrPath=$(realpath ${data}|sed "s/\/$(basename ${data})//g")
# awk -vFS='[\t=;]' -vOFS="\t" '{if ($3=="gene") print $1,$4,$5,$10,$6,$7}' demo.gff3 > refseq.bed  
refseq=/root/genome/oryza_sativa/refseq.bed
# -------------------------------------------------

if [ -d "${SrrPath}/${SrrName%%.*}-result" ];then
echo "${SrrName%%.*} analysis has been done"
exit 1
else
echo "${SrrName%%.*} analysis start"
fi


if [ -a ${SrrPath}/${SrrName} ];then
# fasterq-dump.3.0.0

fasterq-dump -p --include-technical -S -e 30 -O ${SrrPath} ${SrrPath}/${SrrName}
echo "($(date)) Step1: ${SrrName%%.*} split end"

# ------------------------------------------------------------------
if [ "$(ls ${SrrPath}/${SrrName%%.*}_* | wc -l)" == "2" ];then
# fastp version: 0.22.0
fastp --thread 16 \
-i ${SrrPath}/${SrrName%%.*}_1.fastq \
-o ${SrrPath}/${SrrName%%.*}_1.clean.fastq \
-I ${SrrPath}/${SrrName%%.*}_2.fastq \
-O ${SrrPath}/${SrrName%%.*}_2.clean.fastq
echo "($(date)) Step2: ${SrrName%%.*} filter end"

total_reads=$(head -17 fastp.json | tail -1 | cut -d ":" -f2 | sed 's/,//g')
q30_rate=$(echo "scale=2;100*$(head -22 fastp.json | tail -1 | cut -d ":" -f2 | sed 's/,//g')/1"|bc)%
read_mean_length=$(head -23 fastp.json | tail -1 | cut -d ":" -f2 | sed 's/,//g')bp
gc_content=$(echo "scale=2;100*$(head -25 fastp.json | tail -1 | cut -d ":" -f2 | sed 's/,//g')/1"|bc)%
duplication_rate=$(echo "scale=2;100*$(head -36 fastp.json | tail -1 | cut -d ":" -f2 | sed 's/,//g')/1"|bc)%
rm fastp*

if [ "$(ls $Index*|wc -l)" == "5" ];then
# bwa version: 0.7.17-r1188
bwa mem -M -t 30 \
-R "@RG\tID:${SrrName%%.*}\tSM:${SrrName%%.*}\tLB:WXS\tPL:Illumina" \
${Index} \
${SrrPath}/${SrrName%%.*}_1.clean.fastq \
${SrrPath}/${SrrName%%.*}_2.clean.fastq \
| samtools sort -O bam -@ 30 -o ${SrrPath}/${SrrName%%.*}.raw.bam
rm ${SrrPath}/${SrrName%%.*}_1* ${SrrPath}/${SrrName%%.*}_2*
echo "($(date)) Step3: ${SrrName%%.*} alignment end"

alignment_rate=$(sambamba flagstat ${SrrPath}/${SrrName%%.*}.raw.bam |head -5|tail -1|cut -d "(" -f2|cut -d ":" -f1)

else
echo "Error: can't find Index, please specify full Index pathway in this script"
exit 1
fi

# sambamba 0.6.6
sambamba markdup -t 30 -r \
${SrrPath}/${SrrName%%.*}.raw.bam ${SrrPath}/${SrrName%%.*}.temp.bam
samtools view -h -f 2 -q 30 -@ 30 ${SrrPath}/${SrrName%%.*}.temp.bam \
| samtools sort -O bam -@ 30 -o ${SrrPath}/${SrrName%%.*}.bam
rm ${SrrPath}/${SrrName%%.*}.temp.bam* ${SrrPath}/${SrrName%%.*}.raw.bam
sambamba index ${SrrPath}/${SrrName%%.*}.bam
echo "($(date)) Step4: ${SrrName%%.*} remove duplicates end"

final_reads=$(sambamba flagstat ${SrrPath}/${SrrName%%.*}.bam |head -5|tail -1|cut -d " " -f1)

picard CollectInsertSizeMetrics \
I=${SrrPath}/${SrrName%%.*}.bam \
O=${SrrPath}/${SrrName%%.*}.insertsize \
H=${SrrPath}/${SrrName%%.*}.hist.pdf
rm ${SrrPath}/${SrrName%%.*}.insertsize
echo "($(date)) Step5: ${SrrName%%.*} collect insertSize end"

# ------------------------------------------------------------------
elif  [ "$(ls ${SrrName%%.*}.* | wc -l)" == "1" ];then

# fastp version: 0.22.0
fastp --thread 16 \
-i ${SrrPath}/${SrrName%%.*}.fastq \
-o ${SrrPath}/${SrrName%%.*}.clean.fastq
echo "($(date)) Step2: ${SrrName%%.*} filter end"

total_reads=$(head -16 fastp.json | tail -1 | cut -d ":" -f2 | sed 's/,//g')
q30_rate=$(echo "scale=2;100*$(head -21 fastp.json | tail -1 | cut -d ":" -f2 | sed 's/,//g')/1"|bc)%
read_mean_length=$(head -22 fastp.json | tail -1 | cut -d ":" -f2 | sed 's/,//g')bp
gc_content=$(echo "scale=2;100*$(head -23 fastp.json | tail -1 | cut -d ":" -f2 | sed 's/,//g')/1"|bc)%
duplication_rate=$(echo "scale=2;100*$(head -34 fastp.json | tail -1 | cut -d ":" -f2 | sed 's/,//g')/1"|bc)%
rm fastp*

if [ "$(ls $Index*|wc -l)" == "5" ];then

bwa mem -M -t 30 \
-R "@RG\tID:${SrrName%%.*}\tSM:${SrrName%%.*}\tLB:WXS\tPL:Illumina" \
${Index} \
${SrrPath}/${SrrName%%.*}.clean.fastq \
| samtools sort -O bam -@ 30 -o ${SrrPath}/${SrrName%%.*}.raw.bam
rm ${SrrPath}/${SrrName%%.*}.clean.fastq ${SrrPath}/${SrrName%%.*}.fastq
echo "($(date)) Step3: ${SrrName%%.*} alignment end"

alignment_rate=$(sambamba flagstat ${SrrPath}/${SrrName%%.*}.raw.bam |head -5|tail -1|cut -d "(" -f2|cut -d ":" -f1)

else
echo "Error: can't find Index, please specify full Index pathway in this script"
exit 1
fi

# sambamba 0.6.6
sambamba markdup -t 30 -r \
${SrrPath}/${SrrName%%.*}.raw.bam ${SrrPath}/${SrrName%%.*}.temp.bam
samtools view -h -q 30 -@ 30 ${SrrPath}/${SrrName%%.*}.temp.bam \
| samtools sort -O bam -@ 30 -o ${SrrPath}/${SrrName%%.*}.bam
rm ${SrrPath}/${SrrName%%.*}.temp.bam* ${SrrPath}/${SrrName%%.*}.raw.bam
sambamba index ${SrrPath}/${SrrName%%.*}.bam
echo "($(date)) Step4: ${SrrName%%.*} remove duplicates end"

final_reads=$(sambamba flagstat ${SrrPath}/${SrrName%%.*}.bam |head -5|tail -1|cut -d " " -f1)

echo "($(date)) Step5: ${SrrName%%.*} is single library, not have insertSize"

else
echo "Error: can't recognize ${SrrName%%.*} splite file, maybe it is single cell sample"
exit 1
fi

# ------------------------------------------------------------------
# bedtools v2.30.0
bedtools bamtobed -i ${SrrPath}/${SrrName%%.*}.bam > ${SrrPath}/${SrrName%%.*}.bed
echo "($(date)) Step6: ${SrrName%%.*} bam to bed end"

# macs2 v2.2.7.1
macs2 callpeak \
-t ${SrrPath}/${SrrName%%.*}.bed \
-g ${effective_genome_size} \
-n ${SrrName%%.*} \
--nomodel --shift -100 --extsize 200 \
-q 0.01 --outdir ${SrrPath}
echo "($(date)) Step7: ${SrrName%%.*} call peak end"

# bc
narrowPeakReads=$(bedtools intersect -a ${SrrPath}/${SrrName%%.*}.bed -b ${SrrPath}/${SrrName%%.*}_peaks.narrowPeak |wc -l|awk '{print $1}')
FRiP_score=$(echo "scale=2;100*${narrowPeakReads}/${final_reads}"|bc)%
rm ${SrrPath}/${SrrName%%.*}.bed

bamCoverage -p 30 --binSize 10 \
--normalizeTo1x ${effective_genome_size} \
-b ${SrrPath}/${SrrName%%.*}.bam \
-o ${SrrPath}/${SrrName%%.*}.bw

echo "($(date)) Step8: ${SrrName%%.*}  bam to bw end"

computeMatrix scale-regions -p 30 --skipZeros \
-b 3000 -a 3000 \
-R ${refseq} \
-S ${SrrPath}/${SrrName%%.*}.bw \
-o ${SrrPath}/${SrrName%%.*}.gz

plotProfile -m ${SrrPath}/${SrrName%%.*}.gz \
--plotWidth 7 --plotHeight 7 \
-out ${SrrPath}/${SrrName%%.*}.pdf --plotFileFormat pdf
rm ${SrrPath}/${SrrName%%.*}.gz
echo "($(date)) Step9: ${SrrName%%.*} plot Profile end"

else
echo "Error: can't find ${SrrName} file, please retry"
exit 1
fi

mkdir ${SrrPath}/${SrrName%%.*}-result
mv ${SrrPath}/${SrrName%%.*}.* ${SrrPath}/${SrrName%%.*}-result
mv ${SrrPath}/${SrrName%%.*}_* ${SrrPath}/${SrrName%%.*}-result

echo -e "${SrrName%%.*}\t${total_reads}\t${q30_rate}\t${read_mean_length}\t${gc_content}\t${duplication_rate}\t${alignment_rate}\t${final_reads}\t${FRiP_score}" > ${SrrPath}/${SrrName%%.*}-result/${SrrName%%.*}.Statistics
echo -e "${SrrName%%.*}\t${total_reads}\t${q30_rate}\t${read_mean_length}\t${gc_content}\t${duplication_rate}\t${alignment_rate}\t${final_reads}\t${FRiP_score}"
