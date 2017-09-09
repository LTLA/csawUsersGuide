###############################################################
# This script is design to obtain the BAM files from SRA files.
###############################################################

set -e
set -u

# <<ASSUMPTION>>: the relevant files have already been obtained from GEO:
#	Tiwari et al. = GSE25532 
#	Revila-i-Domingo et al. = GSE38046
#	Zhang et al. = GSE31233
#	Pal et al. = GSE43212

tiwari=(SRR074398.sra SRR074399.sra SRR074417.sra SRR074418.sra SRR074401.sra)
domingo=(SRR499714.sra SRR499715.sra SRR499716.sra SRR499717.sra SRR499729.sra SRR499730.sra SRR499731.sra SRR499732.sra SRR499733.sra)
zhang=(SRR330784.sra SRR330785.sra SRR330786.sra SRR330800.sra SRR330801.sra SRR330814.sra SRR330815.sra SRR330816.sra)
pal=(SRR642390.sra)

# <<ASSUMPTION>>: subread is installed.

subcmd=subread-align

# <<ASSUMPTION>>: mm10 indices have been built.
#	This can be done by making a FASTA file from a BSGenome object:
#
#	> bs <- BSGenome.Mmusculus.UCSC.mm10
#	> outfile <- "mm10.fa"
#	> for (chr in seqnames(bs)) {
#	+	  y <- getSeq(bs, names = chr, start = 1, end = length(bs[[chr]]))
#	+	  y <- DNAStringSet(y)
#	+	  names(y) <- chr
#	+	  writeXStringSet(filepath = outfile, y, append = TRUE)
#	+ }
#
#	... and then running "subread-buildindex mm10.fa -o mm10_index/mm10" in the shell.
#	Alternatively, you can use your own FASTA file.

genome=mm10_index/mm10

# <<ASSUMPTION>>: FixMateInformation and MarkDuplicates are available from the Picard suite.

markcmd=MarkDuplicates

# <<ASSUMPTION>>: samtools has been installed.

samcmd=samtools

# <<ASSUMPTION>>: fastq-dump (from NCBI's SRA toolkit) has been installed.

fqdcmd=fastq-dump

###############################################################
# Generating temporaries.

vtmp=`mktemp -d --tmpdir=.`
tmpfile=blah.txt

function renameBAM {
	mv $1.bam $2.bam
	mv $1.bam.bai $2.bam.bai
}

###############################################################
# Running all files through the pipeline.

for i in {1..4}
do 
	if [[ $i -eq 1 ]]; then
		allfiles=(${tiwari[@]})
		curphred=3
		conthresh=2
		pet=0
	elif [[ $i -eq 2 ]]; then
		allfiles=(${domingo[@]})
		curphred=3
		conthresh=2
		pet=0
	elif [[ $i -eq 3 ]]; then
		allfiles=(${zhang[@]})
		curphred=3
		conthresh=2
		pet=0
	elif [[ $i -eq 4 ]]; then
		allfiles=(${pal[@]})
		curphred=3
		conthresh=2
		pet=1
	fi
	
	# Aligning.
	aligncmd="${subcmd} -u -H -T 4 -P ${curphred} -m ${conthresh} --BAMoutput -i ${genome}"
	for sra in "${allfiles[@]}"
	do
		rawbam=temp.bam
		if [[ $pet -eq 0 ]]; then
			fastq-dump $sra
			prefix=`echo $sra | sed "s/\.sra//g"`
			fastq=`ls | egrep "$prefix\.(fastq|fq)"`
			${aligncmd} -r $fastq -o $rawbam
		else 
			fastq-dump --split-files $sra
			prefix=`echo $sra | sed "s/\.sra//g"`
			fastq1=`ls | egrep "${prefix}_1\.(fastq|fq)"`
			fastq2=`ls | egrep "${prefix}_2\.(fastq|fq)"`
			${aligncmd} -r $fastq1 -R $fastq2 -o $rawbam
		fi

		${samcmd} sort $rawbam sorted
		mv sorted.bam $rawbam
		finalbam=${prefix}.bam
		${markcmd} I=$rawbam O=$finalbam M=$tmpfile TMP_DIR=$vtmp AS=true REMOVE_DUPLICATES=false VALIDATION_STRINGENCY=SILENT
		${samcmd} index $finalbam

		rm $rawbam
		rm ${rawbam}.indel
		if [[ $pet -eq 0 ]]; then
			rm $fastq
		else
			rm $fastq1 $fastq2
		fi
	done

	# Merging technical replicates and renaming things.
	if [[ $i -eq 1 ]]; then
		renameBAM SRR074398 es_1
		renameBAM SRR074399 es_2
		renameBAM SRR074417 tn_1
		renameBAM SRR074418 tn_2
		renameBAM SRR074401 input
	elif [[ $i -eq 2 ]]; then
		${samcmd} merge -f h3k4me3_pro.bam SRR499716.bam SRR499717.bam
		${samcmd} index h3k4me3_pro.bam
		${samcmd} merge -f h3k4me3_pro_2.bam SRR499714.bam SRR499715.bam
		${samcmd} index h3k4me3_pro_2.bam
		${samcmd} merge -f h3k4me3_mat.bam SRR499732.bam SRR499733.bam
		${samcmd} index h3k4me3_mat.bam
		${samcmd} merge -f h3k4me3_mat_2.bam SRR499729.bam SRR499730.bam SRR499731.bam
		${samcmd} index h3k4me3_mat_2.bam
	elif [[ $i -eq 3 ]]; then
		${samcmd} merge -f h3ac.bam SRR330784.bam SRR330785.bam
		${samcmd} index h3ac.bam
		renameBAM SRR330786 h3ac_2
		${samcmd} merge -f h3k4me2.bam SRR330800.bam SRR330801.bam
		${samcmd} index h3k4me2.bam
		${samcmd} merge -f h3k27me3.bam SRR330814.bam SRR330815.bam SRR330816.bam
		${samcmd} index h3k27me3.bam
	elif [[ $i -eq 4 ]]; then
		renameBAM SRR642390 example-pet
	fi
done

###############################################################
# Mopping up.

rm $tmpfile
rm -r $vtmp

# Printing out a log file with version numbers.
ticket=success.log
if [[ -e $ticket ]]; then
        rm $ticket
fi

set +e
stored=`subread-align -v 2>&1 | sed "s/.*v/v/" | sed "s/\n//g"`
printf "$subcmd (" >> $ticket
printf $stored >> $ticket
printf ")\n" >> $ticket
stored=`$samcmd 2>&1 | grep -i "Version:"`
printf "Samtools $stored\n" >> $ticket
stored=`$markcmd --version 2>&1`
printf "MarkDuplicates version $stored\n" >> $ticket
$fqdcmd --version | grep "." >> $ticket

###############################################################
