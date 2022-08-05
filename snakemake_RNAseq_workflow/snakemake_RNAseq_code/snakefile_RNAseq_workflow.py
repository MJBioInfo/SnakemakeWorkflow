import glob
import os

#genome = ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_32/GRCh38.primary_assembly.genome.fa.gz
#gtf = ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_32/gencode.v32.annotation.gtf.gz

IDS,REPS = glob_wildcards("RawReads/{id}_{rep}.fastq.gz")


rule all:
	input:
		expand("QC/rawfastqc/{id}_{rep}_fastqc.{format}", id=IDS, rep=[1,2], format=["html","zip"]),
		expand("QC/trimming/{id}.log", id=IDS),
		expand("Trimmed/{id}_{direction}_{paired}.fastq.gz", id=IDS, direction=["forward","reverse"], paired=["paired","unpaired"]),
		expand("QC/trimmedqc/{id}_{direction}_paired_fastqc.{format}", id=IDS, direction=["forward","reverse"], format=["html","zip"]),
		expand("starOut/{id}{extension}", id=IDS, extension=["Unmapped.out.mate1","Unmapped.out.mate2","Aligned.sortedByCoord.out.bam"]),
		expand("rRNAContam/{id}_{ext}", id=IDS, ext=["rrna.bam","rrna.out","rrna.sam","unmapped.bam"]),
		expand("QC/star/{id}Log.final.out", id=IDS),
		expand("QC/rRNA/{id}_rrna.out", id=IDS),
        expand("QC/markduplicate_QC/{id}AlignedSortedMarkedMetrices.txt", id=IDS),
        expand("QC/markduplicate_QC/{id}.log", id=IDS),
        expand("QC/Rseqc/infer/{id}", id=IDS),
        expand("starAlignedMarkedSortedByName/{id}{extn}",  id=IDS, extn=["AlignedSortedMarkedByName.bam"]),
        expand("Featurecounts/{id}.txt", id=IDS)



rule fastqc:
	input:
		"RawReads/{id}_{rep}.fastq.gz"
	output:
		zip="QC/rawfastqc/{id}_{rep}_fastqc.zip",
		html="QC/rawfastqc/{id}_{rep}_fastqc.html"
	params:
		path="QC/rawfastqc/"
	threads: 10
	shell:
		"fastqc {input} --threads {threads} -O {params}"

rule trimmomatic:
	input:
		read1="RawReads/{id}_1.fastq.gz",
		read2="RawReads/{id}_2.fastq.gz"
	output:
		fp="Trimmed/{id}_forward_paired.fastq.gz",
		fu="Trimmed/{id}_forward_unpaired.fastq.gz",
		rp="Trimmed/{id}_reverse_paired.fastq.gz",
		ru="Trimmed/{id}_reverse_unpaired.fastq.gz"
	log:
		"QC/trimming/{id}.log"
	threads: 4
	shell:
		"trimmomatic PE -threads {threads} {input.read1} {input.read2} {output.fp} {output.fu} {output.rp} {output.ru} ILLUMINACLIP:adapters.fa:2:30:10:2:keepBothReads SLIDINGWINDOW:4:20 TRAILING:3 MINLEN:36 2>{log}"

rule trimmedQC:
	input:
		"Trimmed/{id}_{direction}_paired.fastq.gz"
	output:
		zip="QC/trimmedqc/{id}_{direction}_paired_fastqc.zip",
		html="QC/trimmedqc/{id}_{direction}_paired_fastqc.html"
	params:
		path="QC/trimmedqc/"
	threads: 10
	shell:
		"fastqc {input} --threads {threads} -O {params}"

rule rRNAContamination:
	input:
		trimmed1="Trimmed/{id}_forward_paired.fastq.gz",
		trimmed2="Trimmed/{id}_reverse_paired.fastq.gz"
	output:
		rnaSam="rRNAContam/{id}_rrna.sam"
	params:
		rna="RNAindex/human_rRNA_db.fasta"
	threads: 10
	shell:
		"""
		bwa mem -t {threads} -P {params} {input.trimmed1} {input.trimmed2} > {output.rnaSam}
		"""

rule rRNAContaminationConversion:
	input:
		rnaSam="rRNAContam/{id}_rrna.sam"
	output:
		rnaBam="rRNAContam/{id}_rrna.bam",
		rnaOut="rRNAContam/{id}_rrna.out",
		rnaUnm="rRNAContam/{id}_unmapped.bam"
	threads: 10
	shell:
		"""
		samtools view -@ {threads} -bS -o {output.rnaBam} {input.rnaSam}
		samtools flagstat -@ {threads} {output.rnaBam} > {output.rnaOut}
		samtools view -@ {threads} -u -f 12 -F 256 {output.rnaBam} > {output.rnaUnm}
		"""

rule rRNAfreeFastQ:
	input:
		rnaUnm="rRNAContam/{id}_unmapped.bam"
	output:
		fwd="rRNAfreeTrimmed/{id}_forward.fastq",
		rvs="rRNAfreeTrimmed/{id}_reverse.fastq"
	threads: 10
	shell:
		"""
		bamToFastq -i {input} -fq {output.fwd} -fq2 {output.rvs}
		"""

rule starAlignment:
	input:
		trimmed1="rRNAfreeTrimmed/{id}_forward.fastq",
		trimmed2="rRNAfreeTrimmed/{id}_reverse.fastq"
	output:
		"starOut/{id}Unmapped.out.mate1",
		"starOut/{id}Unmapped.out.mate2",
		"starOut/{id}Aligned.sortedByCoord.out.bam",
		"starOut/{id}Log.final.out"
	params:
		prefix="starOut/{id}"
	threads: 20
	shell:
		"""
		STAR --runThreadN {threads} --genomeLoad LoadAndKeep --genomeDir genomeIndex --readFilesIn {input.trimmed1} {input.trimmed2} --outFilterIntronMotifs RemoveNoncanonical --outFileNamePrefix {params.prefix} --limitBAMsortRAM 15000000000 --quantMode GeneCounts --outSAMtype BAM SortedByCoordinate  --outReadsUnmapped Fastx



        """


rule mark_duplicates:
    input:
        starbam="starOut/{id}Aligned.sortedByCoord.out.bam"


    output:
        markedbam="markDuplicate/{id}AlignedSortedMarked.bam",
        metrics="QC/markduplicate_QC/{id}AlignedSortedMarkedMetrices.txt"

    threads:
        10

    log:
        "QC/markduplicate_QC/{id}.log"

    shell:
        "picard MarkDuplicates -I {input.starbam} -O {output.markedbam}  -M {output.metrics}"



rule starAlignedMarkedSortedByName:
    input:
        markedbam_input="markDuplicate/{id}AlignedSortedMarked.bam"

    output:
        "starAlignedMarkedSortedByName/{id}AlignedSortedMarkedByName.bam"

    shell:
        """
        samtools sort  -O bam {input.markedbam_input} > {output}
        samtools index {output}


        """




rule rseqcstrandness:
    input:
        markedbam="markDuplicate/{id}AlignedSortedMarked.bam"
    output:
        "QC/Rseqc/infer/{id}"
    shell:
        "infer_experiment.py -i {input.markedbam} -r genome/gencode.v38.primary_assembly.annotation.bed > {output}"




rule featureCounts:
    input:
        bamForCount="starAlignedMarkedSortedByName/{id}AlignedSortedMarkedByName.bam",
        gtf="genome/gencode.v38.primary_assembly.annotation.gtf"

    output:
        "Featurecounts/{id}.txt"

    threads:
        20
    shell:
        " featureCounts -p -B -T {threads} -a {input.gtf} -s 2 -o {output}  {input.bamForCount}"




rule copyToQC:
	input:
		star="starOut/{id}Log.final.out",
		rrna="rRNAContam/{id}_rrna.out"


	output:
		starout="QC/star/{id}Log.final.out",
		rrnaout="QC/rRNA/{id}_rrna.out"


	shell:
		"""
		cp {input.star} {output.starout}
		cp {input.rrna} {output.rrnaout}


		"""


rule multiqc:
    input:
        logfiles="QC/"

    output:
        "multiqc_report.{zip}"

    shell:
        "multiqc {input.logfiles}."
