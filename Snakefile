configfile: "config/config.json"

shell.prefix("set +o pipefail; ")
shell.prefix("set -x;")
import pandas as pd

barcodes = pd.read_csv("s0RawData/barcode", sep="\t") #?
barcodes["Sample"] = barcodes["Sample"].astype(str)
samples = barcodes.Sample.tolist()
read_order = [1,2]

rule all:
	input:
		"s6TidyData/ASV_samplCounts.csv",
		"s6TidyData/ASV_taxanomy.csv"

rule Demultiplex_fwd:
	input:
		barcode_file = config["barcode"],
		first_fqfile = expand("{fastq}", fastq=config["raw_fastq"])
	output:
		expand("s1Demultiplex/first.{sample}.r1.fastq", sample=samples),
		expand("s1Demultiplex/first.{sample}.r2.fastq", sample=samples),
		expand("s1Demultiplex/first.unmatched.r{read_order}.fastq", read_order=read_order),
		expand("s1Demultiplex/first.Sample.r{read_order}.fastq", read_order=read_order)
	params:
		first_r1 = "s1Demultiplex/first.%.r1.fastq",
		first_r2 = "s1Demultiplex/first.%.r2.fastq"
	shell:
		"""
		fastq-multx -m 0 -x -b \\
			-B {input.barcode_file} \\
			{input.first_fqfile} \\
			-o {params.first_r1} \\
			-o {params.first_r2}
		"""

rule Demultiplex_rev:
	input:
		barcode_file = config["barcode"],
		second_r2 = "s1Demultiplex/first.unmatched.r2.fastq",
		second_r1 = "s1Demultiplex/first.unmatched.r1.fastq"
	output:
		expand("s1Demultiplex/second.{sample}.r1.fastq", sample=samples),
		expand("s1Demultiplex/second.{sample}.r2.fastq", sample=samples),
		expand("s1Demultiplex/second.unmatched.r{read_order}.fastq", read_order=read_order),
		expand("s1Demultiplex/second.Sample.r{read_order}.fastq", read_order=read_order)
	params:
		second_r2 = "s1Demultiplex/second.%.r2.fastq",
		second_r1 = "s1Demultiplex/second.%.r1.fastq"
	shell:
		"""
		fastq-multx -m 0 -x -b \\
			-B {input.barcode_file} \\
			{input.second_r2} \\
			{input.second_r1} \\
			-o {params.second_r2} \\
			-o {params.second_r1}
		"""

def get_barcode(wildcards):
	return barcodes.iloc[samples.index(wildcards.sample)]["BarcodeSequence"]

rule CutAdapt:
	input:
		barcode_file = config["barcode"],
		first_r1 = "s1Demultiplex/first.{sample}.r1.fastq",
		first_r2 = "s1Demultiplex/first.{sample}.r2.fastq",
		second_r2 = "s1Demultiplex/second.{sample}.r2.fastq",
		second_r1 = "s1Demultiplex/second.{sample}.r1.fastq"
	output:
		first_r1 = "s2CutAdapt/first.{sample}.r1.fastq",
		first_r2 = "s2CutAdapt/first.{sample}.r2.fastq", 
		second_r2 = "s2CutAdapt/second.{sample}.r2.fastq", 
		second_r1 = "s2CutAdapt/second.{sample}.r1.fastq"
	params:
		barcode = get_barcode,
		primer_fwd = config["primers"]["fwd"],
		primer_rev = config["primers"]["rev"]
	threads:
		config["cutadapt"]["threads"]
	shell:
		"""
		adapterFwd={params.barcode}{params.primer_fwd}
		adapterRev={params.primer_rev}
		cutadapt -e 0.15 \\
			-g $adapterFwd \\
			-G $adapterRev \\
			-o {output.first_r1} --discard-untrimmed \\
			-p {output.first_r2} --discard-untrimmed \\
			{input.first_r1} \\
			{input.first_r2} \\
			-j {threads}
		cutadapt -e 0.15 \\
			-g $adapterFwd -G $adapterRev \\
			-o {output.second_r2} --discard-untrimmed \\
			-p {output.second_r1} --discard-untrimmed \\
			{input.second_r2} \\
			{input.second_r1} \\
			-j {threads}
		"""

rule concatenate:
	input:
		first_r1 = "s2CutAdapt/first.{sample}.r1.fastq",
		second_r1 = "s2CutAdapt/second.{sample}.r1.fastq",
		first_r2 = "s2CutAdapt/first.{sample}.r2.fastq",
		second_r2 = "s2CutAdapt/second.{sample}.r2.fastq"
	output:
		r1 = "s3Combine/{sample}.r1.fastq",
		r2 = "s3Combine/{sample}.r2.fastq"
	shell:
		"""
		# set -x
		cat {input.first_r1} {input.second_r1} > {output.r1}
		cat {input.first_r2} {input.second_r2} > {output.r2}
		"""

rule fastQC:
	input:
		r1 = expand("s3Combine/{sample}.r1.fastq", sample=samples),
		r2 = expand("s3Combine/{sample}.r2.fastq", sample=samples)
	output:
		# directory("s4FastQC"),
		zipfiles = expand("s4FastQC/r{read_order}_fastqc.zip", read_order=read_order),
		r1 = temp("s4FastQC/r1.fastq.gz"),
		r2 = temp("s4FastQC/r2.fastq.gz")
	shell:
		"""
		cat {input.r1}|gzip -c > s4FastQC/r1.fastq.gz
		cat {input.r2}|gzip -c > s4FastQC/r2.fastq.gz
		fastqc -o s4FastQC s4FastQC/r1.fastq.gz s4FastQC/r2.fastq.gz
		"""

rule fastqcScore:
	input:
		zipfiles = expand("s4FastQC/r{read_order}_fastqc.zip", read_order=read_order)
	output:
		directory("s4FastQC/r1_fastqc"),
		directory("s4FastQC/r2_fastqc"),
		temp('cutLens.csv')
	script:
		"scripts/fastqcscore.py"

rule DownloadRefDB:
	input:
		"urls"
	output:
		config["reference_genome"]["genome_file"]
	params:
		path = config["reference_genome"]["genome_file_path"]
	shell:
		"""
		wget --directory-prefix={params} -i {input} 
		"""

rule DADA2:
	input: 
		cutLens = "cutLens.csv"
	output:
		temp("s5DADA2/DADA2_seqtab_nochim.rda"),
		"s5DADA2/img/ploterrF.png",
		"s5DADA2/img/ploterrR.png"
	params:
		path = "s3Combine"
	script:
		"scripts/dada2.R"

rule taxonomy:
	input:
		"s5DADA2/DADA2_seqtab_nochim.rda",
		config["reference_genome"]["genome_file"]
	output:
		temp("s5DADA2/DADA2_taxaid.rda")
	script:
		"scripts/taxonomy.R"

rule dada2Cleaning:
	input:
		"s5DADA2/DADA2_seqtab_nochim.rda",
		"s5DADA2/DADA2_taxaid.rda"
	output:
		"s6TidyData/ASV_samplCounts.csv",
		"s6TidyData/ASV_taxanomy.csv"
	script:
		"scripts/dada2Cleaning.R"
