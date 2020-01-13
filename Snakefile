configfile: "config/config.yml"

shell.prefix("set +o pipefail; ")
shell.prefix("set -x;")

import pandas as pd
import json
import os

# extract sample names, path, and barcode
barcodes = pd.read_csv(config["barcode"], sep="\t") #?
barcodes["Sample"] = barcodes["Sample"].astype(str)
samples = barcodes.Sample.tolist()
# 1 for R1 and 2 for R2
read_order = [1,2]
# save truncLen to disk if exist
truncLen = config.get("truncLen")
if not os.path.isdir("output"):
	os.mkdir("output")
if truncLen != None:
	with open("output/truncLens.json", "wt") as fh:
		json.dump(config["truncLen"], fh)

rule all:
	input:
		"output/s6TidyData/ASV_samplCounts.csv",
		"output/s6TidyData/ASV_taxanomy.csv"

rule Demultiplex_fwd:
	input:
		barcode_file = config["barcode"],
		first_fqfile = expand("{fastq}", fastq=config["raw_fastq"])
	output:
		expand("output/s1Demultiplex/first.{sample}.r1.fastq", sample=samples),
		expand("output/s1Demultiplex/first.{sample}.r2.fastq", sample=samples),
		expand("output/s1Demultiplex/first.unmatched.r{read_order}.fastq", read_order=read_order),
		expand("output/s1Demultiplex/first.Sample.r{read_order}.fastq", read_order=read_order)
	params:
		first_r1 = "output/s1Demultiplex/first.%.r1.fastq",
		first_r2 = "output/s1Demultiplex/first.%.r2.fastq"
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
		second_r2 = "output/s1Demultiplex/first.unmatched.r2.fastq",
		second_r1 = "output/s1Demultiplex/first.unmatched.r1.fastq"
	output:
		expand("output/s1Demultiplex/second.{sample}.r1.fastq", sample=samples),
		expand("output/s1Demultiplex/second.{sample}.r2.fastq", sample=samples),
		expand("output/s1Demultiplex/second.unmatched.r{read_order}.fastq", read_order=read_order),
		expand("output/s1Demultiplex/second.Sample.r{read_order}.fastq", read_order=read_order)
	params:
		second_r2 = "output/s1Demultiplex/second.%.r2.fastq",
		second_r1 = "output/s1Demultiplex/second.%.r1.fastq"
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
		first_r1 = "output/s1Demultiplex/first.{sample}.r1.fastq",
		first_r2 = "output/s1Demultiplex/first.{sample}.r2.fastq",
		second_r2 = "output/s1Demultiplex/second.{sample}.r2.fastq",
		second_r1 = "output/s1Demultiplex/second.{sample}.r1.fastq"
	output:
		first_r1 = "output/s2CutAdapt/first.{sample}.r1.fastq",
		first_r2 = "output/s2CutAdapt/first.{sample}.r2.fastq", 
		second_r2 = "output/s2CutAdapt/second.{sample}.r2.fastq", 
		second_r1 = "output/s2CutAdapt/second.{sample}.r1.fastq"
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
		first_r1 = "output/s2CutAdapt/first.{sample}.r1.fastq",
		second_r1 = "output/s2CutAdapt/second.{sample}.r1.fastq",
		first_r2 = "output/s2CutAdapt/first.{sample}.r2.fastq",
		second_r2 = "output/s2CutAdapt/second.{sample}.r2.fastq"
	output:
		r1 = "output/s3Combine/{sample}.r1.fastq",
		r2 = "output/s3Combine/{sample}.r2.fastq"
	shell:
		"""
		# set -x
		cat {input.first_r1} {input.second_r1} > {output.r1}
		cat {input.first_r2} {input.second_r2} > {output.r2}
		"""

rule fastQC:
	input:
		r1 = expand("output/s3Combine/{sample}.r1.fastq", sample=samples),
		r2 = expand("output/s3Combine/{sample}.r2.fastq", sample=samples)
	output:
		# directory("s4FastQC"),
		zipfiles = expand("output/s4FastQC/r{read_order}_fastqc.zip", read_order=read_order),
		r1 = temp("output/s4FastQC/r1.fastq.gz"),
		r2 = temp("output/s4FastQC/r2.fastq.gz")
	shell:
		"""
		cat {input.r1}|gzip -c > output/s4FastQC/r1.fastq.gz
		cat {input.r2}|gzip -c > output/s4FastQC/r2.fastq.gz
		fastqc -o output/s4FastQC output/s4FastQC/r1.fastq.gz output/s4FastQC/r2.fastq.gz
		"""

rule fastqcScore:
	input:
		zipfiles = expand("output/s4FastQC/r{read_order}_fastqc.zip", read_order=read_order)
	output:
		directory("output/s4FastQC/r1_fastqc"),
		directory("output/s4FastQC/r2_fastqc"),
		"output/truncLens.json"
	script:
		"scripts/fastqcscore.py"

rule DownloadRefDB:
	output:
		"database/SILVA_SSU.RData"
	params:
		url = config["reference_db"]["dada2_silva_url"]
	shell:
		"""
		wget {params.url} --output-document=database/SILVA_SSU.RData -q
		"""

def dada2_input():
	res = [expand("output/s4FastQC/r{read_order}_fastqc.zip", read_order=read_order)]
	if not config.get("truncLen"):
		res.append("output/truncLens.json")
	return res
		

rule DADA2:
	input: 
		dada2_input()
	output:
		temp("output/s5DADA2/DADA2_seqtab_nochim.rda"),
		"output/s5DADA2/img/ploterrF.png",
		"output/s5DADA2/img/ploterrR.png"
	params:
		path = "output/s3Combine"
	script:
		"scripts/dada2.R"

rule taxonomy:
	input:
		"output/s5DADA2/DADA2_seqtab_nochim.rda",
		"database/SILVA_SSU.RData"
	output:
		temp("output/s5DADA2/DADA2_taxaid.rda")
	script:
		"scripts/taxonomy.R"

rule dada2Cleaning:
	input:
		"output/s5DADA2/DADA2_seqtab_nochim.rda",
		"output/s5DADA2/DADA2_taxaid.rda"
	output:
		"output/s6TidyData/ASV_samplCounts.csv",
		"output/s6TidyData/ASV_taxanomy.csv"
	script:
		"scripts/dada2Cleaning.R"
